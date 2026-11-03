import logging
from collections import defaultdict
from pathlib import Path

logger = logging.getLogger(__name__)

gene_columns = ["query_id", "subject_id", "num_ident", "pct_ident", "aln_length", "mismatches", "gaps", "gap_opens", "q_start", "q_end", "s_start", "s_end", "evalue", "bitscore", "qlength", "slength"]
primer_columns = gene_columns + ["qseq", "sseq"]


def parse_contig_id(subject_id, contig_id_map):
    if "_ctg" in subject_id:
        isolate_id, contig_id = subject_id.rsplit("_ctg", 1)
        contig_local = f"ctg{contig_id}"
    else:
        isolate_id = "unknown"
        contig_local = subject_id

    return {
        "isolate_id": isolate_id,
        "contig_id": subject_id,
        "contig_local": contig_local,
        "contig_orig": contig_id_map.get(subject_id, subject_id)
    }


def read_blast_rows(blast_path, columns):
    #read a BLAST output file and check its column count
    if blast_path is None:
        return

    with Path(blast_path).open("r", encoding="utf-8") as blast_file:
        for line_number, line in enumerate(blast_file, start=1):
            if not line.strip() or line.startswith("#"):
                continue

            fields = line.rstrip("\r\n").split("\t")
            if len(fields) != len(columns):
                raise ValueError(
                    f"{blast_path}, line {line_number}: expected "
                    f"{len(columns)} columns, found {len(fields)}"
                )

            yield dict(zip(columns, fields))


def create_hit(row, contig_id_map):
    #convert the shared BLAST fields into a hit dictionary
    query_length = int(row["qlength"])
    query_start = int(row["q_start"])
    query_end = int(row["q_end"])
    subject_start = int(row["s_start"])
    subject_end = int(row["s_end"])

    if query_length <= 0:
        raise ValueError(f"Invalid query length for {row['query_id']}")

    return {
        "query_id": row["query_id"],
        **parse_contig_id(row["subject_id"], contig_id_map),
        "identity": float(row["pct_ident"]),
        "coverage": (abs(query_end - query_start) + 1) / query_length * 100,
        "evalue": float(row["evalue"]),
        "bitscore": float(row["bitscore"]),
        "aln_length": int(row["aln_length"]),
        "q_length": query_length,
        "q_start": query_start,
        "q_end": query_end,
        "s_length": int(row["slength"]),
        "s_start": subject_start,
        "s_end": subject_end,
        "strand": "forward" if (query_end - query_start) * (subject_end - subject_start) >= 0 else "reverse",
        "mismatches": int(row["mismatches"]),
        "gaps": int(row["gaps"]),
        "gap_opens": int(row["gap_opens"])
    }


def parse_gene_blast(blast_path, config, contig_id_map, query_type="gene", stats=None):
    #parse gene, flank, gene-start or gene-end results
    #if given, stats collects raw alignment rows and rows passing thresholds for reporting
    if stats is not None:
        stats.update({"raw_rows": 0, "passing_rows": 0, "raw_isolates": set(), "passing_isolates": set(), "raw_points": []})
    params = getattr(config, f"{query_type}_blast_params")
    queries = getattr(config, f"{query_type}_queries")
    hits = []

    high_identity = params.get("high_identity", params["min_identity"])
    high_coverage = params.get("high_coverage", params["min_coverage"])

    for row in read_blast_rows(blast_path, gene_columns):
        query_id = row["query_id"]
        if query_id not in queries:
            raise ValueError(f"Unknown {query_type} query in BLAST results: {query_id}")

        hit = create_hit(row, contig_id_map)
        passed = hit["identity"] >= params["min_identity"] and hit["coverage"] >= params["min_coverage"]
        if stats is not None:
            stats["raw_rows"] += 1
            stats["raw_isolates"].add(hit["isolate_id"])
            stats["raw_points"].append({"identity": hit["identity"], "coverage": hit["coverage"], "passed": passed})
            if passed:
                stats["passing_rows"] += 1
                stats["passing_isolates"].add(hit["isolate_id"])
        if not passed:
            continue

        hit["query_type"] = query_type
        hit["high_quality"] = hit["identity"] >= high_identity and hit["coverage"] >= high_coverage

        if query_type == "flank":
            query_data = queries[query_id]
            position = query_data["position"]
            expected_distance = query_data["distance"]

            if position not in ("upstream", "downstream"):
                raise ValueError(f"Invalid flank position for {query_id}: {position}")
            #negative distances mean the flank overlaps the gene by that many bp
            if not isinstance(expected_distance, (int, float)):
                raise ValueError(f"Flank distance must be a number for {query_id}")

            hit["position"] = position
            hit["expected_distance"] = expected_distance

        hits.append(hit)

    logger.info("Parsed %d %s hits", len(hits), query_type)
    return hits


def get_primer_mismatch_positions(query_sequence, subject_sequence, query_start, query_length):
    #report substitutions, subject gaps and unaligned primer bases
    query_position = query_start - 1
    mismatch_positions = []

    for query_base, subject_base in zip(query_sequence.upper(), subject_sequence.upper()):
        if query_base == "-":
            continue

        query_position += 1
        if subject_base == "-" or query_base != subject_base:
            mismatch_positions.append(query_position)

    mismatch_positions.extend(range(1, query_start))
    mismatch_positions.extend(range(query_position + 1, query_length + 1))
    mismatch_positions = sorted(set(mismatch_positions))
    positions_from_3prime = sorted(query_length - pos + 1 for pos in mismatch_positions)

    return mismatch_positions, positions_from_3prime


def parse_primer_blast(blast_path, config, contig_id_map):
    params = config.primer_blast_params
    forward_hits = []
    reverse_hits = []
    complement = str.maketrans("ACGTRYSWKMBDHVN", "TGCAYRSWMKVHDBN")

    for row in read_blast_rows(blast_path, primer_columns):
        query_id = row["query_id"]
        if query_id not in config.primer_queries:
            raise ValueError(f"Unknown primer query in BLAST results: {query_id}")

        position = config.primer_queries[query_id]["position"]
        if position not in ("forward", "reverse"):
            raise ValueError(f"Invalid primer position for {query_id}: {position}")

        hit = create_hit(row, contig_id_map)
        aligned_bases = abs(hit["q_end"] - hit["q_start"]) + 1
        unaligned_bases = hit["q_length"] - aligned_bases
        total_mismatches = hit["mismatches"] + unaligned_bases + hit["gaps"]
        min_coverage = max(params["min_coverage"], params["min_alignment_pct"])

        if hit["identity"] < params["min_identity"] or hit["coverage"] < min_coverage:
            continue
        if total_mismatches > params["max_mismatches"]:
            continue

        query_sequence = row["qseq"].upper()
        subject_sequence = row["sseq"].upper()
        position_start = hit["q_start"]

        if hit["q_start"] > hit["q_end"]:
            query_sequence = query_sequence.translate(complement)[::-1]
            subject_sequence = subject_sequence.translate(complement)[::-1]
            position_start = hit["q_end"]

        mismatch_positions, positions_from_3prime = get_primer_mismatch_positions(
            query_sequence, subject_sequence, position_start, hit["q_length"]
        )

        hit.update({
            "position": position,
            "qseq": row["qseq"],
            "sseq": row["sseq"],
            "total_mismatches": total_mismatches,
            "blast_mismatches": hit["mismatches"],
            "unaligned_bases": unaligned_bases,
            "mismatch_positions": ";".join(map(str, mismatch_positions)),
            "positions_from_3prime": ";".join(map(str, positions_from_3prime))
        })

        if position == "forward":
            forward_hits.append(hit)
        else:
            reverse_hits.append(hit)

    logger.info("Parsed %d forward hits and %d reverse hits", len(forward_hits), len(reverse_hits))
    return forward_hits, reverse_hits


def hit_priority(hit):
    return (
        int(hit.get("high_quality", False)),
        hit["bitscore"],
        hit["coverage"],
        hit["identity"],
        -hit["evalue"]
    )


def merge_gene_variants(hits, gene_name, overlap_fraction=0.80):
    #keep one representative hit per overlapping candidate locus
    def same_locus(hit_a, hit_b):
        if (
            hit_a["isolate_id"] != hit_b["isolate_id"]
            or hit_a["contig_id"] != hit_b["contig_id"]
            or hit_a["strand"] != hit_b["strand"]
        ):
            return False

        a_start, a_end = sorted((hit_a["s_start"], hit_a["s_end"]))
        b_start, b_end = sorted((hit_b["s_start"], hit_b["s_end"]))
        overlap = max(0, min(a_end, b_end) - max(a_start, b_start) + 1)
        shorter_length = min(a_end - a_start + 1, b_end - b_start + 1)
        return overlap > 0 and overlap / shorter_length >= overlap_fraction

    merged_by_contig = defaultdict(list)

    for original_hit in sorted(hits, key=hit_priority, reverse=True):
        candidate = original_hit.copy()
        candidate["blast_query"] = candidate["query_id"]
        candidate["supporting_queries"] = candidate["query_id"]
        candidate["n_query_variants"] = 1
        candidate["query_hit_details"] = [original_hit.copy()]

        key = (candidate["isolate_id"], candidate["contig_id"], candidate["strand"])
        matching_hit = next(
            (hit for hit in merged_by_contig[key] if same_locus(candidate, hit)),
            None
        )

        if matching_hit is None:
            merged_by_contig[key].append(candidate)
            continue

        matching_hit["query_hit_details"].extend(candidate["query_hit_details"])
        supporting_queries = sorted({
            detail["query_id"] for detail in matching_hit["query_hit_details"]
        })
        matching_hit["supporting_queries"] = ";".join(supporting_queries)
        matching_hit["n_query_variants"] = len(supporting_queries)

    merged_hits = [hit for group in merged_by_contig.values() for hit in group]
    logger.info("Merged %d BLAST hits into %d candidate %s loci", len(hits), len(merged_hits), gene_name)
    return merged_hits


def validate_genomic_context(gene_hits, flank_hits, config):
    #check flank position and deviation from its expected distance
    if config.flank_buffer < 0:
        raise ValueError("flank_buffer must be non-negative")

    flanks_by_contig = defaultdict(list)
    for flank in flank_hits:
        flanks_by_contig[(flank["isolate_id"], flank["contig_id"])].append(flank)

    for gene_hit in gene_hits:
        for position in ("upstream", "downstream"):
            gene_hit[f"{position}_found"] = False
            gene_hit[f"{position}_distance"] = None
            gene_hit[f"{position}_flank_coords"] = None
            gene_hit[f"{position}_flank_query"] = None

        gene_hit["genomic_context"] = "no_flanks" if config.flank_queries else "not_checked"
        gene_start, gene_end = sorted((gene_hit["s_start"], gene_hit["s_end"]))
        best_matches = {}
        nearby_flanks = flanks_by_contig.get(
            (gene_hit["isolate_id"], gene_hit["contig_id"]), []
        )

        for flank in nearby_flanks:
            flank_start, flank_end = sorted((flank["s_start"], flank["s_end"]))
            position = flank["position"]
            is_left = (
                position == "upstream" and gene_hit["strand"] == "forward"
            ) or (
                position == "downstream" and gene_hit["strand"] == "reverse"
            )

            if is_left:
                #overlap is allowed (negative distance) but the flank must extend beyond the gene
                if flank_start >= gene_start:
                    continue
                distance = gene_start - flank_end - 1
            else:
                if flank_end <= gene_end:
                    continue
                distance = flank_start - gene_end - 1

            deviation = abs(distance - flank["expected_distance"])
            if deviation > config.flank_buffer:
                continue

            rank = (deviation, -flank["bitscore"], -flank["identity"], -flank["coverage"])
            if position in best_matches and rank >= best_matches[position]:
                continue

            best_matches[position] = rank
            gene_hit[f"{position}_found"] = True
            gene_hit[f"{position}_distance"] = distance
            gene_hit[f"{position}_flank_coords"] = (flank_start - 1, flank_end)
            gene_hit[f"{position}_flank_query"] = flank["query_id"]

        if gene_hit["upstream_found"] and gene_hit["downstream_found"]:
            gene_hit["genomic_context"] = "both_flanks"
        elif gene_hit["upstream_found"]:
            gene_hit["genomic_context"] = "upstream_only"
        elif gene_hit["downstream_found"]:
            gene_hit["genomic_context"] = "downstream_only"

    context_counts = defaultdict(int)
    for hit in gene_hits:
        context_counts[hit["genomic_context"]] += 1

    logger.info("Genomic context counts: %s", dict(context_counts))
    return gene_hits


def filter_top_hits_per_isolate(hits, config, query_type="gene", per_query=False, assign_copies=False):
    #limit gene loci per isolate, or supporting hits per isolate and query
    params = getattr(config, f"{query_type}_blast_params")
    max_hits = params["max_hits_per_isolate"]
    if not isinstance(max_hits, int) or isinstance(max_hits, bool) or max_hits < 1:
        raise ValueError(f"{query_type} max_hits_per_isolate must be a positive integer")

    grouped_hits = defaultdict(list)
    for hit in hits:
        key = (hit["isolate_id"], hit["query_id"]) if per_query else hit["isolate_id"]
        grouped_hits[key].append(hit)

    filtered_hits = []
    for group_hits in grouped_hits.values():
        group_hits.sort(key=lambda hit: (-hit["identity"], -hit["coverage"], hit["evalue"], -hit["bitscore"]))

        for copy_number, original_hit in enumerate(group_hits[:max_hits], start=1):
            hit = original_hit.copy()
            if assign_copies:
                hit["gene_copy"] = f"copy_{copy_number}"
            filtered_hits.append(hit)

    logger.info("Retained %d %s hits after applying hit limits", len(filtered_hits), query_type)
    return filtered_hits


def parse_blast_outputs(blast_outputs, config, contig_id_map):
    #parse all available searches
    gene_stats = {}
    gene_hits = parse_gene_blast(
        blast_outputs.get("gene"), config, contig_id_map, query_type="gene", stats=gene_stats
    )
    flank_hits = parse_gene_blast(
        blast_outputs.get("flank"), config, contig_id_map, query_type="flank"
    )
    gene_start_hits = parse_gene_blast(
        blast_outputs.get("gene_start"), config, contig_id_map, query_type="gene_start"
    )
    gene_end_hits = parse_gene_blast(
        blast_outputs.get("gene_end"), config, contig_id_map, query_type="gene_end"
    )
    forward_hits, reverse_hits = parse_primer_blast(
        blast_outputs.get("primers"), config, contig_id_map
    )

    flank_hits = filter_top_hits_per_isolate(
        flank_hits, config, query_type="flank", per_query=True
    )
    gene_start_hits = filter_top_hits_per_isolate(
        gene_start_hits, config, query_type="gene_start", per_query=True
    )
    gene_end_hits = filter_top_hits_per_isolate(
        gene_end_hits, config, query_type="gene_end", per_query=True
    )
    forward_hits = filter_top_hits_per_isolate(
        forward_hits, config, query_type="primer", per_query=True
    )
    reverse_hits = filter_top_hits_per_isolate(
        reverse_hits, config, query_type="primer", per_query=True
    )

    gene_hits = merge_gene_variants(gene_hits, config.gene_name)
    validated_gene_hits = validate_genomic_context(gene_hits, flank_hits, config)
    final_gene_hits = filter_top_hits_per_isolate(
        validated_gene_hits, config, query_type="gene", assign_copies=True
    )

    return {
        "gene_hits": final_gene_hits,
        "flank_hits": flank_hits,
        "upstream_flank_hits": [hit for hit in flank_hits if hit["position"] == "upstream"],
        "downstream_flank_hits": [hit for hit in flank_hits if hit["position"] == "downstream"],
        "gene_hits_before_filter": validated_gene_hits,
        "gene_start_hits": gene_start_hits,
        "gene_end_hits": gene_end_hits,
        "forward_primer_hits": forward_hits,
        "reverse_primer_hits": reverse_hits,
        "gene_stats": gene_stats,
    }
