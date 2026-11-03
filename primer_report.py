"""
Primer report - are the primers binding where we expect around the gene hits, and how well do they match?
Only primer hits that passed the parser cutoffs (max_mismatches, min_alignment_pct) are in here.
"""

import logging
from collections import defaultdict
from pathlib import Path
import pandas as pd
from blast_prep import active_queries
from reporting import format_workbook, isolate_sort_key, locus_index

logger = logging.getLogger(__name__)

last_bases = 5 #a difference this close to the 3' end can stop extension


def sorted_range(hit):
    return tuple(sorted((hit["s_start"], hit["s_end"])))


def locate_primer(hit, locus):
    #where the primer is compared to the gene locus, in gene orientation
    p_low, p_high = sorted_range(hit)
    l_low, l_high = sorted_range(locus)
    gap = max(0, p_low - l_high, l_low - p_high)

    if gap == 0:
        return {"distance": 0, "relation": "inside gene", "faces_gene": "overlaps gene"}

    left_of_locus = p_high < l_low
    gene_is_forward = locus["strand"] == "forward"
    relation = "upstream of gene" if left_of_locus == gene_is_forward else "downstream of gene"

    primer_points_right = hit["strand"] == "forward"
    faces_gene = "toward gene" if primer_points_right == left_of_locus else "away from gene"

    return {"distance": gap, "relation": relation, "faces_gene": faces_gene}


def nearest_locus(hit, index):
    loci = index.get((hit["isolate_id"], hit["contig_id"]), [])
    if not loci:
        return None, None

    best = min(loci, key = lambda item: locate_primer(hit, item[2])["distance"])
    return best[2], locate_primer(hit, best[2])


def binding_class(locus, location, buffer):
    if locus is None:
        return "off-target: no gene locus on contig"
    if location["distance"] > buffer:
        return f"off-target: >{buffer} bp from gene"
    if location["relation"] == "inside gene":
        return "on gene"
    if location["faces_gene"] == "toward gene":
        return "flanking, facing gene"
    return "flanking, facing away"


def make_primer_hits_table(hits, index, buffer):
    rows = []
    for hit in sorted(hits, key = lambda h: (isolate_sort_key(h["isolate_id"]), h["query_id"], h["total_mismatches"])):
        locus, location = nearest_locus(hit, index)
        from_3prime = [int(p) for p in hit["positions_from_3prime"].split(";") if p]

        rows.append({
            "id": hit["isolate_id"],
            "primer": hit["query_id"],
            "role": hit["position"],
            "contig": hit.get("contig_orig"),
            "strand": hit["strand"],
            "subject_start": hit["s_start"],
            "subject_end": hit["s_end"],
            "total_differences": hit["total_mismatches"],
            "blast_mismatches": hit["blast_mismatches"],
            "unaligned_bases": hit["unaligned_bases"],
            "gaps": hit["gaps"],
            "mismatch_positions_from_5prime": hit["mismatch_positions"],
            "mismatch_positions_from_3prime": hit["positions_from_3prime"],
            f"difference_in_last_{last_bases}_bases": any(p <= last_bases for p in from_3prime),
            "binding": binding_class(locus, location, buffer),
            "nearest_locus_query": locus["query_id"] if locus else None,
            "locus_strand": locus["strand"] if locus else None,
            "distance_to_locus_bp": location["distance"] if location else None,
            "position_vs_gene": location["relation"] if location else None,
            "primer_direction": location["faces_gene"] if location else None,
            "qseq": hit["qseq"],
            "sseq": hit["sseq"],
        })

    return pd.DataFrame(rows)


def hits_for_primer(hits_df, name):
    if hits_df.empty:
        return hits_df
    return hits_df[hits_df["primer"] == name]


def make_primer_summary(hits_df, config, isolate_ids):
    last_col = f"difference_in_last_{last_bases}_bases"
    rows = []

    for name, data in active_queries(config.primer_queries).items():
        df = hits_for_primer(hits_df, name)
        n = len(df)

        row = {
            "primer": name,
            "role": data["position"],
            "sequence": data["query_sequence"],
            "length_bp": len(data["query_sequence"]),
            "max_differences_allowed": config.primer_blast_params["max_mismatches"],
            "isolates_searched": len(set(isolate_ids)),
            "isolates_with_hit": df["id"].nunique() if n else 0,
            "passing_hits": n,
            "perfect_hits": int((df["total_differences"] == 0).sum()) if n else 0,
            f"hits_with_difference_in_last_{last_bases}_bases": int(df[last_col].sum()) if n else 0,
            "mean_differences": round(df["total_differences"].mean(), 2) if n else None,
            "max_differences_seen": int(df["total_differences"].max()) if n else None,
        }

        for label in ("on gene", "flanking, facing gene", "flanking, facing away"):
            row[f"hits_{label}"] = int((df["binding"] == label).sum()) if n else 0
        row["hits_off-target"] = int(df["binding"].str.startswith("off-target").sum()) if n else 0
        row["hits_per_isolate_max"] = int(df.groupby("id").size().max()) if n else 0

        rows.append(row)

    return pd.DataFrame(rows)


def make_mismatch_distribution(hits_df, config):
    rows = []

    for name in active_queries(config.primer_queries):
        df = hits_for_primer(hits_df, name)
        if df.empty:
            continue

        for differences, count in df["total_differences"].value_counts().sort_index().items():
            rows.append({
                "primer": name,
                "total_differences": int(differences),
                "hits": int(count),
                "percent_of_hits": round(count / len(df) * 100, 2)
            })

    return pd.DataFrame(rows)


def make_position_table(hits_df, config):
    #how often each primer base differs, position 1 = 5' end
    rows = []

    for name, data in active_queries(config.primer_queries).items():
        df = hits_for_primer(hits_df, name)
        length = len(data["query_sequence"])

        counts = defaultdict(int)
        for value in df["mismatch_positions_from_5prime"] if len(df) else []:
            for position in str(value).split(";"):
                if position:
                    counts[int(position)] += 1

        for position in range(1, length + 1):
            rows.append({
                "primer": name,
                "position_from_5prime": position,
                "position_from_3prime": length - position + 1,
                "base": data["query_sequence"][position - 1],
                "hits_with_difference": counts[position],
                "percent_of_hits": round(counts[position] / len(df) * 100, 2) if len(df) else 0,
            })

    return pd.DataFrame(rows)


def make_locus_table(gene_loci, primer_hits, config):
    #for each gene locus, the best forward and reverse primer hit nearby and whether the pair faces each other
    buffer = config.flank_buffer
    max_fragment = config.primer_blast_params.get("max_pcr_fragment_length")

    by_contig = defaultdict(list)
    for hit in primer_hits:
        by_contig[(hit["isolate_id"], hit["contig_id"])].append(hit)

    rows = []
    for locus in sorted(gene_loci, key = lambda h: (isolate_sort_key(h["isolate_id"]), h["contig_id"], min(h["s_start"], h["s_end"]))):
        l_low, l_high = sorted_range(locus)
        row = {
            "id": locus["isolate_id"],
            "contig": locus.get("contig_orig"),
            "strand": locus["strand"],
            "representative_query": locus.get("query_id"),
            "genomic_context": locus.get("genomic_context"),
            "locus_start": l_low,
            "locus_end": l_high,
        }

        best = {}
        for role in ("forward", "reverse"):
            candidates = []
            for hit in by_contig.get((locus["isolate_id"], locus["contig_id"]), []):
                if hit["position"] != role:
                    continue
                location = locate_primer(hit, locus)
                if location["distance"] <= buffer:
                    candidates.append((hit["total_mismatches"], location["distance"], hit, location))

            if candidates:
                _, _, hit, location = min(candidates, key = lambda c: (c[0], c[1]))
                best[role] = hit
                row[f"{role}_primer_differences"] = hit["total_mismatches"]
                row[f"{role}_primer_position_vs_gene"] = location["relation"]
                row[f"{role}_primer_direction"] = location["faces_gene"]
                row[f"{role}_primer_distance_bp"] = location["distance"]
            else:
                row[f"{role}_primer_differences"] = None
                row[f"{role}_primer_position_vs_gene"] = f"no hit within {buffer} bp"
                row[f"{role}_primer_direction"] = None
                row[f"{role}_primer_distance_bp"] = None

        row["primers_face_each_other"] = None
        row["primer_span_bp"] = None
        row["span_covers_aligned_span"] = None
        row["span_within_max_pcr_fragment_length"] = None

        #only checks the geometry - no PCR product is inferred
        if len(best) == 2:
            left, right = sorted(best.values(), key = lambda h: sorted_range(h)[0])
            if left["strand"] == "forward" and right["strand"] == "reverse" and sorted_range(left)[1] < sorted_range(right)[0]:
                span = sorted_range(right)[1] - sorted_range(left)[0] + 1
                row["primers_face_each_other"] = True
                row["primer_span_bp"] = span
                row["span_covers_aligned_span"] = sorted_range(left)[0] <= l_low and sorted_range(right)[1] >= l_high
                row["span_within_max_pcr_fragment_length"] = span <= max_fragment if max_fragment else None
            else:
                row["primers_face_each_other"] = False

        rows.append(row)

    return pd.DataFrame(rows)


def create_primer_report(blast_results, config, isolate_ids, output_dir):
    forward = blast_results.get("forward_primer_hits", [])
    reverse = blast_results.get("reverse_primer_hits", [])
    primer_hits = forward + reverse
    buffer = config.flank_buffer
    gene_loci = blast_results.get("gene_hits_before_filter", [])

    hits_df = make_primer_hits_table(primer_hits, locus_index(gene_loci), buffer)
    summary_df = make_primer_summary(hits_df, config, isolate_ids)

    notes = pd.DataFrame({"Notes": [
        "Only primer hits passing the parser cutoffs are listed; hits with more differences than max_mismatches are not counted.",
        "total_differences = BLAST mismatches + gaps + unaligned primer bases.",
        f"Positions are numbered from the primer 5' end (position_from_5prime) or 3' end; the last {last_bases} bases are flagged because differences there can block extension.",
        "Gene locus = merged BLAST-aligned gene hit (may be a partial gene); 'inside gene' means the primer overlaps that span.",
        f"Primers more than flank_buffer ({buffer} bp) from every locus, or on a contig with no locus, are called off-target.",
        "position_vs_gene is in gene orientation; primer_direction says whether the primer extends toward the gene.",
        "Per-locus span columns describe primer geometry only; no PCR product is inferred or extracted.",
    ]})

    output_path = Path(output_dir) / f"{config.gene_name}_primer_report.xlsx"
    with pd.ExcelWriter(output_path, engine = "openpyxl") as writer:
        summary_df.to_excel(writer, sheet_name = "Primer summary", index = False)
        make_mismatch_distribution(hits_df, config).to_excel(writer, sheet_name = "Differences per hit", index = False)
        make_position_table(hits_df, config).to_excel(writer, sheet_name = "Differences by position", index = False)
        make_locus_table(gene_loci, primer_hits, config).to_excel(writer, sheet_name = "Primers per locus", index = False)
        hits_df.to_excel(writer, sheet_name = "Primer hits", index = False)
        notes.to_excel(writer, sheet_name = "Notes", index = False)
        format_workbook(writer)
        writer.book["Notes"].auto_filter.ref = None

    logger.info("Primer report saved to %s", output_path)
    return output_path
