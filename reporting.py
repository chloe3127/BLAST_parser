import logging
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any
import pandas as pd
from blast_parser import gene_columns, primer_columns
from blast_prep import active_queries

logger = logging.getLogger(__name__)

def write_df_with_spacing(df, writer, sheet_name, current_row, blank_rows_after = 2):
    #helper function to have spaces between tables on summary sheets
    df.to_excel(writer, sheet_name = sheet_name, index = False, startrow = current_row)

    return current_row + len(df) + 1 + blank_rows_after

def format_workbook(writer):
    for worksheet in writer.book.worksheets:
        worksheet.freeze_panes = "A2"
        worksheet.auto_filter.ref = worksheet.dimensions

        for column_cells in worksheet.columns:
            values = [str(cell.value) for cell in column_cells if cell.value is not None]
            width = min(max((len(value) for value in values), default = 10) + 2, 60)
            worksheet.column_dimensions[column_cells[0].column_letter].width = width

    
# 1 - BLAST reports

identity_bins = [
        ("<50%", lambda value: value < 50),
        ("50-60%", lambda value: 50 <= value < 60),
        ("60-70%", lambda value: 60 <= value < 70),
        ("70-80%", lambda value: 70 <= value < 80),
        ("80-90%", lambda value: 80 <= value < 90),
        ("90-95%", lambda value: 90 <= value < 95),
        ("95-98%", lambda value: 95 <= value < 98),
        ("98-100%", lambda value: 98 <= value <= 100),
        ("100%", lambda value: value == 100)
    ]
coverage_bins = [
        ("<50%", lambda value: value < 50),
        ("50-60%", lambda value: 50 <= value < 60),
        ("60-70%", lambda value: 60 <= value < 70),
        ("70-80%", lambda value: 70 <= value < 80),
        ("80-90%", lambda value: 80 <= value < 90),
        ("90-95%", lambda value: 90 <= value < 95),
        ("95-98%", lambda value: 95 <= value < 98),
        ("98-100%", lambda value: 98 <= value <= 100),
        ("100%", lambda value: value == 100)
    ]

def count_bins(values, bins):
        return[f"{label}: {sum(predicate(value) for value in values)}"
               for label, predicate in bins]
               
def read_blast_rows(blast_path, columns = None):
    #raw BLAST rows with identity and query coverage
    rows = []
    with Path(blast_path).open("r", encoding = "utf-8") as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\r\n").split("\t")
            columns = columns or gene_columns
            if len(fields) != len(columns):
                continue
            row = dict(zip(columns, fields))
            row["identity"] = float(row["pct_ident"])
            row["coverage"] = (abs(int(row["q_end"]) - int(row["q_start"])) + 1) / int(row["qlength"]) * 100
            rows.append(row)
    return rows

def write_text_report(report_lines, output_file):
    for line in report_lines: logger.info(line)
    output_file = Path(output_file)
    output_file.parent.mkdir(parents = True, exist_ok = True)
    output_file.write_text("\n".join(report_lines), encoding = "utf-8")

def report_gene_hit_stats(blast_path, output_file, queries):
    #summarise each hit result from raw gene BLAST output, split by reference query
    rows = read_blast_rows(blast_path)
    identities = [row["identity"] for row in rows]
    coverages = [row["coverage"] for row in rows]
    isolates = {row["subject_id"].split("_ctg", 1)[0] for row in rows}
    query_counts = Counter(row["query_id"] for row in rows)

    report_lines = [
        "Raw gene BLAST hit statistics",
        f"Source file: {blast_path}",
        f"Raw total hits: {len(rows)}",
        f"Isolates: {len(isolates)}",
        "",
        "Identity breakdown:",
        *count_bins(identities, identity_bins),
        "",
        "Coverage breakdown:",
        *count_bins(coverages, coverage_bins),
        "",
        "Query breakdown:",
        *[f"{query}: {query_counts.get(query, 0)}" for query in queries]
    ]
    write_text_report(report_lines, output_file)

def report_flank_hit_stats(blast_path, output_file, flank_queries):
    #same as report_gene_hit_stats but for the flanking gene hits - splits upstream and downstream
    flank_rows = {"upstream": [], "downstream": []}
    for row in read_blast_rows(blast_path):
        query = flank_queries.get(row["query_id"])
        if query is not None:
            flank_rows[query["position"]].append(row)

    report_lines = [
        "Raw flank BLAST hit statistics",
        f"Source file: {blast_path}",
        f"Raw total flank hits: {sum(len(rows) for rows in flank_rows.values())}",
        ""
    ]
    for position, rows in flank_rows.items():
        isolates = {row["subject_id"].split("_ctg", 1)[0] for row in rows}
        report_lines.extend([
            f"{position.title()} flank hits: {len(rows)}",
            f"{position.title()} flank isolates: {len(isolates)}",
            "Identity breakdown:",
            *count_bins([row["identity"] for row in rows], identity_bins),
            "Coverage breakdown:",
            *count_bins([row["coverage"] for row in rows], coverage_bins),
            ""
        ])
    write_text_report(report_lines, output_file)

def isolate_sort_key(isolate_id):
    isolate_id = str(isolate_id)
    return (0, int(isolate_id)) if isolate_id.isdigit() else (1, isolate_id)
    
def sort_report_df_by_id(df):
    #sorts reports by numeric id, if id has text they come after (some ids won't split properly and needs further debugging)
    id_column = next((column for column in ("id", "isolate_id") if column in df.columns), None)

    if id_column is None or df.empty:
        return df

    order = sorted(range(len(df)), key = lambda index: isolate_sort_key(df.iloc[index][id_column]))
    return df.iloc[order].reset_index(drop=True)

def clean_sequence_for_report(sequence):
    return "".join(str(sequence).split()).upper() if sequence else ""

def make_query_table(config):
    rows = []
    groups = [
        ("Gene reference", config.gene_queries),
        ("Gene start", config.gene_start_queries),
        ("Gene end", config.gene_end_queries),
        ("Flank", config.flank_queries),
        ("Primer", config.primer_queries)
    ]
    for group, queries in groups:
        for name, data in active_queries(queries).items():
            sequence = clean_sequence_for_report(data["query_sequence"])
            label = data.get("position")
            rows.append({
                "Query type": f"{group}: {name}" + (f" ({label})" if label else ""),
                "Sequence length (bp)": len(sequence),
                "Sequence": sequence
            })

    return pd.DataFrame(rows)

def make_parameter_table(config):
    rows = []
    searches = [
        ("Gene BLAST", config.gene_blast_params, config.gene_queries),
        ("Flank BLAST", config.flank_blast_params, config.flank_queries),
        ("Gene start BLAST", config.gene_start_blast_params, config.gene_start_queries),
        ("Gene end BLAST", config.gene_end_blast_params, config.gene_end_queries),
        ("Primer BLAST", config.primer_blast_params, config.primer_queries)
    ]
    for search, params, queries in searches:
        if not active_queries(queries):
            continue
        for name, value in params.items():
            rows.append({"Search": search, "Parameter": name, "Value": value})

    for name, data in active_queries(config.flank_queries).items():
        rows.append({"Search": "Genomic context", "Parameter": f"{name} expected distance ({data['position']})", "Value": data.get("distance")})
    rows.append({"Search": "Genomic context", "Parameter": "Flank buffer", "Value": config.flank_buffer})

    return pd.DataFrame(rows)


def make_isolate_summary(blast_results, isolate_ids, config, extracted_sequences = None, before_N_filtering = None):
    searches_run = blast_results.get("searches_run", set())
    gene_hits = blast_results.get("gene_hits", [])
    gene_start_hits = blast_results.get("gene_start_hits", [])
    gene_end_hits = blast_results.get("gene_end_hits", [])
    forward_primer_hits = blast_results.get("forward_primer_hits", [])
    reverse_primer_hits = blast_results.get("reverse_primer_hits", [])
    contig_counts = blast_results.get("contig_counts", {})

    def by_isolate(hits):
        grouped = defaultdict(list)
        for hit in hits:
            grouped[str(hit["isolate_id"])].append(hit)
        return grouped

    candidates = by_isolate(gene_hits)
    extracted_before = by_isolate(before_N_filtering or [])
    retained = by_isolate(extracted_sequences or [])
    gene_start = by_isolate(gene_start_hits)
    gene_end = by_isolate(gene_end_hits)
    locus_hits = locus_index(gene_hits)
    forward = by_isolate(forward_primer_hits)
    reverse = by_isolate(reverse_primer_hits)

    #skipped searches say "not searched" instead of 0 hits
    def count_or_skipped(grouped, isolate_id, searched):
        return len(grouped.get(isolate_id, [])) if searched else "not searched"

    rows = []
    for isolate_id in sorted({str(i) for i in isolate_ids}, key = isolate_sort_key):
        isolate_hits = candidates.get(isolate_id, [])
        contexts = sorted({h.get("genomic_context", "unknown") for h in isolate_hits})
        n_candidates = len(isolate_hits)
        n_extracted = len(extracted_before.get(isolate_id, []))
        n_retained = len(retained.get(isolate_id, []))

        if not n_candidates:
            outcome = "no candidate locus"
        elif not n_retained:
            outcome = "all candidate loci lost in extraction/QC"
        elif n_retained < n_candidates:
            outcome = "some candidate loci lost in extraction/QC"
        else:
            outcome = "all candidate loci retained"

        row = {
            "id": isolate_id,
            "n_contigs": contig_counts.get(isolate_id),
            "gene_found": bool(n_candidates),
            "n_candidate_loci": n_candidates,
            "n_extracted_before_N_filter": n_extracted,
            "n_retained_loci": n_retained,
            "locus_outcome": outcome,
            "n_high_quality": sum(bool(h.get("high_quality")) for h in isolate_hits),
            "genomic_context": ";".join(contexts) if contexts else "unknown"
        }
        if active_queries(config.gene_start_queries):
            row["n_gene_start_hits"] = count_or_skipped(gene_start, isolate_id, "gene_start" in searches_run)
        if active_queries(config.gene_end_queries):
            row["n_gene_end_hits"] = count_or_skipped(gene_end, isolate_id, "gene_end" in searches_run)
        row["n_forward_primer_hits"] = count_or_skipped(forward, isolate_id, "primers" in searches_run)
        row["n_reverse_primer_hits"] = count_or_skipped(reverse, isolate_id, "primers" in searches_run)
        if "primers" in searches_run:
            #primer hits that overlap a candidate gene locus
            for label, grouped in (("forward", forward), ("reverse", reverse)):
                row[f"n_{label}_primer_hits_on_locus"] = sum(distance_to_locus(h, locus_hits) == 0 for h in grouped.get(isolate_id, []))
        else:
            row["n_forward_primer_hits_on_locus"] = row["n_reverse_primer_hits_on_locus"] = "not searched"
        rows.append(row)

    return pd.DataFrame(rows)

search_types = [
    ("gene", "Gene", "gene_queries", "gene_blast_params"),
    ("flank", "Flank", "flank_queries", "flank_blast_params"),
    ("gene_start", "Gene start", "gene_start_queries", "gene_start_blast_params"),
    ("gene_end", "Gene end", "gene_end_queries", "gene_end_blast_params"),
    ("primers", "Primer", "primer_queries", "primer_blast_params")
]

def passes_thresholds(row, search, params):
    #same cutoffs as the parser (primers use coverage and mismatches)
    if search == "primers":
        aligned = abs(int(row["q_end"]) - int(row["q_start"])) + 1
        total = int(row["mismatches"]) + int(row["gaps"]) + int(row["qlength"]) - aligned
        return row["coverage"] >= params["min_alignment_pct"] and total <= params["max_mismatches"]
    return row["identity"] >= params["min_identity"] and row["coverage"] >= params["min_coverage"]

def read_search_rows(blast_outputs, config):
    #raw rows for every search that was run, with whether they pass the thresholds
    rows = []
    for search, label, queries_attr, params_attr in search_types:
        path = blast_outputs.get(search)
        if path is None:
            continue
        params = getattr(config, params_attr)
        columns = primer_columns if search == "primers" else gene_columns
        for row in read_blast_rows(path, columns):
            row["search"] = label
            row["passed"] = passes_thresholds(row, search, params)
            row["isolate_id"] = row["subject_id"].split("_ctg", 1)[0]
            rows.append(row)
    return rows

def make_query_summary(search_rows, blast_results, config):
    #one row per query: raw and passing rows, and isolates
    representative = Counter(h.get("blast_query") for h in blast_results.get("gene_hits_before_filter", []))
    retained = Counter(h.get("blast_query") for h in blast_results.get("gene_hits", []))
    rows = []
    for search, label, queries_attr, _ in search_types:
        for name in active_queries(getattr(config, queries_attr)):
            raw = [r for r in search_rows if r["search"] == label and r["query_id"] == name]
            passing = [r for r in raw if r["passed"]]
            row = {
                "Search": label,
                "Query": name,
                "Raw alignment rows": len(raw),
                "Raw isolates": len({r["isolate_id"] for r in raw}),
                "Rows passing thresholds": len(passing),
                "Passing isolates": len({r["isolate_id"] for r in passing}),
                "Median identity (passing)": round(pd.Series([r["identity"] for r in passing]).median(), 2) if passing else None,
                "Median coverage (passing)": round(pd.Series([r["coverage"] for r in passing]).median(), 2) if passing else None
            }
            if label == "Gene":
                row["Loci where query is representative"] = representative.get(name, 0)
                row["Representative loci retained after cap"] = retained.get(name, 0)
            rows.append(row)
    return pd.DataFrame(rows)

def make_hits_by_query_table(search_rows):
    #every raw row so each query can be filtered on its own
    rows = [{
        "Search": r["search"],
        "Query": r["query_id"],
        "id": r["isolate_id"],
        "contig": r["subject_id"],
        "identity": round(r["identity"], 2),
        "coverage": round(r["coverage"], 2),
        "aln_length": int(r["aln_length"]),
        "s_start": int(r["s_start"]),
        "s_end": int(r["s_end"]),
        "strand": "forward" if int(r["s_start"]) <= int(r["s_end"]) else "reverse",
        "evalue": r["evalue"],
        "bitscore": r["bitscore"],
        "passed_thresholds": "Yes" if r["passed"] else "No"
    } for r in search_rows if r["search"] != "Primer"]
    return sort_report_df_by_id(pd.DataFrame(rows))

def make_isolate_query_matrix(search_rows, isolate_ids, config):
    #isolates x gene queries, number of passing rows
    queries = list(active_queries(config.gene_queries))
    passing = defaultdict(lambda: defaultdict(list))
    for r in search_rows:
        if r["search"] == "Gene" and r["passed"]:
            passing[r["isolate_id"]][r["query_id"]].append(r)

    rows = []
    for isolate_id in sorted({str(i) for i in isolate_ids}, key = isolate_sort_key):
        per_query = passing.get(isolate_id, {})
        row = {"id": isolate_id}
        for q in queries:
            row[q] = len(per_query.get(q, []))
        best = max(
            (r for hits in per_query.values() for r in hits),
            key = lambda r: (r["identity"], r["coverage"]), default = None
        )
        row["best_query"] = best["query_id"] if best else "no passing hit"
        row["best_identity"] = round(best["identity"], 2) if best else None
        row["best_coverage"] = round(best["coverage"], 2) if best else None
        rows.append(row)
    return pd.DataFrame(rows)

def flank_key(isolate_id, contig_id, query_id, coords):
    return (str(isolate_id), contig_id, query_id, tuple(coords) if coords else None)

def make_gene_flank_table(gene_hits_before_filter, retained_hits, config, forward_primer_hits = (), reverse_primer_hits = ()):
    #every merged gene locus with its flank hits and nearby primer hits
    retained = {(h["isolate_id"], h["contig_id"], h["s_start"], h["s_end"], h["strand"]): h.get("gene_copy") for h in retained_hits}
    rows = []
    for hit in gene_hits_before_filter:
        row = {
            "id": hit["isolate_id"],
            "gene_copy": retained.get((hit["isolate_id"], hit["contig_id"], hit["s_start"], hit["s_end"], hit["strand"])),
            "kept_after_cap": "Yes" if (hit["isolate_id"], hit["contig_id"], hit["s_start"], hit["s_end"], hit["strand"]) in retained else "No",
            "contig": hit.get("contig_orig"),
            "strand": hit.get("strand"),
            "s_start": hit["s_start"],
            "s_end": hit["s_end"],
            "representative_query": hit.get("blast_query"),
            "supporting_queries": (hit.get("supporting_queries") or "").replace(";", "; "),
            "identity": round(hit["identity"], 2),
            "coverage": round(hit["coverage"], 2),
            "reference_query_length": hit.get("q_length"),
            "gene_length": abs(hit["s_end"] - hit["s_start"]) + 1,
            "genomic_context": hit.get("genomic_context")
        }
        row["span_fraction_of_reference"] = round(row["gene_length"] / row["reference_query_length"], 3) if row["reference_query_length"] else None
        if forward_primer_hits or reverse_primer_hits:
            counts = primer_counts_for_locus(hit, forward_primer_hits, reverse_primer_hits, config.flank_buffer)
            row.update({
                "forward_primer_hits_on_locus": counts["forward"][0],
                "reverse_primer_hits_on_locus": counts["reverse"][0],
                f"forward_primer_hits_within_{config.flank_buffer}bp": counts["forward"][1],
                f"reverse_primer_hits_within_{config.flank_buffer}bp": counts["reverse"][1]
            })
        for position in ("upstream", "downstream"):
            query = hit.get(f"{position}_flank_query")
            expected = config.flank_queries.get(query, {}).get("distance") if query else None
            distance = hit.get(f"{position}_distance")
            coords = hit.get(f"{position}_flank_coords")
            row.update({
                f"{position}_flank_query": query,
                f"{position}_distance": distance,
                f"{position}_expected_distance": expected,
                f"{position}_deviation": abs(distance - expected) if distance is not None and expected is not None else None,
                f"{position}_flank_start": coords[0] if coords else None,
                f"{position}_flank_end": coords[1] if coords else None
            })
        rows.append(row)
    return sort_report_df_by_id(pd.DataFrame(rows))

def make_totals_by_query(search_rows, blast_results, config):
    #raw rows, passing rows and hits kept after the per-isolate cap, per query
    kept = {
        "Flank": blast_results.get("flank_hits", []),
        "Gene start": blast_results.get("gene_start_hits", []),
        "Gene end": blast_results.get("gene_end_hits", []),
        "Primer": blast_results.get("forward_primer_hits", []) + blast_results.get("reverse_primer_hits", [])
    }
    rows = []
    for search, label, queries_attr, _ in search_types:
        if search not in blast_results.get("searches_run", {search}):
            continue
        for name in active_queries(getattr(config, queries_attr)):
            raw = [r for r in search_rows if r["search"] == label and r["query_id"] == name]
            rows.append({
                "Search": label,
                "Query": name,
                "Raw alignment rows": len(raw),
                "Rows passing thresholds": sum(r["passed"] for r in raw),
                "Hits kept after per-isolate cap": sum(h["query_id"] == name for h in kept[label]) if label in kept else None
            })
    df = pd.DataFrame(rows)
    for label, group in df.groupby("Search", sort = False):
        total = group["Rows passing thresholds"].sum()
        df.loc[group.index, "% of search's passing rows"] = (group["Rows passing thresholds"] / total * 100).round(2) if total else 0
    summary = [
        ("Gene", "Merged candidate loci (all queries)", None, None, len(blast_results.get("gene_hits_before_filter", []))),
        ("Gene", "Loci after per-isolate cap", None, None, len(blast_results.get("gene_hits", [])))
    ]
    extra = pd.DataFrame(summary, columns = ["Search", "Query", "Raw alignment rows", "Rows passing thresholds", "Hits kept after per-isolate cap"])
    return pd.concat([df, extra], ignore_index = True)

def locus_index(gene_hits):
    #gene loci grouped by isolate and contig
    index = defaultdict(list)
    for hit in gene_hits:
        low, high = sorted((hit["s_start"], hit["s_end"]))
        index[(hit["isolate_id"], hit["contig_id"])].append((low, high, hit))
    return index

def distance_to_locus(primer_hit, index):
    #0 if the primer overlaps a locus, None if the contig has no locus
    loci = index.get((primer_hit["isolate_id"], primer_hit["contig_id"]), [])
    if not loci:
        return None
    low, high = sorted((primer_hit["s_start"], primer_hit["s_end"]))
    return min(max(0, low - locus_high, locus_low - high) for locus_low, locus_high, _ in loci)

def primer_counts_for_locus(locus, forward_hits, reverse_hits, buffer):
    #forward/reverse primer hits on the locus and within flank_buffer of it
    low, high = sorted((locus["s_start"], locus["s_end"]))
    counts = {}
    for label, hits in (("forward", forward_hits), ("reverse", reverse_hits)):
        on = near = 0
        for hit in hits:
            if hit["isolate_id"] != locus["isolate_id"] or hit["contig_id"] != locus["contig_id"]:
                continue
            p_low, p_high = sorted((hit["s_start"], hit["s_end"]))
            gap = max(0, p_low - high, low - p_high)
            on += gap == 0
            near += gap <= buffer
        counts[label] = (on, near)
    return counts

def flank_unassigned_reason(flank, index, buffer):
    #why a flank hit wasn't assigned to a locus
    loci = index.get((flank["isolate_id"], flank["contig_id"]), [])
    if not loci:
        return "no gene locus on this contig"
    f_low, f_high = sorted((flank["s_start"], flank["s_end"]))
    reasons = set()
    for g_low, g_high, gene in loci:
        is_left = (flank["position"] == "upstream") == (gene["strand"] == "forward")
        if is_left:
            if f_low >= g_low:
                reasons.add("not beyond gene on expected side")
                continue
            distance = g_low - f_high - 1
        else:
            if f_high <= g_high:
                reasons.add("not beyond gene on expected side")
                continue
            distance = f_low - g_high - 1
        if abs(distance - flank["expected_distance"]) > buffer:
            reasons.add("outside flank_buffer of expected distance")
        else:
            reasons.add("a better-ranked flank was assigned")
    for reason in ("a better-ranked flank was assigned", "outside flank_buffer of expected distance", "not beyond gene on expected side"):
        if reason in reasons:
            return reason
    return "not assigned"

def make_stage_summary(search_rows, blast_results, isolate_ids, extracted_sequences, before_N_filtering):
    #counts at each stage, with the unit since isolates, rows and loci differ
    gene_rows_all = [r for r in search_rows if r["search"] == "Gene"]
    passing = [r for r in gene_rows_all if r["passed"]]
    merged = blast_results.get("gene_hits_before_filter", [])
    capped = blast_results.get("gene_hits", [])
    n_before = len(before_N_filtering or [])
    n_retained = len(extracted_sequences or [])

    def isolates(items, key = "isolate_id"):
        return len({str(item[key]) for item in items})

    rows = [
        ("Isolates searched", "isolates", len(set(map(str, isolate_ids))), len(set(map(str, isolate_ids)))),
        ("Raw gene BLAST alignment rows (all reference queries)", "alignment rows", len(gene_rows_all), isolates(gene_rows_all)),
        ("Rows passing parser thresholds", "alignment rows", len(passing), isolates(passing)),
        ("Merged candidate loci", "loci", len(merged), isolates(merged)),
        ("Loci after per-isolate cap", "loci", len(capped), isolates(capped)),
        ("Extracted before N filtering", "loci", n_before, isolates(before_N_filtering or [])),
        ("Excluded by N filter", "loci (exclusion)", n_before - n_retained, isolates(before_N_filtering or []) - isolates(extracted_sequences or [])),
        ("Retained after N filtering", "loci", n_retained, isolates(extracted_sequences or []))
    ]
    df = pd.DataFrame(rows, columns = ["Stage", "Count unit", "Count", "Unique isolates"])
    not_extracted = len(capped) - n_before
    df.loc[len(df)] = ["Loci not extracted (contig not found)", "loci (exclusion)", not_extracted, None]
    return df

def make_copy_number_table(isolate_df):
    #how many isolates have 0, 1, 2... loci
    candidates = isolate_df["Candidate Loci"]
    retained = isolate_df["Retained Loci"]
    rows = [{
        "Loci per isolate": n,
        "Isolates (candidate loci)": int((candidates == n).sum()),
        "Isolates (retained loci)": int((retained == n).sum())
    } for n in range(0, int(max(candidates.max(), retained.max(), 1)) + 1)]
    return pd.DataFrame(rows)

def make_gene_hits_table(gene_hits):
    rows = []

    for hit in sorted(gene_hits, key = lambda h: (isolate_sort_key(h["isolate_id"]), h.get("gene_copy", ""))):
        rows.append({
            "id": hit["isolate_id"],
            "gene_copy": hit.get("gene_copy"),
            "n_gene_copies": hit.get("n_gene_copies"),
            "contig": hit.get("contig_orig"),
            "identity": round(hit.get("identity"), 2),
            "coverage": round(hit.get("coverage"), 2),
            "high_quality": bool(hit.get("high_quality")),
            "reference": hit.get("blast_query"),
            "genomic_context": hit.get("genomic_context"),
            "upstream_found": bool(hit.get("upstream_found")),
            "downstream_found": bool(hit.get("downstream_found")),
            "upstream_distance": hit.get("upstream_distance"),
            "downstream_distance": hit.get("downstream_distance"),
            "subject_start": hit.get("s_start"),
            "subject_end": hit.get("s_end"),
            "strand": hit.get("strand"),
            "aln_length": hit.get("aln_length"),
            "mismatches": hit.get("mismatches"),
            "gaps": hit.get("gaps"),
            "gap_opens": hit.get("gap_opens"),
            "evalue": hit.get("evalue"),
            "bitscore": hit.get("bitscore")
        })

    return pd.DataFrame(rows)

def make_gene_start_table(gene_start_hits):
    rows = []

    for hit in sorted(gene_start_hits, key = lambda h: (isolate_sort_key(h["isolate_id"]), h.get("gene_copy", ""))):
        aligned_query_bases = abs(hit["q_end"] - hit["q_start"]) + 1
        unaligned_bases = hit["q_length"] - aligned_query_bases
        total_differences = hit["mismatches"] + hit["gaps"] + hit["gap_opens"] + unaligned_bases

        rows.append({
            "id": hit["isolate_id"],
            "contig": hit.get("contig_orig"),
            "identity": round(hit.get("identity"), 2),
            "coverage": round(hit.get("coverage"), 2),
            "high_quality": bool(hit.get("high_quality")),
            "aln_length": hit.get("aln_length"),
            "blast_mismatches": hit.get("mismatches"),
            "gaps": hit.get("gaps"),
            "gap_opens": hit.get("gap_opens"),
            "unaligned_bases": unaligned_bases,
            "total_differences": total_differences,
            "perfect_match": total_differences == 0,
            "subject_start": hit.get("s_start"),
            "subject_end": hit.get("s_end"),
            "strand": hit.get("strand"),
            "evalue": hit.get("evalue"),
            "bitscore": hit.get("bitscore")
        })

    return pd.DataFrame(rows)

def make_gene_assembly_table(blast_results, config):
    rows = []

    gene_hits = blast_results.get("gene_hits_before_filter", blast_results.get("gene_hits", []))
    gene_start_hits = blast_results.get("gene_start_hits", [])
    gene_end_hits = blast_results.get("gene_end_hits", [])

    def add_hit(isolate_id, region, hit, detail = None):
        detail = detail or hit
        rows.append({
            "id": str(isolate_id),
            "ref_region": region,
            "query": detail.get("query_id", region),
            "contig": detail.get("contig_orig"),
            "subject_start": detail.get("s_start"),
            "subject_end": detail.get("s_end"),
            "identity": round(detail["identity"], 2) if detail.get("identity") is not None else None,
            "coverage": round(detail["coverage"], 2) if detail.get("coverage") is not None else None,
            "aln_length": detail.get("aln_length"),
            "strand": detail.get("strand"),
            "evalue": detail.get("evalue"),
            "bitscore": detail.get("bitscore")
        })

    for hit in gene_hits:
        details = hit.get("query_hit_details") or [hit]
        for detail in details:
            add_hit(hit["isolate_id"], "Gene", hit, detail)

    for hit in gene_start_hits:
        add_hit(hit["isolate_id"], "Gene start", hit)

    for hit in gene_end_hits:
        add_hit(hit["isolate_id"], "Gene end", hit)

    assembly_df = pd.DataFrame(rows)
    if assembly_df.empty:
        return assembly_df

    contigs_by_isolate = assembly_df.groupby("id")["contig"].transform(
        lambda values: ";".join(sorted({str(value) for value in values if pd.notna(value)}))
    )

    assembly_df["contigs_by_isolate"] = contigs_by_isolate
    assembly_df["contig_status"] = assembly_df.groupby("id")["contig"].transform(
        lambda values: "single" if values.dropna().nunique() == 1 else "multiple"
    )

    return sort_report_df_by_id(assembly_df)

def make_flank_table(flank_hits, gene_hits = (), buffer = 0):
    rows = []
    index = locus_index(gene_hits)
    assignments = defaultdict(list)
    for g in gene_hits:
        for position in ("upstream", "downstream"):
            if g.get(f"{position}_found"):
                label = f"{g.get('contig_orig')}:{g['s_start']}-{g['s_end']} ({position})"
                assignments[flank_key(g["isolate_id"], g["contig_id"], g[f"{position}_flank_query"], g[f"{position}_flank_coords"])].append(label)

    for hit in sorted(flank_hits, key = lambda h: (isolate_sort_key(h["isolate_id"]), h["query_id"])):
        rows.append({
            "id": hit["isolate_id"],
            "flank_type": hit.get("position"),
            "query": hit["query_id"],
            "contig": hit.get("contig_orig"),
            "s_start": hit.get("s_start"),
            "s_end": hit.get("s_end"),
            "identity": round(hit["identity"], 2) if hit.get("identity") is not None else None,
            "coverage": round(hit["coverage"], 2) if hit.get("coverage") is not None else None,
            "aln_length": hit.get("aln_length"),
            "mismatches": hit.get("mismatches"),
            "gaps": hit.get("gaps"),
            "gap_opens": hit.get("gap_opens"),
            "strand": hit.get("strand"),
            "evalue": hit.get("evalue"),
            "bitscore": hit.get("bitscore"),
            "assigned_to": "; ".join(assignments.get(flank_key(hit["isolate_id"], hit["contig_id"], hit["query_id"], (min(hit["s_start"], hit["s_end"]) - 1, max(hit["s_start"], hit["s_end"]))), [])) or "not assigned",
            "unassigned_reason": ""
        })
        if rows[-1]["assigned_to"] == "not assigned":
            rows[-1]["unassigned_reason"] = flank_unassigned_reason(hit, index, buffer)

    return pd.DataFrame(rows)

def make_primer_table(primer_hits):
    rows = []

    for hit in sorted(primer_hits, key = lambda h: (isolate_sort_key(h["isolate_id"]), h["total_mismatches"], -h["bitscore"])):
        rows.append({
            "id": hit["isolate_id"],
            "contig": hit.get("contig_orig"),
            "identity": round(hit["identity"], 2) if hit.get("identity") is not None else None,
            "coverage": round(hit["coverage"], 2) if hit.get("coverage") is not None else None,
            "total_differences": hit.get("total_mismatches"),
            "blast_mismatches": hit.get("blast_mismatches"),
            "unaligned_bases": hit.get("unaligned_bases"),
            "gaps": hit.get("gaps"),
            "gap_opens": hit.get("gap_opens"),
            "mismatch_positions": hit.get("mismatch_positions"),
            "mismatch_positions_from_3prime": hit.get("positions_from_3prime"),
            "perfect_match": hit.get("total_mismatches") == 0,
            "subject_start": hit.get("s_start"),
            "subject_end": hit.get("s_end"),
            "strand": hit.get("strand"),
            "evalue": hit.get("evalue"),
            "bitscore": hit.get("bitscore")
        })

    return pd.DataFrame(rows)


def make_totals_table(blast_results, isolate_ids):
    gene_hits = blast_results.get("gene_hits", [])
    gene_hits_before_filter = blast_results.get("gene_hits_before_filter", [])
    gene_start_hits = blast_results.get("gene_start_hits", [])
    gene_end_hits = blast_results.get("gene_end_hits", [])
    upstream_flank_hits = blast_results.get("upstream_flank_hits", [])
    downstream_flank_hits = blast_results.get("downstream_flank_hits", [])
    forward_primer_hits = blast_results.get("forward_primer_hits", [])
    reverse_primer_hits = blast_results.get("reverse_primer_hits", [])

    totals = [
        ("input isolates", len(set(isolate_ids))),
        ("gene hits before filter", len(gene_hits_before_filter)),
        ("gene hits", len(gene_hits)),
        ("high quality gene hits", sum(bool(h.get("high_quality")) for h in gene_hits)),
        ("gene start hits", len(gene_start_hits)),
        ("gene end hits", len(gene_end_hits)),
        ("upstream flank hits", len(upstream_flank_hits)),
        ("downstream flank hits", len(downstream_flank_hits)),
        ("gene hits with both flanks", sum(h.get("genomic_context") == "both_flanks" for h in gene_hits)),
        ("gene hits with upstream flank", sum(h.get("genomic_context") == "upstream_only" for h in gene_hits)),
        ("gene hits with downstream flank", sum(h.get("genomic_context") == "downstream_only" for h in gene_hits)),
        ("forward primer hits", len(forward_primer_hits)),
        ("reverse primer hits", len(reverse_primer_hits)),
        ("perfect forward primer hits", sum(h.get("total_mismatches") == 0 for h in forward_primer_hits)),
        ("perfect reverse primer hits", sum(h.get("total_mismatches") == 0 for h in reverse_primer_hits))
    ]

    return pd.DataFrame(totals, columns=["Metric", "Count"])

def make_primer_mismatch_table(forward_primer_hits, reverse_primer_hits):
    rows = []
    for primer, hits in [("forward", forward_primer_hits), ("reverse", reverse_primer_hits)]:
        mismatch_counts = sorted({hit["total_mismatches"] for hit in hits})
        for mismatch_count in mismatch_counts:
            count = sum(hit["total_mismatches"] == mismatch_count for hit in hits)
            percentage = count / len(hits) * 100 if hits else 0

            rows.append({
                "Primer": primer,
                "Mismatches": mismatch_count,
                "Count": count,
                "Percentage": round(percentage, 2)
            })

    return pd.DataFrame(rows)

def create_blast_summary_report(blast_results, config, isolate_ids, output_dir, extracted_sequences = None, before_N_filtering = None, blast_outputs = None):
    gene_name = config.gene_name
    output_dir = Path(output_dir)
    output_dir.mkdir(parents = True, exist_ok = True)
    output_path = output_dir / f"{gene_name}_blast_summary.xlsx"

    gene_hits = blast_results.get("gene_hits", [])
    gene_hits_before_filter = blast_results.get("gene_hits_before_filter", [])
    gene_start_hits = blast_results.get("gene_start_hits", [])
    gene_end_hits = blast_results.get("gene_end_hits", [])
    upstream_flank_hits = blast_results.get("upstream_flank_hits", [])
    downstream_flank_hits = blast_results.get("downstream_flank_hits", [])
    flank_hits = blast_results.get("flank_hits", [])
    forward_primer_hits = blast_results.get("forward_primer_hits", [])
    reverse_primer_hits = blast_results.get("reverse_primer_hits", [])

    logger.info("Creating BLAST summary report for gene: %s", gene_name)

    query_df = make_query_table(config)
    parameter_df = make_parameter_table(config)

    blast_queries = pd.concat([query_df, parameter_df.rename(columns={"Search": "Query Type", "Parameter": "Parameter"})], ignore_index = True, sort = False)

    isolate_df = make_isolate_summary(blast_results, isolate_ids, config, extracted_sequences, before_N_filtering).rename(columns = {
        "gene_found": "Gene Found",
        "n_candidate_loci": "Candidate Loci",
        "n_extracted_before_N_filter": "Extracted Before N Filter",
        "n_retained_loci": "Retained Loci",
        "locus_outcome": "Locus Outcome",
        "genomic_context": "Genomic Context",
        "n_forward_primer_hits": "Forward Primer Hits",
        "n_reverse_primer_hits": "Reverse Primer Hits"
    })

    isolate_df["Gene Found"] = isolate_df["Gene Found"].map({True: "Yes", False: "No"})

    def gene_rows(hits):
        return sort_report_df_by_id(pd.DataFrame([{
            "id": hit["isolate_id"],
            "gene copy": hit.get("gene_copy", "NA"),
            "contig": hit.get("contig_orig"),
            "identity": round(hit["identity"], 2),
            "coverage": round(hit["coverage"], 2),
            "high quality": "Yes" if hit.get("high_quality", False) else "No",
            "Query": hit.get("blast_query", "NA"),
            "Genomic context": hit.get("genomic_context", "NA"),
            "Upstream flank": hit.get("upstream_found", "NA"),
            "Downstream flank": hit.get("downstream_found", "NA"),
            "Subject start": hit.get("s_start"),
            "Subject end": hit.get("s_end"),
            "Strand": hit.get("strand", "NA"),
            "Aln_length": hit.get("aln_length"),
            "evalue": hit.get("evalue"),
            "bitscore": hit.get("bitscore"),
        } for hit in hits]))

    def primer_rows(hits):
        return sort_report_df_by_id(pd.DataFrame([{
            "id": hit["isolate_id"],
            "primer": hit.get("query_id", "NA"),
            "contig": hit.get("contig_orig"),
            "identity": round(hit["identity"], 2),
            "coverage": round(hit["coverage"], 2),
            "aln_length": hit.get("aln_length"),
            "query_start": hit.get("q_start"),
            "query_end": hit.get("q_end"),
            "subject_start": hit.get("s_start"),
            "subject_end": hit.get("s_end"),
            "total_differences": hit.get("total_mismatches"),
            "blast_mismatches": hit.get("blast_mismatches"),
            "unaligned_bases": hit.get("unaligned_bases"),
            "gaps": hit.get("gaps"),
            "gap_opens": hit.get("gap_opens"),
            "mismatch_positions": hit.get("mismatch_positions"),
            "mismatch_positions_from_3prime": hit.get("positions_from_3prime"),
            "strand": hit.get("strand", "NA"),
            "evalue": hit.get("evalue"),
            "bitscore": hit.get("bitscore"),
            "distance_to_nearest_locus": distance_to_locus(hit, primer_locus_index),
            "on_gene_locus": {True: "Yes", False: "No", None: "no locus on contig"}[None if distance_to_locus(hit, primer_locus_index) is None else distance_to_locus(hit, primer_locus_index) == 0],
        } for hit in hits]))

    primer_locus_index = locus_index(gene_hits_before_filter)
    flank_df = make_flank_table(flank_hits, gene_hits_before_filter, config.flank_buffer)
    assembly_df = make_gene_assembly_table(blast_results, config)

    search_rows = read_search_rows(blast_outputs or {}, config)
    query_summary_df = make_query_summary(search_rows, blast_results, config)
    isolate_query_df = make_isolate_query_matrix(search_rows, isolate_ids, config)
    hits_by_query_df = make_hits_by_query_table(search_rows)

    totals_df = make_totals_by_query(search_rows, blast_results, config)
    stage_df = make_stage_summary(search_rows, blast_results, isolate_ids, extracted_sequences, before_N_filtering)
    copy_number_df = make_copy_number_table(isolate_df)
    gene_flank_df = make_gene_flank_table(gene_hits_before_filter, gene_hits, config, forward_primer_hits, reverse_primer_hits)

    with pd.ExcelWriter(output_path, engine = "openpyxl") as writer:
        blast_queries.to_excel(writer, sheet_name = "Blast Queries", index = False)
        stage_df.to_excel(writer, sheet_name = "Stage summary", index = False)
        query_summary_df.to_excel(writer, sheet_name = "Query summary", index = False)
        isolate_query_df.to_excel(writer, sheet_name = "Isolate x Query", index = False)
        hits_by_query_df.to_excel(writer, sheet_name = "Hits by Query", index = False)
        totals_df.to_excel(writer, sheet_name = "Hit totals", index = False)
        isolate_df.to_excel(writer, sheet_name = "Isolate summary", index = False)
        copy_number_df.to_excel(writer, sheet_name = "Copy number", index = False)
        gene_rows(gene_hits).to_excel(writer, sheet_name = "Gene Hits", index = False)
        if active_queries(config.gene_start_queries):
            gene_rows(gene_start_hits).to_excel(writer, sheet_name = "Gene Start Hits", index = False)
        if active_queries(config.gene_end_queries):
            gene_rows(gene_end_hits).to_excel(writer, sheet_name = "Gene End Hits", index = False)
        gene_flank_df.to_excel(writer, sheet_name = "Gene Hits + Flanks", index = False)
        flank_df.to_excel(writer, sheet_name = "Flank Hits", index = False)
        assembly_df.to_excel(writer, sheet_name = "Gene Assembly", index = False)
        primer_rows(forward_primer_hits).to_excel(writer, sheet_name = "Forward Primer Hits", index = False)
        primer_rows(reverse_primer_hits).to_excel(writer, sheet_name = "Reverse Primer Hits", index = False)
        make_primer_mismatch_table(forward_primer_hits, reverse_primer_hits).to_excel(writer, sheet_name = "Hit totals", index = False, startrow = len(totals_df) + 3)
        format_workbook(writer)

    logger.info("BLAST summary report saved to %s", output_path)
    return output_path

#2 - Sequence extraction report
def make_extraction_table(sequences, kept_sequences, config):
    #gene sequences are the BLAST-aligned region so may be partial genes
    columns = [
        "id", "gene_copy", "contig", "outcome",
        "N_filter", "N_in_gene", "N_in_flank_only", "gene_N_count", "boundary_status", "blast_query",
        "genomic_context", "contig_length", "gene_start", "gene_end", "gene_length",
        "reference_query_length", "span_fraction_of_reference", "gene_with_flank_start", "gene_with_flank_end", "gene_with_flank_length",
        "dist_to_contig_start", "dist_to_contig_end", "strand", "flags"
    ]

    rows = []
    for s in sequences:
        key = (s["isolate_id"], s["gene_copy"], s["contig_id"])
        with_flank = s.get("with_flank")
        rows.append({
            "id": s["isolate_id"],
            "gene_copy": s["gene_copy"],
            "contig": s.get("original_contig_id", s["contig_id"]),
            "outcome": s.get("outcome"),
            "N_filter": "Kept" if key in kept_sequences else "Removed",
            "N_in_gene": bool(s.get("n_in_gene")),
            "N_in_flank_only": bool(s.get("n_in_flank_only")),
            "gene_N_count": s["sequence"].upper().count("N"),
            "boundary_status": s["boundary_status"],
            "blast_query": s.get("blast_query", "unknown"),
            "genomic_context": s.get("genomic_context", "unknown"),
            "contig_length": s["contig_length"],
            "gene_start": s["gene_start"],
            "gene_end": s["gene_end"],
            "gene_length": len(s["sequence"]),
            "reference_query_length": s.get("query_length"),
            "span_fraction_of_reference": round(len(s["sequence"]) / s["query_length"], 3) if s.get("query_length") else None,
            "gene_with_flank_start": with_flank["extract_start"] if with_flank else None,
            "gene_with_flank_end": with_flank["extract_end"] if with_flank else None,
            "gene_with_flank_length": len(with_flank["sequence"]) if with_flank else None,
            "dist_to_contig_start": s["dist_to_contig_start"],
            "dist_to_contig_end": s["dist_to_contig_end"],
            "strand": s["strand"],
            "flags": "; ".join(s.get("flags") or [])
        })

    return sort_report_df_by_id(pd.DataFrame(rows, columns = columns))

def create_sequence_extraction_report(before_N_filtering, extracted_sequences, config, output_dir):
    gene_name = config.gene_name
    params = config.extraction_params
    output_dir = Path(output_dir)
    output_dir.mkdir(parents = True, exist_ok = True)
    output_path = output_dir / f"{gene_name}_sequence_extraction.xlsx"

    logger.info("Creating sequence extraction report for gene: %s", gene_name)

    kept_sequences = {
        (s["isolate_id"], s["gene_copy"], s["contig_id"])
        for s in extracted_sequences
    }
    before_df = make_extraction_table(before_N_filtering, kept_sequences, config)
    after_df = make_extraction_table(extracted_sequences, kept_sequences, config)
    removed_df = before_df[before_df["N_filter"] == "Removed"]

    parameter_rows = [
        ("Gene name", gene_name),
        ("Flank extracted", bool(params.get("extract_flank"))),
        ("Flanking bp used", params["flanking_bp"] if params.get("extract_flank") else 0),
        ("Critical boundary distance (bp)", params["boundary_critical_distance"]),
        ("Warning boundary distance (bp)", params["boundary_warning_distance"]),
        ("Coordinate system", "0-based, end-exclusive"),
        ("Boundary assessment", "BLAST hit boundaries"),
        ("Sequence type", "Gene = BLAST-aligned region (may be a partial gene)"),
        ("N filter scope", "Gene plus flank (if extracted)")
    ]
    parameter_df = pd.DataFrame(parameter_rows, columns = ["Parameter", "Value"])

    n_before = len(before_df)
    n_removed = len(removed_df)
    totals_df = pd.DataFrame([
        ("Loci before N filtering", n_before),
        ("Loci after N filtering", len(after_df)),
        ("Loci removed by N filtering", n_removed),
        ("Percentage removed", round(n_removed / n_before * 100, 2) if n_before else 0)
    ], columns = ["Metric", "Value"])

    boundary_df = pd.DataFrame([
        {
            "Boundary status": status,
            "Before N filtering": sum(s["boundary_status"] == status for s in before_N_filtering),
            "After N filtering": sum(s["boundary_status"] == status for s in extracted_sequences)
        }
        for status in ["OK", "WARNING", "CRITICAL"]
    ])

    with pd.ExcelWriter(output_path, engine = "openpyxl") as writer:
        parameter_df.to_excel(writer, sheet_name = "Extraction Parameters", index = False)
        current_row = write_df_with_spacing(totals_df, writer, "Extraction Summary", 0)
        write_df_with_spacing(boundary_df, writer, "Extraction Summary", current_row)
        before_df.to_excel(writer, sheet_name = "Before N Filtering", index = False)
        after_df.to_excel(writer, sheet_name = "After N Filtering", index = False)
        removed_df.to_excel(writer, sheet_name = "Removed by N Filtering", index = False)
        format_workbook(writer)

        #summary contains two separate tables, so remove the worksheet-wide filter
        writer.book["Extraction Summary"].auto_filter.ref = None

    logger.info("Sequence extraction report saved to %s", output_path)
    return output_path
#one function to make all reports 

def generate_all_reports(config, blast_results, output_dir, isolate_ids, blast_outputs, extracted_sequences, before_N_filtering):
    output_dir = Path(output_dir)
    gene_name = config.gene_name
    reports_generated = {
        "gene_hits": None,
        "flank_hits": None
    }
    if blast_outputs.get("gene") is not None:
        gene_stats_path = output_dir / f"{gene_name}_gene_hit_stats.txt"
        report_gene_hit_stats(blast_outputs["gene"], gene_stats_path, active_queries(config.gene_queries))
        reports_generated["gene_hits"] = gene_stats_path
    if blast_outputs.get("flank") is not None:
        flank_stats_path = output_dir / f"{gene_name}_flank_hit_stats.txt"
        report_flank_hit_stats(blast_outputs["flank"], flank_stats_path, active_queries(config.flank_queries))
        reports_generated["flank_hits"] = flank_stats_path
    reports_generated["blast_summary"] = create_blast_summary_report(blast_results, config, isolate_ids, output_dir, extracted_sequences, before_N_filtering, blast_outputs)
    reports_generated["sequence_extraction"] = create_sequence_extraction_report(before_N_filtering, extracted_sequences, config, output_dir)
    if blast_outputs.get("primers") is not None:
        from primer_report import create_primer_report
        reports_generated["primer_report"] = create_primer_report(blast_results, config, isolate_ids, output_dir)

    return reports_generated
