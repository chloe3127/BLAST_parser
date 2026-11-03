import logging
from pathlib import Path
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.backends.backend_agg import FigureCanvasAgg
import numpy as np
import pandas as pd
import seaborn as sns
from collections import defaultdict, Counter
from blast_parser import gene_columns, primer_columns
from blast_prep import active_queries
from reporting import read_search_rows, isolate_sort_key

logger = logging.getLogger(__name__)

#1 - plot the RAW BLAST results
def read_raw_blast(blast_path, is_primer = False):
    hits = defaultdict(list)
    columns = primer_columns if is_primer else gene_columns
    with Path(blast_path).open("r", encoding = "utf-8") as f:
        for line_number, line in enumerate(f, start = 1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) != len(columns):
                raise ValueError(f"{blast_path}:{line_number} expected {len(columns)} BLAST columns, got {len(fields)}")
            row = dict(zip(columns, fields))
            identity = float(row["pct_ident"])
            query_length = int(row["qlength"])
            query_start, query_end = int(row["q_start"]), int(row["q_end"])
            if query_length <= 0 or not 1 <= query_start <= query_length or not 1 <= query_end <= query_length:
                raise ValueError(f"{blast_path}:{line_number} has invalid query coordinates: qlength={query_length}, q_start={query_start}, q_end={query_end}")
            aligned_bases = abs(query_end - query_start) + 1
            coverage = aligned_bases / query_length * 100
            if not np.isfinite(identity) or not 0 <= identity <= 100:
                raise ValueError(f"{blast_path}:{line_number} has invalid identity: {identity}")
            if is_primer:
                unaligned_bases = query_length - aligned_bases
                total_mismatches = int(row["mismatches"]) + unaligned_bases + int(row["gaps"])
                hits[row["query_id"]].append((identity, coverage, total_mismatches))
            else:
                hits[row["query_id"]].append((identity, coverage))
    return hits

def plot_blast_query(hits, title, params, output_path, is_primer = False):
    fig = Figure(figsize = (15, 4.8), layout = "constrained")
    FigureCanvasAgg(fig)
    axes = fig.subplots(1, 3)
    n_columns = 3 if is_primer else 2
    values = np.asarray(hits, dtype = float).reshape(-1, n_columns)
    identity, coverage = values[:, 0], values[:, 1]

    if is_primer:
        min_identity = None
        high_identity = None
        min_coverage = params["min_alignment_pct"]
        high_coverage = None
        total_mismatches = values[:, 2]
        passing = int(np.sum((coverage >= min_coverage) & (total_mismatches <= params["max_mismatches"])))
    else:
        min_identity = params["min_identity"]
        min_coverage = params["min_coverage"]
        high_identity = params.get("high_quality_identity", params.get("high_identity", min_identity))
        high_coverage = params.get("high_coverage", min_coverage)
        passing = int(np.sum((identity >= min_identity) & (coverage >= min_coverage)))

    thresholds = [(min_identity, high_identity), (min_coverage, high_coverage)]
    labels = ["BLAST identity (%)", "Query coverage (%)"]

    for ax, data, label, (minimum, high) in zip(axes[:2], [identity, coverage], labels, thresholds):
        ax.hist(data, bins = np.arange(0, 101, 1), color = "#9da6cf", edgecolor = "white", linewidth = 0.3)
        if minimum is not None:
            ax.axvline(minimum, color = "#d18c20", linestyle = "--", label = f"Minimum: {minimum:g}%")
        if high is not None and high != minimum:
            ax.axvline(high, color = "#d94c4c", linestyle = ":", label = f"High quality: {high:g}%")
        ax.set(xlabel = label, ylabel = "BLAST hits", xlim = (0, 102))
        if minimum is not None or high is not None:
            ax.legend(fontsize = 8)
        ax.grid(axis = "y", alpha = 0.2)

    ax = axes[2]
    if len(values):
        counts = ax.hexbin(identity, coverage, gridsize = 45, extent = (0, 100, 0, 100), mincnt = 1, bins = "log", cmap = "viridis")
        fig.colorbar(counts, ax = ax, label = "BLAST hits (log scale)")
    if min_identity is not None:
        ax.axvline(min_identity, color = "#d18c20", linestyle = "--")
    ax.axhline(min_coverage, color = "#d18c20", linestyle = "--")
    if high_identity is not None and high_identity != min_identity:
        ax.axvline(high_identity, color = "#d94c4c", linestyle = ":")
    if high_coverage is not None and high_coverage != min_coverage:
        ax.axhline(high_coverage, color = "#d94c4c", linestyle = ":")
    ax.set(xlabel = labels[0], ylabel = labels[1], xlim = (-2, 102), ylim = (-2, 102))

    if is_primer:
        ax.set_title(
            f"BLAST hits passing filters: {passing:,}/{len(values):,}\n"
            f"Coverage ≥ {min_coverage:g}%; total mismatches ≤ {params['max_mismatches']}",
            fontsize = 10
        )
    else:
        ax.set_title(f"BLAST hits passing both thresholds: {passing:,}/{len(values):,}", fontsize = 10)

    if not len(values):
        for ax in axes:
            ax.text(0.5, 0.5, "No BLAST hits", transform = ax.transAxes, ha = "center")
    fig.suptitle(f"{title} - raw BLAST hits (n={len(values):,})", fontsize = 13)
    output_path = Path(output_path)
    output_path.parent.mkdir(parents = True, exist_ok = True)
    fig.savefig(output_path, dpi = 300)
    fig.clear()
    logger.info("Saved raw BLAST plot: %s", output_path)
    return output_path

def plot_raw_blast_results(blast_outputs, config, output_dir):
    #keep every reported HSP and plot each reference query separately
    gene_name = config.gene_name
    searches = [
        ("gene", config.gene_queries, config.gene_blast_params),
        ("flank", config.flank_queries, config.flank_blast_params),
        ("gene_start", config.gene_start_queries, config.gene_start_blast_params),
        ("gene_end", config.gene_end_queries, config.gene_end_blast_params),
        ("primers", config.primer_queries, config.primer_blast_params)
    ]
    outputs = []
    for search, queries, params in searches:
        blast_path = blast_outputs.get(search)
        if blast_path is None:
            continue
        is_primer = search == "primers"
        raw_hits = read_raw_blast(blast_path, is_primer = is_primer)
        for query in active_queries(queries):
            output_path = Path(output_dir) / f"{gene_name}_{search}_{query}_raw_blast.png"
            outputs.append(plot_blast_query(
                raw_hits.get(query, []), f"{gene_name}: {search} - {query}", params,
                output_path, is_primer = is_primer
            ))
    return outputs

#2 - Plot primer mismatches
def plot_primer_mismatch(blast_outputs, config, output_dir):
    blast_path = blast_outputs.get("primers")
    if blast_path is None:
        return None
    raw_hits = read_raw_blast(blast_path, is_primer = True)
    params = config.primer_blast_params
    min_coverage = params["min_alignment_pct"]
    max_mismatches = params["max_mismatches"]
    #pool queries by primer direction
    primer_queries = active_queries(config.primer_queries)
    primer_names = ["forward", "reverse"]
    primer_data = {}

    for primer in primer_names:
        rows = [hit for name, data in primer_queries.items() if data["position"] == primer for hit in raw_hits.get(name, [])]
        values = np.asarray(rows, dtype = float).reshape(-1, 3)
        mismatches = values[values[:, 1] >= min_coverage, 2].astype(int)
        primer_data[primer] = (len(values), mismatches)

    largest_mismatch = max(
        [max_mismatches + 1] +
        [int(mismatches.max()) for _, mismatches in primer_data.values() if len(mismatches)]
    )
    mismatch_counts = np.arange(largest_mismatch + 1)
    fig = Figure(figsize = (12, 4.8), layout = "constrained")
    FigureCanvasAgg(fig)
    axes = fig.subplots(1, 2)

    for ax, primer in zip(axes, primer_names):
        raw_count, mismatches = primer_data[primer]
        counts = np.bincount(mismatches, minlength = largest_mismatch + 1)
        passing = int(np.sum(mismatches <= max_mismatches))
        colours = ["#9da6cf" if count <= max_mismatches else "#dddddd" for count in mismatch_counts]
        ax.bar(mismatch_counts, counts, color = colours, edgecolor = "white", width = 0.85)

        for mismatch, count in zip(mismatch_counts, counts):
            if count:
                ax.annotate(
                    f"{count:,}", (mismatch, count),
                    xytext = (0, 4), textcoords = "offset points",
                    ha = "center", fontsize = 9
                )

        ax.axvline(
            max_mismatches + 0.5, color = "#d18c20", linestyle = "--",
            label = f"Maximum total mismatches: {max_mismatches}"
        )
        ax.set(
            xlabel = "Total mismatches",
            ylabel = "BLAST hits meeting coverage threshold",
            xlim = (-0.5, largest_mismatch + 0.5)
        )
        ax.set_xticks(mismatch_counts)
        ax.set_title(
            f"{primer.capitalize()} primer\n"
            f"Coverage retained: {len(mismatches):,}/{raw_count:,}; passing mismatch filter: {passing:,}",
            fontsize = 10
        )
        if not len(mismatches):
            ax.text(
                0.5, 0.5, "No hits meet coverage threshold",
                transform = ax.transAxes, ha = "center"
            )
        ax.grid(axis = "y", alpha = 0.2)
        ax.set_axisbelow(True)
        ax.legend(fontsize = 8)

    upper_y = max(ax.get_ylim()[1] for ax in axes)
    for ax in axes:
        ax.set_ylim(0, max(1, upper_y * 1.15))

    fig.suptitle(
        f"{config.gene_name} - primer mismatches among hits with coverage ≥ {min_coverage:g}%",
        fontsize = 13
    )
    output_path = Path(output_dir) / f"{config.gene_name}_primer_mismatches.png"
    output_path.parent.mkdir(parents = True, exist_ok = True)
    fig.savefig(output_path, dpi = 300)
    fig.clear()
    logger.info("Saved primer mismatch plot: %s", output_path)
    return output_path

#3 - Plot sequence extraction QC before and after N filtering
def plot_sequence_extraction(before_sequences, after_sequences, config, output_dir):
    if not before_sequences:
        logger.info("No extracted sequences available for plotting")
        return None

    gene_name = config.gene_name
    params = config.extraction_params
    #flank is only extracted when extract_flank is on
    flank_on = bool(params.get("extract_flank"))
    requested = int(params["flanking_bp"]) if flank_on else 0

    fig = Figure(figsize = (15, 9), layout = "constrained")
    FigureCanvasAgg(fig)
    axes = fig.subplots(2, 3, sharex = "col", sharey = "col")

    #use identical bins for both rows
    sequence_lengths = [len(seq["sequence"]) for seq in before_sequences]
    length_bins = np.histogram_bin_edges(sequence_lengths, bins = 30)
    context_bins = np.linspace(-0.5, max(requested, 1) + 0.5, 31)
    boundary_max = 1

    for row, (stage, sequences) in enumerate([
        ("Before N filtering", before_sequences),
        ("After N filtering", after_sequences)
    ]):
        boundary_ax, context_ax, length_ax = axes[row]
        n_sequences = len(sequences)

        #boundary status based on BLAST hit coordinates
        statuses = ["OK", "WARNING", "CRITICAL"]
        counts = [sum(seq["boundary_status"] == status for seq in sequences) for status in statuses]
        boundary_max = max(boundary_max, max(counts))
        bars = boundary_ax.bar(
            statuses, counts,
            color = ["#9da6cf", "#d18c20", "#d94c4c"],
            edgecolor = "white"
        )

        for bar, count in zip(bars, counts):
            percentage = count / n_sequences * 100 if n_sequences else 0
            boundary_ax.annotate(
                f"{count:,}\n({percentage:.1f}%)",
                (bar.get_x() + bar.get_width() / 2, count),
                xytext = (0, 4), textcoords = "offset points",
                ha = "center", va = "bottom", fontsize = 9
            )

        boundary_ax.set_title(f"{stage}: {n_sequences:,} loci\nContig boundary status", fontsize = 10)
        boundary_ax.set_ylabel("Loci")
        boundary_ax.set_xlabel("Based on BLAST hit boundaries")

        if not sequences:
            for ax in (context_ax, length_ax):
                ax.text(
                    0.5, 0.5, "No sequences retained",
                    transform = ax.transAxes, ha = "center", va = "center"
                )
            context_ax.set_title("Flanks recovered around gene", fontsize = 10)
            length_ax.set_title("Gene length", fontsize = 10)
            continue

        #flank recovered either side of the gene, clipped at contig edges
        if not flank_on:
            context_ax.text(
                0.5, 0.5, "Flank not extracted\n(extract_flank is off)",
                transform = context_ax.transAxes, ha = "center", va = "center"
            )
            context_ax.set_title("Flanks recovered around gene", fontsize = 10)
        else:
            upstream_lengths = []
            downstream_lengths = []

            for seq in sequences:
                with_flank = seq.get("with_flank")
                if not with_flank:
                    continue
                left_flank = seq["gene_start"] - with_flank["extract_start"]
                right_flank = with_flank["extract_end"] - seq["gene_end"]

                if seq["strand"] == "forward":
                    upstream_lengths.append(left_flank)
                    downstream_lengths.append(right_flank)
                elif seq["strand"] == "reverse":
                    upstream_lengths.append(right_flank)
                    downstream_lengths.append(left_flank)
                else:
                    raise ValueError(f"Invalid strand: {seq['strand']}")

            upstream_complete = sum(length == requested for length in upstream_lengths)
            downstream_complete = sum(length == requested for length in downstream_lengths)

            context_ax.hist(
                [upstream_lengths, downstream_lengths],
                bins = context_bins,
                color = ["#9da6cf", "#70b6a1"],
                edgecolor = "white", linewidth = 0.3,
                label = [
                    f"Upstream: {upstream_complete:,} full length",
                    f"Downstream: {downstream_complete:,} full length"
                ]
            )
            context_ax.set_title("Flanks recovered around gene", fontsize = 10)
            context_ax.axvline(
                requested, color = "#d18c20", linestyle = "--",
                label = f"Requested: {requested:,} bp"
            )
            context_ax.legend(fontsize = 8)

        #gene length
        lengths = [len(seq["sequence"]) for seq in sequences]
        median_length = np.median(lengths)

        length_ax.hist(
            lengths, bins = length_bins,
            color = "#9da6cf", edgecolor = "white", linewidth = 0.3
        )
        length_ax.axvline(
            median_length, color = "#d18c20", linestyle = "--",
            label = f"Median: {median_length:,.0f} bp"
        )
        length_ax.set_title("Gene length", fontsize = 10)
        length_ax.legend(fontsize = 8)

    for row_axes in axes:
        row_axes[0].set_ylim(0, boundary_max * 1.35)
        row_axes[1].set_xlabel("Flanking sequence length (bp)")
        row_axes[1].set_ylabel("Loci per direction")
        row_axes[2].set_xlabel("Gene length (bp)")
        row_axes[2].set_ylabel("Loci")

        for ax in row_axes:
            ax.grid(axis = "y", alpha = 0.2)
            ax.set_axisbelow(True)
            ax.tick_params(axis = "x", labelbottom = True)

    removed = len(before_sequences) - len(after_sequences)
    fig.suptitle(
        f"{gene_name} - sequence extraction QC\n"
        f"N filtering removed {removed:,}/{len(before_sequences):,} loci "
        f"({removed / len(before_sequences):.1%}); "
        f"boundary thresholds: critical ≤ {params['boundary_critical_distance']} bp, "
        f"warning ≤ {params['boundary_warning_distance']} bp",
        fontsize = 13
    )
    output_path = Path(output_dir) / f"{gene_name}_sequence_extraction_QC.png"
    output_path.parent.mkdir(parents = True, exist_ok = True)
    fig.savefig(output_path, dpi = 300)
    fig.clear()
    logger.info("Saved sequence extraction QC plot: %s", output_path)
    return output_path

#4 - Plot BLAST hit positions relative to contig boundaries
def plot_contig_boundaries(sequences, config, output_dir):
    if not sequences:
        logger.info("No extracted sequences available for contig boundary plotting")
        return None

    critical = config.extraction_params["boundary_critical_distance"]
    warning = config.extraction_params["boundary_warning_distance"]
    if not 0 <= critical < warning:
        raise ValueError("Boundary thresholds must satisfy 0 <= critical < warning")

    start_distances = np.asarray([s["dist_to_contig_start"] for s in sequences])
    end_distances = np.asarray([s["dist_to_contig_end"] for s in sequences])
    nearest = np.minimum(start_distances, end_distances)

    fig = Figure(figsize = (15, 5.5), layout = "constrained")
    FigureCanvasAgg(fig)
    axes = fig.subplots(1, 2)

    #scatter of distances to both contig edges
    groups = [
        ("OK", nearest > warning, "#70b6a1", f"OK: >{warning:g} bp"),
        ("WARNING", (nearest > critical) & (nearest <= warning), "#e9c46a", f"Warning: >{critical:g}–{warning:g} bp"),
        ("CRITICAL", nearest <= critical, "#d94c4c", f"Critical: ≤{critical:g} bp")
    ]
    for status, mask, colour, label in groups:
        axes[0].scatter(
            start_distances[mask], end_distances[mask],
            s = 30, color = colour, edgecolors = "#555555",
            linewidths = 0.4, alpha = 0.6,
            label = f"{label} (n={int(mask.sum()):,})"
        )

    for threshold, colour in [(critical, "#d94c4c"), (warning, "#d18c20")]:
        axes[0].axvline(threshold, color = colour, linestyle = "--", linewidth = 1)
        axes[0].axhline(threshold, color = colour, linestyle = "--", linewidth = 1)

    axes[0].set_xscale("symlog", linthresh = 1)
    axes[0].set_yscale("symlog", linthresh = 1)
    axes[0].set(
        xlabel = "Distance from BLAST hit to contig start (bp)",
        ylabel = "Distance from BLAST hit to contig end (bp)",
        title = "BLAST hit position relative to contig boundaries"
    )
    upper_distance = max(float(start_distances.max()), float(end_distances.max()), warning, 1) * 1.3
    axes[0].set_xlim(-0.2, upper_distance)
    axes[0].set_ylim(-0.2, upper_distance)
    axes[0].legend(fontsize = 8)

    #non-overlapping distance bins, adjusted to configured thresholds
    upper_bounds = [critical, warning]
    upper_bounds.extend(bound for bound in [500, 1000, 5000, 10000] if bound > warning)
    counts = []
    labels = []
    lower = None

    for upper in upper_bounds:
        if lower is None:
            mask = nearest <= upper
            label = f"0–{upper:g}\n(Critical)"
        else:
            mask = (nearest > lower) & (nearest <= upper)
            label = f">{lower:g}–{upper:g}"
            if upper == warning:
                label += "\n(Warning)"
        counts.append(int(mask.sum()))
        labels.append(label)
        lower = upper

    counts.append(int(np.sum(nearest > upper_bounds[-1])))
    labels.append(f">{upper_bounds[-1]:g}")
    colours = ["#d94c4c", "#e9c46a"] + ["#9da6cf"] * (len(counts) - 2)
    bars = axes[1].bar(
        np.arange(len(counts)), counts,
        color = colours, edgecolor = "white"
    )
    for bar, count in zip(bars, counts):
        axes[1].annotate(
            f"{count:,}", (bar.get_x() + bar.get_width() / 2, count),
            xytext = (0, 4), textcoords = "offset points",
            ha = "center", fontsize = 9
        )

    axes[1].set_xticks(np.arange(len(labels)))
    axes[1].set_xticklabels(labels, fontsize = 8)
    axes[1].set(
        xlabel = "Distance to nearest contig edge (bp)",
        ylabel = "Loci",
        title = "Distance to nearest contig edge",
        ylim = (0, max(1, max(counts) * 1.2))
    )

    for ax in axes:
        ax.grid(alpha = 0.2)
        ax.set_axisbelow(True)

    fig.suptitle(
        f"{config.gene_name} - contig boundary QC (n={len(sequences):,} loci)",
        fontsize = 13
    )
    output_path = Path(output_dir) / f"{config.gene_name}_contig_boundary_QC.png"
    output_path.parent.mkdir(parents = True, exist_ok = True)
    fig.savefig(output_path, dpi = 300)
    fig.clear()
    logger.info("Saved contig boundary QC plot: %s", output_path)
    return output_path

#5 - Compare the alternative gene reference queries
def plot_query_comparison(blast_outputs, config, output_dir):
    blast_path = blast_outputs.get("gene")
    queries = list(active_queries(config.gene_queries))
    if blast_path is None or len(queries) < 2:
        return None

    params = config.gene_blast_params
    raw_hits = read_raw_blast(blast_path)
    raw_counts, passing_counts, passing_identity = [], [], []
    for query in queries:
        values = np.asarray(raw_hits.get(query, []), dtype = float).reshape(-1, 2)
        passed = (values[:, 0] >= params["min_identity"]) & (values[:, 1] >= params["min_coverage"])
        raw_counts.append(len(values))
        passing_counts.append(int(passed.sum()))
        passing_identity.append(values[passed, 0])

    fig = Figure(figsize = (15, 5.5), layout = "constrained")
    FigureCanvasAgg(fig)
    axes = fig.subplots(1, 2)

    positions = np.arange(len(queries))
    width = 0.4
    raw_bars = axes[0].bar(positions - width / 2, raw_counts, width, color = "#dddddd", edgecolor = "white", label = "Raw alignment rows")
    pass_bars = axes[0].bar(positions + width / 2, passing_counts, width, color = "#9da6cf", edgecolor = "white", label = "Passing thresholds")
    for bar in list(raw_bars) + list(pass_bars):
        axes[0].annotate(
            f"{int(bar.get_height()):,}", (bar.get_x() + bar.get_width() / 2, bar.get_height()),
            xytext = (0, 4), textcoords = "offset points", ha = "center", fontsize = 9
        )
    axes[0].set_xticks(positions)
    axes[0].set_xticklabels(queries, rotation = 30, ha = "right", fontsize = 8)
    axes[0].set(ylabel = "BLAST alignment rows", title = "BLAST alignment rows per reference query")
    axes[0].legend(fontsize = 8)

    box = axes[1].boxplot(
        [v if len(v) else [np.nan] for v in passing_identity], positions = positions,
        patch_artist = True, widths = 0.5, medianprops = {"color": "#d18c20"}
    )
    for patch in box["boxes"]:
        patch.set(facecolor = "#9da6cf", edgecolor = "#555555")
    axes[1].axhline(params["min_identity"], color = "#d18c20", linestyle = "--", label = f"Minimum: {params['min_identity']:g}%")
    axes[1].set_xticks(positions)
    axes[1].set_xticklabels(queries, rotation = 30, ha = "right", fontsize = 8)
    axes[1].set(ylabel = "BLAST identity (%)", title = "Identity of passing hits per reference query")
    axes[1].legend(fontsize = 8)

    for ax in axes:
        ax.grid(axis = "y", alpha = 0.2)
        ax.set_axisbelow(True)

    fig.suptitle(f"{config.gene_name} - alternative reference query comparison", fontsize = 13)
    output_path = Path(output_dir) / f"{config.gene_name}_query_comparison.png"
    output_path.parent.mkdir(parents = True, exist_ok = True)
    fig.savefig(output_path, dpi = 300)
    fig.clear()
    logger.info("Saved query comparison plot: %s", output_path)
    return output_path

#6 - Show how flanks were assigned to gene hits
def plot_flank_assignment(blast_results, config, output_dir):
    gene_hits = (blast_results or {}).get("gene_hits", [])
    flank_queries = active_queries(config.flank_queries)
    if not gene_hits or not flank_queries:
        return None

    buffer = config.flank_buffer
    fig = Figure(figsize = (15, 5.5), layout = "constrained")
    FigureCanvasAgg(fig)
    axes = fig.subplots(1, 2)

    #observed distance between each gene hit and its assigned flank, per flank query
    colours = ["#9da6cf", "#70b6a1", "#e9c46a", "#d94c4c"]
    all_distances = []
    expected_values = [d["distance"] for d in flank_queries.values()]
    observed = [h[f"{d['position']}_distance"] for d in flank_queries.values() for h in gene_hits if h.get(f"{d['position']}_flank_query")]
    low = int(min(observed + expected_values, default = 0)) - 10
    high = int(max(observed + expected_values, default = 0)) + 10
    bins = np.arange(low - 0.5, high + 1.5, max(1, (high - low) // 60))
    for index, (name, data) in enumerate(flank_queries.items()):
        position = data["position"]
        distances = [h[f"{position}_distance"] for h in gene_hits if h.get(f"{position}_flank_query") == name]
        all_distances.extend(distances)
        colour = colours[index % len(colours)]
        axes[0].hist(distances, bins = bins, color = colour, edgecolor = "white", linewidth = 0.3, alpha = 0.8,
                     label = f"{name} ({position}): {len(distances):,}/{len(gene_hits):,} loci")
        axes[0].axvline(data["distance"], color = colour, linestyle = "--")
    axes[0].set(
        xlabel = "Distance from gene hit to assigned flank (bp; negative = overlap)",
        ylabel = "Loci", title = f"Assigned flank distance (dashed = expected; accepted within ±{buffer:g} bp)",
        xlim = (low, high)
    )
    axes[0].legend(fontsize = 8)

    #genomic context by strand
    contexts = ["both_flanks", "upstream_only", "downstream_only", "no_flanks"]
    strands = ["forward", "reverse"]
    bottoms = np.zeros(len(strands))
    for context, colour in zip(contexts, ["#70b6a1", "#9da6cf", "#e9c46a", "#d94c4c"]):
        counts = np.array([sum(h["genomic_context"] == context and h["strand"] == strand for h in gene_hits) for strand in strands])
        bars = axes[1].bar(strands, counts, bottom = bottoms, color = colour, edgecolor = "white", label = f"{context}: {int(counts.sum()):,}")
        for bar, count, bottom in zip(bars, counts, bottoms):
            if count:
                axes[1].text(bar.get_x() + bar.get_width() / 2, bottom + count / 2, f"{count:,}", ha = "center", va = "center", fontsize = 9)
        bottoms += counts
    axes[1].set_ylim(0, max(1, bottoms.max() * 1.15))
    axes[1].set(xlabel = "Gene hit strand", ylabel = "Loci", title = "Genomic context assigned to each locus")
    axes[1].legend(fontsize = 8)

    for ax in axes:
        ax.grid(axis = "y", alpha = 0.2)
        ax.set_axisbelow(True)

    fig.suptitle(f"{config.gene_name} - flank assignment (n={len(gene_hits):,} loci)", fontsize = 13)
    output_path = Path(output_dir) / f"{config.gene_name}_flank_assignment.png"
    output_path.parent.mkdir(parents = True, exist_ok = True)
    fig.savefig(output_path, dpi = 300)
    fig.clear()
    logger.info("Saved flank assignment plot: %s", output_path)
    return output_path

status_colours = {"OK": "#70b6a1", "WARNING": "#e9c46a", "CRITICAL": "#d94c4c"}

def save_figure(fig, output_path, label):
    output_path = Path(output_path)
    output_path.parent.mkdir(parents = True, exist_ok = True)
    fig.savefig(output_path, dpi = 300)
    fig.clear()
    logger.info("Saved %s plot: %s", label, output_path)
    return output_path

#7 - gene length compared with the reference query
def plot_span_vs_reference(before_sequences, config, output_dir):
    sequences = [s for s in before_sequences if s.get("query_length")]
    if not sequences:
        return None

    spans = np.asarray([len(s["sequence"]) for s in sequences], dtype = float)
    references = np.asarray([s["query_length"] for s in sequences], dtype = float)
    fractions = spans / references
    statuses = np.asarray([s["boundary_status"] for s in sequences])

    fig = Figure(figsize = (15, 5.5), layout = "constrained")
    FigureCanvasAgg(fig)
    axes = fig.subplots(1, 2)

    low = min(float(fractions.min()), 1.0) - 0.05
    upper = max(1.05, float(fractions.max()) + 0.02)
    bins = np.linspace(low, upper, 40)
    axes[0].hist(fractions, bins = bins, color = "#9da6cf", edgecolor = "white", linewidth = 0.3)
    axes[0].axvline(1, color = "#d18c20", linestyle = "--", label = "Span = reference length")
    axes[0].axvline(float(np.median(fractions)), color = "#d94c4c", linestyle = ":", label = f"Median: {np.median(fractions):.2f}")
    axes[0].set(xlabel = "Gene length / reference query length", ylabel = "Loci", title = "Fraction of reference covered by the gene")
    axes[0].legend(fontsize = 8)

    for status, colour in status_colours.items():
        mask = statuses == status
        axes[1].scatter(references[mask], spans[mask], s = 30, color = colour, edgecolors = "#555555", linewidths = 0.4, alpha = 0.7, label = f"{status} (n={int(mask.sum()):,})")
    values_low = min(float(references.min()), float(spans.min()))
    values_high = max(float(references.max()), float(spans.max()))
    pad = max(5.0, (values_high - values_low) * 0.05)
    limit = [values_low - pad, values_high + pad]
    axes[1].plot(limit, limit, color = "#d18c20", linestyle = "--", linewidth = 1)
    axes[1].set(xlabel = "Reference query length (bp)", ylabel = "Gene length (bp)", title = "Gene length vs reference (colour = contig boundary status)", xlim = limit, ylim = limit)
    axes[1].legend(fontsize = 8)

    for ax in axes:
        ax.grid(alpha = 0.2)
        ax.set_axisbelow(True)
    fig.suptitle(f"{config.gene_name} - gene length vs reference length (n={len(sequences):,} loci, before N filtering)", fontsize = 13)
    return save_figure(fig, Path(output_dir) / f"{config.gene_name}_span_vs_reference.png", "span vs reference")

#8 - Loci per isolate
def plot_copy_number(blast_results, extracted_sequences, config, output_dir):
    isolate_ids = {str(i) for i in (blast_results or {}).get("isolate_ids", [])}
    if not isolate_ids:
        return None
    candidates = Counter(str(h["isolate_id"]) for h in blast_results.get("gene_hits", []))
    retained = Counter(str(s["isolate_id"]) for s in extracted_sequences)
    largest = max([1] + list(candidates.values()) + list(retained.values()))
    copies = np.arange(largest + 1)
    candidate_counts = [sum(candidates.get(i, 0) == n for i in isolate_ids) for n in copies]
    retained_counts = [sum(retained.get(i, 0) == n for i in isolate_ids) for n in copies]

    fig = Figure(figsize = (9, 5.5), layout = "constrained")
    FigureCanvasAgg(fig)
    ax = fig.subplots()
    width = 0.4
    for offset, counts, colour, label in [(-width / 2, candidate_counts, "#9da6cf", "Candidate loci"), (width / 2, retained_counts, "#70b6a1", "Retained after N filtering")]:
        bars = ax.bar(copies + offset, counts, width, color = colour, edgecolor = "white", label = label)
        for bar, count in zip(bars, counts):
            ax.annotate(f"{count:,}", (bar.get_x() + bar.get_width() / 2, count), xytext = (0, 4), textcoords = "offset points", ha = "center", fontsize = 9)
    ax.set_xticks(copies)
    ax.set(xlabel = "Loci per isolate", ylabel = "Isolates", title = f"Gene copy number per isolate (n={len(isolate_ids):,} isolates)")
    ax.set_ylim(0, max(1, max(candidate_counts + retained_counts) * 1.2))
    ax.legend(fontsize = 8)
    ax.grid(axis = "y", alpha = 0.2)
    ax.set_axisbelow(True)
    fig.suptitle(f"{config.gene_name} - copy number", fontsize = 13)
    return save_figure(fig, Path(output_dir) / f"{config.gene_name}_copy_number.png", "copy number")

#9 - Identity and coverage of the retained loci
def plot_retained_quality(extracted_sequences, config, output_dir):
    if not extracted_sequences:
        return None
    params = config.gene_blast_params
    identity = np.asarray([s["blast_identity"] for s in extracted_sequences], dtype = float)
    coverage = np.asarray([s["blast_coverage"] for s in extracted_sequences], dtype = float)

    fig = Figure(figsize = (15, 4.8), layout = "constrained")
    FigureCanvasAgg(fig)
    axes = fig.subplots(1, 2)
    for ax, data, label, minimum, high in [
        (axes[0], identity, "BLAST identity (%)", params["min_identity"], params.get("high_identity", params["min_identity"])),
        (axes[1], coverage, "Query coverage (%)", params["min_coverage"], params.get("high_coverage", params["min_coverage"]))
    ]:
        ax.hist(data, bins = np.arange(0, 101, 1), color = "#9da6cf", edgecolor = "white", linewidth = 0.3)
        ax.axvline(minimum, color = "#d18c20", linestyle = "--", label = f"Minimum: {minimum:g}%")
        if high != minimum:
            ax.axvline(high, color = "#d94c4c", linestyle = ":", label = f"High quality: {high:g}%")
        ax.set(xlabel = label, ylabel = "Retained loci", xlim = (min(minimum, float(data.min())) - 2, 102))
        ax.legend(fontsize = 8)
        ax.grid(axis = "y", alpha = 0.2)
        ax.set_axisbelow(True)
    fig.suptitle(f"{config.gene_name} - BLAST identity and coverage of retained loci (n={len(extracted_sequences):,})", fontsize = 13)
    return save_figure(fig, Path(output_dir) / f"{config.gene_name}_retained_loci_quality.png", "retained loci quality")

#10 - Best identity for each isolate against each alternative reference query
def plot_query_heatmap(blast_outputs, config, isolate_ids, output_dir):
    queries = list(active_queries(config.gene_queries))
    if blast_outputs.get("gene") is None or len(queries) < 2 or not isolate_ids:
        return None
    best = defaultdict(dict)
    for row in read_search_rows({"gene": blast_outputs["gene"]}, config):
        if row["passed"]:
            current = best[row["isolate_id"]].get(row["query_id"], -1)
            best[row["isolate_id"]][row["query_id"]] = max(current, row["identity"])
    isolates = sorted({str(i) for i in isolate_ids}, key = isolate_sort_key)
    matrix = np.full((len(isolates), len(queries)), np.nan)
    for i, isolate in enumerate(isolates):
        for j, query in enumerate(queries):
            matrix[i, j] = best.get(isolate, {}).get(query, np.nan)

    fig = Figure(figsize = (max(7, 1.2 * len(queries) + 3), min(40, max(5, 0.2 * len(isolates) + 2))), layout = "constrained")
    FigureCanvasAgg(fig)
    ax = fig.subplots()
    cmap = plt.get_cmap("viridis").copy()
    cmap.set_bad("#eeeeee")
    image = ax.imshow(np.ma.masked_invalid(matrix), aspect = "auto", cmap = cmap, vmin = config.gene_blast_params["min_identity"], vmax = 100)
    fig.colorbar(image, ax = ax, label = "Best passing BLAST identity (%)")
    ax.set_xticks(range(len(queries)))
    ax.set_xticklabels(queries, rotation = 30, ha = "right", fontsize = 8)
    if len(isolates) <= 80:
        ax.set_yticks(range(len(isolates)))
        ax.set_yticklabels(isolates, fontsize = 7)
    ax.set(xlabel = "Reference query", ylabel = "Isolate", title = "Best passing hit per isolate (grey = no passing hit)")
    fig.suptitle(f"{config.gene_name} - isolates vs alternative reference queries", fontsize = 13)
    return save_figure(fig, Path(output_dir) / f"{config.gene_name}_query_heatmap.png", "query heatmap")

#11 - Do the primers hit the gene locus?
def plot_primer_locus(blast_results, config, output_dir):
    #aggregate view so it still works with lots of isolates, positions are in gene orientation (0 = start of gene)
    blast_results = blast_results or {}
    if "primers" not in blast_results.get("searches_run", set()):
        return None
    loci = blast_results.get("gene_hits_before_filter", [])
    primer_hits = blast_results.get("forward_primer_hits", []) + blast_results.get("reverse_primer_hits", [])
    if not loci:
        return None
    buffer = config.flank_buffer
    names = list(active_queries(config.primer_queries))
    palette = ["#9da6cf", "#d18c20", "#70b6a1", "#e9c46a"]
    colours = {name: palette[i % len(palette)] for i, name in enumerate(names)}

    by_contig = defaultdict(list)
    for hit in primer_hits:
        by_contig[(hit["isolate_id"], hit["contig_id"])].append(hit)

    positions = defaultdict(list) #primer -> (start, end, difference near 3' end)
    locus_status = {name: Counter() for name in names}
    outside = Counter()
    lengths = []
    for locus in loci:
        low, high = sorted((locus["s_start"], locus["s_end"]))
        lengths.append(high - low + 1)
        nearby = defaultdict(list)
        for hit in by_contig.get((locus["isolate_id"], locus["contig_id"]), []):
            p_low, p_high = sorted((hit["s_start"], hit["s_end"]))
            if max(0, p_low - high, low - p_high) > buffer:
                outside[hit["query_id"]] += 1
                continue
            if locus["strand"] == "forward":
                start, end = p_low - low, p_high - low + 1
            else:
                start, end = high - p_high, high - p_low + 1
            near_3prime = any(int(v) <= 5 for v in hit["positions_from_3prime"].split(";") if v)
            positions[hit["query_id"]].append((start, end, near_3prime))
            nearby[hit["query_id"]].append(hit["total_mismatches"])
        for name in names:
            if not nearby[name]:
                locus_status[name]["No hit within window"] += 1
            elif min(nearby[name]) == 0:
                locus_status[name]["Perfect match"] += 1
            else:
                locus_status[name]["Best hit has differences"] += 1

    median_length = float(np.median(lengths))
    fig = Figure(figsize = (14, 6), layout = "constrained")
    FigureCanvasAgg(fig)
    ax, bar_ax = fig.subplots(1, 2, gridspec_kw = {"width_ratios": [2.2, 1]})
    bins = np.arange(-buffer, median_length + buffer + 10, 10)
    top = 1
    for name in names:
        if not positions[name]:
            continue
        starts = [p[0] for p in positions[name]]
        counts, _, _ = ax.hist(starts, bins = bins, histtype = "stepfilled", alpha = 0.6, color = colours[name],
                               label = f"{name} ({config.primer_queries[name]['position']}, n={len(starts):,} hits)")
        top = max(top, counts.max())
        flagged = [p[0] for p in positions[name] if p[2]]
        if flagged:
            flagged_counts, _ = np.histogram(flagged, bins = bins)
            keep = flagged_counts > 0
            ax.bar(bins[:-1][keep], flagged_counts[keep], width = 10, align = "edge", facecolor = "none", edgecolor = "#d94c4c", linewidth = 1.2)
    ax.add_patch(plt.Rectangle((0, -top * 0.12), median_length, top * 0.09, color = "#9da6cf"))
    ax.text(median_length / 2, -top * 0.075, f"gene span (median {median_length:,.0f} bp)", ha = "center", va = "center", fontsize = 9)
    ax.plot([], [], color = "#d94c4c", label = "hits with a difference in last 5 bases (outline)")
    ax.axvline(0, color = "#555555", linewidth = 0.8)
    ax.axvline(median_length, color = "#555555", linewidth = 0.8)
    ax.set(xlabel = "Primer start in gene orientation (bp from start of gene, 10 bp bins)", ylabel = "Primer hits",
           title = f"Where primers bind around the gene (within {buffer:g} bp)", ylim = (-top * 0.15, top * 1.15))
    ax.set_yticks([t for t in ax.get_yticks() if t >= 0])
    ax.legend(fontsize = 8)
    ax.grid(axis = "y", alpha = 0.2)
    ax.set_axisbelow(True)

    statuses = [("Perfect match", "#70b6a1"), ("Best hit has differences", "#e9c46a"), ("No hit within window", "#d94c4c")]
    left = np.zeros(len(names))
    for status, colour in statuses:
        values = np.array([locus_status[n][status] / len(loci) * 100 for n in names])
        bar_ax.barh(names, values, left = left, color = colour, edgecolor = "white", label = status)
        for y, (value, base) in enumerate(zip(values, left)):
            if value >= 4:
                bar_ax.text(base + value / 2, y, f"{value:.0f}%", ha = "center", va = "center", fontsize = 9)
        left += values
    bar_ax.set(xlim = (0, 100), xlabel = f"% of candidate loci (n={len(loci):,})", title = "Primer coverage per locus")
    bar_ax.legend(loc = "upper center", bbox_to_anchor = (0.5, -0.14), ncol = 1, fontsize = 8, frameon = False)
    note = ", ".join(f"{n}: {c:,} hit(s) >{buffer:g} bp from a locus not shown" for n, c in sorted(outside.items()))
    fig.suptitle(f"{config.gene_name} - primer positions around candidate loci" + (f"\n{note}" if note else ""), fontsize = 12)
    return save_figure(fig, Path(output_dir) / f"{config.gene_name}_primer_locus.png", "primer locus")

#9 - in silico PCR fragment lengths
def plot_pcr_fragment_lengths(blast_results, config, output_dir):
    fragments = (blast_results or {}).get("pcr_fragments", [])
    if not fragments:
        return None

    lengths = np.asarray([f["length"] for f in fragments], dtype = float)
    differences = np.asarray([f["forward_primer_differences"] + f["reverse_primer_differences"] for f in fragments])
    relations = np.asarray([f["gene_relation"] for f in fragments])
    bins = np.arange(lengths.min() - 5, lengths.max() + 15, 10)

    fig = Figure(figsize = (15, 5.5), layout = "constrained")
    FigureCanvasAgg(fig)
    axes = fig.subplots(1, 2)

    relation_names = sorted(set(relations))
    palette = ["#9da6cf", "#d18c20", "#70b6a1", "#e9c46a"]
    axes[0].hist([lengths[relations == r] for r in relation_names], bins = bins, stacked = True, color = palette[:len(relation_names)],
                 edgecolor = "white", linewidth = 0.3, label = [f"{r} (n={int((relations == r).sum()):,})" for r in relation_names])
    axes[0].set(xlabel = "Fragment length (bp)", ylabel = "Fragments", title = "Fragment length by overlap with the gene")

    groups = [("No differences", differences == 0, "#70b6a1"), ("1-2 differences", (differences >= 1) & (differences <= 2), "#e9c46a"), ("3+ differences", differences >= 3, "#d94c4c")]
    axes[1].hist([lengths[mask] for _, mask, _ in groups], bins = bins, stacked = True, color = [c for _, _, c in groups],
                 edgecolor = "white", linewidth = 0.3, label = [f"{name} (n={int(mask.sum()):,})" for name, mask, _ in groups])
    axes[1].set(xlabel = "Fragment length (bp)", ylabel = "Fragments", title = "Fragment length by total primer differences")

    for ax in axes:
        ax.axvline(float(np.median(lengths)), color = "#d94c4c", linestyle = ":", label = f"Median: {np.median(lengths):.0f} bp")
        max_length = config.primer_blast_params.get("max_pcr_fragment_length")
        if max_length and max_length <= lengths.max() + 50:
            ax.axvline(max_length, color = "#d18c20", linestyle = "--", label = f"max_pcr_fragment_length ({max_length})")
        ax.legend(fontsize = 8)
        ax.grid(alpha = 0.2)
        ax.set_axisbelow(True)

    isolates = len({f["isolate_id"] for f in fragments})
    fig.suptitle(f"{config.gene_name} - in silico PCR fragment lengths (n={len(fragments):,} fragments, {isolates:,} isolates)", fontsize = 13)
    return save_figure(fig, Path(output_dir) / f"{config.gene_name}_pcr_fragment_lengths.png", "PCR fragment length")

#one function to make all plots

def generate_all_plots(config, blast_outputs, extracted_sequences, before_N_filtering, output_dir, blast_results = None):
    plots_generated = {
        "raw_blast": plot_raw_blast_results(blast_outputs, config, output_dir),
        "query_comparison": plot_query_comparison(blast_outputs, config, output_dir),
        "flank_assignment": plot_flank_assignment(blast_results, config, output_dir),
        "query_heatmap": plot_query_heatmap(blast_outputs, config, (blast_results or {}).get("isolate_ids", []), output_dir),
        "copy_number": plot_copy_number(blast_results, extracted_sequences, config, output_dir),
        "primer_locus": plot_primer_locus(blast_results, config, output_dir),
        "primer_mismatch": plot_primer_mismatch(blast_outputs, config, output_dir),
        "pcr_fragment_lengths": plot_pcr_fragment_lengths(blast_results, config, output_dir),
        "sequence_extraction": plot_sequence_extraction(before_N_filtering, extracted_sequences, config, output_dir),
        "span_vs_reference": plot_span_vs_reference(before_N_filtering, config, output_dir),
        "retained_loci_quality": plot_retained_quality(extracted_sequences, config, output_dir),
        "contig_boundary_QC": plot_contig_boundaries(extracted_sequences, config, output_dir)
    }

    return plots_generated
