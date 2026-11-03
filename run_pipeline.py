import os
import sys
import logging
import shutil
import argparse
import importlib.util
from collections import Counter
from pathlib import Path

from blast_prep import create_combined_genome_fasta, make_blast_db, create_reference_queries
from blast_runner import run_all_blast_searches
from blast_parser import parse_blast_outputs
from extract_sequence import extract_sequences
from reporting import generate_all_reports
from plots import generate_all_plots

def create_output_dirs(config, base_output_dir):
    gene_name = config.gene_name
    identifier = config.identifier

    output_label = identifier if identifier else gene_name

    base_output_path = Path(base_output_dir) / output_label

    dirs = {
        "base": base_output_path,
        "blast_temp": Path(base_output_dir) / "blast_temp",
        "blast": base_output_path / "blast",
        "sequences": base_output_path / "sequences",
        "reports": base_output_path / "reports",
        "plots": base_output_path / "plots"
    }

    dirs["multicopy"] = dirs["sequences"] / "multicopy"
    dirs["pcr_fragments"] = dirs["sequences"] / "pcr_fragments"

    if config.extraction_params.get("extract_igr", False):
        dirs["igr"] = dirs["sequences"] / "igr"

    for dir_path in dirs.values():
        os.makedirs(dir_path, exist_ok=True)

    return dirs



def setup_logging(output_dir, gene_name):
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    log_file = output_dir / f"{gene_name}_pipeline.log"

    logging.basicConfig(
        level = logging.INFO,
        format = "%(asctime)s [%(levelname)s] %(message)s",
        handlers = [logging.FileHandler(log_file), logging.StreamHandler(sys.stdout)]
    )

    return logging.getLogger("pipeline_logger")

logger = logging.getLogger(__name__)





def load_config(config_path):
    config_path = Path(config_path)

    if not config_path.exists():
        raise FileNotFoundError(f"Config file not found: {config_path}")

    spec = importlib.util.spec_from_file_location(config_path.stem, config_path)
    if spec is None or spec.loader is None:
        raise ImportError(f"Could not load config file from {config_path}")

    config = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(config)
    return config


def get_configs_to_run(gene_args, config_dir="configs"):
    #one config per gene so several genes can share the same combined FASTA and BLAST database
    config_dir = Path(config_dir)

    if "all" in gene_args:
        return sorted(
            p for p in config_dir.glob("config_*.py")
            if p.name not in ("config_template.py", "config_example.py")
        )

    config_paths = []
    for gene in gene_args:
        config_path = config_dir / f"config_{gene}.py"
        if not config_path.exists():
            raise FileNotFoundError(f"Config file not found: {config_path}")
        config_paths.append(config_path)
    return config_paths


def run_gene(config, dirs, combined_fasta, contig_id_map, isolate_ids, db_path):
    pipeline_logger = logging.getLogger("pipeline_logger")
    pipeline_logger.info("Running pipeline for gene: %s", config.gene_name)

    #query groups; optional groups configured as {} are skipped
    query_groups = {
        "gene": config.gene_queries,
        "flank": config.flank_queries,
        "gene_start": config.gene_start_queries,
        "gene_end": config.gene_end_queries,
        "primer": config.primer_queries,
    }
    if not config.gene_queries:
        raise ValueError("gene_queries must contain at least one reference sequence")

    #optional groups: skip template placeholders with no sequence; gene_queries stays strict
    for group, queries in query_groups.items():
        if group == "gene":
            continue
        blank = [name for name, data in queries.items() if not data["query_sequence"].strip()]
        if blank:
            pipeline_logger.warning("Skipping %s queries with empty sequence: %s", group, ", ".join(blank))
            query_groups[group] = {name: data for name, data in queries.items() if name not in blank}

    query_fastas = {}
    for group, queries in query_groups.items():
        _, query_fastas[group] = create_reference_queries(queries, dirs["blast"], f"{config.gene_name}_{group}_queries.fasta")

    blast_outputs = run_all_blast_searches(
        config, query_fastas["gene"], query_fastas["flank"], query_fastas["gene_start"],
        query_fastas["gene_end"], query_fastas["primer"], db_path, dirs["blast"]
    )

    parsed = parse_blast_outputs(blast_outputs, config, contig_id_map)
    extracted, igrs, extraction_info = extract_sequences(parsed["gene_hits"], combined_fasta, dirs["sequences"], config)

    pipeline_logger.info(
        "Isolates searched: %d; isolates with a retained locus: %d; loci retained: %d",
        len(isolate_ids), len({seq["isolate_id"] for seq in extracted}), len(extracted)
    )
    before_n_filtering = extraction_info["extracted_before_n_filter"]
    parsed["contig_counts"] = Counter(cid.rsplit("_ctg", 1)[0] for cid in contig_id_map)
    parsed["searches_run"] = {search for search, path in blast_outputs.items() if path is not None}
    parsed["isolate_ids"] = list(isolate_ids)

    report_paths = generate_all_reports(config, parsed, dirs["reports"], isolate_ids, blast_outputs, extracted, before_n_filtering)
    plot_paths = generate_all_plots(config, blast_outputs, extracted, before_n_filtering, dirs["plots"], parsed)
    return {"parsed": parsed, "extracted": extracted, "igrs": igrs, "isolate_ids": isolate_ids,
            "extraction_info": extraction_info, "report_paths": report_paths, "plot_paths": plot_paths, "dirs": dirs}


def parse_args():
    parser = argparse.ArgumentParser(description="BLAST, parse, extract and report for one or more genes")
    parser.add_argument("--gene", nargs="+", required=True,
                        help="gene name(s) matching configs/config_<gene>.py; use 'all' for every config")
    parser.add_argument("--input", required=True, help="genome FASTA folder or .zip")
    parser.add_argument("--output", required=True, help="output folder")
    parser.add_argument("--keep_blast", action="store_true", help="keep the combined FASTA and BLAST database")
    return parser.parse_args()


def main():
    args = parse_args()

    configs = [load_config(p) for p in get_configs_to_run(args.gene)]
    if not configs:
        raise FileNotFoundError("No gene config files found to run")

    base_output_dir = Path(args.output)
    setup_logging(base_output_dir, "run")
    pipeline_logger = logging.getLogger("pipeline_logger")
    pipeline_logger.info("Starting pipeline for %s", ", ".join(c.gene_name for c in configs))

    output_dirs = [create_output_dirs(config, base_output_dir) for config in configs]

    blast_temp = base_output_dir / "blast_temp"
    combined_fasta, contig_id_map, isolate_ids = create_combined_genome_fasta(args.input, blast_temp)
    db_path = make_blast_db(combined_fasta)

    all_succeeded = False
    try:
        for config, dirs in zip(configs, output_dirs):
            run_gene(config, dirs, combined_fasta, contig_id_map, isolate_ids, db_path)
        all_succeeded = True
    finally:
        #keep the database after a failure so the run can be inspected
        if all_succeeded and not args.keep_blast:
            shutil.rmtree(blast_temp)
            pipeline_logger.info("Removed temporary BLAST directory: %s", blast_temp)


if __name__ == "__main__":
    main()
