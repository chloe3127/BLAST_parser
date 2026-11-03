import logging
import subprocess
from pathlib import Path

logger = logging.getLogger(__name__)

gene_outfmt = "6 qseqid sseqid nident pident length mismatch gaps gapopen qstart qend sstart send evalue bitscore qlen slen"
primer_outfmt = "6 qseqid sseqid nident pident length mismatch gaps gapopen qstart qend sstart send evalue bitscore qlen slen qseq sseq"

def run_blast(cmd, search_name, output_file):
    #run the blast cmd and return output path
    logger.info("Running BLAST search for %s", search_name)
    logger.info("Command: %s:", " ".join(cmd))

    try:
        subprocess.run(cmd, check = True, capture_output = True, text = True)
    except FileNotFoundError:
        logger.error("blastn was not found. Check that BLAST+ is installed and available in this environment")
        raise
    except subprocess.CalledProcessError as e:
        logger.error("BLAST search for %s failed: %s", search_name, e.stderr)
        raise

    logger.info("%s BLAST search completed successfully. Raw results saved to %s", search_name, output_file)

    return Path(output_file)

def run_gene_blast(config, queries, db_path, output_file, query_type="gene"):
    params = getattr(config, f"{query_type}_blast_params")
    search_name = f"{query_type}_blast"

    cmd = [
        "blastn",
        "-task", str(params["task"]),
        "-query", str(queries),
        "-db", str(db_path),
        "-out", str(output_file),
        "-outfmt", gene_outfmt,
        "-max_target_seqs", str(params["max_target_seqs"]),
        "-max_hsps", str(params["max_hsps"]),
        "-word_size", str(params["word_size"]),
        "-reward", str(params["reward"]),
        "-penalty", str(params["penalty"]),
        "-gapopen", str(params["gapopen"]),
        "-gapextend", str(params["gapextend"]),
        "-evalue", str(params["evalue"])
    ]

    return run_blast(cmd, search_name, output_file)

def run_gene_start_blast(config, queries, db_path, output_file):
    params = config.gene_start_blast_params
    search_name = "gene_start_blast"

    cmd = [
        "blastn",
        "-task", str(params["task"]),
        "-query", str(queries),
        "-db", str(db_path),
        "-out", str(output_file),
        "-outfmt", gene_outfmt,
        "-max_target_seqs", str(params["max_target_seqs"]),
        "-max_hsps", str(params["max_hsps"]),
        "-word_size", str(params["word_size"]),
        "-reward", str(params["reward"]),
        "-penalty", str(params["penalty"]),
        "-gapopen", str(params["gapopen"]),
        "-gapextend", str(params["gapextend"]),
        "-evalue", str(params["evalue"])
    ]

    return run_blast(cmd, search_name, output_file)

def run_gene_end_blast(config, queries, db_path, output_file):
    params = config.gene_end_blast_params
    search_name = "gene_end_blast"

    cmd = [
        "blastn",
        "-task", str(params["task"]),
        "-query", str(queries),
        "-db", str(db_path),
        "-out", str(output_file),
        "-outfmt", gene_outfmt,
        "-max_target_seqs", str(params["max_target_seqs"]),
        "-max_hsps", str(params["max_hsps"]),
        "-word_size", str(params["word_size"]),
        "-reward", str(params["reward"]),
        "-penalty", str(params["penalty"]),
        "-gapopen", str(params["gapopen"]),
        "-gapextend", str(params["gapextend"]),
        "-evalue", str(params["evalue"])
    ]

    return run_blast(cmd, search_name, output_file)

def run_primer_blast(config, queries, db_path, output_file):    
    params = config.primer_blast_params
    search_name = "primer_blast"

    cmd = [
        "blastn",
        "-task", str(params["task"]),
        "-query", str(queries),
        "-db", str(db_path),
        "-out", str(output_file),
        "-outfmt", primer_outfmt,
        "-max_target_seqs", str(params["max_target_seqs"]),
        "-max_hsps", str(params["max_hsps"]),
        "-word_size", str(params["word_size"]),
        "-reward", str(params["reward"]),
        "-penalty", str(params["penalty"]),
        "-gapopen", str(params["gapopen"]),
        "-gapextend", str(params["gapextend"]),
        "-evalue", str(params["evalue"])
    ]

    return run_blast(cmd, search_name, output_file)

def run_all_blast_searches(config, gene_ref_fasta, flank_ref_fasta, gene_start_ref_fasta, gene_end_ref_fasta, primer_ref_fasta, db_path, output_dir):
    #run all configured BLAST searches for gene
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    gene_name = config.gene_name

    blast_outputs = {
        "gene": None,
        "flank": None,
        "gene_start": None,
        "gene_end": None,
        "primers": None,
    }

    if gene_ref_fasta is not None:
        blast_outputs["gene"] = run_gene_blast(
            config, gene_ref_fasta, db_path, output_dir / f"{gene_name}_gene_blast.tsv"
        )
    if flank_ref_fasta is not None:
        blast_outputs["flank"] = run_gene_blast(
            config, flank_ref_fasta, db_path, output_dir / f"{gene_name}_flank_blast.tsv", query_type="flank"
        )
    if gene_start_ref_fasta is not None:
        blast_outputs["gene_start"] = run_gene_start_blast(
            config, gene_start_ref_fasta, db_path, output_dir / f"{gene_name}_gene_start_blast.tsv"
        )
    if gene_end_ref_fasta is not None:
        blast_outputs["gene_end"] = run_gene_end_blast(
            config, gene_end_ref_fasta, db_path, output_dir / f"{gene_name}_gene_end_blast.tsv"
        )
    if primer_ref_fasta is not None:
        blast_outputs["primers"] = run_primer_blast(
            config, primer_ref_fasta, db_path, output_dir / f"{gene_name}_primer_blast.tsv"
        )

    logger.info("Completed all BLAST searches: %s", blast_outputs)
    return blast_outputs