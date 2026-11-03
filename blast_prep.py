import io
import logging
import zipfile
import subprocess
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from pathlib import Path

logger = logging.getLogger(__name__)

fasta_extensions = (".fasta", ".fa", ".fna", ".fas")

def create_combined_genome_fasta(genome_dir, output_dir):
    #combine all genomes into single fasta file for blastdb
    combined_fasta = Path(output_dir) / "combined_genomes.fasta"
    genome_dir = Path(genome_dir)

    contig_id_map = {}
    isolate_ids = set()
    n_contigs = 0
    n_isolates = 0

    with open(combined_fasta, "w") as combined_handle:
        if genome_dir.suffix == ".zip":
            with zipfile.ZipFile(genome_dir, "r") as zf:
                genome_entries = [
                    entry
                    for entry in zf.namelist()
                    if entry.lower().endswith(fasta_extensions)
                    and not entry.rsplit("/", 1)[-1].startswith("._")
                ]

                logger.info("Combining %d genome entries from %s", len(genome_entries), genome_dir)

                for entry in genome_entries:
                    isolate_id = Path(entry).stem
                    isolate_ids.add(isolate_id)
                    n_isolates += 1

                    with zf.open(entry) as f:
                        text_handle = io.TextIOWrapper(f, encoding = "utf-8")
                        for contig_no, record in enumerate(SeqIO.parse(text_handle, "fasta"), start = 1): 
                            original_id = record.id
                            internal_id = f"{isolate_id}_ctg{contig_no:06d}"
                            contig_id_map[internal_id] = original_id
                            record.id = internal_id
                            record.name = internal_id
                            record.description = ""
                            SeqIO.write(record, combined_handle, "fasta")
                            n_contigs += 1

        elif Path(genome_dir).is_dir():
            genome_files = sorted(
                path
                for path in Path(genome_dir).iterdir()
                if path.suffix.lower() in fasta_extensions
            )

            logger.info("Combining %d genome files from %s", len(genome_files), genome_dir)

            for genome_file in genome_files:
                isolate_id = genome_file.stem
                isolate_ids.add(isolate_id)
                n_isolates += 1

                with open(genome_file, "r", encoding="utf-8") as f:
                    for contig_no, record in enumerate(SeqIO.parse(f, "fasta"), start=1):
                        original_id = record.id
                        internal_id = f"{isolate_id}_ctg{contig_no:06d}"
                        contig_id_map[internal_id] = original_id
                        record.id = internal_id
                        record.name = internal_id
                        record.description = ""
                        SeqIO.write(record, combined_handle, "fasta")
                        n_contigs += 1

            logger.info("Combined %d contigs from %d isolates into %s", n_contigs, n_isolates, combined_fasta)

    if n_contigs == 0:
        raise FileNotFoundError(f"No genome FASTA files found in {genome_dir}")

    return combined_fasta, contig_id_map, sorted(isolate_ids)

def make_blast_db(fasta_path):
    logger.info("Creating BLAST database for %s", fasta_path)

    db_prefix = Path(fasta_path).parent / "blast_db"

    cmd = [
        "makeblastdb",
        "-in", str(fasta_path),
        "-dbtype", "nucl",
        "-out", str(db_prefix)
    ]

    try:
        result = subprocess.run(cmd, check = True, capture_output = True, text = True)
        logger.info("BLAST database created successfully")
        return db_prefix
    except subprocess.CalledProcessError as e:
        logger.error("Failed to create BLAST database: %s", e.stderr)
        raise

def clean_sequence(sequence, sequence_name):
    #remove whitespace, convert to uppercase and check for invalid characters
    sequence = "".join(sequence.split()).upper()
    valid_bases = set("ACGTRYSWKMBDHVN") #all IUPAC DNA coddes
    invalid_bases = set(sequence) - valid_bases
    if not sequence:
        raise ValueError(f"Sequence for {sequence_name} is empty")

    if invalid_bases:
        raise ValueError(f"Sequence for {sequence_name} contains invalid bases: {', '.join(invalid_bases)}")

    return sequence


def active_queries(queries):
    #queries with a sequence; empty template placeholders are ignored
    return {name: data for name, data in queries.items() if data.get("query_sequence", "").strip()}


def create_reference_queries(queries, output_dir, filename):
    reference_queries = []
    for query_name, query_data in queries.items():
        query_sequence = clean_sequence(query_data["query_sequence"], query_name)
        reference_queries.append(SeqRecord(Seq(query_sequence), id=query_name, description=""))

    if not reference_queries:
        return [], None

    ref_fasta_path = Path(output_dir) / filename
    SeqIO.write(reference_queries, ref_fasta_path, "fasta")
    return reference_queries, ref_fasta_path