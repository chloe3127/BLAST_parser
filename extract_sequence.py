import logging
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

logger = logging.getLogger(__name__)

def load_genome_sequences(combined_fasta):
    #load the combined genome FASTA used to build the BLAST database, keyed by the BLAST subject ID
    logger.info(f"Loading genome sequences from {combined_fasta}")

    combined_path = Path(combined_fasta)
    if not combined_path.is_file():
        raise FileNotFoundError(f"Combined genome FASTA not found: {combined_fasta}")

    return {record.id: record for record in SeqIO.parse(combined_path, "fasta")}

def blast_to_zero_based(s_start, s_end):
    #BLAST 1-based inclusive coordinates -> 0-based, end-exclusive (left_0, right_0)
    return min(s_start, s_end) - 1, max(s_start, s_end)

def extract_sequence_no_flank(contig_seq, left_0, right_0, strand):
    #extract the aligned span; left_0/right_0 are 0-based, end-exclusive
    contig_length = len(contig_seq.seq)

    extract_start = max(0, left_0)
    extract_end = min(contig_length, right_0)

    extracted_seq = contig_seq.seq[extract_start:extract_end]

    if strand == "reverse":
        extracted_seq = extracted_seq.reverse_complement()

    dist_to_contig_start = left_0
    dist_to_contig_end = contig_length - right_0

    return {
        "sequence": str(extracted_seq),
        "contig_length": contig_length,
        "strand": strand,
        "extract_start": extract_start,
        "extract_end": extract_end,
        "gene_left_0": left_0,
        "gene_right_0": right_0,
        "dist_to_contig_start": dist_to_contig_start,
        "dist_to_contig_end": dist_to_contig_end
    }

def extract_sequence_with_flank(contig_seq, left_0, right_0, strand, flanking_bp):
    #extract the aligned span plus flanking_bp each side, clipped at contig edges; left_0/right_0 are 0-based, end-exclusive
    contig_length = len(contig_seq.seq)

    extract_start = max(0, left_0 - flanking_bp)
    extract_end = min(contig_length, right_0 + flanking_bp)

    extracted_seq = contig_seq.seq[extract_start:extract_end]

    if strand == "reverse":
        extracted_seq = extracted_seq.reverse_complement()

    dist_to_contig_start = left_0
    dist_to_contig_end = contig_length - right_0

    return {
        "sequence": str(extracted_seq),
        "contig_length": contig_length,
        "extract_start": extract_start,
        "extract_end": extract_end,
        "gene_left_0": left_0,
        "gene_right_0": right_0,
        "strand": strand,
        "dist_to_contig_start": dist_to_contig_start,
        "dist_to_contig_end": dist_to_contig_end
    }


def extract_igr_from_coords(contig_seq, gene_left_0, gene_right_0, flank_coords, igr_type, strand):
    #extract intergenic region between gene and flanking gene
    #gene and flank coordinates are both 0-based, end-exclusive (the parser already converts flank coordinates)
    flank_left_0, flank_right_0 = flank_coords
    contig_length = len(contig_seq.seq)

    if strand not in ("forward", "reverse"):
        raise ValueError(f"Invalid strand: {strand}")
    if igr_type not in ("upstream", "downstream"):
        raise ValueError(f"Invalid IGR type: {igr_type}")
    if not (0 <= gene_left_0 < gene_right_0 <= contig_length):
        raise ValueError("Gene coordinates are outside the contig or invalid")
    if not (0 <= flank_left_0 < flank_right_0 <= contig_length):
        raise ValueError("Flank coordinates are outside the contig or invalid")

    if strand == "forward":
        if igr_type == "upstream":
            start, end = flank_right_0, gene_left_0
        else:
            start, end = gene_right_0, flank_left_0
    else:
        if igr_type == "upstream":
            start, end = gene_right_0, flank_left_0
        else:
            start, end = flank_right_0, gene_left_0

    if start >= end:
        return None, 0

    igr_seq = contig_seq.seq[start:end]
    if strand == "reverse":
        igr_seq = igr_seq.reverse_complement()

    return str(igr_seq), len(igr_seq)

def write_igr_fastas(upstream_igrs, downstream_igrs, output_dir, config):
    gene_name = config.gene_name
    Path(output_dir).mkdir(parents=True, exist_ok=True)

    for igr_type, igrs in (("upstream", upstream_igrs), ("downstream", downstream_igrs)):
        if not igrs:
            continue

        records = []
        for igr in igrs:
            header_parts = [
                f"isolate_id={igr['isolate_id']}",
                f"gene={gene_name}",
                f"region={igr_type}_IGR",
                f"contig={igr['contig_id']}",
                f"length={len(igr['sequence'])}bp",
                f"strand={igr['strand']}"
            ]
            record = SeqRecord(
                Seq(igr["sequence"]),
                id = f"{igr['isolate_id']}_{igr['gene_copy']}",
                description = " ".join(header_parts)
            )
            records.append(record)

        path = Path(output_dir) / f"{gene_name}_{igr_type}_IGR.fasta"
        SeqIO.write(records, path, "fasta")
        logger.info(f"Wrote {len(records)} {igr_type} IGR sequences to {path}")


def check_sequence_quality(seq_data, boundary_critical, boundary_warning):
    #check if sequence is too close to contig boundaries
    flags = []
    boundary_status = "OK"

    if (seq_data["dist_to_contig_start"] <= boundary_critical or seq_data["dist_to_contig_end"] <= boundary_critical):
        boundary_status = "CRITICAL"
        flags.append(f"Boundary CRITICAL: Gene within {boundary_critical}bp of contig edge")

    elif (seq_data["dist_to_contig_start"] <= boundary_warning or seq_data["dist_to_contig_end"] <= boundary_warning):
        boundary_status = "WARNING"
        flags.append(f"Boundary WARNING: Gene within {boundary_warning}bp of contig edge")

    return {
        "boundary_status": boundary_status,
        "flags": flags
    }

def remove_ambiguous_bases(extracted_seqs):
    #drop a locus if any sequence written for it (aligned span or padded) contains an N, so all outputs keep the same loci
    total = len(extracted_seqs)
    filtered = []
    for s in extracted_seqs:
        s["n_in_aligned_span"] = 'N' in s["sequence"].upper()
        s["n_in_padding_only"] = (
            not s["n_in_aligned_span"] and "padded" in s and 'N' in s["padded"]["sequence"].upper()
        )
        if s["n_in_aligned_span"]:
            s["outcome"] = "excluded_N_in_aligned_span"
        elif s["n_in_padding_only"]:
            s["outcome"] = "excluded_N_in_padding_only"
        else:
            s["outcome"] = "retained"
            filtered.append(s)
    n_removed = total - len(filtered)

    logger.info(f"Sequences before filtering for ambiguous bases: {total}")
    logger.info(f"Sequences removed due to ambiguous bases: {n_removed}")
    logger.info(f"Sequences after filtering for ambiguous bases: {len(filtered)}")
    return filtered

def write_extracted_fasta(sequences, output_dir, config, padded=False):
    #write the BLAST-aligned span (may be only part of the gene) or the fixed-length padded sequence
    gene_name = config.gene_name
    if padded:
        label = f"{config.extraction_params['flanking_bp']}bp_padded"
        region = "padded"
    else:
        label = "aligned_span"
        region = "aligned_span"

    all_records = []
    for seq in sequences:
        data = seq["padded"] if padded else seq
        header_parts = [
            f"isolate={seq['isolate_id']}",
            f"copy={seq['gene_copy']}",
            f"contig={seq['contig_id']}",
            f"region={region}",
            f"coords={data['extract_start']}-{data['extract_end']}",
            "coordinate_system=0-based_end-exclusive",
            f"strand={seq['strand']}",
            f"query={seq['blast_query']}",
            f"identity={seq['blast_identity']:.1f}%",
            f"query_coverage={seq['blast_coverage']:.1f}%",
            f"genomic_context={seq['genomic_context']}"
        ]

        record = SeqRecord(
            Seq(data["sequence"]),
            id = f"{seq['isolate_id']}_{seq['gene_copy']}",
            description = " ".join(header_parts)
        )
        all_records.append(record)

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok = True)
    all_fasta_path = output_dir / f"{gene_name}_{label}_all_sequences.fasta"
    SeqIO.write(all_records, all_fasta_path, "fasta")
    logger.info(f"{len(all_records)} {region} sequences written to {all_fasta_path}")

    #boundary issues
    for status in ("CRITICAL", "WARNING"):
        boundary_records = [r for r, s in zip(all_records, sequences) if s.get("boundary_status") == status]
        if boundary_records:
            boundary_path = output_dir / f"{gene_name}_{label}_boundary_{status.lower()}.fasta"
            SeqIO.write(boundary_records, boundary_path, "fasta")
            logger.info(f"{status} boundary sequences written to {boundary_path}")
    return all_fasta_path

def write_multicopy_fastas(sequences, output_dir, config):
    #group retained loci by isolate
    sequences_by_isolate = {}
    for seq in sequences:
        sequences_by_isolate.setdefault(seq["isolate_id"], []).append(seq)

    multicopy_isolates = {
        isolate_id: copies
        for isolate_id, copies in sequences_by_isolate.items()
        if len(copies) > 1
    }
    if not multicopy_isolates:
        logger.info("No isolates with multiple retained gene copies")
        return None

    output_dir = Path(output_dir)
    output_dir.mkdir(parents = True, exist_ok = True)
    all_records = []

    for isolate_id, copies in multicopy_isolates.items():
        records = []
        for seq in copies:
            header_parts = [
                f"isolate_id={seq['isolate_id']}",
                f"copy={seq['gene_copy']}",
                f"contig={seq['contig_id']}",
                f"coords={seq['extract_start']}-{seq['extract_end']}",
                "coordinate_system=0-based_end-exclusive",
                f"strand={seq['strand']}",
                f"boundary_status={seq['boundary_status']}"
            ]
            records.append(SeqRecord(
                Seq(seq["sequence"]),
                id = f"{seq['isolate_id']}_{seq['gene_copy']}",
                description = " ".join(header_parts)
            ))

        path = output_dir / f"{config.gene_name}_{isolate_id}_multicopy.fasta"
        SeqIO.write(records, path, "fasta")
        all_records.extend(records)

    combined_path = output_dir / f"{config.gene_name}_multicopy_all_sequences.fasta"
    SeqIO.write(all_records, combined_path, "fasta")
    logger.info(
        "Saved %d retained loci from %d multicopy isolates to %s",
        len(all_records), len(multicopy_isolates), output_dir
    )
    return combined_path

def extract_sequences(blast_results, combined_fasta, output_dir, config):
    #blast_results is the parsed gene_hits list; combined_fasta is the FASTA the BLAST database was built from
    extraction_params = config.extraction_params
    extract_flank = extraction_params["extract_flank"]
    flanking_bp = extraction_params["flanking_bp"]
    logger.info("Extracting sequence for %d blast results", len(blast_results))

    genome_seqs = load_genome_sequences(combined_fasta)
    extracted = []
    skipped = []

    for hit in blast_results:
        contig_id = hit["contig_id"]

        if contig_id not in genome_seqs:
            logger.warning("Contig %s not found in genome sequences, skipping", contig_id)
            skipped.append({
                "isolate_id": hit["isolate_id"], "gene_copy": hit["gene_copy"], "contig_id": contig_id,
                "strand": hit["strand"], "blast_query": hit.get("blast_query", "unknown"),
                "supporting_queries": hit.get("supporting_queries", ""),
                "blast_identity": hit["identity"], "blast_coverage": hit["coverage"],
                "query_length": hit["q_length"], "genomic_context": hit.get("genomic_context", "unknown"),
                "outcome": "skipped_contig_not_found"
            })
            continue

        left_0, right_0 = blast_to_zero_based(hit["s_start"], hit["s_end"])
        seq_no_flank = extract_sequence_no_flank(genome_seqs[contig_id], left_0, right_0, hit["strand"])

        no_flank_quality = check_sequence_quality(seq_no_flank, boundary_critical = extraction_params["boundary_critical_distance"], boundary_warning = extraction_params["boundary_warning_distance"])

        extracted_seq = {
            "isolate_id": hit["isolate_id"],
            "gene_copy": hit["gene_copy"],
            "contig_id": contig_id,
            "original_contig_id": hit.get("contig_orig", contig_id),
            "blast_identity": hit["identity"],
            "blast_coverage": hit["coverage"],
            "query_length": hit["q_length"],
            "sequence": seq_no_flank["sequence"],
            "strand": seq_no_flank["strand"],
            "boundary_status": no_flank_quality["boundary_status"],
            "flags": no_flank_quality["flags"],
            "contig_length": seq_no_flank["contig_length"],
            "gene_start": seq_no_flank["gene_left_0"],
            "gene_end": seq_no_flank["gene_right_0"],
            "extract_start": seq_no_flank["extract_start"],
            "extract_end": seq_no_flank["extract_end"],
            "dist_to_contig_start": seq_no_flank["dist_to_contig_start"],
            "dist_to_contig_end": seq_no_flank["dist_to_contig_end"],
            "genomic_context": hit.get("genomic_context", "unknown"),
            "blast_query": hit.get("blast_query", "unknown"),
            "supporting_queries": hit.get("supporting_queries", ""),
            "upstream_found": hit.get("upstream_found", False),
            "downstream_found": hit.get("downstream_found", False),
            "upstream_flank_coords": hit.get("upstream_flank_coords", None),
            "downstream_flank_coords": hit.get("downstream_flank_coords", None)
        }

        if extract_flank:
            seq_with_flank = extract_sequence_with_flank(genome_seqs[contig_id], left_0, right_0, hit["strand"], flanking_bp)
            extracted_seq["padded"] = {
                "sequence": seq_with_flank["sequence"],
                "extract_start": seq_with_flank["extract_start"],
                "extract_end": seq_with_flank["extract_end"]
            }

        extracted.append(extracted_seq)

    extracted_before_n_filter = list(extracted)
    extracted = remove_ambiguous_bases(extracted)

    if extracted:
        write_extracted_fasta(extracted, Path(output_dir), config)
        if extract_flank:
            write_extracted_fasta(extracted, Path(output_dir), config, padded=True)
        write_multicopy_fastas(extracted, Path(output_dir) / "multicopy", config)
    else:
        logger.warning("No sequences retained for extraction; no FASTA files written")

    upstream_igrs = []
    downstream_igrs = []

    if extraction_params.get("extract_igr", False):
        logger.info("Extracting intergenic regions (IGRs)...")
        for seq in extracted:
            contig_seq = genome_seqs[seq["contig_id"]]

            for igr_type, igrs in (("upstream", upstream_igrs), ("downstream", downstream_igrs)):
                flank_coords = seq.get(f"{igr_type}_flank_coords")
                if not seq[f"{igr_type}_found"] or flank_coords is None:
                    continue

                igr_seq, igr_len = extract_igr_from_coords(
                    contig_seq = contig_seq,
                    gene_left_0 = seq["gene_start"],
                    gene_right_0 = seq["gene_end"],
                    flank_coords = flank_coords,
                    igr_type = igr_type,
                    strand = seq["strand"]
                )
                if igr_seq:
                    igrs.append({
                        "isolate_id": seq["isolate_id"],
                        "gene_copy": seq["gene_copy"],
                        "contig_id": seq["contig_id"],
                        "sequence": igr_seq,
                        "length": igr_len,
                        "strand": seq["strand"]
                    })

        write_igr_fastas(upstream_igrs, downstream_igrs, Path(output_dir) / "igr", config)
        logger.info(f"IGR extraction complete: {len(upstream_igrs)} upstream, {len(downstream_igrs)} downstream IGRs")

    extraction_info = {
        "extracted_before_n_filter": extracted_before_n_filter,
        "skipped": skipped,
        "boundary_critical": extraction_params["boundary_critical_distance"],
        "boundary_warning": extraction_params["boundary_warning_distance"],
    }
    return extracted, {"upstream": upstream_igrs, "downstream": downstream_igrs}, extraction_info
