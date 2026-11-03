"""
In silico PCR fragments from the parsed primer hits.

A fragment is a forward primer hit and a reverse primer hit on the same contig, on opposite strands,
facing each other and no further apart than primer_blast_params["max_pcr_fragment_length"].
The fragment includes both primer binding sites. Only primer hits that passed the parser cutoffs are used.
"""

import logging
from collections import defaultdict
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from reporting import isolate_sort_key

logger = logging.getLogger(__name__)


def hit_range(hit):
    return tuple(sorted((hit["s_start"], hit["s_end"])))


def gene_relation(start, end, loci):
    #start/end are 1-based inclusive, compared with the gene loci on the same contig
    best, best_overlap = None, 0
    for locus in loci:
        l_low, l_high = hit_range(locus)
        overlap = min(end, l_high) - max(start, l_low) + 1
        if overlap > best_overlap:
            best, best_overlap = locus, overlap

    if best is None:
        return "no gene overlap", None, None

    l_low, l_high = hit_range(best)
    if start <= l_low and end >= l_high:
        relation = "contains gene"
    elif start >= l_low and end <= l_high:
        relation = "inside gene"
    else:
        relation = "overlaps part of gene"

    return relation, best.get("query_id"), best.get("genomic_context")


def find_pcr_fragments(blast_results, combined_fasta, config):
    forward_hits = blast_results.get("forward_primer_hits", [])
    reverse_hits = blast_results.get("reverse_primer_hits", [])
    if not forward_hits or not reverse_hits:
        return []

    max_length = config.primer_blast_params.get("max_pcr_fragment_length")

    primers_by_contig = defaultdict(list)
    for hit in forward_hits + reverse_hits:
        primers_by_contig[(hit["isolate_id"], hit["contig_id"])].append(hit)

    loci_by_contig = defaultdict(list)
    for locus in blast_results.get("gene_hits_before_filter", []):
        loci_by_contig[(locus["isolate_id"], locus["contig_id"])].append(locus)

    genome = SeqIO.index(str(combined_fasta), "fasta")
    fragments = []
    counts = defaultdict(int)

    try:
        for (isolate_id, contig_id), hits in sorted(primers_by_contig.items(), key = lambda item: (isolate_sort_key(item[0][0]), item[0][1])):
            #the left primer is on the + strand and the right one on the - strand, one of each primer
            lefts = [h for h in hits if h["strand"] == "forward"]
            rights = [h for h in hits if h["strand"] == "reverse"]

            for left in sorted(lefts, key = lambda h: hit_range(h)[0]):
                for right in sorted(rights, key = lambda h: hit_range(h)[0]):
                    if left["position"] == right["position"]:
                        continue

                    start, end = hit_range(left)[0], hit_range(right)[1]
                    if hit_range(left)[1] >= hit_range(right)[0]:
                        continue
                    length = end - start + 1
                    if max_length and length > max_length:
                        continue

                    if contig_id not in genome:
                        logger.warning("Contig %s not found in genome sequences, skipping PCR fragment", contig_id)
                        continue

                    sequence = genome[contig_id].seq[start - 1:end]

                    #write it from the forward primer, so flip it if the forward primer is the right-hand one
                    forward_primer, reverse_primer = (left, right) if left["position"] == "forward" else (right, left)
                    orientation = "forward" if left["position"] == "forward" else "reverse"
                    if orientation == "reverse":
                        sequence = sequence.reverse_complement()
                    sequence = str(sequence)

                    relation, gene_query, context = gene_relation(start, end, loci_by_contig.get((isolate_id, contig_id), []))

                    counts[isolate_id] += 1
                    fragments.append({
                        "isolate_id": isolate_id,
                        "fragment_id": f"{isolate_id}_pcr_{counts[isolate_id]}",
                        "contig_id": contig_id,
                        "contig_orig": left.get("contig_orig"),
                        "start": start - 1,
                        "end": end,
                        "length": length,
                        "orientation": orientation,
                        "forward_primer": forward_primer["query_id"],
                        "reverse_primer": reverse_primer["query_id"],
                        "forward_primer_differences": forward_primer["total_mismatches"],
                        "reverse_primer_differences": reverse_primer["total_mismatches"],
                        "gene_relation": relation,
                        "gene_query": gene_query,
                        "genomic_context": context,
                        "n_count": sequence.upper().count("N"),
                        "sequence": sequence,
                    })
    finally:
        genome.close()

    logger.info("Found %d in silico PCR fragments in %d isolates", len(fragments), len(counts))
    return fragments


def write_pcr_fragments(fragments, output_dir, config):
    if not fragments:
        return None

    records = []
    for f in fragments:
        description = " ".join([
            f"isolate={f['isolate_id']}",
            f"contig={f['contig_id']}",
            f"coords={f['start']}-{f['end']}",
            "coordinate_system=0-based_end-exclusive",
            f"length={f['length']}",
            f"orientation={f['orientation']}",
            f"primers={f['forward_primer']}+{f['reverse_primer']}",
            f"differences={f['forward_primer_differences']}+{f['reverse_primer_differences']}",
            f"gene={f['gene_relation'].replace(' ', '_')}",
        ])
        records.append(SeqRecord(Seq(f["sequence"]), id = f["fragment_id"], description = description))

    output_dir = Path(output_dir)
    output_dir.mkdir(parents = True, exist_ok = True)
    output_path = output_dir / f"{config.gene_name}_pcr_fragments.fasta"
    SeqIO.write(records, output_path, "fasta")
    logger.info("%d PCR fragments written to %s", len(records), output_path)

    return output_path
