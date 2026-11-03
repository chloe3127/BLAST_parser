"""
Define what BLAST queries you want to run and the different gene/primer BLAST parameters
"""

gene_name = "" 
identifier = "" #if dont want output folders with long gene names can have a unique identifier instead

gene_queries = {
    "query_name": {
        "query_sequence": "",
    },
} #can extend this with as many gene sequences you want

flank_queries = {
    "upstream_flank": {
        "query_sequence": "",
        "position": "upstream",
        "distance": 0
    },
    "downstream_flank": {
        "query_sequence": "",
        "position": "downstream",
        "distance": 0
    }
}

flank_buffer = 500 #buffer to check flanking distances

gene_start_queries = {
    "query_name": {
        "query_sequence": "",
    }
}

gene_end_queries = {
    "query_name": {
        "query_sequence": "",
    }
}

primer_queries = {
    "f_primer": {
        "query_sequence": "",
        "position": "forward",
    },
    "r_primer": {
        "query_sequence": "",
        "position": "reverse",
    },
}

for queries in (gene_queries, flank_queries, gene_start_queries, gene_end_queries, primer_queries):
    for query_data in queries.values():
        query_data["query_sequence"] = "".join(query_data["query_sequence"].split()).upper()


#BLAST search parameters
gene_blast_params = {
    "task": "blastn",
    "min_identity": 80,
    "high_identity": 90,
    "min_coverage": 70,
    "high_coverage": 90,
    "max_hits_per_isolate": 10,
    "max_target_seqs": 1000000,
    "max_hsps": 10,
    "word_size": 11,
    "reward": 2,
    "penalty": -3,
    "gapopen": 5,
    "gapextend": 2,
    "evalue": 1e-10
}

flank_blast_params = {
    "task": "blastn",
    "min_identity": 80,
    "high_identity": 90,
    "min_coverage": 70,
    "high_coverage": 90,
    "max_hits_per_isolate": 10,
    "max_target_seqs": 1000000,
    "max_hsps": 10,
    "word_size": 11,
    "reward": 2,
    "penalty": -3,
    "gapopen": 5,
    "gapextend": 2,
    "evalue": 1e-10
}

gene_start_blast_params = {
    "task": "blastn",
    "min_identity": 80,
    "high_identity": 95,
    "min_coverage": 70,
    "word_size": 7,
    "evalue": 1000,
    "max_hits_per_isolate": 5,
    "max_target_seqs": 1000000,
    "max_hsps": 10,
    "reward": 2,
    "penalty": -3,
    "gapopen": 5,
    "gapextend": 2,
}

gene_end_blast_params = {
    "task": "blastn",
    "min_identity": 80,
    "high_identity": 90,
    "min_coverage": 70,
    "high_coverage": 90,
    "max_hits_per_isolate": 10,
    "max_target_seqs": 1000000,
    "max_hsps": 10,
    "word_size": 11,
    "reward": 2,
    "penalty": -3,
    "gapopen": 5,
    "gapextend": 2,
    "evalue": 1e-10
}

primer_blast_params = {
    "task": "blastn-short",
    "min_identity": 80,
    "min_coverage": 70,
    "max_mismatches": 4,
    "word_size": 7,
    "evalue": 1000,
    "min_alignment_pct": 90,
    "max_hits_per_isolate": 5,
    "max_target_seqs": 1000000,
    "max_pcr_fragment_length": 1000, #in silico PCR and gene dupes mean you get some 'fragments' between the genes which wouldn't happen in vivo so filter them out,
    "max_hsps": 10,
    "reward": 2,
    "penalty": -3,
    "gapopen": 5,
    "gapextend": 2,
}

extraction_params = {
    "extract_flank": True, #save one fasta without flank and save a second with specified flanking length
    "flanking_bp": 1000,
    "boundary_critical_distance": 25,
    "boundary_warning_distance": 50,
    "extract_igr": False
}