"""
Generate tiny example FASTA and FASTQ files for testing the kallisto pipeline.
Two short transcripts (GeneA, GeneB) with 4 samples (2 control, 2 treatment).
"""

import os
import random
import csv
from pathlib import Path


# Two short "transcript" sequences (~300 bp each)
# These use pseudo-random but diverse nucleotide sequences so kallisto
# can build a useful k-mer index (default k=31 requires diverse 31-mers).
GENE_A_SEQ = (
    "ATGGCTCAGTTCAAGGATCTGAACGTCATCGACACCTTCAGTGCGAATATCGTCAACGAGG"
    "TCGTCAAAGCGATCATCGGCTTTCCTGTCAGCAAAGCCTATGGCACTGTCGGTTTTCATACG"
    "GCATCTTCACAGAGTTCAGGATCTTCAGCGATCAACGGGCCGAATATCTGGAAGCGCTCAGC"
    "AATGATCCAGATGCCTATTCGGAAATCAAGCAAGGCTATCAGTCGCTTAAGGCTTCAAATAC"
    "TGGCCAATGTGGCATCGTTCAGCACTTCAATGGTCGACTTCACATCCAAGGATCCTTCAGA"
)  # 299 bp

GENE_B_SEQ = (
    "ATGAACGGTCTGCGCAAAGACGTCAAGATCATCGACTCCTACTTCGTCTGCAACTACATCG"
    "ACGGCGTCAAATCCGATCTCGGCAACCCTGATAGCAAAGCGTTTGGCACGGTCGGCTTCCAC"
    "ACCGCATCGTCACAGAACTCGGGATCTTCAGCAATAAACGGACCCAACATCTGGAAACGATC"
    "AGCCATGACCCTGATGCATACTCAGAGATCAAGCAGGGTTACCAGAGCCTTAAAGCATCGAA"
    "CACCGGTCAATGCGGCATCGTCCAACACTTCAACGGTCGGCTTCACATCCAAAGACCCTGA"
)  # 299 bp

READ_LENGTH = 75


def _random_quality(length: int) -> str:
    """Generate a fake quality string (all high quality)."""
    return 'I' * length  # Phred 40


def _make_reads_from_seq(seq: str, n_reads: int, read_len: int = READ_LENGTH) -> list:
    """Extract random subsequences from a transcript to simulate reads."""
    reads = []
    max_start = len(seq) - read_len
    if max_start < 0:
        max_start = 0
        read_len = len(seq)
    for _ in range(n_reads):
        start = random.randint(0, max_start)
        read_seq = seq[start:start + read_len]
        reads.append(read_seq)
    return reads


def create_example_fasta(output_path: str) -> str:
    """Write a minimal 2-transcript FASTA file."""
    os.makedirs(os.path.dirname(output_path) or '.', exist_ok=True)
    with open(output_path, 'w') as f:
        f.write('>GeneA_transcript\n')
        # Wrap at 80 chars
        for i in range(0, len(GENE_A_SEQ), 80):
            f.write(GENE_A_SEQ[i:i+80] + '\n')
        f.write('>GeneB_transcript\n')
        for i in range(0, len(GENE_B_SEQ), 80):
            f.write(GENE_B_SEQ[i:i+80] + '\n')
    return output_path


def create_example_fastqs(
    output_dir: str,
    n_samples: int = 4,
    seed: int = 42,
) -> dict:
    """Generate example single-end FASTQ files and a metadata CSV.

    Returns dict with keys:
        - samples: list of {name, condition, fastq_path}
        - metadata_path: path to metadata.csv
        - fasta_path: path to generated reference FASTA
    """
    random.seed(seed)
    os.makedirs(output_dir, exist_ok=True)

    # Also create the FASTA in the same directory
    fasta_path = os.path.join(output_dir, 'example_transcripts.fa')
    create_example_fasta(fasta_path)

    # Sample definitions
    # Control: roughly equal reads from both genes
    # Treatment: more GeneA, fewer GeneB (simulating DE)
    sample_defs = [
        ('Control_1', 'control', 50, 50),   # (name, condition, geneA_reads, geneB_reads)
        ('Control_2', 'control', 55, 45),
        ('Treatment_1', 'treatment', 90, 10),
        ('Treatment_2', 'treatment', 85, 15),
    ]

    if n_samples > 4:
        # Add more samples if requested
        for i in range(4, n_samples):
            if i % 2 == 0:
                sample_defs.append((f'Control_{i//2+1}', 'control', 50 + random.randint(-5, 5), 50 + random.randint(-5, 5)))
            else:
                sample_defs.append((f'Treatment_{i//2+1}', 'treatment', 85 + random.randint(-5, 5), 15 + random.randint(-5, 5)))

    samples = []
    for name, condition, n_a, n_b in sample_defs[:n_samples]:
        reads_a = _make_reads_from_seq(GENE_A_SEQ, n_a)
        reads_b = _make_reads_from_seq(GENE_B_SEQ, n_b)
        all_reads = reads_a + reads_b
        random.shuffle(all_reads)

        fastq_path = os.path.join(output_dir, f'{name}.fastq')
        with open(fastq_path, 'w') as f:
            for i, seq in enumerate(all_reads):
                f.write(f'@{name}_read_{i}\n')
                f.write(seq + '\n')
                f.write('+\n')
                f.write(_random_quality(len(seq)) + '\n')

        samples.append({
            'name': name,
            'condition': condition,
            'fastq_path': fastq_path,
        })

    # Write metadata CSV
    metadata_path = os.path.join(output_dir, 'metadata.csv')
    with open(metadata_path, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['sample', 'condition'])
        for s in samples:
            writer.writerow([s['name'], s['condition']])

    return {
        'samples': samples,
        'metadata_path': metadata_path,
        'fasta_path': fasta_path,
    }
