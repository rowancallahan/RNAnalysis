"""
Generate synthetic RNA-seq count data for testing.
Creates realistic-looking data with known differentially expressed genes.
"""

import numpy as np
import pandas as pd
from typing import Tuple, List


# Realistic gene names for DE genes
REALISTIC_GENE_NAMES = [
    # Oncogenes and tumor suppressors
    'TP53', 'BRCA1', 'BRCA2', 'EGFR', 'MYC', 'KRAS', 'PIK3CA', 'PTEN', 'RB1',
    'APC', 'CDKN2A', 'VEGFA', 'ERBB2', 'BRAF', 'NRAS', 'IDH1', 'IDH2',
    # Hematological markers
    'FLT3', 'NPM1', 'DNMT3A', 'TET2', 'JAK2', 'CALR', 'MPL', 'BCL2',
    'BCL6', 'MYD88', 'SF3B1', 'ASXL1', 'RUNX1', 'CEBPA',
    # Stem cell markers
    'TP63', 'SOX2', 'KLF4', 'NANOG', 'POU5F1', 'LIN28A', 'PROM1', 'ALDH1A1',
    # Cytokines and immune
    'IL6', 'TNF', 'IFNG', 'IL1B', 'IL10', 'TGFB1', 'IL2', 'IL4', 'IL17A',
    'CXCL8', 'CCL2', 'CXCL12', 'CD274', 'PDCD1', 'CTLA4', 'LAG3',
    # Signaling
    'STAT3', 'NFKB1', 'MAPK1', 'MAPK3', 'AKT1', 'MTOR', 'HIF1A',
    'NOTCH1', 'WNT1', 'SHH', 'PTCH1', 'SMO', 'GLI1', 'CTNNB1',
    # Cell cycle
    'CCND1', 'CCNE1', 'CDK4', 'CDK6', 'CDKN1A', 'CDKN1B', 'E2F1',
    # Apoptosis
    'BAX', 'BAK1', 'BID', 'BCL2L1', 'MCL1', 'XIAP', 'CASP3', 'CASP9',
    # DNA repair
    'ATM', 'ATR', 'CHEK1', 'CHEK2', 'RAD51', 'XRCC1', 'PARP1',
    # Metabolism
    'LDHA', 'PKM', 'HK2', 'GLUT1', 'ACLY', 'FASN', 'SCD1',
    # EMT markers
    'CDH1', 'CDH2', 'VIM', 'SNAI1', 'SNAI2', 'ZEB1', 'TWIST1',
    # Housekeeping (for reference)
    'GAPDH', 'ACTB', 'TUBB', 'RPL13A', 'B2M', 'HPRT1', 'TBP'
]


def generate_example_data(
    n_samples: int = 20,
    n_genes: int = 5000,
    n_de_genes: int = 500,
    effect_size: float = 2.0,
    seed: int = 42
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Generate synthetic RNA-seq count data with known differential expression.

    Parameters
    ----------
    n_samples : int
        Total number of samples (split evenly between conditions)
    n_genes : int
        Total number of genes
    n_de_genes : int
        Number of differentially expressed genes
    effect_size : float
        Log2 fold change for DE genes
    seed : int
        Random seed for reproducibility

    Returns
    -------
    counts : pd.DataFrame
        Count matrix (genes x samples)
    metadata : pd.DataFrame
        Sample metadata with condition assignments
    """
    np.random.seed(seed)

    n_per_group = n_samples // 2

    # Generate gene names
    gene_names = []

    # Use realistic names for DE genes
    n_realistic = min(len(REALISTIC_GENE_NAMES), n_de_genes)
    gene_names.extend(REALISTIC_GENE_NAMES[:n_realistic])

    # Fill remaining DE genes with synthetic names
    for i in range(n_realistic, n_de_genes):
        gene_names.append(f"DE_GENE_{i:04d}")

    # Add non-DE genes
    for i in range(n_de_genes, n_genes):
        gene_names.append(f"GENE_{i:05d}")

    # Generate sample names
    control_samples = [f"Control_{i+1:02d}" for i in range(n_per_group)]
    treatment_samples = [f"Treatment_{i+1:02d}" for i in range(n_per_group)]
    sample_names = control_samples + treatment_samples

    # Generate base expression levels (log-normal distribution)
    # Most genes lowly expressed, few highly expressed
    base_means = np.exp(np.random.normal(5, 2, n_genes))
    base_means = np.clip(base_means, 10, 100000)

    # Generate dispersion parameters (inversely related to mean)
    dispersions = 0.1 + 0.5 / np.sqrt(base_means)

    # Generate counts matrix
    counts = np.zeros((n_genes, n_samples), dtype=int)

    for i in range(n_genes):
        mean_expr = base_means[i]
        disp = dispersions[i]

        # Control samples
        for j in range(n_per_group):
            # Add sample-specific variation
            sample_mean = mean_expr * np.exp(np.random.normal(0, 0.2))
            # Negative binomial parameters
            r = 1 / disp
            p = r / (r + sample_mean)
            counts[i, j] = np.random.negative_binomial(max(1, int(r)), min(0.99, max(0.01, p)))

        # Treatment samples
        for j in range(n_per_group):
            # Apply fold change for DE genes
            if i < n_de_genes:
                # Alternate up/down regulation
                if i % 2 == 0:
                    fold_change = 2 ** effect_size
                else:
                    fold_change = 2 ** (-effect_size)
                sample_mean = mean_expr * fold_change
            else:
                sample_mean = mean_expr

            # Add sample-specific variation
            sample_mean *= np.exp(np.random.normal(0, 0.2))
            # Negative binomial parameters
            r = 1 / disp
            p = r / (r + sample_mean)
            counts[i, n_per_group + j] = np.random.negative_binomial(max(1, int(r)), min(0.99, max(0.01, p)))

    # Create counts DataFrame (genes as rows, samples as columns)
    counts_df = pd.DataFrame(
        counts,
        index=gene_names,
        columns=sample_names
    )
    counts_df.index.name = 'gene'

    # Create metadata DataFrame
    conditions = ['control'] * n_per_group + ['treatment'] * n_per_group
    batches = np.random.choice(['batch1', 'batch2'], n_samples)
    sexes = np.random.choice(['M', 'F'], n_samples)
    ages = np.random.randint(25, 75, n_samples)
    rna_quality = np.round(np.random.uniform(7.5, 10.0, n_samples), 1)

    metadata_df = pd.DataFrame({
        'condition': conditions,
        'batch': batches,
        'sex': sexes,
        'age': ages,
        'rna_quality': rna_quality
    }, index=sample_names)
    metadata_df.index.name = 'sample'

    return counts_df, metadata_df


def get_expected_de_genes(n_de_genes: int = 500) -> Tuple[List[str], List[str]]:
    """
    Return lists of expected up and down regulated genes.

    Parameters
    ----------
    n_de_genes : int
        Number of DE genes generated

    Returns
    -------
    up_genes : list
        Expected upregulated genes (even indices)
    down_genes : list
        Expected downregulated genes (odd indices)
    """
    n_realistic = min(len(REALISTIC_GENE_NAMES), n_de_genes)

    up_genes = []
    down_genes = []

    for i in range(n_de_genes):
        if i < n_realistic:
            gene = REALISTIC_GENE_NAMES[i]
        else:
            gene = f"DE_GENE_{i:04d}"

        if i % 2 == 0:
            up_genes.append(gene)
        else:
            down_genes.append(gene)

    return up_genes, down_genes


if __name__ == '__main__':
    # Test the generator
    counts, metadata = generate_example_data()
    print(f"Counts shape: {counts.shape}")
    print(f"Metadata shape: {metadata.shape}")
    print(f"\nCounts preview:\n{counts.iloc[:5, :5]}")
    print(f"\nMetadata preview:\n{metadata.head()}")

    up, down = get_expected_de_genes()
    print(f"\nExpected upregulated: {up[:10]}...")
    print(f"Expected downregulated: {down[:10]}...")
