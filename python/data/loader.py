"""
Data loading utilities for counts and metadata files.
Supports CSV, TSV, and Excel formats.
"""

import pandas as pd
from pathlib import Path


def load_counts(file_path: str, index_col: int = 0) -> pd.DataFrame:
    """
    Load count matrix from file.

    Parameters
    ----------
    file_path : str
        Path to counts file (CSV, TSV, or Excel)
    index_col : int
        Column to use as index (gene names)

    Returns
    -------
    pd.DataFrame
        Count matrix with genes as rows, samples as columns
    """
    path = Path(file_path)
    suffix = path.suffix.lower()

    if suffix in ['.xlsx', '.xls']:
        df = pd.read_excel(file_path, index_col=index_col)
    elif suffix == '.tsv':
        df = pd.read_csv(file_path, sep='\t', index_col=index_col)
    else:  # Default to CSV
        df = pd.read_csv(file_path, index_col=index_col)

    # Ensure all values are numeric
    df = df.apply(pd.to_numeric, errors='coerce')

    # Remove rows with all NaN
    df = df.dropna(how='all')

    # Fill remaining NaN with 0
    df = df.fillna(0)

    # Convert to integers (counts should be integers)
    df = df.astype(int)

    df.index.name = 'gene'

    return df


def load_metadata(file_path: str, index_col: int = 0) -> pd.DataFrame:
    """
    Load sample metadata from file.

    Parameters
    ----------
    file_path : str
        Path to metadata file (CSV, TSV, or Excel)
    index_col : int
        Column to use as index (sample names)

    Returns
    -------
    pd.DataFrame
        Metadata with samples as rows, annotations as columns
    """
    path = Path(file_path)
    suffix = path.suffix.lower()

    if suffix in ['.xlsx', '.xls']:
        df = pd.read_excel(file_path, index_col=index_col)
    elif suffix == '.tsv':
        df = pd.read_csv(file_path, sep='\t', index_col=index_col)
    else:  # Default to CSV
        df = pd.read_csv(file_path, index_col=index_col)

    df.index.name = 'sample'

    return df
