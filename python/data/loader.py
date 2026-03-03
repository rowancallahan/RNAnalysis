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


def validate_data(counts: pd.DataFrame, metadata: pd.DataFrame) -> dict:
    """
    Validate that counts and metadata are compatible.

    Parameters
    ----------
    counts : pd.DataFrame
        Count matrix (genes x samples)
    metadata : pd.DataFrame
        Sample metadata

    Returns
    -------
    dict
        Validation results with 'valid' bool and 'messages' list
    """
    messages = []
    valid = True

    # Check if sample names match
    count_samples = set(counts.columns)
    metadata_samples = set(metadata.index)

    missing_in_metadata = count_samples - metadata_samples
    missing_in_counts = metadata_samples - count_samples

    if missing_in_metadata:
        messages.append(f"Samples in counts but not metadata: {list(missing_in_metadata)[:5]}")
        valid = False

    if missing_in_counts:
        messages.append(f"Samples in metadata but not counts: {list(missing_in_counts)[:5]}")
        valid = False

    # Check for negative values
    if (counts < 0).any().any():
        messages.append("Warning: Negative values found in counts")
        valid = False

    # Check for numeric columns in counts
    non_numeric = counts.select_dtypes(exclude=['number']).columns
    if len(non_numeric) > 0:
        messages.append(f"Non-numeric columns in counts: {list(non_numeric)}")
        valid = False

    return {
        'valid': valid,
        'messages': messages,
        'matching_samples': len(count_samples & metadata_samples),
        'total_genes': len(counts),
        'total_samples': len(count_samples)
    }
