"""
Count normalization utilities.
"""

import numpy as np
import pandas as pd
from typing import Optional


def normalize_counts(
    counts: pd.DataFrame,
    size_factors: Optional[pd.Series] = None
) -> pd.DataFrame:
    """
    Normalize counts using size factors (DESeq2-style normalization).

    If size factors are not provided, calculates median-of-ratios normalization.

    Parameters
    ----------
    counts : pd.DataFrame
        Raw count matrix (genes x samples or samples x genes)
    size_factors : pd.Series, optional
        Pre-computed size factors indexed by sample name

    Returns
    -------
    pd.DataFrame
        Normalized counts
    """
    # Ensure counts are numeric
    counts = counts.apply(pd.to_numeric, errors='coerce').fillna(0)

    if size_factors is None:
        size_factors = calculate_size_factors(counts)

    # Normalize
    normalized = counts.div(size_factors, axis='columns')

    return normalized


def calculate_size_factors(counts: pd.DataFrame) -> pd.Series:
    """
    Calculate size factors using median-of-ratios method (DESeq2 approach).

    Parameters
    ----------
    counts : pd.DataFrame
        Raw count matrix (genes x samples)

    Returns
    -------
    pd.Series
        Size factors indexed by sample name
    """
    # Add pseudocount to avoid division by zero
    counts_pseudo = counts + 1

    # Calculate geometric mean per gene
    log_counts = np.log(counts_pseudo)
    geo_means = np.exp(log_counts.mean(axis=1))

    # Calculate ratios to geometric mean
    ratios = counts_pseudo.div(geo_means, axis=0)

    # Size factors are median of ratios per sample
    size_factors = ratios.median(axis=0)

    # Normalize size factors to have geometric mean of 1
    size_factors = size_factors / np.exp(np.log(size_factors).mean())

    return size_factors


def log2_transform(
    counts: pd.DataFrame,
    pseudocount: float = 1.0
) -> pd.DataFrame:
    """
    Log2 transform counts.

    Parameters
    ----------
    counts : pd.DataFrame
        Count matrix (can be raw or normalized)
    pseudocount : float
        Pseudocount to add before log transformation

    Returns
    -------
    pd.DataFrame
        Log2 transformed counts
    """
    return np.log2(counts + pseudocount)
