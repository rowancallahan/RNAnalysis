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


def vst_transform(
    counts: pd.DataFrame,
    size_factors: Optional[pd.Series] = None
) -> pd.DataFrame:
    """
    Variance stabilizing transformation (simplified version).

    For the full VST, use DESeq2's implementation.
    This is a simplified approximation using asinh transformation.

    Parameters
    ----------
    counts : pd.DataFrame
        Raw count matrix
    size_factors : pd.Series, optional
        Size factors for normalization

    Returns
    -------
    pd.DataFrame
        Variance-stabilized counts
    """
    # Normalize first
    if size_factors is not None:
        normalized = counts.div(size_factors, axis='columns')
    else:
        normalized = normalize_counts(counts)

    # Apply asinh transformation (approximates VST for moderate counts)
    transformed = np.arcsinh(normalized)

    return transformed


def rlog_transform(
    counts: pd.DataFrame,
    size_factors: Optional[pd.Series] = None
) -> pd.DataFrame:
    """
    Regularized log transformation (simplified version).

    For the full rlog, use DESeq2's implementation.
    This is a simplified approximation.

    Parameters
    ----------
    counts : pd.DataFrame
        Raw count matrix
    size_factors : pd.Series, optional
        Size factors for normalization

    Returns
    -------
    pd.DataFrame
        Regularized log transformed counts
    """
    # Normalize first
    if size_factors is not None:
        normalized = counts.div(size_factors, axis='columns')
    else:
        normalized = normalize_counts(counts)

    # Apply log transformation with regularization
    # Use a dispersion-like shrinkage factor
    dispersion = 0.1
    transformed = np.log2((normalized + 1) / (1 + dispersion * normalized.mean(axis=1).values[:, None]))

    return transformed


def zscore_normalize(
    data: pd.DataFrame,
    axis: int = 0
) -> pd.DataFrame:
    """
    Z-score normalize data.

    Parameters
    ----------
    data : pd.DataFrame
        Data matrix
    axis : int
        0 for row-wise (per gene), 1 for column-wise (per sample)

    Returns
    -------
    pd.DataFrame
        Z-score normalized data
    """
    if axis == 0:
        # Row-wise (per gene)
        means = data.mean(axis=1)
        stds = data.std(axis=1)
        stds = stds.replace(0, 1)  # Avoid division by zero
        return data.sub(means, axis=0).div(stds, axis=0)
    else:
        # Column-wise (per sample)
        means = data.mean(axis=0)
        stds = data.std(axis=0)
        stds = stds.replace(0, 1)
        return data.sub(means, axis=1).div(stds, axis=1)
