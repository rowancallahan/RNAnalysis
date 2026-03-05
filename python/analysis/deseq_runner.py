"""
DESeq2 analysis runner using PyDESeq2.
Wraps the PyDESeq2 pipeline and provides convenient methods for the GUI.
"""

import sys
import os
import numpy as np
import pandas as pd
from typing import Optional
from contextlib import contextmanager


@contextmanager
def suppress_stdout():
    """Context manager to suppress stdout (used to prevent pydeseq2 from corrupting JSON output)."""
    # Save the real stdout
    old_stdout = sys.stdout
    # Redirect stdout to devnull
    sys.stdout = open(os.devnull, 'w')
    try:
        yield
    finally:
        sys.stdout.close()
        sys.stdout = old_stdout


class DESeqRunner:
    """
    Wrapper for PyDESeq2 differential expression analysis.
    """

    def __init__(self):
        self.dds = None
        self.stats = None
        self.results_df = None
        self._normalized_counts = None

    def run_analysis(
        self,
        counts: pd.DataFrame,
        metadata: pd.DataFrame,
        condition_column: str = 'condition',
        reference_level: str = 'control',
        refit_cooks: bool = True,
        n_cpus: int = 4
    ) -> pd.DataFrame:
        """
        Run complete DESeq2 analysis pipeline.

        Parameters
        ----------
        counts : pd.DataFrame
            Raw count matrix (genes x samples)
        metadata : pd.DataFrame
            Sample metadata with condition column
        condition_column : str
            Column name for the condition factor
        reference_level : str
            Reference level for the condition factor
        refit_cooks : bool
            Whether to refit outliers
        n_cpus : int
            Number of CPUs for parallel processing

        Returns
        -------
        pd.DataFrame
            Results with columns: baseMean, log2FoldChange, lfcSE, stat, pvalue, padj
        """
        from pydeseq2.dds import DeseqDataSet
        from pydeseq2.ds import DeseqStats
        from pydeseq2.default_inference import DefaultInference

        # Ensure samples match between counts and metadata
        common_samples = list(set(counts.columns) & set(metadata.index))
        counts = counts[common_samples]
        metadata = metadata.loc[common_samples]

        # PyDESeq2 expects counts as samples x genes
        counts_t = counts.T

        # Create inference object
        inference = DefaultInference(n_cpus=n_cpus)

        # Create design formula
        design_formula = f"~ {condition_column}"

        # Create DeseqDataSet - suppress stdout to prevent JSON corruption
        with suppress_stdout():
            self.dds = DeseqDataSet(
                counts=counts_t,
                metadata=metadata,
                design_factors=condition_column,
                refit_cooks=refit_cooks,
                inference=inference,
                ref_level=[condition_column, reference_level]
            )

            # Run DESeq2 pipeline
            self.dds.deseq2()

            # Get the non-reference level for contrast
            levels = metadata[condition_column].unique().tolist()
            non_ref_level = [l for l in levels if l != reference_level][0]

            # Statistical testing - need to specify contrast
            self.stats = DeseqStats(
                self.dds,
                contrast=[condition_column, non_ref_level, reference_level],
                inference=inference
            )
            self.stats.summary()

        # Get results DataFrame
        self.results_df = self.stats.results_df.copy()

        # Store normalized counts
        self._normalized_counts = pd.DataFrame(
            self.dds.layers['normed_counts'],
            index=self.dds.obs_names,
            columns=self.dds.var_names
        )

        return self.results_df

    def get_normalized_counts(self) -> pd.DataFrame:
        """
        Get DESeq2 normalized counts.

        Returns
        -------
        pd.DataFrame
            Normalized count matrix (samples x genes)
        """
        if self._normalized_counts is None:
            raise ValueError("Run analysis first")
        return self._normalized_counts

