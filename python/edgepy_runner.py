"""
edgePython runner — wraps the edgepython package for differential expression analysis.
Provides the same interface shape as the DESeqRunner so the frontend can switch tools seamlessly.
"""

import numpy as np
import pandas as pd


class EdgePythonRunner:
    """Run edgePython differential expression analysis."""

    def __init__(self):
        self.dgelist = None
        self.fit = None
        self.results_df = None
        self.counts = None
        self.metadata = None
        self.norm_factors = None

    def run_analysis(
        self,
        counts: pd.DataFrame,
        metadata: pd.DataFrame,
        condition_column: str = 'condition',
        reference_level: str = None,
        normalization: str = 'TMM',
        test_method: str = 'qlf',
        dispersion: str = 'trended',
        robust: bool = True,
        min_count: int = 10,
    ) -> pd.DataFrame:
        """
        Run edgePython differential expression.

        Parameters
        ----------
        counts : DataFrame
            Gene x sample count matrix (genes as rows, samples as columns).
        metadata : DataFrame
            Sample metadata with index matching count columns.
        condition_column : str
            Column in metadata with group labels.
        reference_level : str
            Reference group for comparison.
        normalization : str
            Normalization method: TMM, TMMwsp, RLE, upperquartile, none.
        test_method : str
            Test method: 'qlf' (quasi-likelihood F-test) or 'exact' (exact test).
        dispersion : str
            Dispersion type: trended, common, tagwise.
        robust : bool
            Use robust fitting for GLM.
        min_count : int
            Minimum total count to keep a gene.

        Returns
        -------
        DataFrame with columns: baseMean, log2FoldChange, lfcSE, pvalue, padj
        """
        import edgepython as ep

        self.counts = counts.copy()
        self.metadata = metadata.copy()

        # Ensure sample order matches
        shared_samples = [s for s in counts.columns if s in metadata.index]
        counts_mat = counts[shared_samples]
        meta = metadata.loc[shared_samples]

        # Filter low-count genes
        if min_count > 0:
            keep = counts_mat.sum(axis=1) >= min_count
            counts_mat = counts_mat[keep]

        # Get group labels
        groups = meta[condition_column].values
        unique_groups = list(pd.unique(groups))

        # Set reference level
        if reference_level and reference_level in unique_groups:
            unique_groups.remove(reference_level)
            unique_groups.insert(0, reference_level)

        # Encode groups numerically for design matrix
        group_numeric = np.array([unique_groups.index(g) for g in groups], dtype=float)

        # Convert to numpy
        counts_np = counts_mat.values.astype(float)
        gene_names = list(counts_mat.index)

        # Create DGEList
        y = ep.make_dgelist(counts=counts_np, group=group_numeric)

        # Normalize
        if normalization != 'none':
            y = ep.calc_norm_factors(y, method=normalization)

        # Store norm factors
        self.dgelist = y

        # Estimate dispersion
        y = ep.estimate_disp(y)

        n_samples = len(shared_samples)

        if test_method == 'exact':
            # Exact test (two-group comparison)
            results = ep.exact_test(y)
            top = ep.top_tags(results, n=len(gene_names))
        else:
            # Quasi-likelihood F-test with design matrix
            # Design: intercept + condition indicator
            design = np.column_stack([
                np.ones(n_samples),
                group_numeric,
            ])

            fit = ep.glm_ql_fit(y, design, robust=robust)
            self.fit = fit

            results = ep.glm_ql_ftest(fit, coef=1)
            top = ep.top_tags(results, n=len(gene_names))

        # Extract results table
        res_table = top['table']

        # Map back to gene names using original index order
        if len(res_table) <= len(gene_names):
            res_table.index = [gene_names[i] if i < len(gene_names) else f'gene_{i}'
                               for i in res_table.index]

        # Standardize column names to match DESeq2 output format
        result_df = pd.DataFrame(index=res_table.index)

        # baseMean — average count across samples (like DESeq2)
        if 'logCPM' in res_table.columns:
            # Convert logCPM back to approximate mean: 2^logCPM * lib_size / 1e6
            result_df['baseMean'] = counts_mat.loc[result_df.index].mean(axis=1) if all(
                g in counts_mat.index for g in result_df.index
            ) else np.power(2, res_table['logCPM'].values)
        else:
            result_df['baseMean'] = counts_mat.loc[result_df.index].mean(axis=1)

        # log2FoldChange
        if 'logFC' in res_table.columns:
            result_df['log2FoldChange'] = res_table['logFC'].values
        else:
            result_df['log2FoldChange'] = 0.0

        # lfcSE (not directly available from edgeR, estimate from logFC confidence)
        result_df['lfcSE'] = np.abs(result_df['log2FoldChange']) / 2.0  # rough estimate

        # p-value
        if 'PValue' in res_table.columns:
            result_df['pvalue'] = res_table['PValue'].values
        else:
            result_df['pvalue'] = 1.0

        # adjusted p-value (FDR)
        if 'FDR' in res_table.columns:
            result_df['padj'] = res_table['FDR'].values
        else:
            result_df['padj'] = result_df['pvalue']

        result_df.index.name = 'gene'
        self.results_df = result_df

        return result_df

    def get_normalized_counts(self) -> pd.DataFrame:
        """
        Return normalized counts (CPM-style using norm factors).
        """
        if self.counts is None:
            raise ValueError('Run analysis first')

        # Simple CPM normalization using library sizes
        lib_sizes = self.counts.sum(axis=0)
        cpm = self.counts.div(lib_sizes, axis=1) * 1e6
        return cpm
