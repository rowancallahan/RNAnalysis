"""
GSEA pathway analysis runner using GSEApy.
Supports both prerank GSEA and over-representation analysis via Enrichr.
"""

import numpy as np
import pandas as pd
from typing import List, Optional, Dict, Union


# Available gene set libraries
GENE_SET_LIBRARIES = [
    'KEGG_2021_Human',
    'GO_Biological_Process_2021',
    'GO_Molecular_Function_2021',
    'GO_Cellular_Component_2021',
    'Reactome_2022',
    'MSigDB_Hallmark_2020',
    'WikiPathway_2021_Human',
    'Panther_2016',
    'BioCarta_2016',
    'CORUM',
    'Transcription_Factor_PPIs',
    'miRTarBase_2017',
    'ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X'
]


class GSEARunner:
    """
    Wrapper for GSEApy pathway enrichment analysis.
    """

    def __init__(self):
        self.prerank_results = None
        self.enrichr_results = None

    @staticmethod
    def create_ranking_metric(deseq_results: pd.DataFrame) -> pd.Series:
        """
        Create ranking metric from DESeq2 results.
        Uses -log10(pvalue) * sign(log2FoldChange).

        Parameters
        ----------
        deseq_results : pd.DataFrame
            DESeq2 results with pvalue and log2FoldChange columns

        Returns
        -------
        pd.Series
            Ranking metric indexed by gene name
        """
        df = deseq_results.dropna(subset=['pvalue', 'log2FoldChange']).copy()

        # Handle zero p-values
        pvals = df['pvalue'].replace(0, 1e-300)

        # Create ranking: -log10(p) * sign(log2FC)
        ranking = -np.log10(pvals) * np.sign(df['log2FoldChange'])
        ranking = ranking.sort_values(ascending=False)

        return ranking

    def run_prerank(
        self,
        ranking: pd.Series,
        gene_sets: str = 'KEGG_2021_Human',
        min_size: int = 15,
        max_size: int = 500,
        permutation_num: int = 1000,
        seed: int = 42,
        outdir: str = None
    ) -> pd.DataFrame:
        """
        Run prerank GSEA on a ranked gene list.

        Parameters
        ----------
        ranking : pd.Series
            Gene names as index, ranking metric as values
        gene_sets : str
            Gene set library name
        min_size, max_size : int
            Gene set size filters
        permutation_num : int
            Number of permutations
        seed : int
            Random seed
        outdir : str
            Output directory (None to skip file output)

        Returns
        -------
        pd.DataFrame
            GSEA results with NES, p-value, FDR
        """
        import gseapy as gp

        # Create ranking DataFrame
        rnk = ranking.reset_index()
        rnk.columns = ['gene', 'rank']

        # Convert gene names to uppercase for compatibility
        rnk['gene'] = rnk['gene'].str.upper()

        try:
            self.prerank_results = gp.prerank(
                rnk=rnk,
                gene_sets=gene_sets,
                min_size=min_size,
                max_size=max_size,
                permutation_num=permutation_num,
                seed=seed,
                outdir=outdir,
                verbose=False
            )
            return self.prerank_results.res2d
        except Exception as e:
            # Return empty DataFrame if GSEA fails
            return pd.DataFrame(columns=[
                'Term', 'ES', 'NES', 'NOM p-val', 'FDR q-val', 'FWER p-val',
                'Tag %', 'Gene %', 'Lead_genes'
            ])

    def run_enrichr(
        self,
        gene_list: List[str],
        gene_sets: Union[str, List[str]] = 'KEGG_2021_Human',
        organism: str = 'human'
    ) -> pd.DataFrame:
        """
        Run Enrichr over-representation analysis.

        Parameters
        ----------
        gene_list : list
            List of gene symbols
        gene_sets : str or list
            Gene set library(s) to query
        organism : str
            Organism for gene sets

        Returns
        -------
        pd.DataFrame
            Enrichr results with columns:
            - Term, Overlap, P-value, Adjusted P-value, Odds Ratio, Combined Score, Genes
        """
        import gseapy as gp

        # Convert to list if single string
        if isinstance(gene_sets, str):
            gene_sets = [gene_sets]

        # Convert gene names to uppercase
        gene_list_upper = [g.upper() for g in gene_list if g]

        if not gene_list_upper:
            return pd.DataFrame()

        try:
            enr = gp.enrichr(
                gene_list=gene_list_upper,
                gene_sets=gene_sets,
                organism=organism,
                outdir=None,
                no_plot=True
            )
            self.enrichr_results = enr.results
            return self.enrichr_results
        except Exception as e:
            # Return empty DataFrame if Enrichr fails
            return pd.DataFrame(columns=[
                'Gene_set', 'Term', 'Overlap', 'P-value', 'Adjusted P-value',
                'Odds Ratio', 'Combined Score', 'Genes'
            ])

    def get_leading_edge_genes(self, term: str) -> List[str]:
        """
        Get leading edge genes for a specific pathway term.

        Parameters
        ----------
        term : str
            Pathway term name

        Returns
        -------
        list
            Leading edge genes for the term
        """
        if self.prerank_results is None:
            return []

        res = self.prerank_results.res2d
        if term in res['Term'].values:
            row = res[res['Term'] == term].iloc[0]
            genes = row.get('Lead_genes', '')
            if isinstance(genes, str):
                return genes.split(';') if genes else []
        return []

    @staticmethod
    def get_available_libraries() -> List[str]:
        """Get list of available gene set libraries."""
        return GENE_SET_LIBRARIES.copy()


def run_pathway_analysis(
    deseq_results: pd.DataFrame,
    method: str = 'prerank',
    gene_sets: str = 'KEGG_2021_Human',
    lfc_threshold: float = 1.0,
    padj_threshold: float = 0.05
) -> pd.DataFrame:
    """
    Convenience function to run pathway analysis.

    Parameters
    ----------
    deseq_results : pd.DataFrame
        DESeq2 results
    method : str
        'prerank' for GSEA prerank, 'enrichr' for ORA
    gene_sets : str
        Gene set library
    lfc_threshold : float
        Log2 fold change threshold (for enrichr)
    padj_threshold : float
        Adjusted p-value threshold (for enrichr)

    Returns
    -------
    pd.DataFrame
        Pathway analysis results
    """
    runner = GSEARunner()

    if method == 'prerank':
        ranking = runner.create_ranking_metric(deseq_results)
        return runner.run_prerank(ranking=ranking, gene_sets=gene_sets)
    else:
        # Get significant upregulated genes for enrichr
        sig_genes = deseq_results[
            (deseq_results['padj'] < padj_threshold) &
            (deseq_results['log2FoldChange'] > lfc_threshold)
        ].index.tolist()

        if not sig_genes:
            return pd.DataFrame()

        return runner.run_enrichr(gene_list=sig_genes, gene_sets=gene_sets)
