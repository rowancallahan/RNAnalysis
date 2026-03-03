"""
Plot generator for creating publication-quality matplotlib figures.
Generates SVG plots that can be displayed in the Electron app.
"""

import os
import io
import base64
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')  # Non-interactive backend
import matplotlib.pyplot as plt
import seaborn as sns
from typing import Dict, List, Optional, Tuple


class PlotGenerator:
    """Generate matplotlib plots and return as SVG or base64 PNG."""

    def __init__(self, output_dir: Optional[str] = None):
        self.output_dir = output_dir or '/tmp/deseq_plots'
        os.makedirs(self.output_dir, exist_ok=True)

        # Set default style - use a clean white background
        plt.style.use('seaborn-v0_8-white')
        plt.rcParams.update({
            'figure.dpi': 150,
            'savefig.dpi': 300,
            'figure.facecolor': 'white',
            'axes.facecolor': 'white',
            'savefig.facecolor': 'white',
            # Use DejaVu Sans which has full Unicode support including subscripts
            'font.family': 'sans-serif',
            'font.sans-serif': ['DejaVu Sans', 'Helvetica', 'Arial'],
            'font.size': 12,
            'axes.titlesize': 16,
            'axes.labelsize': 14,
            'xtick.labelsize': 12,
            'ytick.labelsize': 12,
            'legend.fontsize': 11,
            'axes.linewidth': 1.5,
            'xtick.major.width': 1.5,
            'ytick.major.width': 1.5,
            'xtick.major.size': 6,
            'ytick.major.size': 6,
        })

    def _fig_to_svg(self, fig: plt.Figure) -> str:
        """Convert matplotlib figure to SVG string."""
        buf = io.BytesIO()
        fig.savefig(buf, format='svg', bbox_inches='tight', facecolor='white', edgecolor='none')
        buf.seek(0)
        svg_string = buf.getvalue().decode('utf-8')
        plt.close(fig)
        return svg_string

    def _fig_to_base64_png(self, fig: plt.Figure) -> str:
        """Convert matplotlib figure to base64 PNG."""
        buf = io.BytesIO()
        fig.savefig(buf, format='png', bbox_inches='tight', dpi=150, facecolor='white', edgecolor='none')
        buf.seek(0)
        img_base64 = base64.b64encode(buf.getvalue()).decode('utf-8')
        plt.close(fig)
        return f"data:image/png;base64,{img_base64}"

    def create_volcano_plot(
        self,
        results: pd.DataFrame,
        lfc_threshold: float = 1.0,
        padj_threshold: float = 0.05,
        highlight_genes: Optional[List[str]] = None,
        output_format: str = 'svg'
    ) -> Dict:
        """
        Create publication-quality volcano plot.

        Parameters
        ----------
        results : pd.DataFrame
            DESeq2 results with log2FoldChange and padj columns
        lfc_threshold : float
            Log2 fold change threshold for significance
        padj_threshold : float
            Adjusted p-value threshold
        highlight_genes : list, optional
            Genes to label on the plot
        output_format : str
            'svg' or 'png'

        Returns
        -------
        dict with 'image' (SVG string or base64 PNG) and 'stats'
        """
        fig, ax = plt.subplots(figsize=(10, 8))

        # Clean white background
        ax.set_facecolor('white')
        fig.patch.set_facecolor('white')

        # Prepare data
        df = results.dropna(subset=['log2FoldChange', 'padj']).copy()
        df['neg_log10_padj'] = -np.log10(df['padj'].replace(0, 1e-300))

        # Categorize genes
        up_mask = (df['log2FoldChange'] >= lfc_threshold) & (df['padj'] < padj_threshold)
        down_mask = (df['log2FoldChange'] <= -lfc_threshold) & (df['padj'] < padj_threshold)
        ns_mask = ~(up_mask | down_mask)

        # Plot non-significant genes first (background) - with dark edge
        ax.scatter(
            df.loc[ns_mask, 'log2FoldChange'],
            df.loc[ns_mask, 'neg_log10_padj'],
            c='#CCCCCC', alpha=0.5, s=30, label='NS',
            edgecolors='#666666', linewidths=0.5, rasterized=True
        )

        # Plot significant genes with dark borders - using muted pastel colors
        ax.scatter(
            df.loc[up_mask, 'log2FoldChange'],
            df.loc[up_mask, 'neg_log10_padj'],
            c='#C97B7B', alpha=0.8, s=50, label=f'Up ({up_mask.sum()})',
            edgecolors='#8B4545', linewidths=1.0
        )
        ax.scatter(
            df.loc[down_mask, 'log2FoldChange'],
            df.loc[down_mask, 'neg_log10_padj'],
            c='#7BA3C9', alpha=0.8, s=50, label=f'Down ({down_mask.sum()})',
            edgecolors='#4A6D8C', linewidths=1.0
        )

        # Add threshold lines - thicker, darker dashes
        ax.axhline(-np.log10(padj_threshold), color='#333333', linestyle='--',
                   alpha=0.8, linewidth=2, dashes=(8, 4))
        ax.axvline(lfc_threshold, color='#333333', linestyle='--',
                   alpha=0.8, linewidth=2, dashes=(8, 4))
        ax.axvline(-lfc_threshold, color='#333333', linestyle='--',
                   alpha=0.8, linewidth=2, dashes=(8, 4))

        # Highlight specific genes with labels and lines
        if highlight_genes:
            for gene in highlight_genes:
                if gene in df.index:
                    x = df.loc[gene, 'log2FoldChange']
                    y = df.loc[gene, 'neg_log10_padj']
                    # Plot highlighted point
                    ax.scatter([x], [y], c='#FFD700', s=150,
                              edgecolors='#000000', linewidths=2.5, zorder=10)
                    # Add annotation with arrow/line
                    ax.annotate(
                        gene, xy=(x, y), xytext=(25, 25),
                        textcoords='offset points',
                        fontsize=11, fontweight='bold',
                        bbox=dict(boxstyle='round,pad=0.4', facecolor='#FFFF99',
                                 edgecolor='#000000', linewidth=1.5),
                        arrowprops=dict(arrowstyle='->', color='#000000',
                                       linewidth=2, connectionstyle='arc3,rad=0.2')
                    )

        # Styling - use plain text labels to avoid font issues
        ax.set_xlabel('Log2 Fold Change', fontsize=14, fontweight='bold')
        ax.set_ylabel('-Log10 Adjusted P-value', fontsize=14, fontweight='bold')
        ax.set_title('Volcano Plot', fontsize=16, fontweight='bold', pad=15)

        # Legend with frame
        legend = ax.legend(loc='upper right', frameon=True, fancybox=False,
                          edgecolor='#333333', fontsize=12)
        legend.get_frame().set_linewidth(1.5)

        # Clean up spines - only left and bottom, make them thicker
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_linewidth(2)
        ax.spines['bottom'].set_linewidth(2)
        ax.spines['left'].set_color('#333333')
        ax.spines['bottom'].set_color('#333333')

        # Tick styling
        ax.tick_params(axis='both', which='major', labelsize=12, width=2, length=6)

        plt.tight_layout()

        # Convert to output format
        if output_format == 'svg':
            image = self._fig_to_svg(fig)
        else:
            image = self._fig_to_base64_png(fig)

        return {
            'image': image,
            'format': output_format,
            'stats': {
                'total': len(df),
                'up': int(up_mask.sum()),
                'down': int(down_mask.sum()),
                'ns': int(ns_mask.sum())
            }
        }

    def create_boxplot(
        self,
        normalized_counts: pd.DataFrame,
        metadata: pd.DataFrame,
        genes: List[str],
        group_column: str = 'condition',
        results: Optional[pd.DataFrame] = None,
        normalization: str = 'log2',
        output_format: str = 'svg'
    ) -> Dict:
        """
        Create grouped boxplots with significance stars.

        Parameters
        ----------
        normalized_counts : pd.DataFrame
            Normalized count matrix (samples x genes)
        metadata : pd.DataFrame
            Sample metadata
        genes : list
            Genes to plot
        group_column : str
            Column for grouping
        results : pd.DataFrame, optional
            DESeq2 results for significance
        normalization : str
            'log2' or 'raw'
        output_format : str
            'svg' or 'png'

        Returns
        -------
        dict with 'image' and metadata
        """
        if not genes:
            return {'error': 'No genes selected'}

        # Filter to valid genes
        valid_genes = [g for g in genes if g in normalized_counts.columns]
        if not valid_genes:
            return {'error': 'None of the selected genes found in data'}

        # Prepare expression data
        expr = normalized_counts[valid_genes].copy()

        if normalization == 'log2':
            expr = np.log2(expr + 1)
            ylabel = 'Log2(Normalized Counts + 1)'
        else:
            ylabel = 'Normalized Counts'

        # Melt for seaborn - reset index to avoid ambiguity
        expr = expr.reset_index(drop=True)
        expr['sample_id'] = normalized_counts.index.tolist()

        # Add group column from metadata
        sample_to_group = metadata[group_column].to_dict()
        expr[group_column] = expr['sample_id'].map(sample_to_group)

        plot_data = expr.melt(
            id_vars=['sample_id', group_column],
            var_name='gene',
            value_name='expression'
        )

        # Create figure
        n_genes = len(valid_genes)
        fig_width = max(8, n_genes * 2.5)
        fig, ax = plt.subplots(figsize=(fig_width, 7))

        # Clean white background
        ax.set_facecolor('white')
        fig.patch.set_facecolor('white')

        # Color palette - muted pastels
        unique_groups = metadata[group_column].unique()
        if 'control' in unique_groups and 'treatment' in unique_groups:
            palette = {'control': '#7BA3C9', 'treatment': '#C97B7B'}
        else:
            palette = dict(zip(unique_groups, sns.color_palette('pastel', len(unique_groups))))

        # Boxplot with thicker lines
        box_props = dict(linewidth=2)
        whisker_props = dict(linewidth=2)
        cap_props = dict(linewidth=2)
        median_props = dict(linewidth=2.5, color='black')

        sns.boxplot(
            data=plot_data, x='gene', y='expression', hue=group_column,
            ax=ax, palette=palette, width=0.6, linewidth=2,
            boxprops=box_props, whiskerprops=whisker_props,
            capprops=cap_props, medianprops=median_props
        )
        sns.stripplot(
            data=plot_data, x='gene', y='expression', hue=group_column,
            ax=ax, dodge=True, alpha=0.7, size=6, legend=False,
            edgecolor='#333333', linewidth=0.5
        )

        # Add significance stars
        def pvalue_to_stars(pval):
            if pd.isna(pval):
                return ''
            if pval < 0.001:
                return '***'
            if pval < 0.01:
                return '**'
            if pval < 0.05:
                return '*'
            return 'ns'

        if results is not None:
            ymax = ax.get_ylim()[1]
            for i, gene in enumerate(valid_genes):
                if gene in results.index:
                    padj = results.loc[gene, 'padj']
                    stars = pvalue_to_stars(padj)
                    if stars:
                        ax.text(i, ymax * 0.97, stars, ha='center', va='top',
                               fontsize=14, fontweight='bold')

        # Styling
        ax.set_xlabel('Gene', fontsize=14, fontweight='bold')
        ax.set_ylabel(ylabel, fontsize=14, fontweight='bold')
        ax.set_title('Gene Expression by Condition', fontsize=16, fontweight='bold', pad=15)

        # Legend
        legend = ax.legend(title=group_column.capitalize(), loc='upper right',
                          frameon=True, fancybox=False, edgecolor='#333333', fontsize=11)
        legend.get_frame().set_linewidth(1.5)
        legend.get_title().set_fontweight('bold')

        plt.xticks(rotation=45, ha='right', fontsize=12)

        # Clean up spines
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_linewidth(2)
        ax.spines['bottom'].set_linewidth(2)
        ax.spines['left'].set_color('#333333')
        ax.spines['bottom'].set_color('#333333')

        ax.tick_params(axis='both', which='major', labelsize=12, width=2, length=6)

        plt.tight_layout()

        # Convert to output format
        if output_format == 'svg':
            image = self._fig_to_svg(fig)
        else:
            image = self._fig_to_base64_png(fig)

        return {
            'image': image,
            'format': output_format,
            'genes': valid_genes
        }

    def create_heatmap(
        self,
        normalized_counts: pd.DataFrame,
        metadata: pd.DataFrame,
        deseq_results: pd.DataFrame = None,
        genes: Optional[List[str]] = None,
        top_n: int = 50,
        gene_selection: str = 'variance',
        annotation_columns: List[str] = None,
        scale: str = 'row',
        cluster_cols: bool = True,
        cluster_rows: bool = True,
        show_dendrogram: bool = True,
        output_format: str = 'svg'
    ) -> Dict:
        """
        Create clustered heatmap with annotations and dendrograms.

        Parameters
        ----------
        normalized_counts : pd.DataFrame
            Normalized count matrix (samples x genes)
        metadata : pd.DataFrame
            Sample metadata
        deseq_results : pd.DataFrame, optional
            DESeq2 results for gene selection
        genes : list, optional
            Specific genes to include
        top_n : int
            Number of top genes if genes not specified
        gene_selection : str
            'variance', 'up', 'down', or 'significant'
        annotation_columns : list
            Metadata columns for annotation
        scale : str
            'row', 'column', or 'none'
        cluster_cols : bool
            Whether to cluster columns (samples)
        cluster_rows : bool
            Whether to cluster rows (genes)
        show_dendrogram : bool
            Whether to show dendrograms
        output_format : str
            'svg' or 'png'

        Returns
        -------
        dict with 'image' and metadata
        """
        annotation_columns = annotation_columns or ['condition']

        # Select genes based on selection method
        if genes:
            valid_genes = [g for g in genes if g in normalized_counts.columns]
        elif deseq_results is not None and gene_selection != 'variance':
            df = deseq_results.dropna(subset=['log2FoldChange', 'padj']).copy()
            if gene_selection == 'up':
                df = df[df['log2FoldChange'] > 0].nlargest(top_n, 'log2FoldChange')
            elif gene_selection == 'down':
                df = df[df['log2FoldChange'] < 0].nsmallest(top_n, 'log2FoldChange')
            elif gene_selection == 'significant':
                df['abs_lfc'] = df['log2FoldChange'].abs()
                df = df.nsmallest(top_n, 'padj')
            valid_genes = [g for g in df.index if g in normalized_counts.columns]
        else:
            # Use top variable genes
            gene_var = normalized_counts.var()
            valid_genes = list(gene_var.nlargest(top_n).index)

        if len(valid_genes) == 0:
            return {'error': 'No valid genes found for heatmap'}

        # Get expression matrix
        expr = normalized_counts[valid_genes].copy()

        # Log transform
        expr = np.log2(expr + 1)

        # Scale
        if scale == 'row':
            expr = expr.apply(lambda x: (x - x.mean()) / x.std() if x.std() > 0 else x, axis=0)
        elif scale == 'column':
            expr = expr.apply(lambda x: (x - x.mean()) / x.std() if x.std() > 0 else x, axis=1)

        # Prepare annotation colors - muted pastel palette
        color_maps = {
            'condition': {'control': '#7BA3C9', 'treatment': '#C97B7B'},
            'batch': {'batch1': '#8FB996', 'batch2': '#B39DDB', 'batch3': '#E8C87D'},
            'sex': {'M': '#7BA3C9', 'F': '#D4A5A5', 'male': '#7BA3C9', 'female': '#D4A5A5'}
        }

        col_colors = pd.DataFrame(index=expr.index)
        for col in annotation_columns:
            if col in metadata.columns:
                colors = color_maps.get(col, {})
                if not colors:
                    unique_vals = metadata[col].unique()
                    palette = sns.color_palette('pastel', len(unique_vals))
                    colors = dict(zip(unique_vals, [plt.matplotlib.colors.rgb2hex(c) for c in palette]))
                col_colors[col] = metadata.loc[expr.index, col].map(colors)

        # Dendrogram ratio based on show_dendrogram setting
        dendro_ratio = 0.15 if show_dendrogram else 0.001

        # Create clustermap with improved styling
        g = sns.clustermap(
            expr.T,
            cmap='RdBu_r',
            figsize=(14, max(10, len(valid_genes) * 0.25)),
            row_cluster=cluster_rows,
            col_cluster=cluster_cols,
            col_colors=col_colors if len(col_colors.columns) > 0 else None,
            linewidths=0.3,
            linecolor='white',
            xticklabels=True,
            yticklabels=True if len(valid_genes) <= 100 else False,
            center=0,
            cbar_pos=(0.02, 0.8, 0.03, 0.15),
            dendrogram_ratio=(dendro_ratio, dendro_ratio),
            tree_kws={'linewidths': 1.5, 'colors': '#5A6B7A'} if show_dendrogram else {}
        )

        # Hide dendrograms if not showing
        if not show_dendrogram:
            g.ax_row_dendrogram.set_visible(False)
            g.ax_col_dendrogram.set_visible(False)

        g.ax_heatmap.set_xlabel('Samples', fontsize=14, fontweight='bold')
        g.ax_heatmap.set_ylabel('Genes', fontsize=14, fontweight='bold')
        plt.setp(g.ax_heatmap.xaxis.get_majorticklabels(), rotation=45, ha='right', fontsize=10)
        plt.setp(g.ax_heatmap.yaxis.get_majorticklabels(), fontsize=9)

        # Style the colorbar
        cbar = g.ax_cbar
        cbar.tick_params(labelsize=10, width=1.5, length=4)
        for spine in cbar.spines.values():
            spine.set_linewidth(1.5)

        # Convert to output format
        if output_format == 'svg':
            buf = io.BytesIO()
            g.savefig(buf, format='svg', bbox_inches='tight', facecolor='white', edgecolor='none')
            buf.seek(0)
            image = buf.getvalue().decode('utf-8')
        else:
            buf = io.BytesIO()
            g.savefig(buf, format='png', bbox_inches='tight', dpi=150, facecolor='white', edgecolor='none')
            buf.seek(0)
            image = f"data:image/png;base64,{base64.b64encode(buf.getvalue()).decode('utf-8')}"

        plt.close('all')

        return {
            'image': image,
            'format': output_format,
            'genes': valid_genes,
            'samples': list(expr.index)
        }

    def create_gsea_dotplot(
        self,
        gsea_results: pd.DataFrame,
        top_n: int = 15,
        output_format: str = 'svg'
    ) -> Dict:
        """
        Create GSEA enrichment dot plot.

        Parameters
        ----------
        gsea_results : pd.DataFrame
            GSEA results
        top_n : int
            Number of top pathways to show
        output_format : str
            'svg' or 'png'

        Returns
        -------
        dict with 'image'
        """
        if gsea_results is None or len(gsea_results) == 0:
            return {'error': 'No GSEA results available'}

        df = gsea_results.head(top_n).copy()

        # Get term and NES columns (handle different naming conventions)
        term_col = 'Term' if 'Term' in df.columns else 'term'
        nes_col = 'NES' if 'NES' in df.columns else 'nes'
        fdr_col = 'FDR q-val' if 'FDR q-val' in df.columns else ('fdr' if 'fdr' in df.columns else 'Adjusted P-value')

        fig, ax = plt.subplots(figsize=(12, max(8, len(df) * 0.5)))

        # Clean white background
        ax.set_facecolor('white')
        fig.patch.set_facecolor('white')

        # Get values
        terms = df[term_col].tolist()
        nes_values = df[nes_col].tolist()
        fdr_values = df[fdr_col].replace(0, 1e-10).tolist()

        # Color by NES sign - muted pastels
        colors = ['#C97B7B' if x > 0 else '#7BA3C9' for x in nes_values]

        # Size by -log10(FDR)
        sizes = [-np.log10(x) * 40 for x in fdr_values]
        sizes = [max(80, min(s, 400)) for s in sizes]  # Clamp sizes

        # Truncate long term names
        terms = [t[:45] + '...' if len(t) > 45 else t for t in terms]

        y_positions = range(len(terms))

        scatter = ax.scatter(nes_values, y_positions, c=colors, s=sizes, alpha=0.8,
                            edgecolors='#333333', linewidths=1.5)

        ax.set_yticks(y_positions)
        ax.set_yticklabels(terms, fontsize=11)
        ax.set_xlabel('Normalized Enrichment Score (NES)', fontsize=14, fontweight='bold')
        ax.set_title('Pathway Enrichment Analysis', fontsize=16, fontweight='bold', pad=15)
        ax.axvline(0, color='#333333', linestyle='--', alpha=0.8, linewidth=2, dashes=(8, 4))

        # Clean up spines
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_linewidth(2)
        ax.spines['bottom'].set_linewidth(2)
        ax.spines['left'].set_color('#333333')
        ax.spines['bottom'].set_color('#333333')

        ax.tick_params(axis='both', which='major', labelsize=12, width=2, length=6)
        ax.invert_yaxis()

        plt.tight_layout()

        # Convert to output format
        if output_format == 'svg':
            image = self._fig_to_svg(fig)
        else:
            image = self._fig_to_base64_png(fig)

        return {
            'image': image,
            'format': output_format,
            'n_pathways': len(df)
        }

    def create_pca_plot(
        self,
        normalized_counts: pd.DataFrame,
        metadata: pd.DataFrame,
        color_by: str = 'condition',
        normalization: str = 'log2',
        center: bool = True,
        scale: bool = True,
        show_labels: bool = True,
        n_components: int = 2,
        output_format: str = 'svg'
    ) -> Dict:
        """
        Create PCA scatter plot with scree plot inset, plus per-sample hover data.
        """
        from sklearn.decomposition import PCA
        from sklearn.preprocessing import StandardScaler

        # Transform data based on normalization choice
        if normalization == 'log2':
            data = np.log2(normalized_counts + 1)
        elif normalization == 'raw':
            data = normalized_counts.copy()
        else:
            data = normalized_counts.copy()

        # Center and/or scale
        if center and scale:
            scaler = StandardScaler()
            scaled_counts = scaler.fit_transform(data)
        elif center:
            scaled_counts = data - data.mean()
        elif scale:
            scaled_counts = data / data.std()
        else:
            scaled_counts = data.values if hasattr(data, 'values') else data

        # Run PCA with more components for scree plot
        max_comps = min(10, min(scaled_counts.shape))
        pca_full = PCA(n_components=max_comps)
        pca_full_result = pca_full.fit_transform(scaled_counts)
        all_var = pca_full.explained_variance_ratio_ * 100

        # Use first 2 for the scatter
        pca_result = pca_full_result[:, :2]
        var_explained = all_var[:2]

        # Create DataFrame with PCA results
        pca_df = pd.DataFrame(
            pca_result,
            index=normalized_counts.index,
            columns=['PC1', 'PC2']
        )

        # Add metadata
        for col in metadata.columns:
            pca_df[col] = metadata.loc[pca_df.index, col].values

        # --- Build per-sample hover data for the frontend ---
        sample_data = []
        for idx, row in pca_df.iterrows():
            entry = {'name': str(idx), 'pc1': float(row['PC1']), 'pc2': float(row['PC2'])}
            for col in metadata.columns:
                entry[col] = str(row[col])
            sample_data.append(entry)

        # Create figure: PCA scatter (top) + scree plot (bottom)
        fig, (ax, ax_scree) = plt.subplots(2, 1, figsize=(10, 11),
                                            gridspec_kw={'height_ratios': [3, 1.2], 'hspace': 0.35})
        ax.set_facecolor('white')
        ax_scree.set_facecolor('white')
        fig.patch.set_facecolor('white')

        # Get unique values and create color palette
        if color_by in pca_df.columns:
            unique_vals = pca_df[color_by].unique()

            color_maps = {
                'condition': {'control': '#7BA3C9', 'treatment': '#C97B7B'},
                'batch': {'batch1': '#8FB996', 'batch2': '#B39DDB', 'batch3': '#E8C87D'},
                'sex': {'M': '#7BA3C9', 'F': '#D4A5A5', 'male': '#7BA3C9', 'female': '#D4A5A5'}
            }

            if color_by in color_maps:
                palette = color_maps[color_by]
                for val in unique_vals:
                    if val not in palette:
                        palette[val] = sns.color_palette('husl', len(unique_vals))[list(unique_vals).index(val)]
            else:
                colors = sns.color_palette('husl', len(unique_vals))
                palette = dict(zip(unique_vals, [plt.matplotlib.colors.rgb2hex(c) for c in colors]))

            for val in unique_vals:
                mask = pca_df[color_by] == val
                ax.scatter(
                    pca_df.loc[mask, 'PC1'],
                    pca_df.loc[mask, 'PC2'],
                    c=palette.get(val, '#666666'),
                    s=120,
                    label=str(val),
                    alpha=0.8,
                    edgecolors='#333333',
                    linewidths=1.5
                )

            if show_labels:
                for idx, row in pca_df.iterrows():
                    ax.annotate(
                        idx,
                        (row['PC1'], row['PC2']),
                        xytext=(5, 5),
                        textcoords='offset points',
                        fontsize=8,
                        alpha=0.7
                    )

            legend = ax.legend(title=color_by.capitalize(), loc='best',
                              frameon=True, fancybox=False, edgecolor='#333333', fontsize=11)
            legend.get_frame().set_linewidth(1.5)
            legend.get_title().set_fontweight('bold')
        else:
            ax.scatter(pca_df['PC1'], pca_df['PC2'], c='#3498DB', s=120,
                      alpha=0.8, edgecolors='#333333', linewidths=1.5)

        ax.set_xlabel(f'PC1 ({var_explained[0]:.1f}% variance)', fontsize=14, fontweight='bold')
        ax.set_ylabel(f'PC2 ({var_explained[1]:.1f}% variance)', fontsize=14, fontweight='bold')
        ax.set_title('Principal Component Analysis', fontsize=16, fontweight='bold', pad=15)

        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_linewidth(2)
        ax.spines['bottom'].set_linewidth(2)
        ax.spines['left'].set_color('#333333')
        ax.spines['bottom'].set_color('#333333')

        ax.tick_params(axis='both', which='major', labelsize=12, width=2, length=6)

        ax.axhline(0, color='#CCCCCC', linestyle='--', linewidth=1, zorder=0)
        ax.axvline(0, color='#CCCCCC', linestyle='--', linewidth=1, zorder=0)

        # --- Scree plot (bottom subplot) ---
        n_show = min(len(all_var), 10)
        pcs = list(range(1, n_show + 1))
        cumulative = np.cumsum(all_var[:n_show])

        ax_scree.bar(pcs, all_var[:n_show], color='#c45a3c', alpha=0.75, edgecolor='#333', linewidth=0.8)
        ax_scree.plot(pcs, cumulative, 'o-', color='#333333', markersize=5, linewidth=1.5)

        ax_scree.set_xlabel('Principal Component', fontsize=12, fontweight='bold')
        ax_scree.set_ylabel('% Variance Explained', fontsize=12, fontweight='bold')
        ax_scree.set_title('Scree Plot', fontsize=14, fontweight='bold')
        ax_scree.set_xticks(pcs)
        ax_scree.tick_params(axis='both', labelsize=10, width=1.5, length=5)
        ax_scree.spines['top'].set_visible(False)
        ax_scree.spines['right'].set_visible(False)
        ax_scree.spines['left'].set_linewidth(1.5)
        ax_scree.spines['bottom'].set_linewidth(1.5)
        ax_scree.spines['left'].set_color('#333333')
        ax_scree.spines['bottom'].set_color('#333333')
        ax_scree.set_ylim(0, max(all_var[:n_show]) * 1.15)

        # Right y-axis for cumulative
        ax_scree2 = ax_scree.twinx()
        ax_scree2.set_ylabel('Cumulative %', fontsize=11, color='#555')
        ax_scree2.set_ylim(0, 105)
        ax_scree2.tick_params(axis='y', labelsize=10, colors='#555')
        ax_scree2.spines['top'].set_visible(False)
        ax_scree2.spines['left'].set_linewidth(1.5)

        plt.tight_layout()

        if output_format == 'svg':
            image = self._fig_to_svg(fig)
        else:
            image = self._fig_to_base64_png(fig)

        return {
            'image': image,
            'format': output_format,
            'variance_explained': [float(v) for v in all_var],
            'metadata_columns': list(metadata.columns),
            'sample_data': sample_data,
        }


# Singleton instance
plot_generator = PlotGenerator()
