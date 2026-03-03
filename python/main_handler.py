#!/usr/bin/env python3
"""
Main handler for Python backend.
Receives JSON commands from Electron via stdin and returns JSON results via stdout.
"""

import os
import sys
import json
import traceback
from typing import Any, Dict

# Add parent directory to path for imports
sys.path.insert(0, '.')

from data.example_generator import generate_example_data
from data.loader import load_counts, load_metadata
from analysis.deseq_runner import DESeqRunner
from analysis.gsea_runner import GSEARunner
from analysis.normalization import normalize_counts, log2_transform
from analysis.plot_generator import PlotGenerator
from script_generation.script_builder import ScriptBuilder
from kallisto.runner import KallistoRunner, APP_CONFIG
from kallisto.example_data import create_example_fasta, create_example_fastqs
from edgepy_runner import EdgePythonRunner

# Global state
class AppState:
    def __init__(self):
        self.counts = None
        self.metadata = None
        self.deseq_runner = None
        self.deseq_results = None
        self.normalized_counts = None
        self.gsea_runner = None
        self.gsea_results = None
        self.selected_genes = []
        self.plot_generator = PlotGenerator()
        self.kallisto_runner = KallistoRunner()
        self.kallisto_index_path = None
        self.edgepy_runner = None
        self.active_tool = 'pydeseq2'  # or 'edgepython'
        self.kallisto_sample_results = {}  # sample_name -> output_dir
        self.analysis_params = {
            'condition_column': 'condition',
            'reference_level': 'control',
            'lfc_threshold': 1.0,
            'padj_threshold': 0.05
        }

state = AppState()


def send_response(response_id: int, data: Any = None, error: str = None):
    """Send JSON response to stdout."""
    response = {'id': response_id}
    if error:
        response['error'] = error
    else:
        response['data'] = data
    print(json.dumps(response), flush=True)


def handle_load_example_data(params: Dict) -> Dict:
    """Generate and load example data."""
    n_samples = params.get('n_samples', 20)
    n_genes = params.get('n_genes', 5000)
    n_de_genes = params.get('n_de_genes', 500)

    counts_df, metadata_df = generate_example_data(
        n_samples=n_samples,
        n_genes=n_genes,
        n_de_genes=n_de_genes
    )

    state.counts = counts_df
    state.metadata = metadata_df

    # Return preview data
    return {
        'counts_shape': list(counts_df.shape),
        'counts_preview': counts_df.iloc[:10, :5].to_dict(),
        'counts_columns': list(counts_df.columns),
        'counts_index': list(counts_df.index[:20]),
        'metadata_shape': list(metadata_df.shape),
        'metadata_preview': metadata_df.head(10).to_dict(),
        'metadata_columns': list(metadata_df.columns),
        'metadata_index': list(metadata_df.index)
    }


def handle_load_counts(params: Dict) -> Dict:
    """Load counts from file."""
    file_path = params['file_path']
    state.counts = load_counts(file_path)

    return {
        'counts_shape': list(state.counts.shape),
        'counts_preview': state.counts.iloc[:10, :5].to_dict(),
        'counts_columns': list(state.counts.columns),
        'counts_index': list(state.counts.index[:20])
    }


def handle_load_metadata(params: Dict) -> Dict:
    """Load metadata from file."""
    file_path = params['file_path']
    state.metadata = load_metadata(file_path)

    return {
        'metadata_shape': list(state.metadata.shape),
        'metadata_preview': state.metadata.head(10).to_dict(),
        'metadata_columns': list(state.metadata.columns),
        'metadata_index': list(state.metadata.index)
    }


def handle_list_project_files(params: Dict) -> Dict:
    """List files in a project directory."""
    path = params.get('path', '')
    if not path or not os.path.isdir(path):
        return {'files': []}

    files = []
    valid_extensions = {'.csv', '.tsv', '.txt', '.fastq', '.fq', '.gz', '.fa', '.fasta', '.h5ad'}
    try:
        for name in sorted(os.listdir(path)):
            full_path = os.path.join(path, name)
            if os.path.isfile(full_path):
                ext = os.path.splitext(name)[1].lower()
                # Check double extensions like .fastq.gz
                if ext == '.gz':
                    base_ext = os.path.splitext(os.path.splitext(name)[0])[1].lower()
                    if base_ext in valid_extensions:
                        files.append({'name': name, 'path': full_path})
                elif ext in valid_extensions:
                    files.append({'name': name, 'path': full_path})
    except OSError:
        pass

    return {'files': files}


def handle_set_params(params: Dict) -> Dict:
    """Update analysis parameters."""
    state.analysis_params.update(params)
    return {'params': state.analysis_params}


def handle_run_deseq(params: Dict) -> Dict:
    """Run DESeq2 analysis."""
    if state.counts is None or state.metadata is None:
        raise ValueError("Load data first")

    # Update params if provided
    if params:
        state.analysis_params.update(params)

    condition_col = state.analysis_params['condition_column']
    reference = state.analysis_params['reference_level']
    lfc_thresh = state.analysis_params['lfc_threshold']
    padj_thresh = state.analysis_params['padj_threshold']

    # Run DESeq2
    state.deseq_runner = DESeqRunner()
    state.deseq_results = state.deseq_runner.run_analysis(
        counts=state.counts,
        metadata=state.metadata,
        condition_column=condition_col,
        reference_level=reference
    )

    # Get normalized counts
    state.normalized_counts = state.deseq_runner.get_normalized_counts()

    # Calculate summary statistics
    sig_mask = (state.deseq_results['padj'] < padj_thresh) & \
               (abs(state.deseq_results['log2FoldChange']) > lfc_thresh)
    up_mask = sig_mask & (state.deseq_results['log2FoldChange'] > 0)
    down_mask = sig_mask & (state.deseq_results['log2FoldChange'] < 0)

    # Return results
    results_dict = state.deseq_results.reset_index().to_dict(orient='records')

    return {
        'results': results_dict,
        'summary': {
            'total_genes': len(state.deseq_results),
            'significant': int(sig_mask.sum()),
            'upregulated': int(up_mask.sum()),
            'downregulated': int(down_mask.sum())
        }
    }


def handle_set_active_tool(params: Dict) -> Dict:
    """Set the active differential testing tool."""
    tool = params.get('tool', 'pydeseq2')
    if tool not in ('pydeseq2', 'edgepython'):
        raise ValueError(f'Unknown tool: {tool}')
    state.active_tool = tool
    return {'tool': tool}


def handle_run_edgepython(params: Dict) -> Dict:
    """Run edgePython differential expression analysis."""
    if state.counts is None or state.metadata is None:
        raise ValueError("Load data first")

    # Update params if provided
    if params:
        state.analysis_params.update(params)

    condition_col = state.analysis_params['condition_column']
    reference = state.analysis_params['reference_level']
    lfc_thresh = state.analysis_params['lfc_threshold']
    padj_thresh = state.analysis_params['padj_threshold']

    # edgePython-specific params
    normalization = params.get('normalization', 'TMM')
    test_method = params.get('test_method', 'qlf')
    dispersion = params.get('dispersion', 'trended')
    robust = params.get('robust', True)
    if isinstance(robust, str):
        robust = robust.lower() == 'true'
    min_count = int(params.get('min_count', 10))

    # Run edgePython
    state.edgepy_runner = EdgePythonRunner()
    state.active_tool = 'edgepython'

    state.deseq_results = state.edgepy_runner.run_analysis(
        counts=state.counts,
        metadata=state.metadata,
        condition_column=condition_col,
        reference_level=reference,
        normalization=normalization,
        test_method=test_method,
        dispersion=dispersion,
        robust=robust,
        min_count=min_count,
    )

    # Get normalized counts
    state.normalized_counts = state.edgepy_runner.get_normalized_counts()

    # Calculate summary statistics (same as DESeq2)
    sig_mask = (state.deseq_results['padj'] < padj_thresh) & \
               (abs(state.deseq_results['log2FoldChange']) > lfc_thresh)
    up_mask = sig_mask & (state.deseq_results['log2FoldChange'] > 0)
    down_mask = sig_mask & (state.deseq_results['log2FoldChange'] < 0)

    # Return results in same format
    results_dict = state.deseq_results.reset_index().to_dict(orient='records')

    return {
        'results': results_dict,
        'tool': 'edgepython',
        'summary': {
            'total_genes': len(state.deseq_results),
            'significant': int(sig_mask.sum()),
            'upregulated': int(up_mask.sum()),
            'downregulated': int(down_mask.sum())
        }
    }


def handle_get_volcano_data(params: Dict) -> Dict:
    """Get data for volcano plot."""
    if state.deseq_results is None:
        raise ValueError("Run differential expression analysis first")

    lfc_thresh = params.get('lfc_threshold', state.analysis_params['lfc_threshold'])
    padj_thresh = params.get('padj_threshold', state.analysis_params['padj_threshold'])

    df = state.deseq_results.copy()
    df = df.dropna(subset=['log2FoldChange', 'padj'])

    # Categorize genes
    df['category'] = 'ns'
    df.loc[(df['log2FoldChange'] >= lfc_thresh) & (df['padj'] < padj_thresh), 'category'] = 'up'
    df.loc[(df['log2FoldChange'] <= -lfc_thresh) & (df['padj'] < padj_thresh), 'category'] = 'down'

    # Calculate -log10(padj) with floor for zeros
    import numpy as np
    df['neg_log10_padj'] = -np.log10(df['padj'].replace(0, 1e-300))

    return {
        'genes': list(df.index),
        'log2FoldChange': list(df['log2FoldChange']),
        'neg_log10_padj': list(df['neg_log10_padj']),
        'padj': list(df['padj']),
        'category': list(df['category']),
        'lfc_threshold': lfc_thresh,
        'padj_threshold': padj_thresh
    }


def handle_get_boxplot_data(params: Dict) -> Dict:
    """Get data for boxplots of selected genes."""
    if state.normalized_counts is None:
        raise ValueError("Run DESeq2 analysis first")

    genes = params.get('genes', state.selected_genes)
    normalization = params.get('normalization', 'deseq2')  # 'deseq2' or 'log2'
    group_column = params.get('group_column', state.analysis_params['condition_column'])

    if not genes:
        return {'error': 'No genes selected'}

    # Filter to valid genes
    valid_genes = [g for g in genes if g in state.normalized_counts.columns]
    if not valid_genes:
        return {'error': 'None of the selected genes found in data'}

    import numpy as np

    # Get expression data - make sure samples are aligned with metadata
    samples = list(state.normalized_counts.index)
    expr_data = state.normalized_counts[valid_genes].copy()

    if normalization == 'log2':
        expr_data = np.log2(expr_data + 1)

    # Get group assignments - must be aligned with expr_data's index
    groups = [state.metadata.loc[sample, group_column] for sample in samples]
    unique_groups = list(state.metadata[group_column].unique())

    # Get significance from DESeq2 results
    significance = {}
    for gene in valid_genes:
        if gene in state.deseq_results.index:
            padj = state.deseq_results.loc[gene, 'padj']
            if padj < 0.001:
                significance[gene] = '***'
            elif padj < 0.01:
                significance[gene] = '**'
            elif padj < 0.05:
                significance[gene] = '*'
            else:
                significance[gene] = 'ns'
        else:
            significance[gene] = 'na'

    return {
        'genes': valid_genes,
        'expression': expr_data.to_dict(orient='list'),
        'samples': samples,
        'groups': groups,
        'group_column': group_column,
        'unique_groups': unique_groups,
        'significance': significance,
        'normalization': normalization
    }


def handle_get_heatmap_data(params: Dict) -> Dict:
    """Get data for heatmap."""
    if state.normalized_counts is None:
        raise ValueError("Run DESeq2 analysis first")

    import numpy as np
    from scipy.cluster.hierarchy import linkage, leaves_list
    from scipy.spatial.distance import pdist

    top_n = params.get('top_n', 50)
    genes = params.get('genes', None)
    scale = params.get('scale', 'row')  # 'row', 'column', 'none'
    annotation_columns = params.get('annotation_columns', ['condition'])

    # Select genes
    if genes:
        valid_genes = [g for g in genes if g in state.normalized_counts.columns]
    else:
        # Use top variable genes
        gene_var = state.normalized_counts.var()
        valid_genes = list(gene_var.nlargest(top_n).index)

    # Get expression matrix (samples x genes)
    expr = state.normalized_counts[valid_genes].copy()

    # Log transform
    expr = np.log2(expr + 1)

    # Scale
    if scale == 'row':
        # Z-score per gene (across samples)
        expr = expr.apply(lambda x: (x - x.mean()) / x.std() if x.std() > 0 else x, axis=0)
    elif scale == 'column':
        # Z-score per sample (across genes)
        expr = expr.apply(lambda x: (x - x.mean()) / x.std() if x.std() > 0 else x, axis=1)

    # Cluster samples and genes
    try:
        # Cluster samples (rows)
        sample_dist = pdist(expr.values)
        sample_linkage = linkage(sample_dist, method='average')
        sample_order = leaves_list(sample_linkage)

        # Cluster genes (columns)
        gene_dist = pdist(expr.values.T)
        gene_linkage = linkage(gene_dist, method='average')
        gene_order = leaves_list(gene_linkage)

        # Reorder
        expr = expr.iloc[sample_order, gene_order]
    except Exception:
        # If clustering fails, keep original order
        pass

    # Get annotations
    annotations = {}
    for col in annotation_columns:
        if col in state.metadata.columns:
            annotations[col] = state.metadata.loc[expr.index, col].tolist()

    return {
        'genes': list(expr.columns),
        'samples': list(expr.index),
        'values': expr.values.tolist(),
        'z': expr.values.T.tolist(),  # genes x samples for Plotly heatmap
        'annotations': annotations,
        'scale': scale
    }


def handle_run_gsea(params: Dict) -> Dict:
    """Run GSEA pathway analysis."""
    if state.deseq_results is None:
        raise ValueError("Run DESeq2 analysis first")

    gene_sets = params.get('gene_sets', 'KEGG_2021_Human')
    method = params.get('method', 'prerank')  # 'prerank' or 'enrichr'

    state.gsea_runner = GSEARunner()

    if method == 'prerank':
        # Create ranking metric
        ranking = state.gsea_runner.create_ranking_metric(state.deseq_results)
        state.gsea_results = state.gsea_runner.run_prerank(
            ranking=ranking,
            gene_sets=gene_sets
        )
    else:
        # Use significant genes for enrichr
        lfc_thresh = state.analysis_params['lfc_threshold']
        padj_thresh = state.analysis_params['padj_threshold']

        sig_up = state.deseq_results[
            (state.deseq_results['padj'] < padj_thresh) &
            (state.deseq_results['log2FoldChange'] > lfc_thresh)
        ].index.tolist()

        if sig_up:
            state.gsea_results = state.gsea_runner.run_enrichr(
                gene_list=sig_up,
                gene_sets=gene_sets
            )
        else:
            return {'results': [], 'message': 'No significant upregulated genes found'}

    # Convert results to serializable format
    results_list = state.gsea_results.to_dict(orient='records')

    return {
        'results': results_list,
        'method': method,
        'gene_sets': gene_sets
    }


def handle_set_selected_genes(params: Dict) -> Dict:
    """Update selected genes list."""
    state.selected_genes = params.get('genes', [])
    return {'selected_genes': state.selected_genes}


def handle_generate_script(params: Dict) -> Dict:
    """Generate reproducible Python script adapted to the active tool."""
    include_sections = params.get('include_sections', {
        'data_loading': True,
        'deseq2': True,
        'volcano': True,
        'boxplot': True,
        'heatmap': True,
        'gsea': False
    })

    builder = ScriptBuilder(
        analysis_params=state.analysis_params,
        selected_genes=state.selected_genes,
        active_tool=state.active_tool,
    )

    script = builder.generate_script(include_sections)

    return {'script': script}


def handle_get_gene_info(params: Dict) -> Dict:
    """Get detailed info for a specific gene."""
    gene = params['gene']

    if state.deseq_results is None:
        raise ValueError("Run DESeq2 analysis first")

    if gene not in state.deseq_results.index:
        return {'error': f'Gene {gene} not found'}

    result = state.deseq_results.loc[gene].to_dict()

    # Add expression data if available
    if state.normalized_counts is not None and gene in state.normalized_counts.columns:
        result['expression'] = state.normalized_counts[gene].to_dict()

    return result


def handle_generate_volcano_plot(params: Dict) -> Dict:
    """Generate matplotlib volcano plot as SVG/PNG."""
    if state.deseq_results is None:
        raise ValueError("Run DESeq2 analysis first")

    lfc_thresh = params.get('lfc_threshold', state.analysis_params['lfc_threshold'])
    padj_thresh = params.get('padj_threshold', state.analysis_params['padj_threshold'])
    highlight_genes = params.get('highlight_genes', state.selected_genes)
    output_format = params.get('format', 'svg')

    return state.plot_generator.create_volcano_plot(
        results=state.deseq_results,
        lfc_threshold=lfc_thresh,
        padj_threshold=padj_thresh,
        highlight_genes=highlight_genes,
        output_format=output_format
    )


def handle_generate_boxplot(params: Dict) -> Dict:
    """Generate matplotlib boxplot as SVG/PNG."""
    if state.normalized_counts is None:
        raise ValueError("Run DESeq2 analysis first")

    genes = params.get('genes', state.selected_genes)
    if not genes:
        return {'error': 'No genes selected'}

    group_column = params.get('group_column', state.analysis_params['condition_column'])
    normalization = params.get('normalization', 'log2')
    output_format = params.get('format', 'svg')

    return state.plot_generator.create_boxplot(
        normalized_counts=state.normalized_counts,
        metadata=state.metadata,
        genes=genes,
        group_column=group_column,
        results=state.deseq_results,
        normalization=normalization,
        output_format=output_format
    )


def handle_generate_heatmap_plot(params: Dict) -> Dict:
    """Generate matplotlib heatmap as SVG/PNG."""
    if state.normalized_counts is None:
        raise ValueError("Run DESeq2 analysis first")

    genes = params.get('genes', None)
    top_n = params.get('top_n', 50)
    gene_selection = params.get('gene_selection', 'variance')
    scale = params.get('scale', 'row')
    cluster_cols = params.get('cluster_cols', True)
    cluster_rows = params.get('cluster_rows', True)
    show_dendrogram = params.get('show_dendrogram', True)
    annotation_columns = params.get('annotation_columns', ['condition'])
    output_format = params.get('format', 'svg')

    return state.plot_generator.create_heatmap(
        normalized_counts=state.normalized_counts,
        metadata=state.metadata,
        deseq_results=state.deseq_results,
        genes=genes,
        top_n=top_n,
        gene_selection=gene_selection,
        annotation_columns=annotation_columns,
        scale=scale,
        cluster_cols=cluster_cols,
        cluster_rows=cluster_rows,
        show_dendrogram=show_dendrogram,
        output_format=output_format
    )


def handle_generate_gsea_plot(params: Dict) -> Dict:
    """Generate matplotlib GSEA plot as SVG/PNG."""
    if state.gsea_results is None:
        raise ValueError("Run GSEA analysis first")

    top_n = params.get('top_n', 15)
    output_format = params.get('format', 'svg')

    return state.plot_generator.create_gsea_dotplot(
        gsea_results=state.gsea_results,
        top_n=top_n,
        output_format=output_format
    )


def handle_generate_pca_plot(params: Dict) -> Dict:
    """Generate matplotlib PCA plot as SVG/PNG."""
    if state.normalized_counts is None:
        raise ValueError("Run DESeq2 analysis first")

    color_by = params.get('color_by', 'condition')
    normalization = params.get('normalization', 'log2')
    center = params.get('center', True)
    scale = params.get('scale', True)
    show_labels = params.get('show_labels', True)
    genes = params.get('genes', None)
    output_format = params.get('format', 'svg')

    # Filter to selected genes if provided
    counts = state.normalized_counts
    if genes:
        valid_genes = [g for g in genes if g in counts.columns]
        if valid_genes:
            counts = counts[valid_genes]

    return state.plot_generator.create_pca_plot(
        normalized_counts=counts,
        metadata=state.metadata,
        color_by=color_by,
        normalization=normalization,
        center=center,
        scale=scale,
        show_labels=show_labels,
        output_format=output_format
    )


def handle_get_metadata_columns(params: Dict) -> Dict:
    """Get available metadata columns for dropdown."""
    if state.metadata is None:
        return {'columns': []}

    return {
        'columns': list(state.metadata.columns)
    }


# ---------------------------------------------------------------------------
# Kallisto handlers
# ---------------------------------------------------------------------------

def handle_kallisto_check_index(params: Dict) -> Dict:
    """Check if the prebuilt human index is available."""
    return {
        'has_index': state.kallisto_runner.has_prebuilt_index(),
        'index_path': state.kallisto_runner.get_prebuilt_index_path(),
    }


def handle_kallisto_download_index(params: Dict) -> Dict:
    """Download the prebuilt human transcriptome index."""
    idx_path = state.kallisto_runner.download_human_index()
    state.kallisto_index_path = idx_path
    state.kallisto_runner.load_t2g()
    return {'index_path': idx_path}


def handle_kallisto_load_example(params: Dict) -> Dict:
    """Generate example FASTA + FASTQ files and return their paths."""
    import tempfile
    output_dir = os.path.join(tempfile.gettempdir(), 'kallisto_example')
    result = create_example_fastqs(output_dir)

    return {
        'fasta_path': result['fasta_path'],
        'metadata_path': result['metadata_path'],
        'samples': result['samples'],
    }


def handle_kallisto_build_index(params: Dict) -> Dict:
    """Build a kallisto index from a FASTA file."""
    fasta_path = params.get('fasta_path')
    if not fasta_path:
        raise ValueError('fasta_path is required')

    import tempfile
    output_path = params.get('output_path')
    if not output_path:
        output_path = os.path.join(tempfile.gettempdir(), 'kallisto_custom.idx')

    idx = state.kallisto_runner.build_index(fasta_path, output_path)
    state.kallisto_index_path = idx
    return {'index_path': idx}


def handle_kallisto_quant_sample(params: Dict) -> Dict:
    """Quantify one sample with kallisto."""
    import tempfile

    index_path = params.get('index_path') or state.kallisto_index_path
    if not index_path:
        raise ValueError('No index available. Build or download an index first.')

    sample_name = params['sample_name']
    fastq_files = params['fastq_files']
    single_end = params.get('single_end', True)
    fragment_length = params.get('fragment_length', 200)
    sd = params.get('sd', 20)

    output_dir = os.path.join(tempfile.gettempdir(), 'kallisto_quant', sample_name)

    state.kallisto_runner.quant_sample(
        index_path=index_path,
        output_dir=output_dir,
        fastq_files=fastq_files,
        single_end=single_end,
        fragment_length=fragment_length,
        sd=sd,
    )

    state.kallisto_sample_results[sample_name] = output_dir

    # Return the abundance preview
    ab = state.kallisto_runner.parse_abundance(output_dir)

    # Read run_info.json for QC stats
    import json
    run_info = {}
    run_info_path = os.path.join(output_dir, 'run_info.json')
    if os.path.exists(run_info_path):
        with open(run_info_path) as f:
            run_info = json.load(f)

    return {
        'sample_name': sample_name,
        'output_dir': output_dir,
        'n_targets': len(ab),
        'total_counts': float(ab['est_counts'].sum()),
        'run_info': run_info,
    }


def handle_kallisto_combine_counts(params: Dict) -> Dict:
    """Combine per-sample abundance files into a counts matrix + metadata."""
    sample_dirs = params.get('sample_dirs') or state.kallisto_sample_results
    metadata_path = params.get('metadata_path')
    use_gene_names = params.get('use_gene_names', True)

    if not sample_dirs:
        raise ValueError('No quantified samples available')

    counts_df = state.kallisto_runner.combine_counts(sample_dirs, use_gene_names=use_gene_names)

    # Load or create metadata
    if metadata_path:
        import pandas as pd
        meta_df = pd.read_csv(metadata_path, index_col=0)
    else:
        # Auto-generate minimal metadata from sample names
        import pandas as pd
        meta_df = pd.DataFrame({'sample': list(sample_dirs.keys())})
        meta_df = meta_df.set_index('sample')
        meta_df['condition'] = 'unknown'

    # Store in state (same as load_counts + load_metadata)
    state.counts = counts_df
    state.metadata = meta_df

    return {
        'counts_shape': list(counts_df.shape),
        'counts_preview': counts_df.iloc[:10, :5].to_dict(),
        'counts_columns': list(counts_df.columns),
        'counts_index': list(counts_df.index[:20]),
        'metadata_shape': list(meta_df.shape),
        'metadata_preview': meta_df.head(10).to_dict(),
        'metadata_columns': list(meta_df.columns),
        'metadata_index': list(meta_df.index),
    }


def handle_kallisto_get_qc_stats(params: Dict) -> Dict:
    """Return QC stats (from run_info.json) for all quantified samples."""
    import json as _json
    sample_dirs = state.kallisto_sample_results
    if not sample_dirs:
        return {'samples': []}

    samples = []
    for name, output_dir in sample_dirs.items():
        info = {}
        run_info_path = os.path.join(output_dir, 'run_info.json')
        if os.path.exists(run_info_path):
            with open(run_info_path) as f:
                info = _json.load(f)

        ab = state.kallisto_runner.parse_abundance(output_dir)
        samples.append({
            'name': name,
            'n_processed': info.get('n_processed', 0),
            'n_pseudoaligned': info.get('n_pseudoaligned', 0),
            'n_unique': info.get('n_unique', 0),
            'p_pseudoaligned': info.get('p_pseudoaligned', 0),
            'p_unique': info.get('p_unique', 0),
            'kallisto_version': info.get('kallisto_version', ''),
            'index_version': info.get('index_version', 0),
            'n_targets': info.get('n_targets', len(ab)),
            'n_bootstraps': info.get('n_bootstraps', 0),
            'total_counts': float(ab['est_counts'].sum()),
            'mean_tpm': float(ab['tpm'].mean()) if 'tpm' in ab.columns else 0,
        })

    return {'samples': samples}


def handle_generate_kallisto_qc_plot(params: Dict) -> Dict:
    """Generate a matplotlib QC summary plot for all quantified samples."""
    import json as _json
    import io
    import base64

    sample_dirs = state.kallisto_sample_results
    if not sample_dirs:
        raise ValueError('No quantified samples available')

    # Gather stats
    names = []
    n_processed_list = []
    n_aligned_list = []
    pct_aligned_list = []
    total_counts_list = []

    for name, output_dir in sample_dirs.items():
        info = {}
        run_info_path = os.path.join(output_dir, 'run_info.json')
        if os.path.exists(run_info_path):
            with open(run_info_path) as f:
                info = _json.load(f)

        ab = state.kallisto_runner.parse_abundance(output_dir)
        names.append(name)
        n_proc = info.get('n_processed', 0)
        n_aligned = info.get('n_pseudoaligned', 0)
        n_processed_list.append(n_proc)
        n_aligned_list.append(n_aligned)
        pct_aligned_list.append(info.get('p_pseudoaligned', 0))
        total_counts_list.append(float(ab['est_counts'].sum()))

    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import numpy as np

    fig, axes = plt.subplots(3, 1, figsize=(10, 12))
    fig.suptitle('Kallisto QC Summary', fontsize=14, fontweight='bold')
    x = np.arange(len(names))
    bar_width = 0.6

    # Panel 1: Reads processed vs pseudoaligned (stacked)
    ax = axes[0]
    unaligned = [p - a for p, a in zip(n_processed_list, n_aligned_list)]
    ax.bar(x, n_aligned_list, bar_width, label='Pseudoaligned', color='#2196F3')
    ax.bar(x, unaligned, bar_width, bottom=n_aligned_list, label='Not aligned', color='#BBDEFB')
    ax.set_xticks(x)
    ax.set_xticklabels(names, rotation=45, ha='right', fontsize=8)
    ax.set_ylabel('Reads')
    ax.set_title('Reads Processed')
    ax.legend(fontsize=8)

    # Panel 2: Pseudoalignment rate (%)
    ax = axes[1]
    colors = ['#4CAF50' if p >= 50 else '#FF9800' if p >= 20 else '#F44336' for p in pct_aligned_list]
    ax.bar(x, pct_aligned_list, bar_width, color=colors)
    ax.set_xticks(x)
    ax.set_xticklabels(names, rotation=45, ha='right', fontsize=8)
    ax.set_ylabel('Pseudoalignment Rate (%)')
    ax.set_title('Mapping Rate')
    ax.set_ylim(0, 105)
    ax.axhline(y=50, color='gray', linestyle='--', alpha=0.5, linewidth=0.8)

    # Panel 3: Total estimated counts per sample
    ax = axes[2]
    ax.bar(x, total_counts_list, bar_width, color='#9C27B0')
    ax.set_xticks(x)
    ax.set_xticklabels(names, rotation=45, ha='right', fontsize=8)
    ax.set_ylabel('Estimated Counts')
    ax.set_title('Total Estimated Counts')

    plt.tight_layout()

    fmt = params.get('format', 'svg')
    buf = io.BytesIO()
    fig.savefig(buf, format=fmt, dpi=150, bbox_inches='tight')
    plt.close(fig)
    buf.seek(0)

    if fmt == 'svg':
        return {'plot': buf.read().decode('utf-8'), 'format': 'svg'}
    else:
        return {'plot': base64.b64encode(buf.read()).decode('utf-8'), 'format': fmt}


def handle_get_environment_info(params: Dict) -> Dict:
    """Get Python environment info for reproducibility export."""
    import importlib.metadata
    import platform as plat

    # Collect installed packages
    packages = []
    for dist in importlib.metadata.distributions():
        packages.append(f"{dist.metadata['Name']}=={dist.version}")
    packages.sort(key=lambda x: x.lower())

    # Build requirements.txt content
    requirements_txt = '\n'.join(packages)

    # Build conda-style yaml
    conda_yaml_lines = [
        'name: rnaseq_analysis',
        'channels:',
        '  - defaults',
        '  - conda-forge',
        '  - bioconda',
        'dependencies:',
        f'  - python={plat.python_version()}',
    ]
    # Add key packages as conda deps
    key_packages = [
        'numpy', 'pandas', 'scipy', 'matplotlib', 'seaborn',
        'scikit-learn', 'pydeseq2', 'gseapy', 'edgepy', 'plotly',
    ]
    for pkg_line in packages:
        pkg_name = pkg_line.split('==')[0].lower()
        if pkg_name in key_packages:
            conda_yaml_lines.append(f'  - {pkg_line.lower()}')
    conda_yaml_lines.append('  - pip:')
    for pkg_line in packages:
        pkg_name = pkg_line.split('==')[0].lower()
        if pkg_name not in key_packages and pkg_name != 'python':
            conda_yaml_lines.append(f'    - {pkg_line}')
    conda_yaml = '\n'.join(conda_yaml_lines)

    return {
        'python_version': plat.python_version(),
        'requirements_txt': requirements_txt,
        'conda_yaml': conda_yaml,
        'package_count': len(packages),
    }


def handle_init_config(params: Dict) -> Dict:
    """Receive app configuration from Electron (paths, platform, packaged flag)."""
    APP_CONFIG.update({
        'resources_path': params.get('resources_path'),
        'user_data_path': params.get('user_data_path'),
        'platform': params.get('platform', sys.platform),
        'is_packaged': params.get('is_packaged', False),
    })
    return {'status': 'ok'}


def handle_kallisto_check_platform(params: Dict) -> Dict:
    """Check if kallisto is available on this platform."""
    available = state.kallisto_runner.is_available
    return {
        'available': available,
        'platform': sys.platform,
        'message': (
            'Kallisto FASTQ processing is not available on Windows. '
            'Please use macOS or Linux for raw FASTQ analysis, '
            'or load a pre-computed counts matrix in Step 1.'
        ) if not available else None,
    }


# Command router
COMMANDS = {
    'init_config': handle_init_config,
    'kallisto_check_platform': handle_kallisto_check_platform,
    'load_example_data': handle_load_example_data,
    'load_counts': handle_load_counts,
    'load_metadata': handle_load_metadata,
    'set_params': handle_set_params,
    'run_deseq': handle_run_deseq,
    'get_volcano_data': handle_get_volcano_data,
    'get_boxplot_data': handle_get_boxplot_data,
    'get_heatmap_data': handle_get_heatmap_data,
    'run_gsea': handle_run_gsea,
    'set_selected_genes': handle_set_selected_genes,
    'generate_script': handle_generate_script,
    'get_gene_info': handle_get_gene_info,
    # Matplotlib plot generation
    'generate_volcano_plot': handle_generate_volcano_plot,
    'generate_boxplot': handle_generate_boxplot,
    'generate_heatmap_plot': handle_generate_heatmap_plot,
    'generate_gsea_plot': handle_generate_gsea_plot,
    'generate_pca_plot': handle_generate_pca_plot,
    'get_metadata_columns': handle_get_metadata_columns,
    # Kallisto FASTQ processing
    'list_project_files': handle_list_project_files,
    'set_active_tool': handle_set_active_tool,
    'run_edgepython': handle_run_edgepython,
    'kallisto_check_index': handle_kallisto_check_index,
    'kallisto_download_index': handle_kallisto_download_index,
    'kallisto_load_example': handle_kallisto_load_example,
    'kallisto_build_index': handle_kallisto_build_index,
    'kallisto_quant_sample': handle_kallisto_quant_sample,
    'kallisto_combine_counts': handle_kallisto_combine_counts,
    'kallisto_get_qc_stats': handle_kallisto_get_qc_stats,
    'generate_kallisto_qc_plot': handle_generate_kallisto_qc_plot,
    'get_environment_info': handle_get_environment_info,
}


def main():
    """Main event loop - read commands from stdin, execute, write results to stdout."""
    # Signal that Python is ready
    print(json.dumps({'type': 'ready'}), flush=True)

    for line in sys.stdin:
        line = line.strip()
        if not line:
            continue

        try:
            request = json.loads(line)
            request_id = request.get('id', 0)
            command = request.get('command')
            params = request.get('params', {})

            if command not in COMMANDS:
                send_response(request_id, error=f'Unknown command: {command}')
                continue

            result = COMMANDS[command](params)
            send_response(request_id, data=result)

        except json.JSONDecodeError as e:
            send_response(0, error=f'Invalid JSON: {str(e)}')
        except Exception as e:
            traceback.print_exc(file=sys.stderr)
            request_id = request.get('id', 0) if 'request' in dir() else 0
            send_response(request_id, error=str(e))


if __name__ == '__main__':
    main()
