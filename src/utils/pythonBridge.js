/**
 * Python Bridge
 * Handles communication between JavaScript and Python backend
 */

class PythonBridge {
    constructor() {
        this.ready = false;
        this.readyCallbacks = [];

        // Setup event listeners
        if (window.api) {
            window.api.onPythonReady(() => {
                this.ready = true;
                console.log('Python backend ready');
                this.readyCallbacks.forEach(cb => cb());
                this.readyCallbacks = [];
            });

            window.api.onPythonError((message) => {
                console.error('Python error:', message);
                store.set('error', message);
            });
        }
    }

    /**
     * Wait for Python to be ready
     */
    onReady(callback) {
        if (this.ready) {
            callback();
        } else {
            this.readyCallbacks.push(callback);
        }
    }

    /**
     * Run a Python command
     */
    async run(command, params = {}) {
        if (!this.ready) {
            throw new Error('Python backend not ready');
        }

        try {
            const result = await window.api.runPython(command, params);
            return result;
        } catch (error) {
            console.error(`Python command "${command}" failed:`, error);
            throw error;
        }
    }

    /**
     * Load example data
     */
    async loadExampleData(options = {}) {
        store.set('isLoading', true);
        try {
            const result = await this.run('load_example_data', {
                n_samples: options.nSamples || 20,
                n_genes: options.nGenes || 5000,
                n_de_genes: options.nDeGenes || 500
            });

            store.update({
                counts: result.counts_preview,
                metadata: result.metadata_preview,
                countsShape: result.counts_shape,
                metadataShape: result.metadata_shape,
                metadataColumns: result.metadata_columns,
                countsColumns: result.counts_columns,
                countsIndex: result.counts_index,
                metadataIndex: result.metadata_index
            });

            return result;
        } finally {
            store.set('isLoading', false);
        }
    }

    /**
     * Load counts file
     */
    async loadCounts(filePath) {
        store.set('isLoading', true);
        try {
            let result;
            if (filePath && filePath.startsWith('__web_upload__:')) {
                result = await window.api.uploadFile('load_counts');
            } else {
                result = await this.run('load_counts', { file_path: filePath });
            }
            store.update({
                counts: result.counts_preview,
                countsShape: result.counts_shape,
                countsColumns: result.counts_columns,
                countsIndex: result.counts_index
            });
            return result;
        } finally {
            store.set('isLoading', false);
        }
    }

    /**
     * Load metadata file
     */
    async loadMetadata(filePath) {
        store.set('isLoading', true);
        try {
            let result;
            if (filePath && filePath.startsWith('__web_upload__:')) {
                result = await window.api.uploadFile('load_metadata');
            } else {
                result = await this.run('load_metadata', { file_path: filePath });
            }
            store.update({
                metadata: result.metadata_preview,
                metadataShape: result.metadata_shape,
                metadataColumns: result.metadata_columns,
                metadataIndex: result.metadata_index
            });
            return result;
        } finally {
            store.set('isLoading', false);
        }
    }

    /**
     * Set active differential testing tool
     */
    async setActiveTool(tool) {
        const result = await this.run('set_active_tool', { tool });
        store.set('activeTool', tool);
        return result;
    }

    /**
     * Run DESeq2 analysis
     */
    async runDeseq(params = {}) {
        store.set('isLoading', true);
        try {
            const result = await this.run('run_deseq', params);
            store.update({
                deseqResults: result.results,
                deseqSummary: result.summary,
                activeTool: 'pydeseq2',
            });
            return result;
        } finally {
            store.set('isLoading', false);
        }
    }

    /**
     * Run edgePython analysis
     */
    async runEdgePython(params = {}) {
        store.set('isLoading', true);
        try {
            const result = await this.run('run_edgepython', params);
            store.update({
                deseqResults: result.results,
                deseqSummary: result.summary,
                activeTool: 'edgepython',
            });
            return result;
        } finally {
            store.set('isLoading', false);
        }
    }

    /**
     * Get volcano plot data
     */
    async getVolcanoData(params = {}) {
        return await this.run('get_volcano_data', params);
    }

    /**
     * Get boxplot data
     */
    async getBoxplotData(genes, options = {}) {
        return await this.run('get_boxplot_data', {
            genes,
            normalization: options.normalization || 'log2',
            group_column: options.groupColumn || store.get('params.conditionColumn')
        });
    }

    /**
     * Get heatmap data
     */
    async getHeatmapData(options = {}) {
        return await this.run('get_heatmap_data', {
            top_n: options.topN || 50,
            genes: options.genes,
            scale: options.scale || 'row',
            annotation_columns: options.annotationColumns || ['condition']
        });
    }

    /**
     * Run GSEA analysis
     */
    async runGsea(options = {}) {
        store.set('isLoading', true);
        try {
            const result = await this.run('run_gsea', {
                gene_sets: options.geneSets || 'KEGG_2021_Human',
                method: options.method || 'prerank'
            });
            store.set('gseaResults', result.results);
            return result;
        } finally {
            store.set('isLoading', false);
        }
    }

    /**
     * Generate Python script
     */
    async generateScript(includeSections) {
        const result = await this.run('generate_script', {
            include_sections: {
                data_loading: includeSections.dataLoading,
                deseq2: includeSections.deseq2,
                volcano: includeSections.volcano,
                boxplot: includeSections.boxplot,
                heatmap: includeSections.heatmap,
                gsea: includeSections.gsea
            }
        });
        return result.script;
    }

    /**
     * Generate matplotlib volcano plot (returns SVG or PNG)
     */
    async generateVolcanoPlot(options = {}) {
        return await this.run('generate_volcano_plot', {
            lfc_threshold: options.lfcThreshold || store.get('params.lfcThreshold'),
            padj_threshold: options.padjThreshold || store.get('params.padjThreshold'),
            highlight_genes: options.highlightGenes || store.get('selectedGenes'),
            format: options.format || 'svg'
        });
    }

    /**
     * Generate matplotlib boxplot (returns SVG or PNG)
     */
    async generateBoxplot(genes, options = {}) {
        return await this.run('generate_boxplot', {
            genes: genes || store.get('selectedGenes'),
            group_column: options.groupColumn || store.get('params.conditionColumn'),
            normalization: options.normalization || 'log2',
            format: options.format || 'svg'
        });
    }

    /**
     * Generate matplotlib heatmap (returns SVG or PNG)
     */
    async generateHeatmapPlot(options = {}) {
        return await this.run('generate_heatmap_plot', {
            genes: options.genes,
            top_n: options.topN || 50,
            gene_selection: options.geneSelection || 'variance',
            scale: options.scale || 'row',
            cluster_cols: options.clusterCols !== false,
            cluster_rows: options.clusterRows !== false,
            show_dendrogram: options.showDendrogram !== false,
            annotation_columns: options.annotationColumns || ['condition'],
            format: options.format || 'svg'
        });
    }

    /**
     * Generate matplotlib GSEA plot (returns SVG or PNG)
     */
    async generateGseaPlot(options = {}) {
        return await this.run('generate_gsea_plot', {
            top_n: options.topN || 15,
            format: options.format || 'svg'
        });
    }

    /**
     * Generate matplotlib PCA plot (returns SVG or PNG)
     */
    async generatePcaPlot(options = {}) {
        const params = {
            color_by: options.colorBy || 'condition',
            normalization: options.normalization || 'log2',
            center: options.center !== false,
            scale: options.scale !== false,
            show_labels: options.showLabels !== false,
            show_legend: options.showLegend !== false,
            format: options.format || 'svg'
        };
        if (options.genes && options.genes.length > 0) {
            params.genes = options.genes;
        }
        return await this.run('generate_pca_plot', params);
    }

    // ------------------------------------------------------------------
    // Kallisto / FASTQ processing
    // ------------------------------------------------------------------

    async kallistoCheckIndex() {
        return await this.run('kallisto_check_index', {});
    }

    async kallistoDownloadIndex() {
        return await this.run('kallisto_download_index', {});
    }

    async kallistoLoadExample() {
        return await this.run('kallisto_load_example', {});
    }

    async kallistoBuildIndex(fastaPath) {
        return await this.run('kallisto_build_index', { fasta_path: fastaPath });
    }

    async kallistoQuantSample(sampleName, fastqFiles, options = {}) {
        return await this.run('kallisto_quant_sample', {
            sample_name: sampleName,
            fastq_files: fastqFiles,
            single_end: options.singleEnd !== false,
            fragment_length: options.fragmentLength || 200,
            sd: options.sd || 20,
            index_path: options.indexPath || null,
        });
    }

    async kallistoCombineCounts(options = {}) {
        const result = await this.run('kallisto_combine_counts', {
            metadata_path: options.metadataPath || null,
            use_gene_names: options.useGeneNames !== false,
        });
        // Store in same format as loadCounts + loadMetadata
        store.update({
            counts: result.counts_preview,
            countsShape: result.counts_shape,
            countsColumns: result.counts_columns,
            countsIndex: result.counts_index,
            metadata: result.metadata_preview,
            metadataShape: result.metadata_shape,
            metadataColumns: result.metadata_columns,
            metadataIndex: result.metadata_index,
        });
        return result;
    }

    async kallistoGetQcStats() {
        return await this.run('kallisto_get_qc_stats', {});
    }

    async generateKallistoQcPlot(options = {}) {
        return await this.run('generate_kallisto_qc_plot', {
            format: options.format || 'svg',
        });
    }

    /**
     * Get Python environment info for reproducibility
     */
    async getEnvironmentInfo() {
        return await this.run('get_environment_info', {});
    }

    /**
     * Open file dialog
     */
    async openFileDialog(options = {}) {
        return await window.api.openFileDialog(options);
    }

    /**
     * Save file dialog
     */
    async saveFileDialog(options = {}) {
        return await window.api.saveFileDialog(options);
    }

    /**
     * Save file
     */
    async saveFile(filePath, content) {
        return await window.api.saveFile(filePath, content);
    }
}

// Create global bridge instance
window.python = new PythonBridge();
