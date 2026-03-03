/**
 * BoxPlot Component
 * Grouped boxplots showing gene expression by condition
 * Includes significance stars from DESeq2 results
 */

class BoxPlot {
    constructor(containerId) {
        this.containerId = containerId;
        this.container = document.getElementById(containerId);
        this.data = null;
        this.normalization = 'log2';

        // Color palette for conditions - muted pastels
        this.colors = {
            control: '#7BA3C9',
            treatment: '#C97B7B'
        };

        this.layout = {
            title: {
                text: 'Gene Expression by Condition',
                font: { size: 18, family: 'Arial, sans-serif', weight: 'bold' }
            },
            xaxis: {
                title: { text: 'Gene', font: { size: 14 } },
                tickangle: -45
            },
            yaxis: {
                title: { text: 'Log<sub>2</sub>(Normalized Counts + 1)', font: { size: 14 } },
                gridcolor: '#f0f0f0'
            },
            boxmode: 'group',
            showlegend: true,
            legend: {
                x: 1,
                xanchor: 'right',
                y: 1,
                bgcolor: 'rgba(255,255,255,0.9)',
                bordercolor: '#ccc',
                borderwidth: 1
            },
            margin: { l: 60, r: 40, t: 80, b: 100 },
            paper_bgcolor: 'white',
            plot_bgcolor: 'white',
            annotations: []
        };

        this.config = {
            responsive: true,
            displayModeBar: true,
            toImageButtonOptions: {
                format: 'png',
                filename: 'boxplot',
                height: 600,
                width: 1000,
                scale: 2
            }
        };
    }

    /**
     * Render boxplot for selected genes
     */
    async render(options = {}) {
        if (options.normalization) this.normalization = options.normalization;

        const selectedGenes = store.get('selectedGenes') || [];

        if (selectedGenes.length === 0) {
            this.showPlaceholder();
            return;
        }

        try {
            // Get data from Python
            this.data = await python.getBoxplotData(selectedGenes, {
                normalization: this.normalization
            });

            if (this.data.error) {
                this.showError(this.data.error);
                return;
            }

            // Prepare traces
            const traces = this.prepareTraces();

            // Update layout
            const layout = this.prepareLayout();

            // Render plot
            await Plotly.newPlot(this.containerId, traces, layout, this.config);

        } catch (error) {
            console.error('Error rendering boxplot:', error);
            this.showError(error.message);
        }
    }

    /**
     * Prepare Plotly traces
     */
    prepareTraces() {
        const { genes, expression, samples, groups, unique_groups, significance } = this.data;
        const traces = [];

        // Create a trace for each group
        unique_groups.forEach((group, groupIdx) => {
            // Get indices for this group
            const groupIndices = samples.map((_, i) => groups[i] === group ? i : -1).filter(i => i >= 0);

            genes.forEach((gene, geneIdx) => {
                const geneExpression = expression[gene];
                const groupValues = groupIndices.map(i => geneExpression[i]);

                traces.push({
                    type: 'box',
                    name: group,
                    x: groupValues.map(() => gene),
                    y: groupValues,
                    marker: {
                        color: this.colors[group] || this.getColor(groupIdx),
                        opacity: 0.7
                    },
                    boxpoints: 'all',
                    jitter: 0.3,
                    pointpos: 0,
                    legendgroup: group,
                    showlegend: geneIdx === 0, // Only show legend for first gene
                    hovertemplate:
                        `<b>${group}</b><br>` +
                        `Gene: ${gene}<br>` +
                        `Expression: %{y:.2f}<br>` +
                        '<extra></extra>'
                });
            });
        });

        return traces;
    }

    /**
     * Prepare layout with significance annotations
     */
    prepareLayout() {
        const layout = { ...this.layout };
        const { genes, significance } = this.data;

        // Update y-axis title based on normalization
        if (this.normalization === 'log2') {
            layout.yaxis.title.text = 'Log<sub>2</sub>(Normalized Counts + 1)';
        } else {
            layout.yaxis.title.text = 'Normalized Counts';
        }

        // Add significance stars as annotations
        layout.annotations = genes.map((gene, i) => {
            const sig = significance[gene];
            if (sig && sig !== 'ns' && sig !== 'na') {
                return {
                    x: gene,
                    y: 1.05,
                    xref: 'x',
                    yref: 'paper',
                    text: sig,
                    showarrow: false,
                    font: { size: 16, color: '#333', weight: 'bold' }
                };
            }
            return null;
        }).filter(a => a !== null);

        return layout;
    }

    /**
     * Get color for a group index - muted pastels
     */
    getColor(index) {
        const palette = ['#7BA3C9', '#C97B7B', '#8FB996', '#B39DDB', '#E8C87D', '#D4A5A5'];
        return palette[index % palette.length];
    }

    /**
     * Update normalization and re-render
     */
    async updateNormalization(normalization) {
        this.normalization = normalization;
        await this.render();
    }

    /**
     * Show placeholder
     */
    showPlaceholder() {
        this.container.innerHTML = `
            <div style="display: flex; flex-direction: column; align-items: center; justify-content: center; height: 100%; color: #6c757d; gap: 16px;">
                <svg width="64" height="64" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="1.5">
                    <rect x="3" y="8" width="4" height="12" rx="1"/>
                    <rect x="10" y="5" width="4" height="15" rx="1"/>
                    <rect x="17" y="10" width="4" height="10" rx="1"/>
                </svg>
                <p>Click genes on the volcano plot to add them to the boxplot</p>
            </div>
        `;
    }

    /**
     * Show error message
     */
    showError(message) {
        this.container.innerHTML = `
            <div style="display: flex; align-items: center; justify-content: center; height: 100%; color: #e74c3c;">
                <p>Error: ${message}</p>
            </div>
        `;
    }
}

// Export for use in other files
window.BoxPlot = BoxPlot;
