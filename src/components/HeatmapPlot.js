/**
 * Heatmap Component
 * Clustered heatmap with metadata annotations (pheatmap-style)
 */

class HeatmapPlot {
    constructor(containerId) {
        this.containerId = containerId;
        this.container = document.getElementById(containerId);
        this.data = null;
        this.topN = 50;
        this.scale = 'row';

        // Annotation colors
        this.annotationColors = {
            condition: { control: '#3498db', treatment: '#e74c3c' },
            batch: { batch1: '#2ecc71', batch2: '#9b59b6' },
            sex: { M: '#3498db', F: '#e91e63' }
        };

        this.layout = {
            title: {
                text: 'Expression Heatmap',
                font: { size: 18, family: 'Arial, sans-serif', weight: 'bold' }
            },
            margin: { l: 100, r: 50, t: 80, b: 150 },
            paper_bgcolor: 'white',
            plot_bgcolor: 'white'
        };

        this.config = {
            responsive: true,
            displayModeBar: true,
            toImageButtonOptions: {
                format: 'png',
                filename: 'heatmap',
                height: 800,
                width: 1000,
                scale: 2
            }
        };
    }

    /**
     * Render heatmap
     */
    async render(options = {}) {
        if (options.topN !== undefined) this.topN = options.topN;
        if (options.scale !== undefined) this.scale = options.scale;

        try {
            // Get data from Python
            this.data = await python.getHeatmapData({
                topN: this.topN,
                scale: this.scale,
                annotationColumns: ['condition']
            });

            // Prepare traces
            const traces = this.prepareTraces();

            // Update layout
            const layout = this.prepareLayout();

            // Render plot
            await Plotly.newPlot(this.containerId, traces, layout, this.config);

        } catch (error) {
            console.error('Error rendering heatmap:', error);
            this.showError(error.message);
        }
    }

    /**
     * Prepare Plotly traces
     */
    prepareTraces() {
        const { genes, samples, values, annotations } = this.data;
        const traces = [];

        // Main heatmap trace
        const heatmapTrace = {
            type: 'heatmap',
            z: values,
            x: samples,
            y: genes,
            colorscale: 'RdBu',
            reversescale: true,
            zmid: 0,
            colorbar: {
                title: {
                    text: this.scale === 'row' ? 'Z-score' : 'Expression',
                    side: 'right'
                },
                thickness: 15,
                len: 0.7,
                y: 0.5,
                yanchor: 'middle'
            },
            hovertemplate:
                '<b>Gene:</b> %{y}<br>' +
                '<b>Sample:</b> %{x}<br>' +
                '<b>Value:</b> %{z:.2f}<br>' +
                '<extra></extra>'
        };
        traces.push(heatmapTrace);

        // Add annotation bar for condition if available
        if (annotations && annotations.condition) {
            const conditionColors = annotations.condition.map(c =>
                c === 'control' ? 0 : 1
            );

            traces.push({
                type: 'heatmap',
                z: [conditionColors],
                x: samples,
                y: ['Condition'],
                colorscale: [
                    [0, '#3498db'],
                    [1, '#e74c3c']
                ],
                showscale: false,
                hovertemplate:
                    '<b>Sample:</b> %{x}<br>' +
                    '<b>Condition:</b> ' + annotations.condition.map((c, i) =>
                        `${c}`
                    ).join('') +
                    '<extra></extra>',
                yaxis: 'y2'
            });
        }

        return traces;
    }

    /**
     * Prepare layout
     */
    prepareLayout() {
        const { genes, samples, annotations } = this.data;

        const layout = {
            ...this.layout,
            xaxis: {
                title: { text: 'Samples', font: { size: 14 } },
                tickangle: -45,
                tickfont: { size: 10 },
                side: 'bottom'
            },
            yaxis: {
                title: { text: 'Genes', font: { size: 14 } },
                tickfont: { size: 8 },
                autorange: 'reversed'
            },
            annotations: []
        };

        // Add title annotation for scale
        const scaleText = this.scale === 'row' ? '(Row Z-score)' :
                         this.scale === 'column' ? '(Column Z-score)' : '(Raw values)';
        layout.title.text = `Expression Heatmap ${scaleText}`;

        // Add second y-axis for annotation bar
        if (annotations && Object.keys(annotations).length > 0) {
            layout.yaxis2 = {
                domain: [0.95, 1],
                tickfont: { size: 10 },
                anchor: 'x'
            };
            layout.yaxis.domain = [0, 0.93];

            // Add legend for annotation colors
            layout.annotations.push({
                x: 1.15,
                y: 1,
                xref: 'paper',
                yref: 'paper',
                text: '<b>Condition</b><br>■ Control<br>■ Treatment',
                showarrow: false,
                font: { size: 10 },
                align: 'left'
            });
        }

        return layout;
    }

    /**
     * Update settings and re-render
     */
    async update(options = {}) {
        if (options.topN !== undefined) this.topN = options.topN;
        if (options.scale !== undefined) this.scale = options.scale;
        await this.render(options);
    }

    /**
     * Show placeholder
     */
    showPlaceholder() {
        this.container.innerHTML = `
            <div style="display: flex; flex-direction: column; align-items: center; justify-content: center; height: 100%; color: #6c757d; gap: 16px;">
                <svg width="64" height="64" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="1.5">
                    <rect x="3" y="3" width="18" height="18" rx="2"/>
                    <line x1="3" y1="9" x2="21" y2="9"/>
                    <line x1="3" y1="15" x2="21" y2="15"/>
                    <line x1="9" y1="3" x2="9" y2="21"/>
                    <line x1="15" y1="3" x2="15" y2="21"/>
                </svg>
                <p>Run DESeq2 analysis to generate heatmap</p>
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
window.HeatmapPlot = HeatmapPlot;
