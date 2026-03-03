/**
 * GSEA Plot Component
 * Enrichment dot plot for pathway analysis results
 */

class GSEAPlot {
    constructor(containerId, tableContainerId) {
        this.containerId = containerId;
        this.tableContainerId = tableContainerId;
        this.container = document.getElementById(containerId);
        this.tableContainer = document.getElementById(tableContainerId);
        this.data = null;
        this.topN = 15;

        this.layout = {
            title: {
                text: 'Pathway Enrichment Analysis',
                font: { size: 18, family: 'Arial, sans-serif', weight: 'bold' }
            },
            xaxis: {
                title: { text: 'Normalized Enrichment Score (NES)', font: { size: 14 } },
                zeroline: true,
                zerolinecolor: '#ccc',
                zerolinewidth: 2,
                gridcolor: '#f0f0f0'
            },
            yaxis: {
                title: '',
                automargin: true
            },
            margin: { l: 300, r: 100, t: 60, b: 60 },
            paper_bgcolor: 'white',
            plot_bgcolor: 'white',
            showlegend: true,
            legend: {
                title: { text: '-Log10(FDR)' },
                x: 1.02,
                y: 0.5
            }
        };

        this.config = {
            responsive: true,
            displayModeBar: true,
            toImageButtonOptions: {
                format: 'png',
                filename: 'gsea_dotplot',
                height: 800,
                width: 1200,
                scale: 2
            }
        };
    }

    /**
     * Render GSEA results
     */
    async render(results = null) {
        if (results) {
            this.data = results;
        } else {
            this.data = store.get('gseaResults');
        }

        if (!this.data || this.data.length === 0) {
            this.showPlaceholder();
            return;
        }

        try {
            // Prepare traces
            const traces = this.prepareTraces();

            // Update layout
            const layout = this.prepareLayout();

            // Render plot
            await Plotly.newPlot(this.containerId, traces, layout, this.config);

            // Render table
            this.renderTable();

        } catch (error) {
            console.error('Error rendering GSEA plot:', error);
            this.showError(error.message);
        }
    }

    /**
     * Prepare Plotly traces
     */
    prepareTraces() {
        // Get top N results
        const topResults = this.data.slice(0, this.topN);

        // Extract data
        const terms = topResults.map(r => this.truncateTerm(r.Term || r.term, 50));
        const nes = topResults.map(r => r.NES || r.nes || 0);
        const pvals = topResults.map(r => r['FDR q-val'] || r.fdr || r['Adjusted P-value'] || 0.05);

        // Calculate -log10(FDR)
        const negLogFdr = pvals.map(p => -Math.log10(Math.max(p, 1e-10)));

        // Size based on -log10(FDR)
        const maxLogFdr = Math.max(...negLogFdr);
        const sizes = negLogFdr.map(v => 10 + (v / maxLogFdr) * 20);

        // Colors based on NES (positive = red, negative = blue)
        const colors = nes.map(n => n > 0 ? '#e74c3c' : '#3498db');

        return [{
            type: 'scatter',
            mode: 'markers',
            x: nes,
            y: terms,
            marker: {
                size: sizes,
                color: colors,
                opacity: 0.7,
                line: { color: 'white', width: 1 }
            },
            text: topResults.map((r, i) =>
                `<b>${r.Term || r.term}</b><br>` +
                `NES: ${nes[i].toFixed(2)}<br>` +
                `FDR: ${pvals[i].toExponential(2)}`
            ),
            hoverinfo: 'text'
        }];
    }

    /**
     * Prepare layout
     */
    prepareLayout() {
        const topResults = this.data.slice(0, this.topN);
        const terms = topResults.map(r => this.truncateTerm(r.Term || r.term, 50));

        return {
            ...this.layout,
            yaxis: {
                ...this.layout.yaxis,
                tickvals: terms,
                ticktext: terms,
                autorange: 'reversed'
            },
            height: Math.max(400, terms.length * 30 + 120)
        };
    }

    /**
     * Truncate long term names
     */
    truncateTerm(term, maxLength) {
        if (term.length <= maxLength) return term;
        return term.substring(0, maxLength - 3) + '...';
    }

    /**
     * Render results table
     */
    renderTable() {
        if (!this.tableContainer || !this.data) return;

        const headers = ['Term', 'NES', 'P-value', 'FDR'];
        const rows = this.data.slice(0, 50).map(r => [
            this.truncateTerm(r.Term || r.term || '', 60),
            (r.NES || r.nes || 0).toFixed(3),
            (r['NOM p-val'] || r.pvalue || r['P-value'] || 0).toExponential(2),
            (r['FDR q-val'] || r.fdr || r['Adjusted P-value'] || 0).toExponential(2)
        ]);

        let html = '<table class="data-table"><thead><tr>';
        headers.forEach(h => html += `<th>${h}</th>`);
        html += '</tr></thead><tbody>';

        rows.forEach(row => {
            html += '<tr>';
            row.forEach((cell, i) => {
                const style = i === 1 ?
                    (parseFloat(cell) > 0 ? 'color: #e74c3c' : 'color: #3498db') : '';
                html += `<td style="${style}">${cell}</td>`;
            });
            html += '</tr>';
        });

        html += '</tbody></table>';
        this.tableContainer.innerHTML = html;
    }

    /**
     * Show placeholder
     */
    showPlaceholder() {
        this.container.innerHTML = `
            <div style="display: flex; flex-direction: column; align-items: center; justify-content: center; height: 100%; color: #6c757d; gap: 16px;">
                <svg width="64" height="64" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="1.5">
                    <circle cx="12" cy="12" r="10"/>
                    <path d="M12 6v6l4 2"/>
                </svg>
                <p>Run pathway analysis to see results</p>
            </div>
        `;

        if (this.tableContainer) {
            this.tableContainer.innerHTML = '';
        }
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
window.GSEAPlot = GSEAPlot;
