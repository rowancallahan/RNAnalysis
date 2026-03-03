/**
 * Volcano Plot Component
 * Interactive scatter plot for differential expression visualization
 * Supports click-to-select genes for downstream analysis
 * Features: draggable threshold lines, static rendering for low-LFC genes
 */

class VolcanoPlot {
    constructor(containerId) {
        this.containerId = containerId;
        this.container = document.getElementById(containerId);
        this.data = null;
        this.lfcThreshold = 1.0;
        this.padjThreshold = 0.05;
        this._isDragging = false;
        this._syncPending = false;

        // Plotly layout configuration
        this.layout = {
            title: {
                text: 'Volcano Plot',
                font: { size: 18, family: 'Arial, sans-serif', weight: 'bold' }
            },
            xaxis: {
                title: { text: 'Log<sub>2</sub> Fold Change', font: { size: 14 } },
                zeroline: true,
                zerolinecolor: '#ccc',
                zerolinewidth: 1,
                gridcolor: '#f0f0f0'
            },
            yaxis: {
                title: { text: '-Log<sub>10</sub> Adjusted P-value', font: { size: 14 } },
                gridcolor: '#f0f0f0'
            },
            hovermode: 'closest',
            showlegend: true,
            legend: {
                x: 1,
                xanchor: 'right',
                y: 1,
                bgcolor: 'rgba(255,255,255,0.9)',
                bordercolor: '#ccc',
                borderwidth: 1
            },
            margin: { l: 60, r: 40, t: 60, b: 60 },
            paper_bgcolor: 'white',
            plot_bgcolor: 'white',
            shapes: [] // Will add threshold lines
        };

        this.config = {
            responsive: true,
            displayModeBar: true,
            modeBarButtonsToRemove: ['lasso2d', 'select2d'],
            editable: false,  // We control editability per-shape
            toImageButtonOptions: {
                format: 'png',
                filename: 'volcano_plot',
                height: 800,
                width: 1000,
                scale: 2
            }
        };
    }

    /**
     * Load and render volcano plot data
     */
    async render(options = {}) {
        if (options.lfcThreshold !== undefined) this.lfcThreshold = options.lfcThreshold;
        if (options.padjThreshold !== undefined) this.padjThreshold = options.padjThreshold;

        try {
            // Get data from Python
            this.data = await python.getVolcanoData({
                lfc_threshold: this.lfcThreshold,
                padj_threshold: this.padjThreshold
            });

            // Prepare traces
            const traces = this.prepareTraces();

            // Update layout with threshold lines
            const layout = this.prepareLayout();

            // Render plot
            await Plotly.newPlot(this.containerId, traces, layout, this.config);

            // Setup handlers
            this.setupClickHandler();
            this.setupGrabbers();

        } catch (error) {
            console.error('Error rendering volcano plot:', error);
            this.showError(error.message);
        }
    }

    /**
     * Prepare Plotly traces with performance optimization:
     * - Low-LFC (|LFC| < 0.5) NS genes rendered as static scattergl with no hover
     * - Other genes get full hover interactivity
     */
    prepareTraces() {
        const { genes, log2FoldChange, neg_log10_padj, padj, category } = this.data;
        const selectedGenes = store.get('selectedGenes') || [];

        // Separate data by category
        const nsStaticIndices = [];  // |LFC| < 0.5, no hover (performance)
        const nsInteractiveIndices = [];  // |LFC| >= 0.5, with hover
        const upIndices = [];
        const downIndices = [];
        const selectedIndices = [];

        genes.forEach((gene, i) => {
            if (selectedGenes.includes(gene)) {
                selectedIndices.push(i);
            } else if (category[i] === 'up') {
                upIndices.push(i);
            } else if (category[i] === 'down') {
                downIndices.push(i);
            } else {
                // Split NS into static (low LFC) and interactive
                if (Math.abs(log2FoldChange[i]) < 0.5) {
                    nsStaticIndices.push(i);
                } else {
                    nsInteractiveIndices.push(i);
                }
            }
        });

        const createTrace = (indices, name, color, size = 6, opts = {}) => ({
            type: 'scatter',
            mode: 'markers',
            name: `${name} (${indices.length})`,
            x: indices.map(i => log2FoldChange[i]),
            y: indices.map(i => neg_log10_padj[i]),
            text: indices.map(i => genes[i]),
            customdata: indices.map(i => ({
                gene: genes[i],
                lfc: log2FoldChange[i],
                padj: padj[i]
            })),
            hovertemplate: opts.noHover ? undefined :
                '<b>%{text}</b><br>' +
                'Log2FC: %{x:.2f}<br>' +
                'P-adj: %{customdata.padj:.2e}<br>' +
                '<extra></extra>',
            hoverinfo: opts.noHover ? 'skip' : undefined,
            marker: {
                color: color,
                size: size,
                opacity: opts.noHover ? 0.3 : 0.7,
                line: { color: 'white', width: opts.noHover ? 0 : 0.5 }
            }
        });

        const traces = [
            // Static low-LFC genes: scattergl, no hover, fast
            createTrace(nsStaticIndices, 'NS (low LFC)', '#C8CED4', 3, { useGL: true, noHover: true }),
            // Interactive NS genes with |LFC| >= 0.5
            createTrace(nsInteractiveIndices, 'NS', '#A0AEC0', 5),
            createTrace(downIndices, 'Down', '#7BA3C9', 8),
            createTrace(upIndices, 'Up', '#C97B7B', 8)
        ];

        // Add selected genes trace
        if (selectedIndices.length > 0) {
            const selectedTrace = createTrace(selectedIndices, 'Selected', '#E8C87D', 12);
            selectedTrace.marker.symbol = 'star';
            selectedTrace.marker.line = { color: '#B39DDB', width: 2 };
            traces.push(selectedTrace);

            // Add annotations for selected genes
            this.layout.annotations = selectedIndices.map(i => ({
                x: log2FoldChange[i],
                y: neg_log10_padj[i],
                text: genes[i],
                showarrow: true,
                arrowhead: 2,
                arrowsize: 1,
                arrowwidth: 1,
                arrowcolor: '#333',
                ax: 30,
                ay: -30,
                font: { size: 11, color: '#333', weight: 'bold' },
                bgcolor: 'rgba(255, 248, 220, 0.9)',
                bordercolor: '#E8C87D',
                borderwidth: 1,
                borderpad: 3
            }));
        } else {
            this.layout.annotations = [];
        }

        return traces;
    }

    /**
     * Prepare layout with threshold lines (visual only, non-editable).
     * Draggable grabber handles are HTML overlays created by setupGrabbers().
     */
    prepareLayout() {
        const layout = { ...this.layout };

        const yMax = this.data ? Math.max(...this.data.neg_log10_padj, 10) * 1.1 : 50;
        const xMax = this.data ? Math.max(...this.data.log2FoldChange.map(Math.abs), 5) * 1.2 : 10;
        const padjY = -Math.log10(this.padjThreshold);

        layout.shapes = [
            // Shape 0: padj horizontal dashed line
            {
                type: 'line',
                x0: -xMax, x1: xMax,
                y0: padjY, y1: padjY,
                line: { color: '#c45a3c', width: 2, dash: 'dash' },
                editable: false,
                label: {
                    text: `p-adj = ${this.padjThreshold}`,
                    font: { size: 10, color: '#c45a3c' },
                    padding: 4
                }
            },
            // Shape 1: +LFC vertical dashed line
            {
                type: 'line',
                x0: this.lfcThreshold, x1: this.lfcThreshold,
                y0: 0, y1: yMax,
                line: { color: '#c45a3c', width: 2, dash: 'dash' },
                editable: false,
                label: {
                    text: `LFC = ${this.lfcThreshold.toFixed(1)}`,
                    font: { size: 10, color: '#c45a3c' },
                    padding: 4,
                    textposition: 'end'
                }
            },
            // Shape 2: -LFC vertical dashed line
            {
                type: 'line',
                x0: -this.lfcThreshold, x1: -this.lfcThreshold,
                y0: 0, y1: yMax,
                line: { color: '#c45a3c', width: 2, dash: 'dash' },
                editable: false,
                label: {
                    text: `LFC = -${this.lfcThreshold.toFixed(1)}`,
                    font: { size: 10, color: '#c45a3c' },
                    padding: 4,
                    textposition: 'end'
                }
            }
        ];

        return layout;
    }

    /**
     * Setup click handler for gene selection
     */
    setupClickHandler() {
        // Remove any existing handlers first to prevent duplicates
        if (this.container && this.container.removeAllListeners) {
            this.container.removeAllListeners('plotly_click');
        }

        this.container.on('plotly_click', (data) => {
            // Don't process clicks during drag
            if (this._isDragging) return;

            if (data.points && data.points.length > 0) {
                const point = data.points[0];
                const gene = point.text;

                if (gene) {
                    store.toggleSelectedGene(gene);
                    this.updateSelection();
                    store.emit('geneClicked', gene);
                }
            }
        });
    }

    /**
     * Create HTML overlay grabber handles that are immediately draggable.
     *   - LFC grabbers at the TOP of each vertical line (drag horizontally)
     *   - padj grabber at the RIGHT end of the horizontal line (drag vertically)
     * Uses native mousedown/mousemove/mouseup so no click-to-select is needed.
     */
    setupGrabbers() {
        // Remove previous grabbers
        this.container.querySelectorAll('.volcano-grabber').forEach(el => el.remove());

        // Ensure we can position children absolutely
        if (getComputedStyle(this.container).position === 'static') {
            this.container.style.position = 'relative';
        }

        const yMax = this.data ? Math.max(...this.data.neg_log10_padj, 10) * 1.1 : 50;
        const xMax = this.data ? Math.max(...this.data.log2FoldChange.map(Math.abs), 5) * 1.2 : 10;
        const padjY = -Math.log10(this.padjThreshold);

        // padj grabber — far right on padj line, drag vertically
        this._createGrabber({
            dataX: xMax * 0.85, dataY: padjY,
            axis: 'y', cursor: 'ns-resize', label: '↕',
            onDrag: (dx, dy) => {
                const p = Math.pow(10, -dy);
                if (p > 0 && p < 1) {
                    this.padjThreshold = parseFloat(p.toPrecision(3));
                }
            }
        });

        // +LFC grabber — top of +LFC line, drag horizontally
        this._createGrabber({
            dataX: this.lfcThreshold, dataY: yMax * 0.92,
            axis: 'x', cursor: 'ew-resize', label: '↔',
            onDrag: (dx) => {
                const lfc = Math.abs(dx);
                if (lfc > 0.01) this.lfcThreshold = parseFloat(lfc.toFixed(2));
            }
        });

        // -LFC grabber — top of -LFC line, drag horizontally
        this._createGrabber({
            dataX: -this.lfcThreshold, dataY: yMax * 0.92,
            axis: 'x', cursor: 'ew-resize', label: '↔',
            onDrag: (dx) => {
                const lfc = Math.abs(dx);
                if (lfc > 0.01) this.lfcThreshold = parseFloat(lfc.toFixed(2));
            }
        });
    }

    /** Create a single HTML grabber element and wire up drag events. */
    _createGrabber({ dataX, dataY, axis, cursor, label, onDrag }) {
        const el = document.createElement('div');
        el.className = 'volcano-grabber';
        el.textContent = label;
        el.style.cssText = `
            position: absolute; width: 24px; height: 24px;
            background: #c45a3c; border: 1.5px solid #8b2e1a; border-radius: 4px;
            color: white; font-size: 14px; line-height: 24px; text-align: center;
            cursor: ${cursor}; user-select: none; z-index: 100; opacity: 0.9;
            pointer-events: auto;
        `;
        this.container.appendChild(el);
        this._positionGrabberEl(el, dataX, dataY);

        el.addEventListener('mousedown', (e) => {
            e.preventDefault();
            e.stopPropagation();
            this._isDragging = true;
            el.style.opacity = '1';

            const onMouseMove = (e) => {
                const rect = this.container.getBoundingClientRect();
                const { dataX: newX, dataY: newY } = this._pxToData(
                    e.clientX - rect.left,
                    e.clientY - rect.top
                );

                // Call the drag callback to update threshold
                onDrag(newX, newY);

                // Move ALL grabbers to their correct positions
                this._repositionAllGrabbers();

                // Live-update the dashed lines via Plotly.relayout (throttled)
                this._syncThresholdLines();

                // Live-update input fields
                const lfcInput = document.getElementById('volcano-lfc');
                const padjInput = document.getElementById('volcano-padj');
                if (lfcInput) lfcInput.value = this.lfcThreshold;
                if (padjInput) padjInput.value = this.padjThreshold;
            };

            const onMouseUp = () => {
                this._isDragging = false;
                el.style.opacity = '0.9';
                document.removeEventListener('mousemove', onMouseMove);
                document.removeEventListener('mouseup', onMouseUp);

                // Full re-render to re-categorize genes with new thresholds
                clearTimeout(this._dragTimer);
                this._dragTimer = setTimeout(() => this.render(), 300);
            };

            document.addEventListener('mousemove', onMouseMove);
            document.addEventListener('mouseup', onMouseUp);
        });
    }

    /** Reposition all grabber elements based on current thresholds. */
    _repositionAllGrabbers() {
        const yMax = this.data ? Math.max(...this.data.neg_log10_padj, 10) * 1.1 : 50;
        const xMax = this.data ? Math.max(...this.data.log2FoldChange.map(Math.abs), 5) * 1.2 : 10;
        const padjY = -Math.log10(this.padjThreshold);

        const grabbers = this.container.querySelectorAll('.volcano-grabber');
        // Order matches creation order: [padj, +LFC, -LFC]
        if (grabbers[0]) this._positionGrabberEl(grabbers[0], xMax * 0.85, padjY);
        if (grabbers[1]) this._positionGrabberEl(grabbers[1], this.lfcThreshold, yMax * 0.92);
        if (grabbers[2]) this._positionGrabberEl(grabbers[2], -this.lfcThreshold, yMax * 0.92);
    }

    /** Position an HTML element at data coordinates using Plotly's internal layout. */
    _positionGrabberEl(el, dataX, dataY) {
        const fl = this.container._fullLayout;
        if (!fl || !fl.xaxis || !fl.yaxis) return;
        const xa = fl.xaxis, ya = fl.yaxis, sz = fl._size;
        const px = sz.l + (dataX - xa.range[0]) / (xa.range[1] - xa.range[0]) * sz.w;
        const py = sz.t + (1 - (dataY - ya.range[0]) / (ya.range[1] - ya.range[0])) * sz.h;
        el.style.left = `${px - 12}px`;
        el.style.top = `${py - 12}px`;
    }

    /** Convert container-relative pixel coordinates to data coordinates. */
    _pxToData(px, py) {
        const fl = this.container._fullLayout;
        const xa = fl.xaxis, ya = fl.yaxis, sz = fl._size;
        return {
            dataX: xa.range[0] + (px - sz.l) / sz.w * (xa.range[1] - xa.range[0]),
            dataY: ya.range[0] + (1 - (py - sz.t) / sz.h) * (ya.range[1] - ya.range[0])
        };
    }

    /** Move the Plotly threshold lines to match current thresholds (rAF-throttled). */
    _syncThresholdLines() {
        if (this._syncPending) return;
        this._syncPending = true;
        requestAnimationFrame(() => {
            this._syncPending = false;
            const padjY = -Math.log10(this.padjThreshold);
            Plotly.relayout(this.containerId, {
                'shapes[0].y0': padjY, 'shapes[0].y1': padjY,
                'shapes[0].label.text': `p-adj = ${this.padjThreshold}`,
                'shapes[1].x0': this.lfcThreshold, 'shapes[1].x1': this.lfcThreshold,
                'shapes[1].label.text': `LFC = ${this.lfcThreshold.toFixed(1)}`,
                'shapes[2].x0': -this.lfcThreshold, 'shapes[2].x1': -this.lfcThreshold,
                'shapes[2].label.text': `LFC = -${this.lfcThreshold.toFixed(1)}`
            });
        });
    }

    /**
     * Update plot to reflect current selection
     */
    async updateSelection() {
        const traces = this.prepareTraces();
        const layout = this.prepareLayout();

        await Plotly.react(this.containerId, traces, layout, this.config);
        this.setupClickHandler();
        this.setupGrabbers();

        // Auto-refresh matplotlib plot if visible
        const mode = typeof getPlotMode === 'function' ? getPlotMode('volcano') : 'plotly';
        if ((mode === 'python' || mode === 'both') && typeof generateMplVolcano === 'function') {
            try {
                await generateMplVolcano();
            } catch (e) {
                console.log('Could not auto-refresh matplotlib plot:', e);
            }
        }
    }

    /**
     * Update thresholds and re-render
     */
    async updateThresholds(lfcThreshold, padjThreshold) {
        this.lfcThreshold = lfcThreshold;
        this.padjThreshold = padjThreshold;
        await this.render();
    }

    /**
     * Show error message
     */
    showError(message) {
        this.container.innerHTML = `
            <div style="display: flex; align-items: center; justify-content: center; height: 100%; color: #e74c3c;">
                <p>Error loading volcano plot: ${message}</p>
            </div>
        `;
    }

    /**
     * Show placeholder
     */
    showPlaceholder() {
        this.container.innerHTML = `
            <div style="display: flex; align-items: center; justify-content: center; height: 100%; color: #6c757d;">
                <p>Run DESeq2 analysis to generate volcano plot</p>
            </div>
        `;
    }
}

// Export for use in other files
window.VolcanoPlot = VolcanoPlot;
