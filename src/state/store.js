/**
 * Application State Management
 * Simple reactive store with event-based updates
 */

class Store {
    constructor() {
        this.state = {
            // Data
            counts: null,
            metadata: null,
            countsShape: null,
            metadataShape: null,
            metadataColumns: [],

            // Analysis parameters
            params: {
                conditionColumn: 'condition',
                referenceLevel: 'control',
                lfcThreshold: 1.0,
                padjThreshold: 0.05
            },

            // DESeq2 results
            deseqResults: null,
            deseqSummary: null,
            normalizedCounts: null,

            // Selections
            selectedGenes: [],

            // GSEA results
            gseaResults: null,

            // Include in script flags
            includeInScript: {
                dataLoading: true,
                deseq2: true,
                volcano: false,
                boxplot: false,
                heatmap: false,
                gsea: false
            },

            // UI state
            currentPanel: 'data',
            pythonReady: false,
            isLoading: false,
            error: null
        };

        this.listeners = {};
    }

    /**
     * Get current state value
     */
    get(key) {
        if (key.includes('.')) {
            return key.split('.').reduce((obj, k) => obj?.[k], this.state);
        }
        return this.state[key];
    }

    /**
     * Set state value and notify listeners
     */
    set(key, value) {
        if (key.includes('.')) {
            const keys = key.split('.');
            const lastKey = keys.pop();
            const obj = keys.reduce((o, k) => o[k] = o[k] || {}, this.state);
            obj[lastKey] = value;
        } else {
            this.state[key] = value;
        }
        this.emit(key, value);
        this.emit('change', { key, value });
    }

    /**
     * Update multiple state values
     */
    update(updates) {
        for (const [key, value] of Object.entries(updates)) {
            this.set(key, value);
        }
    }

    /**
     * Add event listener
     */
    on(event, callback) {
        if (!this.listeners[event]) {
            this.listeners[event] = [];
        }
        this.listeners[event].push(callback);
        return () => this.off(event, callback);
    }

    /**
     * Remove event listener
     */
    off(event, callback) {
        if (this.listeners[event]) {
            this.listeners[event] = this.listeners[event].filter(cb => cb !== callback);
        }
    }

    /**
     * Emit event to all listeners
     */
    emit(event, data) {
        if (this.listeners[event]) {
            this.listeners[event].forEach(callback => callback(data));
        }
    }

    /**
     * Add a selected gene
     */
    addSelectedGene(gene) {
        if (!this.state.selectedGenes.includes(gene)) {
            this.state.selectedGenes = [...this.state.selectedGenes, gene];
            this.emit('selectedGenes', this.state.selectedGenes);
            this.emit('change', { key: 'selectedGenes', value: this.state.selectedGenes });
        }
    }

    /**
     * Remove a selected gene
     */
    removeSelectedGene(gene) {
        this.state.selectedGenes = this.state.selectedGenes.filter(g => g !== gene);
        this.emit('selectedGenes', this.state.selectedGenes);
        this.emit('change', { key: 'selectedGenes', value: this.state.selectedGenes });
    }

    /**
     * Toggle a selected gene
     */
    toggleSelectedGene(gene) {
        if (this.state.selectedGenes.includes(gene)) {
            this.removeSelectedGene(gene);
        } else {
            this.addSelectedGene(gene);
        }
    }

    /**
     * Clear all selected genes
     */
    clearSelectedGenes() {
        this.state.selectedGenes = [];
        this.emit('selectedGenes', this.state.selectedGenes);
        this.emit('change', { key: 'selectedGenes', value: this.state.selectedGenes });
    }

    /**
     * Toggle include in script flag
     */
    toggleInclude(key) {
        this.state.includeInScript[key] = !this.state.includeInScript[key];
        this.emit('includeInScript', this.state.includeInScript);
        this.emit('change', { key: 'includeInScript', value: this.state.includeInScript });
    }

    /**
     * Set include in script flag
     */
    setInclude(key, value) {
        this.state.includeInScript[key] = value;
        this.emit('includeInScript', this.state.includeInScript);
        this.emit('change', { key: 'includeInScript', value: this.state.includeInScript });
    }

    /**
     * Check if data is loaded
     */
    hasData() {
        return this.state.counts !== null && this.state.metadata !== null;
    }

    /**
     * Check if DESeq2 results are available
     */
    hasResults() {
        return this.state.deseqResults !== null;
    }
}

// Create global store instance
window.store = new Store();
