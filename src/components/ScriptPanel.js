/**
 * Script Panel Component
 * Generates and displays reproducible Python script
 */

class ScriptPanel {
    constructor() {
        this.codeElement = document.getElementById('generated-script');
        this.copyButton = document.getElementById('btn-copy-script');
        this.saveButton = document.getElementById('btn-save-script');
        this.script = '';

        // Include checkboxes (in script panel - step 5)
        this.checkboxes = {
            dataLoading: document.getElementById('script-include-data-loading'),
            deseq2: document.getElementById('script-include-deseq2'),
            volcano: document.getElementById('script-include-volcano'),
            boxplot: document.getElementById('script-include-boxplot'),
            heatmap: document.getElementById('script-include-heatmap'),
            gsea: document.getElementById('script-include-gsea')
        };

        this.setupEventListeners();
    }

    /**
     * Setup event listeners
     */
    setupEventListeners() {
        // Copy button
        if (this.copyButton) {
            this.copyButton.addEventListener('click', () => this.copyToClipboard());
        }

        // Save button
        if (this.saveButton) {
            this.saveButton.addEventListener('click', () => this.saveScript());
        }

        // Include checkboxes
        Object.entries(this.checkboxes).forEach(([key, checkbox]) => {
            if (checkbox) {
                checkbox.addEventListener('change', () => {
                    store.setInclude(key, checkbox.checked);
                    this.updateScript();
                });
            }
        });

        // Listen for state changes
        store.on('includeInScript', () => this.syncCheckboxes());
        store.on('selectedGenes', () => this.updateScript());
        store.on('deseqResults', () => this.updateScript());
        store.on('activeTool', () => this.updateScript());
    }

    /**
     * Sync checkboxes with store state
     */
    syncCheckboxes() {
        const includes = store.get('includeInScript');
        Object.entries(this.checkboxes).forEach(([key, checkbox]) => {
            if (checkbox && includes[key] !== undefined) {
                checkbox.checked = includes[key];
            }
        });
    }

    /**
     * Update generated script
     */
    async updateScript() {
        if (!python.ready) {
            this.showPlaceholder();
            return;
        }

        try {
            const includes = store.get('includeInScript');
            this.script = await python.generateScript(includes);
            this.renderScript();
        } catch (error) {
            console.error('Error generating script:', error);
            this.script = `# Error generating script: ${error.message}`;
            this.renderScript();
        }
    }

    /**
     * Render script as plain text (safe rendering without HTML injection issues)
     */
    renderScript() {
        if (!this.codeElement) return;

        // Use textContent for safe plain text rendering
        this.codeElement.textContent = this.script || '# No script generated yet';
    }

    /**
     * Copy script to clipboard
     */
    async copyToClipboard() {
        try {
            await navigator.clipboard.writeText(this.script);
            this.showNotification('Script copied to clipboard!');
        } catch (error) {
            console.error('Failed to copy:', error);
            this.showNotification('Failed to copy script', true);
        }
    }

    /**
     * Save script to file
     */
    async saveScript() {
        try {
            const filePath = await python.saveFileDialog({
                defaultPath: 'analysis_script.py',
                filters: [
                    { name: 'Python Script', extensions: ['py'] },
                    { name: 'All Files', extensions: ['*'] }
                ]
            });

            if (filePath) {
                await python.saveFile(filePath, this.script);
                this.showNotification('Script saved successfully!');
            }
        } catch (error) {
            console.error('Failed to save:', error);
            this.showNotification('Failed to save script', true);
        }
    }

    /**
     * Show notification
     */
    showNotification(message, isError = false) {
        // Simple notification
        const originalText = this.copyButton?.textContent;
        if (this.copyButton) {
            this.copyButton.textContent = isError ? 'Error!' : 'Copied!';
            setTimeout(() => {
                this.copyButton.textContent = originalText;
            }, 2000);
        }
    }

    /**
     * Show placeholder
     */
    showPlaceholder() {
        if (this.codeElement) {
            this.codeElement.textContent = '# Load data and run analysis to generate script...';
        }
    }
}

// Export for use in other files
window.ScriptPanel = ScriptPanel;
