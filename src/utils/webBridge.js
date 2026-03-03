/**
 * Web Bridge - Browser-compatible replacement for Electron's preload API.
 * Provides the same window.api interface using HTTP fetch to the dev server.
 */

(function() {
    // Only activate if not running inside Electron
    if (window.api) {
        console.log('[WebBridge] Electron API detected, skipping web bridge.');
        return;
    }

    const API_BASE = window.location.origin;

    window.api = {
        /**
         * Run a Python command via HTTP POST.
         */
        runPython: async function(command, params) {
            const response = await fetch(`${API_BASE}/api/command`, {
                method: 'POST',
                headers: { 'Content-Type': 'application/json' },
                body: JSON.stringify({ command, params: params || {} })
            });

            const result = await response.json();

            if (result.error) {
                throw new Error(result.error);
            }

            return result.data;
        },

        /**
         * Open file dialog - uses browser file input.
         * Returns a File object handle; the caller should use uploadFile() to send it.
         */
        openFileDialog: async function(options) {
            return new Promise((resolve) => {
                const input = document.createElement('input');
                input.type = 'file';

                // Map Electron filter format to accept attribute
                if (options && options.filters) {
                    const extensions = options.filters
                        .flatMap(f => f.extensions.filter(e => e !== '*'))
                        .map(e => '.' + e);
                    if (extensions.length > 0) {
                        input.accept = extensions.join(',');
                    }
                }

                input.onchange = () => {
                    if (input.files && input.files.length > 0) {
                        // Store file reference for later upload
                        const file = input.files[0];
                        window._lastSelectedFile = file;
                        // Return a fake path that indicates web upload
                        resolve('__web_upload__:' + file.name);
                    } else {
                        resolve(null);
                    }
                };

                input.click();
            });
        },

        /**
         * Save file dialog - triggers browser download.
         */
        saveFileDialog: async function(options) {
            // In browser mode, we just return a filename and handle download in saveFile
            const defaultName = (options && options.defaultPath) || 'download.txt';
            return defaultName;
        },

        /**
         * Save file - triggers a browser download.
         */
        saveFile: async function(filePath, content) {
            const blob = new Blob([content], { type: 'text/plain' });
            const url = URL.createObjectURL(blob);
            const a = document.createElement('a');
            a.href = url;
            a.download = filePath.split('/').pop() || 'download.txt';
            document.body.appendChild(a);
            a.click();
            document.body.removeChild(a);
            URL.revokeObjectURL(url);
            return true;
        },

        /**
         * Read file - not directly available in browser, but files are uploaded.
         */
        readFile: async function(filePath) {
            if (window._lastSelectedFile) {
                return new Promise((resolve, reject) => {
                    const reader = new FileReader();
                    reader.onload = () => resolve(reader.result);
                    reader.onerror = () => reject(new Error('Failed to read file'));
                    reader.readAsText(window._lastSelectedFile);
                });
            }
            throw new Error('File reading not available in browser mode');
        },

        /**
         * Python ready event - poll the server status.
         */
        onPythonReady: function(callback) {
            // Poll until the server API is ready
            async function checkReady() {
                try {
                    const response = await fetch(`${API_BASE}/api/status`);
                    const data = await response.json();
                    if (data.ready) {
                        callback();
                        return;
                    }
                } catch (e) {
                    // Server not ready yet
                }
                setTimeout(checkReady, 500);
            }
            checkReady();
            return () => {}; // Return no-op cleanup
        },

        /**
         * Python error event - no-op in browser mode.
         */
        onPythonError: function(callback) {
            // Errors are handled via HTTP responses
            return () => {};
        }
    };

    /**
     * Upload a file to the server for loading.
     * Called instead of passing file paths in browser mode.
     */
    window.api.uploadFile = async function(command) {
        if (!window._lastSelectedFile) {
            throw new Error('No file selected');
        }

        const formData = new FormData();
        formData.append('file', window._lastSelectedFile);
        formData.append('command', command);

        const response = await fetch(`${API_BASE}/api/upload`, {
            method: 'POST',
            body: formData
        });

        const result = await response.json();

        if (result.error) {
            throw new Error(result.error);
        }

        return result.data;
    };

    console.log('[WebBridge] Browser API bridge loaded. Using HTTP backend.');
})();
