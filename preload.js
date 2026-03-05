const { contextBridge, ipcRenderer } = require('electron');

// Expose protected methods to renderer process
contextBridge.exposeInMainWorld('api', {
    // Python communication
    runPython: (command, params) => ipcRenderer.invoke('python-command', command, params),

    // File operations
    openFileDialog: (options) => ipcRenderer.invoke('open-file-dialog', options),
    saveFileDialog: (options) => ipcRenderer.invoke('save-file-dialog', options),
    saveFile: (filePath, content) => ipcRenderer.invoke('save-file', filePath, content),
    readFile: (filePath) => ipcRenderer.invoke('read-file', filePath),

    // Platform & setup
    getPlatformInfo: () => ipcRenderer.invoke('get-platform-info'),
    checkPythonSetup: () => ipcRenderer.invoke('check-python-setup'),
    runPythonSetup: () => ipcRenderer.invoke('run-python-setup'),
    startPythonBackend: () => ipcRenderer.invoke('start-python-backend'),

    // R environment & GEO/recount3
    checkRSetup: () => ipcRenderer.invoke('check-r-setup'),
    runRSetup: () => ipcRenderer.invoke('run-r-setup'),
    runRScript: (scriptName, args) => ipcRenderer.invoke('run-r-script', scriptName, args),
    onRSetupProgress: (callback) => {
        ipcRenderer.on('r-setup-progress', (event, data) => callback(data));
        return () => ipcRenderer.removeListener('r-setup-progress', callback);
    },
    onRScriptProgress: (callback) => {
        ipcRenderer.on('r-script-progress', (event, data) => callback(data));
        return () => ipcRenderer.removeListener('r-script-progress', callback);
    },

    // Event listeners
    onPythonReady: (callback) => {
        ipcRenderer.on('python-ready', callback);
        return () => ipcRenderer.removeListener('python-ready', callback);
    },
    onPythonError: (callback) => {
        ipcRenderer.on('python-error', (event, message) => callback(message));
        return () => ipcRenderer.removeListener('python-error', callback);
    },
    onPythonSetupProgress: (callback) => {
        ipcRenderer.on('python-setup-progress', (event, data) => callback(data));
        return () => ipcRenderer.removeListener('python-setup-progress', callback);
    },
});
