const { app, BrowserWindow, ipcMain, dialog } = require('electron');
const path = require('path');
const { spawn } = require('child_process');
const fs = require('fs');
const { setupPythonEnvironment, isPythonSetUp, getPythonBinPath } = require('./scripts/setup-python');

let mainWindow;
let pythonProcess = null;
let pythonReady = false;
let pendingRequests = new Map();
let requestId = 0;

// Debug logging to file (helps diagnose packaged app issues)
const logFile = path.join(app.getPath('userData'), 'debug.log');
function debugLog(...args) {
    const msg = `[${new Date().toISOString()}] ${args.join(' ')}\n`;
    fs.appendFileSync(logFile, msg);
    console.log(...args);
}

// ---------------------------------------------------------------------------
// Paths
// ---------------------------------------------------------------------------

/** Directory where the first-run Python env is installed. */
function getPythonEnvDir() {
    return path.join(app.getPath('userData'), 'python-env');
}

/** Path to requirements-lock.txt shipped with the app. */
function getRequirementsPath() {
    if (app.isPackaged) {
        return path.join(process.resourcesPath, 'python', 'requirements-lock.txt');
    }
    return path.join(__dirname, 'python', 'requirements-lock.txt');
}

/** Resolve the Python binary for the current platform. */
function getPythonPath() {
    // 1. In development, prefer local venv
    if (!app.isPackaged) {
        const devBin = process.platform === 'win32'
            ? path.join(__dirname, 'venv', 'Scripts', 'python.exe')
            : path.join(__dirname, 'venv', 'bin', 'python');
        if (fs.existsSync(devBin)) {
            return devBin;
        }
    }

    // 2. First-run installed env (userData)
    const envDir = getPythonEnvDir();
    if (isPythonSetUp(envDir)) {
        const binPath = getPythonBinPath(envDir);
        if (fs.existsSync(binPath)) return binPath;
    }

    // 3. Fallback to system Python (should not happen in production)
    if (process.platform === 'win32') return 'python';
    const candidates = [
        '/opt/homebrew/bin/python3',
        '/usr/local/bin/python3',
        '/usr/bin/python3',
        'python3',
    ];
    for (const p of candidates) {
        try { if (p.startsWith('/') && fs.existsSync(p)) return p; } catch (e) { /* next */ }
    }
    return 'python3';
}

function getPythonScriptPath() {
    if (app.isPackaged) {
        return path.join(process.resourcesPath, 'python', 'main_handler.py');
    }
    return path.join(__dirname, 'python', 'main_handler.py');
}

// ---------------------------------------------------------------------------
// Window
// ---------------------------------------------------------------------------

function createWindow() {
    mainWindow = new BrowserWindow({
        width: 1400,
        height: 900,
        minWidth: 1200,
        minHeight: 700,
        webPreferences: {
            preload: path.join(__dirname, 'preload.js'),
            contextIsolation: true,
            nodeIntegration: false
        },
        titleBarStyle: 'hiddenInset',
        show: false
    });

    mainWindow.loadFile('src/index.html');

    mainWindow.once('ready-to-show', () => {
        mainWindow.show();
        // Don't auto-start Python here — renderer will check setup first
    });

    mainWindow.on('closed', () => {
        mainWindow = null;
        if (pythonProcess) {
            pythonProcess.kill();
        }
    });
}

// ---------------------------------------------------------------------------
// Python backend
// ---------------------------------------------------------------------------

function startPythonBackend() {
    const pythonPath = getPythonPath();
    const scriptPath = getPythonScriptPath();

    debugLog('startPythonBackend: pythonPath:', pythonPath, 'scriptPath:', scriptPath);
    debugLog('startPythonBackend: pythonPath exists:', fs.existsSync(pythonPath));
    debugLog('startPythonBackend: scriptPath exists:', fs.existsSync(scriptPath));
    debugLog('startPythonBackend: cwd:', path.dirname(scriptPath));

    pythonProcess = spawn(pythonPath, [scriptPath], {
        stdio: ['pipe', 'pipe', 'pipe'],
        cwd: path.dirname(scriptPath)
    });

    debugLog('startPythonBackend: process spawned, pid:', pythonProcess.pid);

    let buffer = '';

    pythonProcess.stdout.on('data', (data) => {
        const raw = data.toString();
        debugLog('Python stdout:', raw.substring(0, 200));
        buffer += raw;

        const lines = buffer.split('\n');
        buffer = lines.pop();

        for (const line of lines) {
            if (line.trim()) {
                try {
                    const response = JSON.parse(line);

                    if (response.type === 'ready') {
                        pythonReady = true;
                        debugLog('Python backend ready!');

                        // Send config to Python so it knows paths
                        sendToPython('init_config', {
                            resources_path: app.isPackaged ? process.resourcesPath : __dirname,
                            user_data_path: app.getPath('userData'),
                            platform: process.platform,
                            is_packaged: app.isPackaged,
                        }).catch(err => debugLog('init_config failed:', err.message));

                        mainWindow.webContents.send('python-ready');
                    } else if (response.id !== undefined && pendingRequests.has(response.id)) {
                        const { resolve, reject } = pendingRequests.get(response.id);
                        pendingRequests.delete(response.id);

                        if (response.error) {
                            reject(new Error(response.error));
                        } else {
                            resolve(response.data);
                        }
                    }
                } catch (e) {
                    debugLog('Failed to parse Python response:', e.message, 'line:', line.substring(0, 200));
                }
            }
        }
    });

    pythonProcess.stderr.on('data', (data) => {
        debugLog('Python stderr:', data.toString().substring(0, 500));
    });

    pythonProcess.on('close', (code) => {
        debugLog('Python process exited with code:', code);
        pythonReady = false;
        pythonProcess = null;

        for (const [id, { reject }] of pendingRequests) {
            reject(new Error('Python process terminated'));
        }
        pendingRequests.clear();
    });

    pythonProcess.on('error', (err) => {
        debugLog('Failed to start Python process:', err.message);
        mainWindow.webContents.send('python-error', err.message);
    });
}

function sendToPython(command, params = {}) {
    return new Promise((resolve, reject) => {
        if (!pythonProcess || !pythonReady) {
            reject(new Error('Python backend not ready'));
            return;
        }

        const id = requestId++;
        const request = JSON.stringify({ id, command, params }) + '\n';

        pendingRequests.set(id, { resolve, reject });

        setTimeout(() => {
            if (pendingRequests.has(id)) {
                pendingRequests.delete(id);
                reject(new Error('Request timeout'));
            }
        }, 300000);

        pythonProcess.stdin.write(request);
    });
}

// ---------------------------------------------------------------------------
// IPC Handlers
// ---------------------------------------------------------------------------

// Python commands
ipcMain.handle('python-command', async (event, command, params) => {
    try {
        return await sendToPython(command, params);
    } catch (error) {
        throw error;
    }
});

// Platform info
ipcMain.handle('get-platform-info', async () => {
    return {
        platform: process.platform,
        arch: process.arch,
        isPackaged: app.isPackaged,
    };
});

// Check if Python is set up
ipcMain.handle('check-python-setup', async () => {
    debugLog('check-python-setup called, isPackaged:', app.isPackaged);

    // In dev mode with a local venv, always "set up"
    if (!app.isPackaged) {
        const devBin = process.platform === 'win32'
            ? path.join(__dirname, 'venv', 'Scripts', 'python.exe')
            : path.join(__dirname, 'venv', 'bin', 'python');
        debugLog('Dev mode, checking venv at:', devBin, 'exists:', fs.existsSync(devBin));
        if (fs.existsSync(devBin)) {
            return { ready: true, source: 'dev-venv' };
        }
    }

    const envDir = getPythonEnvDir();
    const ready = isPythonSetUp(envDir);
    debugLog('Checking envDir:', envDir, 'isPythonSetUp:', ready);

    if (ready) {
        return { ready: true, source: 'installed' };
    }

    return { ready: false, source: null };
});

// Run first-time Python setup (download + install)
ipcMain.handle('run-python-setup', async () => {
    const envDir = getPythonEnvDir();
    const reqPath = getRequirementsPath();
    debugLog('run-python-setup: envDir:', envDir, 'reqPath:', reqPath, 'exists:', fs.existsSync(reqPath));

    try {
        const pythonBin = await setupPythonEnvironment(envDir, reqPath, (phase, progress, message) => {
            // Stream progress to renderer
            if (mainWindow && !mainWindow.isDestroyed()) {
                mainWindow.webContents.send('python-setup-progress', { phase, progress, message });
            }
        });
        return { success: true, pythonPath: pythonBin };
    } catch (error) {
        return { success: false, error: error.message };
    }
});

// Start Python backend (called by renderer after setup is confirmed)
ipcMain.handle('start-python-backend', async () => {
    debugLog('start-python-backend called, pythonReady:', pythonReady, 'pythonProcess:', !!pythonProcess);
    if (pythonProcess) {
        debugLog('Python process already exists, skipping spawn');
        return pythonReady ? { ready: true } : { started: true, waiting: true };
    }
    const pyPath = getPythonPath();
    debugLog('Will start Python with:', pyPath, 'exists:', fs.existsSync(pyPath));
    startPythonBackend();
    return { started: true };
});

// File dialogs
ipcMain.handle('open-file-dialog', async (event, options) => {
    const result = await dialog.showOpenDialog(mainWindow, {
        properties: ['openFile'],
        filters: options.filters || [
            { name: 'Data Files', extensions: ['csv', 'tsv', 'xlsx', 'xls'] },
            { name: 'All Files', extensions: ['*'] }
        ]
    });

    if (!result.canceled && result.filePaths.length > 0) {
        return result.filePaths[0];
    }
    return null;
});

ipcMain.handle('save-file-dialog', async (event, options) => {
    const result = await dialog.showSaveDialog(mainWindow, {
        defaultPath: options.defaultPath || 'analysis_script.py',
        filters: options.filters || [
            { name: 'Python Script', extensions: ['py'] },
            { name: 'All Files', extensions: ['*'] }
        ]
    });

    if (!result.canceled && result.filePath) {
        return result.filePath;
    }
    return null;
});

ipcMain.handle('save-file', async (event, filePath, content) => {
    try {
        fs.writeFileSync(filePath, content, 'utf8');
        return true;
    } catch (error) {
        throw error;
    }
});

ipcMain.handle('read-file', async (event, filePath) => {
    try {
        return fs.readFileSync(filePath, 'utf8');
    } catch (error) {
        throw error;
    }
});

// ---------------------------------------------------------------------------
// App lifecycle
// ---------------------------------------------------------------------------

app.whenReady().then(createWindow);

app.on('window-all-closed', () => {
    if (process.platform !== 'darwin') {
        app.quit();
    }
});

app.on('activate', () => {
    if (BrowserWindow.getAllWindows().length === 0) {
        createWindow();
    }
});

app.on('will-quit', () => {
    if (pythonProcess) {
        pythonProcess.kill();
    }
});
