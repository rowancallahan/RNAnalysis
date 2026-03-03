/**
 * First-run Python environment setup.
 *
 * Downloads a self-contained Python distribution (python-build-standalone)
 * and installs the locked pip requirements into it.  Everything is stored
 * under the Electron userData directory so it persists across app launches
 * and doesn't require any system-level Python installation.
 */

const { spawn } = require('child_process');
const fs = require('fs');
const path = require('path');
const https = require('https');
const http = require('http');
const { createGunzip } = require('zlib');

// ---------------------------------------------------------------------------
// Configuration
// ---------------------------------------------------------------------------
const PYTHON_VERSION = '3.13.2';
const PBS_RELEASE = '20250205';
const PBS_BASE_URL =
    `https://github.com/astral-sh/python-build-standalone/releases/download/${PBS_RELEASE}`;

function getAssetInfo() {
    const plat = process.platform; // 'darwin' | 'win32' | 'linux'
    const arch = process.arch;     // 'arm64' | 'x64'

    if (plat === 'darwin' && arch === 'arm64') {
        return {
            filename: `cpython-${PYTHON_VERSION}+${PBS_RELEASE}-aarch64-apple-darwin-install_only.tar.gz`,
            pythonBin: path.join('python', 'bin', 'python3'),
            pip: ['-m', 'pip'],
        };
    }
    if (plat === 'darwin' && arch === 'x64') {
        return {
            filename: `cpython-${PYTHON_VERSION}+${PBS_RELEASE}-x86_64-apple-darwin-install_only.tar.gz`,
            pythonBin: path.join('python', 'bin', 'python3'),
            pip: ['-m', 'pip'],
        };
    }
    if (plat === 'win32') {
        return {
            filename: `cpython-${PYTHON_VERSION}+${PBS_RELEASE}-x86_64-pc-windows-msvc-install_only.tar.gz`,
            pythonBin: path.join('python', 'python.exe'),
            pip: ['-m', 'pip'],
        };
    }
    if (plat === 'linux') {
        return {
            filename: `cpython-${PYTHON_VERSION}+${PBS_RELEASE}-x86_64-unknown-linux-gnu-install_only.tar.gz`,
            pythonBin: path.join('python', 'bin', 'python3'),
            pip: ['-m', 'pip'],
        };
    }
    throw new Error(`Unsupported platform: ${plat}-${arch}`);
}

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

/** Follow redirects (GitHub releases redirect to S3). */
function downloadFile(url, destPath, onProgress) {
    return new Promise((resolve, reject) => {
        const file = fs.createWriteStream(destPath);
        const get = url.startsWith('https') ? https.get : http.get;

        get(url, (res) => {
            // Handle redirects
            if (res.statusCode >= 300 && res.statusCode < 400 && res.headers.location) {
                file.close();
                fs.unlinkSync(destPath);
                return downloadFile(res.headers.location, destPath, onProgress).then(resolve, reject);
            }

            if (res.statusCode !== 200) {
                file.close();
                fs.unlinkSync(destPath);
                return reject(new Error(`Download failed: HTTP ${res.statusCode}`));
            }

            const totalBytes = parseInt(res.headers['content-length'], 10) || 0;
            let downloaded = 0;

            res.on('data', (chunk) => {
                downloaded += chunk.length;
                if (onProgress && totalBytes > 0) {
                    onProgress(downloaded / totalBytes, downloaded, totalBytes);
                }
            });

            res.pipe(file);
            file.on('finish', () => file.close(resolve));
            file.on('error', reject);
        }).on('error', reject);
    });
}

/** Extract a .tar.gz to destDir using Node's built-in tar (via child_process). */
function extractTarGz(archivePath, destDir) {
    return new Promise((resolve, reject) => {
        fs.mkdirSync(destDir, { recursive: true });

        if (process.platform === 'win32') {
            // Windows: use tar.exe (available in Windows 10+)
            const proc = spawn('tar', ['-xzf', archivePath, '-C', destDir], {
                stdio: ['ignore', 'pipe', 'pipe'],
            });
            let stderr = '';
            proc.stderr.on('data', (d) => { stderr += d.toString(); });
            proc.on('close', (code) => {
                if (code === 0) resolve();
                else reject(new Error(`tar extraction failed (code ${code}): ${stderr}`));
            });
            proc.on('error', reject);
        } else {
            // macOS / Linux
            const proc = spawn('tar', ['-xzf', archivePath, '-C', destDir], {
                stdio: ['ignore', 'pipe', 'pipe'],
            });
            let stderr = '';
            proc.stderr.on('data', (d) => { stderr += d.toString(); });
            proc.on('close', (code) => {
                if (code === 0) resolve();
                else reject(new Error(`tar extraction failed (code ${code}): ${stderr}`));
            });
            proc.on('error', reject);
        }
    });
}

/** Run pip install with streaming output. */
function pipInstall(pythonBin, requirementsPath, onOutput) {
    return new Promise((resolve, reject) => {
        const proc = spawn(pythonBin, [
            '-m', 'pip', 'install',
            '--no-cache-dir',
            '--disable-pip-version-check',
            '-r', requirementsPath,
        ], {
            stdio: ['ignore', 'pipe', 'pipe'],
            env: { ...process.env, PYTHONDONTWRITEBYTECODE: '1' },
        });

        proc.stdout.on('data', (data) => {
            if (onOutput) onOutput(data.toString());
        });

        proc.stderr.on('data', (data) => {
            if (onOutput) onOutput(data.toString());
        });

        proc.on('close', (code) => {
            if (code === 0) resolve();
            else reject(new Error(`pip install failed with exit code ${code}`));
        });

        proc.on('error', reject);
    });
}

// ---------------------------------------------------------------------------
// Main setup function
// ---------------------------------------------------------------------------

/**
 * Set up the Python environment in the given base directory.
 *
 * @param {string} envDir  - Directory to install Python into (e.g. userData/python-env/)
 * @param {string} requirementsPath - Path to requirements-lock.txt
 * @param {(phase: string, progress: number, message: string) => void} onProgress
 * @returns {Promise<string>} Path to the Python binary
 */
async function setupPythonEnvironment(envDir, requirementsPath, onProgress) {
    const asset = getAssetInfo();
    const pythonBinPath = path.join(envDir, asset.pythonBin);

    // Check if already set up
    if (fs.existsSync(pythonBinPath)) {
        onProgress('complete', 1, 'Python environment already set up');
        return pythonBinPath;
    }

    // Step 1: Download python-build-standalone
    const downloadUrl = `${PBS_BASE_URL}/${asset.filename}`;
    const archivePath = path.join(envDir, asset.filename);

    fs.mkdirSync(envDir, { recursive: true });

    onProgress('download', 0, 'Downloading Python runtime...');

    await downloadFile(downloadUrl, archivePath, (fraction, downloaded, total) => {
        const mb = (downloaded / (1024 * 1024)).toFixed(0);
        const totalMb = (total / (1024 * 1024)).toFixed(0);
        onProgress('download', fraction * 0.4, `Downloading Python... ${mb}/${totalMb} MB`);
    });

    // Step 2: Extract
    onProgress('extract', 0.4, 'Extracting Python runtime...');
    await extractTarGz(archivePath, envDir);

    // Clean up archive
    try { fs.unlinkSync(archivePath); } catch (e) { /* ok */ }

    // Verify extraction
    if (!fs.existsSync(pythonBinPath)) {
        throw new Error(`Python binary not found after extraction at: ${pythonBinPath}`);
    }

    // Make executable on Unix
    if (process.platform !== 'win32') {
        fs.chmodSync(pythonBinPath, 0o755);
    }

    // Step 3: Install pip packages
    // Approximate total download for all packages: ~350 MB
    const ESTIMATED_TOTAL_MB = 350;
    let downloadedMB = 0;
    onProgress('install', 0.45, `Installing packages... 0/${ESTIMATED_TOTAL_MB} MB`);

    await pipInstall(pythonBinPath, requirementsPath, (output) => {
        const lines = output.split('\n').filter(l => l.trim());
        for (const line of lines) {
            // Track download progress from pip output (e.g. "Downloading numpy-2.4.1-...whl (16.5 MB)")
            const dlMatch = line.match(/Downloading\s+\S+\s+\(([\d.]+)\s*([kMG]B)\)/i);
            if (dlMatch) {
                let mb = parseFloat(dlMatch[1]);
                const unit = dlMatch[2].toUpperCase();
                if (unit === 'KB') mb /= 1024;
                else if (unit === 'GB') mb *= 1024;
                downloadedMB += mb;
                const pct = Math.min(downloadedMB / ESTIMATED_TOTAL_MB, 1.0);
                onProgress('install', 0.45 + pct * 0.50, `Installing packages... ${Math.round(downloadedMB)}/${ESTIMATED_TOTAL_MB} MB`);
            } else if (line.startsWith('Installing collected packages')) {
                onProgress('install', 0.93, 'Finalizing installation...');
            } else if (line.startsWith('Successfully installed')) {
                onProgress('install', 0.97, 'Packages installed successfully');
            }
        }
    });

    // Step 4: Write a marker file so we know setup is complete
    const markerPath = path.join(envDir, '.setup-complete');
    fs.writeFileSync(markerPath, JSON.stringify({
        timestamp: new Date().toISOString(),
        pythonVersion: PYTHON_VERSION,
        platform: process.platform,
        arch: process.arch,
    }));

    onProgress('complete', 1, 'Setup complete!');
    return pythonBinPath;
}

/**
 * Check if the Python environment is already set up.
 */
function isPythonSetUp(envDir) {
    const asset = getAssetInfo();
    const pythonBinPath = path.join(envDir, asset.pythonBin);
    const markerPath = path.join(envDir, '.setup-complete');
    return fs.existsSync(pythonBinPath) && fs.existsSync(markerPath);
}

/**
 * Get the path to the Python binary in the environment.
 */
function getPythonBinPath(envDir) {
    const asset = getAssetInfo();
    return path.join(envDir, asset.pythonBin);
}

module.exports = {
    setupPythonEnvironment,
    isPythonSetUp,
    getPythonBinPath,
    getAssetInfo,
};
