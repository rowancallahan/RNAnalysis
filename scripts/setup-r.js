/**
 * On-demand R environment setup.
 *
 * Downloads R from CRAN and installs recount3 + GEOquery packages.
 * Only triggered when user first uses a GEO/recount feature.
 * Everything stored under Electron userData so no system R needed.
 */

const { spawn } = require('child_process');
const fs = require('fs');
const path = require('path');
const https = require('https');
const http = require('http');

// ---------------------------------------------------------------------------
// Configuration
// ---------------------------------------------------------------------------
const R_VERSION = '4.5.2';

function getRDownloadInfo() {
    const plat = process.platform;
    const arch = process.arch;

    if (plat === 'darwin' && arch === 'arm64') {
        return {
            url: `https://cran.r-project.org/bin/macosx/big-sur-arm64/base/R-${R_VERSION}-arm64.pkg`,
            filename: `R-${R_VERSION}-arm64.pkg`,
            type: 'pkg',
        };
    }
    if (plat === 'darwin' && arch === 'x64') {
        return {
            url: `https://cran.r-project.org/bin/macosx/big-sur-x86_64/base/R-${R_VERSION}-x86_64.pkg`,
            filename: `R-${R_VERSION}-x86_64.pkg`,
            type: 'pkg',
        };
    }
    if (plat === 'win32') {
        return {
            url: `https://cran.r-project.org/bin/windows/base/R-${R_VERSION}-win.exe`,
            filename: `R-${R_VERSION}-win.exe`,
            type: 'exe',
        };
    }
    if (plat === 'linux') {
        return {
            url: `https://cran.r-project.org/src/base/R-4/R-${R_VERSION}.tar.gz`,
            filename: `R-${R_VERSION}.tar.gz`,
            type: 'tarball',
        };
    }
    throw new Error(`Unsupported platform: ${plat}-${arch}`);
}

/**
 * Find the actual R version directory inside R.framework/Versions/.
 * CRAN builds may name it "4.5", "4.5-arm64", "4.5-x86_64", etc.
 */
function findRVersionDir(envDir) {
    const versionsDir = path.join(envDir, 'R.framework', 'Versions');
    const majorMinor = R_VERSION.split('.').slice(0, 2).join('.');

    if (!fs.existsSync(versionsDir)) return null;

    // Try exact match first, then architecture-suffixed variants
    const candidates = [majorMinor, `${majorMinor}-arm64`, `${majorMinor}-x86_64`];
    for (const c of candidates) {
        const p = path.join(versionsDir, c);
        if (fs.existsSync(p)) return p;
    }

    // Fallback: use whatever directory exists (skip symlinks like "Current")
    try {
        const entries = fs.readdirSync(versionsDir, { withFileTypes: true });
        for (const e of entries) {
            if (e.isDirectory() && e.name !== 'Current' && e.name.startsWith(majorMinor.split('.')[0])) {
                return path.join(versionsDir, e.name);
            }
        }
    } catch (_) { /* ignore */ }

    return null;
}

/**
 * Get the path to the Rscript binary in the environment.
 */
function getRscriptBinPath(envDir) {
    if (process.platform === 'darwin') {
        const versionDir = findRVersionDir(envDir);
        if (versionDir) {
            return path.join(versionDir, 'Resources', 'bin', 'Rscript');
        }
        // Fallback to simple naming for pre-extraction check
        return path.join(envDir, 'R.framework', 'Versions',
            R_VERSION.split('.').slice(0, 2).join('.'), 'Resources', 'bin', 'Rscript');
    }
    if (process.platform === 'win32') {
        return path.join(envDir, 'R', 'bin', 'Rscript.exe');
    }
    // Linux
    return path.join(envDir, 'R', 'bin', 'Rscript');
}

/**
 * Get R_HOME for environment variables when spawning R.
 */
function getRHome(envDir) {
    if (process.platform === 'darwin') {
        const versionDir = findRVersionDir(envDir);
        if (versionDir) {
            return path.join(versionDir, 'Resources');
        }
        return path.join(envDir, 'R.framework', 'Versions',
            R_VERSION.split('.').slice(0, 2).join('.'), 'Resources');
    }
    if (process.platform === 'win32') {
        return path.join(envDir, 'R');
    }
    return path.join(envDir, 'R');
}

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

/** Follow redirects (CRAN mirrors redirect). */
function downloadFile(url, destPath, onProgress) {
    return new Promise((resolve, reject) => {
        const file = fs.createWriteStream(destPath);
        const get = url.startsWith('https') ? https.get : http.get;

        get(url, (res) => {
            if (res.statusCode >= 300 && res.statusCode < 400 && res.headers.location) {
                file.close();
                try { fs.unlinkSync(destPath); } catch (e) { /* ok */ }
                return downloadFile(res.headers.location, destPath, onProgress).then(resolve, reject);
            }

            if (res.statusCode !== 200) {
                file.close();
                try { fs.unlinkSync(destPath); } catch (e) { /* ok */ }
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

/** Extract macOS .pkg to get R.framework using pkgutil. */
function extractMacPkg(pkgPath, destDir) {
    return new Promise((resolve, reject) => {
        const expandDir = path.join(destDir, '_pkg_expanded');

        // pkgutil fails if target dir exists — clean up any previous attempt
        if (fs.existsSync(expandDir)) {
            fs.rmSync(expandDir, { recursive: true, force: true });
        }

        // Also remove partial R.framework from previous attempt
        const fwDest = path.join(destDir, 'R.framework');
        if (fs.existsSync(fwDest)) {
            fs.rmSync(fwDest, { recursive: true, force: true });
        }

        // pkgutil --expand-full extracts .pkg including Payload contents
        const proc = spawn('pkgutil', ['--expand-full', pkgPath, expandDir], {
            stdio: ['ignore', 'pipe', 'pipe'],
        });
        let stderr = '';
        proc.stderr.on('data', (d) => { stderr += d.toString(); });
        proc.on('close', (code) => {
            if (code !== 0) {
                return reject(new Error(`pkgutil failed (code ${code}): ${stderr}`));
            }

            // R.framework location varies by R version:
            //   Older: R-fw.pkg/Payload/Library/Frameworks/R.framework
            //   Newer (4.5+): R-fw.pkg/Payload/R.framework
            const fwDest = path.join(destDir, 'R.framework');
            const candidates = [
                path.join(expandDir, 'R-fw.pkg', 'Payload', 'R.framework'),
                path.join(expandDir, 'R-fw.pkg', 'Payload', 'Library', 'Frameworks', 'R.framework'),
            ];
            const fwSrc = candidates.find(p => fs.existsSync(p));

            if (!fwSrc) {
                // Debug: list what's actually in the expanded directory
                let contents = '';
                try {
                    const rfwPayload = path.join(expandDir, 'R-fw.pkg', 'Payload');
                    if (fs.existsSync(rfwPayload)) {
                        contents = fs.readdirSync(rfwPayload).join(', ');
                    } else {
                        contents = 'R-fw.pkg/Payload not found. Top-level: ' +
                            fs.readdirSync(expandDir).join(', ');
                    }
                } catch (e) { contents = e.message; }
                return reject(new Error(`R.framework not found in extracted .pkg. Contents: ${contents}`));
            }

            // Move R.framework to destDir
            fs.renameSync(fwSrc, fwDest);

            // Clean up expanded pkg
            fs.rmSync(expandDir, { recursive: true, force: true });

            resolve();
        });
        proc.on('error', reject);
    });
}

/** Run an R script, streaming output. Uses the R shell script (not the Rscript
 *  binary) because the Rscript binary has a hardcoded R_HOME that points to
 *  /Library/Frameworks/... which doesn't exist for our portable install.
 *  The R shell script was patched during extraction to use the correct R_HOME. */
function runRscript(rscriptBin, scriptPath, args, envDir, onOutput) {
    return new Promise((resolve, reject) => {
        const rHome = getRHome(envDir);
        const libraryPath = path.join(envDir, 'library');
        fs.mkdirSync(libraryPath, { recursive: true });

        const env = {
            ...process.env,
            R_HOME: rHome,
            R_LIBS_USER: libraryPath,
            R_LIBS_SITE: libraryPath,
        };

        // On macOS, ensure the R framework's lib is in DYLD path
        if (process.platform === 'darwin') {
            const rLib = path.join(rHome, 'lib');
            env.DYLD_LIBRARY_PATH = rLib + (env.DYLD_LIBRARY_PATH ? ':' + env.DYLD_LIBRARY_PATH : '');
        }

        // Use the R shell script with --no-echo --file= instead of the Rscript binary,
        // because Rscript has a hardcoded R_HOME baked into the binary.
        const rBin = path.join(path.dirname(rscriptBin), 'R');
        const useR = process.platform === 'darwin' && fs.existsSync(rBin);
        const cmd = useR ? rBin : rscriptBin;
        const cmdArgs = useR
            ? ['--no-echo', '--no-restore', `--file=${scriptPath}`, '--args', ...args]
            : [scriptPath, ...args];

        const proc = spawn(cmd, cmdArgs, {
            stdio: ['ignore', 'pipe', 'pipe'],
            env,
        });

        let stdout = '';
        let stderr = '';

        proc.stdout.on('data', (data) => {
            stdout += data.toString();
        });

        proc.stderr.on('data', (data) => {
            stderr += data.toString();
            if (onOutput) onOutput(data.toString());
        });

        proc.on('close', (code) => {
            if (code === 0) resolve({ stdout, stderr });
            else reject(new Error(`Rscript failed (code ${code}): ${stderr}`));
        });

        proc.on('error', reject);
    });
}

// ---------------------------------------------------------------------------
// Main setup function
// ---------------------------------------------------------------------------

/**
 * Set up the R environment: download R, install packages.
 *
 * @param {string} envDir - Directory to install R into (userData/r-env/)
 * @param {string} rScriptsDir - Directory containing R scripts (install_packages.R, etc.)
 * @param {(phase: string, progress: number, message: string) => void} onProgress
 * @returns {Promise<string>} Path to Rscript binary
 */
async function setupREnvironment(envDir, rScriptsDir, onProgress) {
    const rscriptBin = getRscriptBinPath(envDir);

    // Check if already set up
    const markerPath = path.join(envDir, '.r-setup-complete');
    if (fs.existsSync(rscriptBin) && fs.existsSync(markerPath)) {
        onProgress('complete', 1, 'R environment already set up');
        return rscriptBin;
    }

    const dlInfo = getRDownloadInfo();
    const archivePath = path.join(envDir, dlInfo.filename);

    fs.mkdirSync(envDir, { recursive: true });

    // Step 1: Download R (skip if already downloaded)
    if (fs.existsSync(archivePath)) {
        onProgress('download', 0.3, 'R runtime already downloaded, skipping...');
    } else {
        onProgress('download', 0, 'Downloading R runtime...');

        await downloadFile(dlInfo.url, archivePath, (fraction, downloaded, total) => {
            const mb = (downloaded / (1024 * 1024)).toFixed(0);
            const totalMb = (total / (1024 * 1024)).toFixed(0);
            onProgress('download', fraction * 0.3, `Downloading R... ${mb}/${totalMb} MB`);
        });
    }

    // Step 2: Extract
    onProgress('extract', 0.3, 'Extracting R runtime...');

    if (dlInfo.type === 'pkg') {
        await extractMacPkg(archivePath, envDir);
    } else {
        throw new Error(`R extraction not yet implemented for type: ${dlInfo.type}`);
    }

    // Clean up archive
    try { fs.unlinkSync(archivePath); } catch (e) { /* ok */ }

    // Re-resolve path now that R.framework exists on disk
    const rscriptBinResolved = getRscriptBinPath(envDir);

    // Verify extraction
    if (!fs.existsSync(rscriptBinResolved)) {
        // Debug: list what's in Versions/
        let debug = '';
        try {
            const vDir = path.join(envDir, 'R.framework', 'Versions');
            debug = fs.existsSync(vDir) ? fs.readdirSync(vDir).join(', ') : 'Versions/ not found';
        } catch (e) { debug = e.message; }
        throw new Error(`Rscript binary not found after extraction at: ${rscriptBinResolved} (Versions contents: ${debug})`);
    }

    // Make executable on Unix
    if (process.platform !== 'win32') {
        fs.chmodSync(rscriptBinResolved, 0o755);
        // Also chmod the R binary and patch its hardcoded R_HOME_DIR
        const rBin = path.join(path.dirname(rscriptBinResolved), 'R');
        if (fs.existsSync(rBin)) {
            fs.chmodSync(rBin, 0o755);
            // The R shell script has R_HOME_DIR hardcoded to /Library/Frameworks/R.framework/Resources
            // and ignores the R_HOME env var. Patch it to point to our local copy.
            const rHome = getRHome(envDir);
            let rScript = fs.readFileSync(rBin, 'utf8');
            rScript = rScript.replace(
                /^R_HOME_DIR=.*/m,
                `R_HOME_DIR="${rHome}"`
            );
            fs.writeFileSync(rBin, rScript, 'utf8');
            fs.chmodSync(rBin, 0o755);
        }
        // Also chmod the exec/R binary
        const execR = path.join(path.dirname(rscriptBinResolved), 'exec', 'R');
        if (fs.existsSync(execR)) fs.chmodSync(execR, 0o755);
    }

    // Step 3: Install R packages (BiocManager, recount3, GEOquery, jsonlite)
    onProgress('install', 0.4, 'Installing R packages (this may take several minutes)...');

    const libraryPath = path.join(envDir, 'library');
    const installScript = path.join(rScriptsDir, 'install_packages.R');

    await runRscript(rscriptBinResolved, installScript, [libraryPath], envDir, (output) => {
        const lines = output.split('\n').filter(l => l.trim());
        for (const line of lines) {
            // Track package installation progress
            if (line.includes('installing to')) {
                const match = line.match(/installing to.*?'([^']+)'/);
                if (match) {
                    onProgress('install', 0.5, `Installing: ${path.basename(match[1])}...`);
                }
            } else if (line.includes('downloaded')) {
                onProgress('install', 0.6, line.trim().substring(0, 80));
            }
        }
    });

    // Step 4: Marker file
    fs.writeFileSync(markerPath, JSON.stringify({
        timestamp: new Date().toISOString(),
        rVersion: R_VERSION,
        platform: process.platform,
        arch: process.arch,
    }));

    onProgress('complete', 1, 'R environment setup complete!');
    return rscriptBinResolved;
}

/**
 * Check if the R environment is set up.
 */
function isRSetUp(envDir) {
    const rscriptBin = getRscriptBinPath(envDir);
    const markerPath = path.join(envDir, '.r-setup-complete');
    return fs.existsSync(rscriptBin) && fs.existsSync(markerPath);
}

module.exports = {
    setupREnvironment,
    isRSetUp,
    getRscriptBinPath,
    getRHome,
    runRscript,
};
