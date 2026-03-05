/**
 * Main Renderer Process
 * Initializes the application and handles UI interactions
 */

// Global components
let volcanoPlot = null;
let boxPlot = null;
let scriptPanel = null;

/**
 * Initialize application
 */
function initApp() {
    // Initialize components
    volcanoPlot = new VolcanoPlot('volcano-plot');
    boxPlot = new BoxPlot('boxplot-plot');
    scriptPanel = new ScriptPanel();

    // Setup navigation
    setupNavigation();

    // Setup FASTQ processing panel (Step 0)
    setupFastqPanel();

    // Setup data loading
    setupDataLoading();

    // Setup GEO / recount3 fetching
    setupGeoFetching();

    // Setup DESeq2 panel
    setupDeseqPanel();

    // Setup visualization tabs
    setupVisualizationTabs();

    // Setup selected genes display
    setupSelectedGenesDisplay();

    // Setup citations panel
    setupCitationsPanel();

    // Wait for Python backend
    python.onReady(() => {
        updatePythonStatus(true);
        scriptPanel.showPlaceholder();
    });

    // Show placeholders
    volcanoPlot.showPlaceholder();
    boxPlot.showPlaceholder();
}

/**
 * Update Python connection status
 */
function updatePythonStatus(connected) {
    const statusDot = document.getElementById('status-dot');
    const statusText = document.getElementById('status-text');

    if (connected) {
        statusDot.classList.add('connected');
        statusText.textContent = 'Python Ready';
        store.set('pythonReady', true);
    } else {
        statusDot.classList.remove('connected');
        statusText.textContent = 'Connecting to Python...';
    }
}

/**
 * Setup navigation between panels
 */
function setupNavigation() {
    // Tab-based navigation (new horizontal tabs)
    const navTabs = document.querySelectorAll('.nav-tab');
    navTabs.forEach(tab => {
        tab.addEventListener('click', () => {
            const panelId = tab.dataset.panel;

            // Update tabs
            navTabs.forEach(t => t.classList.remove('active'));
            tab.classList.add('active');

            // Show panel
            document.querySelectorAll('.panel').forEach(p => p.classList.add('hidden'));
            document.getElementById(`panel-${panelId}`)?.classList.remove('hidden');

            store.set('currentPanel', panelId);

            // Update script when navigating to script panel
            if (panelId === 'script' && scriptPanel) {
                scriptPanel.updateScript();
            }
        });
    });

    // Tool selector on Analysis panel
    setupToolSelector();

    // Project folder button
    setupProjectFolder();
}

/**
 * Setup tool selector (pyDESeq2 vs edgePython)
 */
function setupToolSelector() {
    const toolOptions = document.querySelectorAll('.tool-option');
    const toolBarName = document.getElementById('tool-bar-name');
    const btnRunLabel = document.getElementById('btn-run-label');
    const edgeBasicParams = document.getElementById('edgepython-basic-params');
    const edgeAdvancedSection = document.getElementById('edgepython-advanced-section');
    const docPydeseq2 = document.getElementById('doc-pydeseq2');
    const docEdgepython = document.getElementById('doc-edgepython');

    // Store active tool globally
    store.set('activeTool', 'pydeseq2');

    toolOptions.forEach(option => {
        option.addEventListener('click', () => {
            const tool = option.dataset.tool;

            // Update UI
            toolOptions.forEach(o => o.classList.remove('active'));
            option.classList.add('active');

            // Update tool bar indicator
            const toolNames = { pydeseq2: 'pyDESeq2', edgepython: 'edgePython' };
            toolBarName.textContent = toolNames[tool];

            // Update run button label
            btnRunLabel.textContent = `Run ${toolNames[tool]} Analysis`;

            // Show/hide tool-specific params
            if (tool === 'edgepython') {
                edgeBasicParams?.classList.remove('hidden');
                edgeAdvancedSection?.classList.remove('hidden');
                docPydeseq2?.classList.add('hidden');
                docEdgepython?.classList.remove('hidden');
            } else {
                edgeBasicParams?.classList.add('hidden');
                edgeAdvancedSection?.classList.add('hidden');
                docPydeseq2?.classList.remove('hidden');
                docEdgepython?.classList.add('hidden');
            }

            store.set('activeTool', tool);

            // Notify Python
            python.setActiveTool(tool).catch(e => console.error('Failed to set tool:', e));
        });
    });
}

/**
 * Setup project folder browser
 */
function setupProjectFolder() {
    const btnOpenProject = document.getElementById('btn-open-project');
    const projectDrawer = document.getElementById('project-drawer');
    const btnCloseDrawer = document.getElementById('btn-close-drawer');
    const projectPathDisplay = document.getElementById('project-path-display');
    const fileListEl = document.getElementById('project-file-list');

    let projectPath = null;

    btnOpenProject?.addEventListener('click', async () => {
        if (projectPath && !projectDrawer.classList.contains('hidden')) {
            // Toggle drawer closed
            projectDrawer.classList.add('hidden');
            return;
        }

        try {
            const dirPath = await python.openFileDialog({
                properties: ['openDirectory'],
            });
            if (!dirPath) return;

            projectPath = dirPath;
            const folderName = dirPath.split('/').pop() || dirPath;
            projectPathDisplay.textContent = folderName;
            store.set('projectPath', dirPath);

            // List files in the directory
            const result = await python.run('list_project_files', { path: dirPath });
            const files = result.files || [];

            fileListEl.innerHTML = files.length === 0
                ? '<p class="text-muted" style="padding:8px;font-size:12px">No files found</p>'
                : files.map(f => {
                    const icon = f.name.endsWith('.csv') || f.name.endsWith('.tsv') ? '📊'
                               : f.name.endsWith('.fastq') || f.name.endsWith('.fq') || f.name.endsWith('.gz') ? '🧬'
                               : f.name.endsWith('.fa') || f.name.endsWith('.fasta') ? '🔬'
                               : '📄';
                    return `<div class="drawer-file" data-path="${f.path}" data-name="${f.name}">
                        <span class="drawer-file-icon">${icon}</span>
                        <span class="drawer-file-name">${f.name}</span>
                    </div>`;
                }).join('');

            // Click to load files
            fileListEl.querySelectorAll('.drawer-file').forEach(el => {
                el.addEventListener('click', async () => {
                    const filePath = el.dataset.path;
                    const fileName = el.dataset.name;

                    if (fileName.match(/count/i) && fileName.match(/\.(csv|tsv|txt)$/i)) {
                        try {
                            await python.loadCounts(filePath);
                            showDataPreview();
                            updateStatusMessage(`Loaded counts: ${fileName}`);
                        } catch (e) {
                            updateStatusMessage('Error loading counts: ' + e.message, true);
                        }
                    } else if (fileName.match(/meta/i) && fileName.match(/\.(csv|tsv|txt)$/i)) {
                        try {
                            await python.loadMetadata(filePath);
                            showDataPreview();
                            updateStatusMessage(`Loaded metadata: ${fileName}`);
                        } catch (e) {
                            updateStatusMessage('Error loading metadata: ' + e.message, true);
                        }
                    } else {
                        updateStatusMessage(`Opened: ${fileName}`);
                    }
                });
            });

            projectDrawer.classList.remove('hidden');
        } catch (e) {
            console.error('Project folder error:', e);
        }
    });

    btnCloseDrawer?.addEventListener('click', () => {
        projectDrawer.classList.add('hidden');
    });
}

/**
 * Setup data loading functionality
 */
function setupDataLoading() {
    const btnLoadExample = document.getElementById('btn-load-example');
    const btnLoadCounts = document.getElementById('btn-load-counts');
    const btnLoadMetadata = document.getElementById('btn-load-metadata');

    // Load example data
    btnLoadExample?.addEventListener('click', async (e) => {
        e.preventDefault();
        const originalText = btnLoadExample.textContent;
        try {
            btnLoadExample.textContent = 'Loading...';
            btnLoadExample.style.pointerEvents = 'none';

            await python.loadExampleData();
            showDataPreview();
            updateStatusMessage('Example data loaded');
            markStepComplete('data');

        } catch (error) {
            console.error('Failed to load example data:', error);
            updateStatusMessage('Error loading data: ' + error.message, true);
        } finally {
            btnLoadExample.textContent = originalText;
            btnLoadExample.style.pointerEvents = '';
        }
    });

    // Load counts file
    btnLoadCounts?.addEventListener('click', async () => {
        try {
            const filePath = await python.openFileDialog();
            if (filePath) {
                await python.loadCounts(filePath);
                showDataPreview();
                updateStatusMessage('Counts loaded from ' + filePath.split('/').pop());
            }
        } catch (error) {
            console.error('Failed to load counts:', error);
            updateStatusMessage('Error loading counts: ' + error.message, true);
        }
    });

    // Load metadata file
    btnLoadMetadata?.addEventListener('click', async () => {
        try {
            const filePath = await python.openFileDialog();
            if (filePath) {
                await python.loadMetadata(filePath);
                showDataPreview();
                updateStatusMessage('Metadata loaded from ' + filePath.split('/').pop());
            }
        } catch (error) {
            console.error('Failed to load metadata:', error);
            updateStatusMessage('Error loading metadata: ' + error.message, true);
        }
    });
}

// ---------------------------------------------------------------------------
// GEO / recount3 fetching
// ---------------------------------------------------------------------------

let rReady = false;

/**
 * Ensure R environment is installed. Shows overlay if needed.
 * Returns true if R is ready, false if install was triggered (caller should wait).
 */
async function ensureRReady() {
    if (rReady) return true;

    try {
        const status = await window.api.checkRSetup();
        if (status.ready) {
            rReady = true;
            return true;
        }
    } catch (e) {
        console.error('check-r-setup failed:', e);
    }

    // Show R setup overlay
    const overlay = document.getElementById('r-setup-overlay');
    overlay.style.display = 'flex';
    return false;
}

/**
 * Setup GEO / recount3 fetch buttons and R install flow.
 */
function setupGeoFetching() {
    const geoInput = document.getElementById('geo-id-input');
    const btnRecount3 = document.getElementById('btn-fetch-recount3');
    const btnGeoMeta = document.getElementById('btn-fetch-geo-meta');
    const btnGeoCounts = document.getElementById('btn-fetch-geo-counts');
    const btnInstallR = document.getElementById('btn-install-r');
    const overlay = document.getElementById('r-setup-overlay');
    const rProgressDiv = document.getElementById('r-setup-progress');
    const rProgressBar = document.getElementById('r-setup-progress-bar');
    const rStatus = document.getElementById('r-setup-status');
    const fetchStatus = document.getElementById('geo-fetch-status');
    const fetchMessage = document.getElementById('geo-fetch-message');

    if (!geoInput || !btnRecount3 || !btnGeoMeta || !btnGeoCounts) return;

    // --- R setup progress listener ---
    window.api.onRSetupProgress((data) => {
        if (rProgressDiv) rProgressDiv.style.display = 'block';
        if (rProgressBar) rProgressBar.style.width = `${(data.progress || 0) * 100}%`;
        if (rStatus) rStatus.textContent = data.message || data.phase || '';

        if (data.phase === 'complete') {
            rReady = true;
            overlay.style.display = 'none';
            updateStatusMessage('R environment ready');
        }
    });

    // --- Install R button ---
    btnInstallR?.addEventListener('click', async () => {
        btnInstallR.disabled = true;
        btnInstallR.textContent = 'Installing...';
        rProgressDiv.style.display = 'block';

        try {
            const result = await window.api.runRSetup();
            if (result.success) {
                rReady = true;
                overlay.style.display = 'none';
                updateStatusMessage('R environment installed successfully');
            } else {
                rStatus.textContent = 'Installation failed: ' + (result.error || 'Unknown error');
                btnInstallR.disabled = false;
                btnInstallR.textContent = 'Retry Install';
            }
        } catch (err) {
            rStatus.textContent = 'Error: ' + err.message;
            btnInstallR.disabled = false;
            btnInstallR.textContent = 'Retry Install';
        }
    });

    // --- Helper: show/hide fetch spinner ---
    function showFetchSpinner(msg) {
        fetchStatus.style.display = 'flex';
        fetchMessage.textContent = msg;
    }
    function hideFetchSpinner() {
        fetchStatus.style.display = 'none';
    }

    // --- R script progress listener (stderr messages from R) ---
    window.api.onRScriptProgress((data) => {
        if (fetchStatus.style.display !== 'none' && data.message) {
            // Show last meaningful line
            const lines = data.message.split('\n').filter(l => l.trim());
            if (lines.length > 0) {
                const lastLine = lines[lines.length - 1].substring(0, 100);
                fetchMessage.textContent = lastLine;
                // Also update the inline status
                const geoStatus = document.getElementById('geo-tool-status');
                if (geoStatus && geoStatus.classList.contains('loading')) {
                    geoStatus.querySelector('.status-text').textContent = lastLine;
                }
            }
        }
    });

    // --- Fetch from recount3 ---
    btnRecount3.addEventListener('click', async () => {
        const geoId = geoInput.value.trim();
        if (!geoId) {
            updateStatusMessage('Please enter a GEO ID (GSE...) or SRA ID (SRP...)', true);
            geoInput.focus();
            return;
        }

        if (!(await ensureRReady())) return; // R not ready, overlay shown

        btnRecount3.disabled = true;
        btnGeoMeta.disabled = true;
        btnGeoCounts.disabled = true;
        showFetchSpinner(`Fetching recount3 data for ${geoId}...`);
        showToolStatus('geo-tool-status', 'loading', `Fetching recount3 data for ${geoId}...`);

        try {
            const result = await window.api.runRScript('fetch_recount3.R', [geoId]);

            if (result.error) {
                throw new Error(result.error);
            }

            // Load counts (always from recount3)
            if (result.counts_path) {
                await python.loadCounts(result.counts_path);
            }
            // Load metadata only if no metadata is already loaded (user may have fetched from GEO separately)
            const hasExistingMeta = !!store.get('metadata');
            if (result.metadata_path && !hasExistingMeta) {
                await python.loadMetadata(result.metadata_path);
            }

            showDataPreview();
            markStepComplete('data');
            let msg = `recount3 counts loaded for ${geoId}: ${result.counts_shape[0]} genes × ${result.counts_shape[1]} samples`;
            if (hasExistingMeta) {
                msg += ' (keeping existing metadata)';
            }
            updateStatusMessage(msg);
            showToolStatus('geo-tool-status', 'success', msg);

        } catch (err) {
            console.error('recount3 fetch failed:', err);
            updateStatusMessage('recount3 fetch failed: ' + err.message, true);
            showToolStatus('geo-tool-status', 'error', 'recount3 fetch failed: ' + err.message);
        } finally {
            btnRecount3.disabled = false;
            btnGeoMeta.disabled = false;
            btnGeoCounts.disabled = false;
            hideFetchSpinner();
        }
    });

    // --- Fetch GEO metadata ---
    btnGeoMeta.addEventListener('click', async () => {
        let geoId = geoInput.value.trim();
        if (!geoId) {
            updateStatusMessage('Please enter a GEO ID (GSE...)', true);
            geoInput.focus();
            return;
        }

        // GEO metadata only works with GSE IDs
        if (!geoId.toUpperCase().startsWith('GSE')) {
            updateStatusMessage('GEO metadata requires a GSE ID (e.g. GSE123456)', true);
            return;
        }

        if (!(await ensureRReady())) return;

        btnRecount3.disabled = true;
        btnGeoMeta.disabled = true;
        btnGeoCounts.disabled = true;
        showFetchSpinner(`Fetching GEO metadata for ${geoId}...`);
        showToolStatus('geo-tool-status', 'loading', `Fetching GEO metadata for ${geoId}...`);

        try {
            const result = await window.api.runRScript('fetch_geo_metadata.R', [geoId]);

            if (result.error) {
                throw new Error(result.error);
            }

            // Check if counts are already loaded (e.g. from recount3 with SRR IDs)
            // If so, merge GEO metadata into existing metadata instead of replacing it
            const existingMeta = store.get('metadata');
            const existingCounts = store.get('counts');
            let msg;

            if (existingMeta && existingCounts && result.metadata_path) {
                // Check if existing metadata has sra.sample_title (recount3 metadata)
                const metaCols = store.get('metadataColumns') || [];
                const hasSampleTitle = metaCols.includes('sra.sample_title');

                if (hasSampleTitle) {
                    // Merge GEO metadata into existing recount3 metadata via title
                    // This expands GEO metadata (1 per GSM) to every SRR technical replicate
                    const mergeResult = await python.run('merge_metadata', {
                        new_meta_path: result.metadata_path,
                        existing_key: 'sra.sample_title',
                        new_key: 'title',
                    });
                    msg = `GEO metadata merged: ${mergeResult.matched} of ${mergeResult.total} samples matched via title`;
                    if (mergeResult.new_columns && mergeResult.new_columns.length > 0) {
                        msg += ` (added ${mergeResult.new_columns.length} columns: ${mergeResult.new_columns.slice(0, 5).join(', ')})`;
                    }
                    // Update store with merged metadata
                    store.update({
                        metadata: mergeResult.metadata_preview,
                        metadataShape: mergeResult.metadata_shape,
                        metadataColumns: mergeResult.metadata_columns,
                        metadataIndex: mergeResult.metadata_index,
                    });
                } else {
                    // No recount3 metadata — just load as normal
                    await python.loadMetadata(result.metadata_path);
                    msg = `GEO metadata loaded for ${geoId}`;
                }
            } else {
                // No existing data, just load metadata directly
                if (result.metadata_path) {
                    await python.loadMetadata(result.metadata_path);
                }
                msg = `GEO metadata loaded for ${geoId}`;
                if (result.expression_available) {
                    msg += ' (expression data also available — use GEO Counts to fetch it separately)';
                } else {
                    msg += ' (use recount3 or GEO Counts for count data)';
                }
            }

            showDataPreview();
            markStepComplete('data');
            updateStatusMessage(msg);
            showToolStatus('geo-tool-status', 'success', msg);

        } catch (err) {
            console.error('GEO fetch failed:', err);
            updateStatusMessage('GEO fetch failed: ' + err.message, true);
            showToolStatus('geo-tool-status', 'error', 'GEO metadata fetch failed: ' + err.message);
        } finally {
            btnRecount3.disabled = false;
            btnGeoMeta.disabled = false;
            btnGeoCounts.disabled = false;
            hideFetchSpinner();
        }
    });

    // --- Fetch GEO Counts ---
    btnGeoCounts.addEventListener('click', async () => {
        let geoId = geoInput.value.trim();
        if (!geoId) {
            updateStatusMessage('Please enter a GEO ID (GSE...)', true);
            geoInput.focus();
            return;
        }

        if (!geoId.toUpperCase().startsWith('GSE')) {
            updateStatusMessage('GEO counts requires a GSE ID (e.g. GSE123456)', true);
            return;
        }

        if (!(await ensureRReady())) return;

        btnRecount3.disabled = true;
        btnGeoMeta.disabled = true;
        btnGeoCounts.disabled = true;
        showFetchSpinner(`Fetching GEO counts for ${geoId}...`);
        showToolStatus('geo-tool-status', 'loading', `Fetching GEO counts for ${geoId}...`);

        try {
            const result = await window.api.runRScript('fetch_geo_counts.R', [geoId]);

            if (result.error) {
                throw new Error(result.error);
            }

            // Load counts only
            let msg = '';
            if (result.counts_path) {
                await python.loadCounts(result.counts_path);
                const src = result.counts_source === 'geo_rnaseq_counts' ? 'GEO RNA-seq counts'
                          : result.counts_source === 'supplementary' ? 'supplementary files'
                          : 'series matrix';
                msg = `GEO counts loaded from ${src}: ${result.counts_shape[0]} genes × ${result.counts_shape[1]} samples`;
            }

            showDataPreview();
            markStepComplete('data');
            updateStatusMessage(msg);
            showToolStatus('geo-tool-status', 'success', msg);

        } catch (err) {
            console.error('GEO counts fetch failed:', err);
            updateStatusMessage('GEO counts fetch failed: ' + err.message, true);
            showToolStatus('geo-tool-status', 'error', 'GEO counts fetch failed: ' + err.message);
        } finally {
            btnRecount3.disabled = false;
            btnGeoMeta.disabled = false;
            btnGeoCounts.disabled = false;
            hideFetchSpinner();
        }
    });
}

/**
 * Show data preview tables using Tabulator
 */
function showDataPreview() {
    const previewSection = document.getElementById('data-preview');
    previewSection.style.display = 'block';

    // Show counts info
    const countsInfo = document.getElementById('counts-info');
    const countsShape = store.get('countsShape');
    if (countsShape) {
        countsInfo.textContent = `${countsShape[0]} genes × ${countsShape[1]} samples`;
    }

    // Show counts preview table with Tabulator
    const counts = store.get('counts');
    if (counts) {
        renderCountsPreviewTable(counts);
    }

    // Show metadata info
    const metadataInfo = document.getElementById('metadata-info');
    const metadataShape = store.get('metadataShape');
    if (metadataShape) {
        metadataInfo.textContent = `${metadataShape[0]} samples × ${metadataShape[1]} columns`;
    }

    // Show metadata preview table with Tabulator
    const metadata = store.get('metadata');
    if (metadata) {
        renderMetadataPreviewTable(metadata);
    }

    // Populate column selectors
    populateColumnSelectors();

    // Check sample matching if both are loaded
    checkSampleMatch();
}

async function checkSampleMatch() {
    const counts = store.get('counts');
    const metadata = store.get('metadata');
    if (!counts || !metadata) {
        hideToolStatus('sample-match-status');
        document.getElementById('sample-match-action').style.display = 'none';
        return;
    }

    try {
        const result = await python.run('check_sample_match', {});
        const actionDiv = document.getElementById('sample-match-action');

        if (result.status === 'match') {
            showToolStatus('sample-match-status', 'success', result.message);
            actionDiv.style.display = 'none';
        } else if (result.status === 'partial') {
            showToolStatus('sample-match-status', 'error', result.message);
            if (result.suggested_column) {
                actionDiv.style.display = 'flex';
                document.getElementById('sample-match-suggestion').textContent =
                    `Column "${result.suggested_column}" matches ${result.suggested_overlap} samples.`;
                document.getElementById('btn-remap-samples').onclick = async () => {
                    await python.run('remap_metadata_index', { column: result.suggested_column });
                    // Refresh metadata display
                    const metaResult = await python.run('check_sample_match', {});
                    showDataPreview();
                    if (metaResult.status === 'match') {
                        showToolStatus('sample-match-status', 'success', metaResult.message);
                        actionDiv.style.display = 'none';
                    }
                };
            } else {
                actionDiv.style.display = 'none';
            }
        } else {
            // mismatch
            showToolStatus('sample-match-status', 'error', result.message);
            if (result.suggested_column) {
                actionDiv.style.display = 'flex';
                document.getElementById('sample-match-suggestion').textContent =
                    `Column "${result.suggested_column}" matches ${result.suggested_overlap} samples.`;
                document.getElementById('btn-remap-samples').onclick = async () => {
                    await python.run('remap_metadata_index', { column: result.suggested_column });
                    showDataPreview();
                };
            } else {
                actionDiv.style.display = 'none';
            }
        }
    } catch (err) {
        console.error('Sample match check failed:', err);
    }
}

/**
 * Render counts preview table using Tabulator
 */
function renderCountsPreviewTable(data) {
    const columns = Object.keys(data);
    if (columns.length === 0) return;

    const rowKeys = Object.keys(data[columns[0]]);

    // Transform data to array of objects
    const tableData = rowKeys.slice(0, 100).map(rowKey => {
        const row = { _gene: rowKey };
        columns.forEach(col => {
            const value = data[col][rowKey];
            row[col] = typeof value === 'number' ?
                (Number.isInteger(value) ? value : parseFloat(value.toFixed(2))) : value;
        });
        return row;
    });

    // Build columns config
    const tabulatorColumns = [
        {
            title: 'Gene',
            field: '_gene',
            frozen: true,
            headerFilter: 'input',
            headerFilterPlaceholder: 'Search...',
            width: 100,
            formatter: function(cell) {
                return `<strong>${cell.getValue()}</strong>`;
            }
        }
    ];

    columns.forEach(col => {
        tabulatorColumns.push({
            title: col,
            field: col,
            sorter: 'number',
            hozAlign: 'right',
            formatter: function(cell) {
                const val = cell.getValue();
                return val != null ? val.toLocaleString() : '-';
            }
        });
    });

    // Destroy existing table
    if (countsPreviewTable) {
        countsPreviewTable.destroy();
    }

    countsPreviewTable = new Tabulator('#counts-preview', {
        data: tableData,
        height: '250px',
        layout: 'fitDataFill',
        movableColumns: true,
        placeholder: 'No data loaded',
        columns: tabulatorColumns
    });
}

/**
 * Render metadata preview table using Tabulator
 */
function renderMetadataPreviewTable(data) {
    const columns = Object.keys(data);
    if (columns.length === 0) return;

    const rowKeys = Object.keys(data[columns[0]]);

    // Transform data to array of objects
    const tableData = rowKeys.map(rowKey => {
        const row = { _sample: rowKey };
        columns.forEach(col => {
            row[col] = data[col][rowKey];
        });
        return row;
    });

    // Build columns config
    const tabulatorColumns = [
        {
            title: 'Sample',
            field: '_sample',
            frozen: true,
            headerFilter: 'input',
            headerFilterPlaceholder: 'Search...',
            width: 100,
            formatter: function(cell) {
                return `<strong>${cell.getValue()}</strong>`;
            }
        }
    ];

    columns.forEach(col => {
        tabulatorColumns.push({
            title: col,
            field: col,
            headerFilter: 'input',
            headerFilterPlaceholder: 'Filter...'
        });
    });

    // Destroy existing table
    if (metadataPreviewTable) {
        metadataPreviewTable.destroy();
    }

    metadataPreviewTable = new Tabulator('#metadata-preview', {
        data: tableData,
        height: '200px',
        layout: 'fitColumns',
        movableColumns: true,
        placeholder: 'No data loaded',
        columns: tabulatorColumns
    });
}

/**
 * Populate condition column and reference level selectors
 */
function populateColumnSelectors() {
    const conditionSelect = document.getElementById('condition-column');
    const referenceSelect = document.getElementById('reference-level');
    const metadataColumns = store.get('metadataColumns') || [];
    const metadata = store.get('metadata');

    // Clear existing options
    conditionSelect.innerHTML = '';
    referenceSelect.innerHTML = '';

    // Add column options
    metadataColumns.forEach(col => {
        const option = document.createElement('option');
        option.value = col;
        option.textContent = col;
        if (col === 'condition') option.selected = true;
        conditionSelect.appendChild(option);
    });

    // Update reference levels when condition changes
    const updateReferenceLevels = () => {
        const selectedColumn = conditionSelect.value;
        referenceSelect.innerHTML = '';

        if (metadata && metadata[selectedColumn]) {
            const uniqueValues = [...new Set(Object.values(metadata[selectedColumn]))];
            uniqueValues.forEach(val => {
                const option = document.createElement('option');
                option.value = val;
                option.textContent = val;
                if (val === 'control') option.selected = true;
                referenceSelect.appendChild(option);
            });
        }

        // Update store
        store.set('params.conditionColumn', selectedColumn);
        store.set('params.referenceLevel', referenceSelect.value);
    };

    conditionSelect.addEventListener('change', updateReferenceLevels);
    referenceSelect.addEventListener('change', () => {
        store.set('params.referenceLevel', referenceSelect.value);
    });

    updateReferenceLevels();
}

/**
 * Setup DESeq2 panel
 */
function setupDeseqPanel() {
    const btnRunDeseq = document.getElementById('btn-run-deseq');
    const lfcInput = document.getElementById('lfc-threshold');
    const padjInput = document.getElementById('padj-threshold');
    const deseqStatus = document.getElementById('deseq-status');
    const deseqSummary = document.getElementById('deseq-summary');
    const resultsSection = document.getElementById('deseq-results-section');

    btnRunDeseq?.addEventListener('click', async () => {
        if (!store.hasData()) {
            updateStatusMessage('Load data first', true);
            return;
        }

        const activeTool = store.get('activeTool') || 'pydeseq2';
        const toolNames = { pydeseq2: 'pyDESeq2', edgepython: 'edgePython' };

        try {
            btnRunDeseq.disabled = true;
            deseqStatus.classList.remove('hidden');
            deseqSummary.classList.add('hidden');
            resultsSection.classList.add('hidden');
            showToolStatus('analysis-tool-status', 'loading', `Running ${toolNames[activeTool]} analysis...`);

            // Shared parameters
            const params = {
                condition_column: store.get('params.conditionColumn'),
                reference_level: store.get('params.referenceLevel'),
                lfc_threshold: parseFloat(lfcInput.value),
                padj_threshold: parseFloat(padjInput.value)
            };

            // Store params
            store.update({
                'params.lfcThreshold': params.lfc_threshold,
                'params.padjThreshold': params.padj_threshold
            });

            // Shared min count filter
            params.min_count = parseInt(document.getElementById('deseq-min-count')?.value || '10');

            let result;

            if (activeTool === 'edgepython') {
                // Gather edgePython-specific params
                params.normalization = document.getElementById('edge-normalization')?.value || 'TMM';
                params.test_method = document.getElementById('edge-test-method')?.value || 'qlf';
                params.dispersion = document.getElementById('edge-dispersion')?.value || 'trended';
                params.robust = document.getElementById('edge-robust')?.value || 'true';
                params.min_count = parseInt(document.getElementById('edge-min-count')?.value || '10');

                result = await python.runEdgePython(params);
            } else {
                result = await python.runDeseq(params);
            }

            // Show summary
            document.getElementById('total-genes').textContent = result.summary.total_genes;
            document.getElementById('sig-genes').textContent = result.summary.significant;
            document.getElementById('up-genes').textContent = result.summary.upregulated;
            document.getElementById('down-genes').textContent = result.summary.downregulated;

            deseqSummary.classList.remove('hidden');
            resultsSection.classList.remove('hidden');

            // Render results table
            renderResultsTable(result.results);

            // Update status
            const filteredMsg = result.summary.genes_filtered ? ` (${result.summary.genes_filtered} low-count genes filtered)` : '';
            const successMsg = `${toolNames[activeTool]} analysis complete: ${result.summary.significant} significant genes${filteredMsg}`;
            updateStatusMessage(successMsg);
            showToolStatus('analysis-tool-status', 'success', successMsg);
            document.getElementById('analysis-status').textContent = `Analysis: Complete (${toolNames[activeTool]})`;
            markStepComplete('deseq');

            // Auto-generate visualizations based on current modes
            renderForMode('volcano').catch(e => console.error('Auto volcano:', e));
            renderForMode('heatmap').catch(e => console.error('Auto heatmap:', e));
            renderForMode('pca').catch(e => console.error('Auto PCA:', e));

            // Update script
            scriptPanel.updateScript();

        } catch (error) {
            console.error(`${toolNames[activeTool]} analysis failed:`, error);
            updateStatusMessage('Analysis failed: ' + error.message, true);
            showToolStatus('analysis-tool-status', 'error', `${toolNames[activeTool]} analysis failed: ${error.message}`);
        } finally {
            btnRunDeseq.disabled = false;
            deseqStatus.classList.add('hidden');
        }
    });

    // Setup results filtering
    setupResultsFiltering();
}

// Global Tabulator instances
let deseqTable = null;
let gseaVizTable = null;
let countsPreviewTable = null;
let metadataPreviewTable = null;

/**
 * Render DESeq2 results table using Tabulator
 */
function renderResultsTable(results) {
    const container = document.getElementById('deseq-results-table');
    const lfcThresh = store.get('params.lfcThreshold');
    const padjThresh = store.get('params.padjThreshold');

    // Format data for Tabulator
    const tableData = results.map(row => ({
        gene: row.index || row.gene,
        baseMean: row.baseMean,
        log2FoldChange: row.log2FoldChange,
        lfcSE: row.lfcSE,
        pvalue: row.pvalue,
        padj: row.padj,
        _isUp: row.log2FoldChange > lfcThresh && row.padj < padjThresh,
        _isDown: row.log2FoldChange < -lfcThresh && row.padj < padjThresh
    }));

    // Destroy existing table if any
    if (deseqTable) {
        deseqTable.destroy();
    }

    // Create Tabulator table
    deseqTable = new Tabulator('#deseq-results-table', {
        data: tableData,
        height: '500px',
        layout: 'fitColumns',
        pagination: true,
        paginationSize: 50,
        paginationSizeSelector: [25, 50, 100, 500],
        movableColumns: true,
        selectable: true,
        placeholder: 'No Data Available',
        columns: [
            {
                title: 'Gene',
                field: 'gene',
                headerFilter: 'input',
                headerFilterPlaceholder: 'Search...',
                frozen: true,
                width: 120,
                formatter: function(cell) {
                    const data = cell.getRow().getData();
                    let color = '#4A5568';
                    if (data._isUp) color = '#C97B7B';
                    if (data._isDown) color = '#7BA3C9';
                    return `<strong style="color: ${color}">${cell.getValue()}</strong>`;
                }
            },
            {
                title: 'baseMean',
                field: 'baseMean',
                headerFilter: 'number',
                headerFilterPlaceholder: 'min',
                headerFilterFunc: '>=',
                sorter: 'number',
                formatter: function(cell) {
                    const val = cell.getValue();
                    return val != null ? val.toFixed(2) : '-';
                }
            },
            {
                title: 'log2FC',
                field: 'log2FoldChange',
                headerFilter: 'number',
                headerFilterPlaceholder: 'min abs',
                headerFilterFunc: function(headerValue, rowValue) {
                    if (!headerValue) return true;
                    return Math.abs(rowValue) >= Math.abs(parseFloat(headerValue));
                },
                sorter: 'number',
                formatter: function(cell) {
                    const val = cell.getValue();
                    if (val == null) return '-';
                    const color = val > 0 ? '#C97B7B' : '#7BA3C9';
                    return `<span style="color: ${color}; font-weight: bold">${val.toFixed(3)}</span>`;
                }
            },
            {
                title: 'lfcSE',
                field: 'lfcSE',
                sorter: 'number',
                formatter: function(cell) {
                    const val = cell.getValue();
                    return val != null ? val.toFixed(3) : '-';
                }
            },
            {
                title: 'pvalue',
                field: 'pvalue',
                headerFilter: 'number',
                headerFilterPlaceholder: 'max',
                headerFilterFunc: '<=',
                sorter: 'number',
                formatter: function(cell) {
                    const val = cell.getValue();
                    return val != null ? val.toExponential(2) : '-';
                }
            },
            {
                title: 'padj',
                field: 'padj',
                headerFilter: 'number',
                headerFilterPlaceholder: 'max',
                headerFilterFunc: '<=',
                sorter: 'number',
                formatter: function(cell) {
                    const val = cell.getValue();
                    if (val == null) return '-';
                    const padjThresh = store.get('params.padjThreshold');
                    const style = val < padjThresh ? 'font-weight: bold; color: #8FB996' : '';
                    return `<span style="${style}">${val.toExponential(2)}</span>`;
                }
            }
        ],
        rowClick: function(e, row) {
            const gene = row.getData().gene;
            if (gene) {
                store.toggleSelectedGene(gene);
                updateStatusMessage(`Toggled gene: ${gene}`);
            }
        },
        rowFormatter: function(row) {
            const data = row.getData();
            if (data._isUp) {
                row.getElement().style.backgroundColor = 'rgba(201, 123, 123, 0.12)';
            } else if (data._isDown) {
                row.getElement().style.backgroundColor = 'rgba(123, 163, 201, 0.12)';
            }
        }
    });
}

/**
 * Setup results table filtering
 */
function setupResultsFiltering() {
    const filterSelect = document.getElementById('filter-significance');
    const exportCsvBtn = document.getElementById('btn-export-csv');
    const exportSelectedBtn = document.getElementById('btn-export-selected');

    const filterResults = () => {
        const results = store.get('deseqResults') || [];
        const filterValue = filterSelect?.value || 'all';
        const lfcThresh = store.get('params.lfcThreshold');
        const padjThresh = store.get('params.padjThreshold');

        let filtered = results.filter(row => {
            if (filterValue === 'significant') {
                return Math.abs(row.log2FoldChange) > lfcThresh && row.padj < padjThresh;
            } else if (filterValue === 'up') {
                return row.log2FoldChange > lfcThresh && row.padj < padjThresh;
            } else if (filterValue === 'down') {
                return row.log2FoldChange < -lfcThresh && row.padj < padjThresh;
            }
            return true;
        });

        renderResultsTable(filtered);
    };

    filterSelect?.addEventListener('change', filterResults);

    // Export all to CSV
    exportCsvBtn?.addEventListener('click', () => {
        if (deseqTable) {
            deseqTable.download('csv', 'deseq2_results.csv');
        }
    });

    // Export selected rows
    exportSelectedBtn?.addEventListener('click', () => {
        if (deseqTable) {
            const selectedData = deseqTable.getSelectedData();
            if (selectedData.length === 0) {
                updateStatusMessage('No rows selected. Click rows to select them.', true);
                return;
            }
            // Create a temporary table with selected data
            const tempDiv = document.createElement('div');
            tempDiv.style.display = 'none';
            document.body.appendChild(tempDiv);
            const tempTable = new Tabulator(tempDiv, {
                data: selectedData,
                autoColumns: true
            });
            tempTable.download('csv', 'deseq2_selected.csv');
            tempTable.destroy();
            document.body.removeChild(tempDiv);
        }
    });
}

/**
 * Get current plot mode for a viz type
 */
function getPlotMode(vizType) {
    const switchEl = document.querySelector(`.plot-mode-switch[data-viz="${vizType}"] .mode-btn.active`);
    return switchEl ? switchEl.dataset.mode : 'plotly';
}

/**
 * Render the appropriate plots based on mode for a given viz type
 */
async function renderForMode(vizType) {
    const mode = getPlotMode(vizType);
    const needsPlotly = (mode === 'plotly' || mode === 'both');
    const needsPython = (mode === 'python' || mode === 'both');

    if (!store.hasResults()) return;

    if (vizType === 'volcano') {
        if (needsPlotly) await volcanoPlot.render().catch(e => console.error('Plotly volcano:', e));
        if (needsPython) await generateMplVolcano().catch(e => console.error('Mpl volcano:', e));
    } else if (vizType === 'boxplot') {
        if (needsPlotly) await boxPlot.render().catch(e => console.error('Plotly boxplot:', e));
        if (needsPython) {
            const sg = store.get('selectedGenes');
            if (sg && sg.length > 0) await generateMplBoxplot().catch(e => console.error('Mpl boxplot:', e));
        }
    } else if (vizType === 'heatmap') {
        if (needsPlotly) await generatePlotlyHeatmap().catch(e => console.error('Plotly heatmap:', e));
        if (needsPython) await generateMplHeatmap().catch(e => console.error('Mpl heatmap:', e));
    } else if (vizType === 'pca') {
        if (needsPlotly) await generatePlotlyPca().catch(e => console.error('Plotly PCA:', e));
        if (needsPython) await generateMplPca().catch(e => console.error('Mpl PCA:', e));
    } else if (vizType === 'gsea') {
        const gseaResults = store.get('gseaResults');
        if (gseaResults && gseaResults.length > 0) {
            document.getElementById('gsea-viz-results')?.classList.remove('hidden');
            const topN = parseInt(document.getElementById('gsea-viz-topn')?.value || 15);
            if (needsPlotly) await generatePlotlyGsea(gseaResults, topN).catch(e => console.error('Plotly GSEA:', e));
            if (needsPython) await generateMplGseaViz(topN).catch(e => console.error('Mpl GSEA:', e));
        }
    }
}

/**
 * Setup plot mode switches for all viz panels
 */
function setupPlotModeSwitches() {
    document.querySelectorAll('.plot-mode-switch').forEach(switchEl => {
        const vizType = switchEl.dataset.viz;
        switchEl.querySelectorAll('.mode-btn').forEach(btn => {
            btn.addEventListener('click', async () => {
                // Update active button
                switchEl.querySelectorAll('.mode-btn').forEach(b => b.classList.remove('active'));
                btn.classList.add('active');

                // Update display container class
                const displayId = {
                    volcano: 'volcano-plot-display',
                    boxplot: 'boxplot-plot-display',
                    heatmap: 'heatmap-plot-display',
                    pca: 'pca-plot-display',
                    gsea: 'gsea-plot-display',
                }[vizType];
                const displayEl = document.getElementById(displayId);
                if (displayEl) {
                    displayEl.classList.remove('show-plotly', 'show-python', 'show-both');
                    displayEl.classList.add(`show-${btn.dataset.mode}`);
                }

                // Generate the appropriate plots
                await renderForMode(vizType);
            });
        });
    });
}

/**
 * Setup visualization tabs
 */
function setupVisualizationTabs() {
    const tabs = document.querySelectorAll('.viz-tab');
    const panels = document.querySelectorAll('.viz-panel');

    // Setup plot mode switches
    setupPlotModeSwitches();

    tabs.forEach(tab => {
        tab.addEventListener('click', async () => {
            const vizType = tab.dataset.viz;

            // Update tabs
            tabs.forEach(t => t.classList.remove('active'));
            tab.classList.add('active');

            // Update panels
            panels.forEach(p => p.classList.add('hidden'));
            document.getElementById(`viz-${vizType}`)?.classList.remove('hidden');

            // Render the right plots for this tab's current mode
            await renderForMode(vizType);
        });
    });

    // Volcano plot controls
    document.getElementById('btn-update-volcano')?.addEventListener('click', async () => {
        const lfc = parseFloat(document.getElementById('volcano-lfc').value);
        const padj = parseFloat(document.getElementById('volcano-padj').value);
        await volcanoPlot.updateThresholds(lfc, padj);
        // Also update mpl if in python/both mode
        const mode = getPlotMode('volcano');
        if (mode === 'python' || mode === 'both') {
            await generateMplVolcano();
        }
    });

    // Boxplot controls
    document.getElementById('btn-update-boxplot')?.addEventListener('click', async () => {
        const norm = document.getElementById('boxplot-norm').value;
        await boxPlot.updateNormalization(norm);
    });

    // Heatmap controls
    document.getElementById('btn-update-heatmap')?.addEventListener('click', async () => {
        await renderForMode('heatmap');
    });

    // PCA plot controls - auto-update on any change
    const pcaColorBy = document.getElementById('pca-color-by');
    const pcaNormalization = document.getElementById('pca-normalization');
    const pcaCenter = document.getElementById('pca-center');
    const pcaScale = document.getElementById('pca-scale');
    const pcaShowLabels = document.getElementById('pca-show-labels');
    const pcaShowLegend = document.getElementById('pca-show-legend');

    const updatePca = async () => {
        if (store.hasResults()) {
            await renderForMode('pca');
        }
    };

    // Auto-update PCA on any option change
    pcaColorBy?.addEventListener('change', updatePca);
    pcaNormalization?.addEventListener('change', updatePca);
    pcaCenter?.addEventListener('change', updatePca);
    pcaScale?.addEventListener('change', updatePca);
    pcaShowLabels?.addEventListener('change', updatePca);
    pcaShowLegend?.addEventListener('change', updatePca);
    document.getElementById('pca-use-selected-genes')?.addEventListener('change', updatePca);

    // Include checkboxes
    document.getElementById('include-volcano')?.addEventListener('change', (e) => {
        store.setInclude('volcano', e.target.checked);
    });
    document.getElementById('include-boxplot')?.addEventListener('change', (e) => {
        store.setInclude('boxplot', e.target.checked);
    });
    document.getElementById('include-heatmap')?.addEventListener('change', (e) => {
        store.setInclude('heatmap', e.target.checked);
    });
    document.getElementById('include-pca')?.addEventListener('change', (e) => {
        store.setInclude('pca', e.target.checked);
    });
    document.getElementById('include-gsea-viz')?.addEventListener('change', (e) => {
        store.setInclude('gsea', e.target.checked);
    });

    // GSEA in visualizations tab
    document.getElementById('btn-run-gsea-viz')?.addEventListener('click', async () => {
        if (!store.hasResults()) {
            updateStatusMessage('Run DESeq2 analysis first', true);
            return;
        }

        const gseaStatus = document.getElementById('gsea-viz-status');
        const gseaResults = document.getElementById('gsea-viz-results');

        try {
            gseaStatus.classList.remove('hidden');
            gseaResults.classList.add('hidden');

            const geneSets = document.getElementById('gsea-gene-set-library').value;
            const method = document.getElementById('gsea-method-select').value;

            const result = await python.runGsea({ geneSets, method });

            gseaResults.classList.remove('hidden');
            store.set('gseaResults', result.results);

            // Generate plots based on current mode
            const topN = parseInt(document.getElementById('gsea-viz-topn')?.value || 15);
            const mode = getPlotMode('gsea');
            if (mode === 'plotly' || mode === 'both') {
                await generatePlotlyGsea(result.results, topN);
            }
            if (mode === 'python' || mode === 'both') {
                await generateMplGseaViz(topN);
            }

            // Render GSEA results table
            renderGseaVizTable(result.results);

            updateStatusMessage('Pathway analysis complete');

        } catch (error) {
            console.error('GSEA failed:', error);
            updateStatusMessage('Pathway analysis failed: ' + error.message, true);
        } finally {
            gseaStatus.classList.add('hidden');
        }
    });

    // Refresh GSEA plot in viz tab
    document.getElementById('btn-refresh-gsea-viz-plot')?.addEventListener('click', async () => {
        const topN = parseInt(document.getElementById('gsea-viz-topn')?.value || 15);
        await renderForMode('gsea');
    });

    // Setup environment export panel
    setupEnvironmentExport();
}

/**
 * Generate matplotlib volcano plot
 */
async function generateMplVolcano() {
    const container = document.getElementById('volcano-plot-mpl');
    if (!container) return;

    try {
        container.innerHTML = '<div class="mpl-placeholder">Generating plot...</div>';
        const lfc = parseFloat(document.getElementById('volcano-lfc').value);
        const padj = parseFloat(document.getElementById('volcano-padj').value);

        const result = await python.generateVolcanoPlot({
            lfcThreshold: lfc,
            padjThreshold: padj,
            highlightGenes: store.get('selectedGenes'),
            format: 'svg'
        });

        if (result.image) {
            container.innerHTML = result.image;
        } else if (result.error) {
            container.innerHTML = `<div class="mpl-placeholder">Error: ${result.error}</div>`;
        }
    } catch (error) {
        console.error('Error generating matplotlib volcano:', error);
        container.innerHTML = `<div class="mpl-placeholder">Error: ${error.message}</div>`;
    }
}

/**
 * Generate matplotlib boxplot
 */
async function generateMplBoxplot() {
    const container = document.getElementById('boxplot-plot-mpl');
    if (!container) return;

    const genes = store.get('selectedGenes');
    if (!genes || genes.length === 0) {
        container.innerHTML = '<div class="mpl-placeholder">Select genes from the volcano plot first</div>';
        return;
    }

    try {
        container.innerHTML = '<div class="mpl-placeholder">Generating plot...</div>';
        const norm = document.getElementById('boxplot-norm').value;

        const result = await python.generateBoxplot(genes, {
            normalization: norm,
            format: 'svg'
        });

        if (result.image) {
            container.innerHTML = result.image;
        } else if (result.error) {
            container.innerHTML = `<div class="mpl-placeholder">Error: ${result.error}</div>`;
        }
    } catch (error) {
        console.error('Error generating matplotlib boxplot:', error);
        container.innerHTML = `<div class="mpl-placeholder">Error: ${error.message}</div>`;
    }
}

/**
 * Generate matplotlib heatmap
 */
async function generateMplHeatmap() {
    const container = document.getElementById('heatmap-plot-mpl');
    if (!container) return;

    try {
        container.innerHTML = '<div class="mpl-placeholder">Generating heatmap...</div>';

        // Get all heatmap options
        const geneSelection = document.getElementById('heatmap-gene-selection')?.value || 'variance';
        const topN = parseInt(document.getElementById('heatmap-topn')?.value || 50);
        const scale = document.getElementById('heatmap-scale')?.value || 'row';
        const clusterCols = document.getElementById('heatmap-cluster-cols')?.checked ?? true;
        const clusterRows = document.getElementById('heatmap-cluster-rows')?.checked ?? true;
        const showDendrogram = document.getElementById('heatmap-show-dendrogram')?.checked ?? true;

        // Use selected genes if that option is chosen
        const selectedGenes = geneSelection === 'selected' ? store.get('selectedGenes') : null;
        if (geneSelection === 'selected' && (!selectedGenes || selectedGenes.length === 0)) {
            container.innerHTML = '<div class="mpl-placeholder">No genes selected. Click genes on the volcano plot first.</div>';
            return;
        }

        const result = await python.generateHeatmapPlot({
            genes: selectedGenes,
            topN: topN,
            scale: scale,
            geneSelection: geneSelection,
            clusterCols: clusterCols,
            clusterRows: clusterRows,
            showDendrogram: showDendrogram,
            format: 'svg'
        });

        if (result.image) {
            container.innerHTML = result.image;
        } else if (result.error) {
            container.innerHTML = `<div class="mpl-placeholder">Error: ${result.error}</div>`;
        }
    } catch (error) {
        console.error('Error generating matplotlib heatmap:', error);
        container.innerHTML = `<div class="mpl-placeholder">Error: ${error.message}</div>`;
    }
}

/**
 * Generate matplotlib PCA plot
 */
async function generateMplPca() {
    const container = document.getElementById('pca-plot-mpl');
    const infoContainer = document.getElementById('pca-variance-info');
    if (!container) return;

    try {
        container.innerHTML = '<div class="mpl-placeholder">Generating PCA plot...</div>';

        // Get all PCA options
        const colorBy = document.getElementById('pca-color-by')?.value || 'condition';
        const normalization = document.getElementById('pca-normalization')?.value || 'log2';
        const center = document.getElementById('pca-center')?.checked ?? true;
        const scale = document.getElementById('pca-scale')?.checked ?? true;
        const showLabels = document.getElementById('pca-show-labels')?.checked ?? true;
        const showLegend = document.getElementById('pca-show-legend')?.checked ?? true;
        const useSelectedGenes = document.getElementById('pca-use-selected-genes')?.checked ?? false;

        const genes = useSelectedGenes ? store.get('selectedGenes') : null;
        if (useSelectedGenes && (!genes || genes.length < 2)) {
            container.innerHTML = '<div class="mpl-placeholder">Select at least 2 genes from the volcano plot to run PCA on selected genes.</div>';
            return;
        }

        const result = await python.generatePcaPlot({
            colorBy: colorBy,
            normalization: normalization,
            center: center,
            scale: scale,
            showLabels: showLabels,
            showLegend: showLegend,
            genes: genes,
            format: 'svg'
        });

        if (result.image) {
            container.innerHTML = result.image;
            // Show variance explained
            if (infoContainer && result.variance_explained) {
                infoContainer.innerHTML = `
                    <span class="variance-info">PC1: ${result.variance_explained[0].toFixed(1)}% variance</span>
                    <span class="variance-info">PC2: ${result.variance_explained[1].toFixed(1)}% variance</span>
                `;
            }
            // Update dropdown with available columns (only first time)
            if (result.metadata_columns) {
                const dropdown = document.getElementById('pca-color-by');
                if (dropdown && dropdown.options.length <= 1) {
                    const currentValue = dropdown.value;
                    dropdown.innerHTML = result.metadata_columns.map(col =>
                        `<option value="${col}" ${col === currentValue ? 'selected' : ''}>${col}</option>`
                    ).join('');
                }
            }
        } else if (result.error) {
            container.innerHTML = `<div class="mpl-placeholder">Error: ${result.error}</div>`;
        }
    } catch (error) {
        console.error('Error generating matplotlib PCA:', error);
        container.innerHTML = `<div class="mpl-placeholder">Error: ${error.message}</div>`;
    }
}

/**
 * Generate matplotlib GSEA plot
 */
async function generateMplGsea(topN) {
    const container = document.getElementById('gsea-plot-mpl');
    if (!container) return;

    try {
        container.innerHTML = '<div class="mpl-placeholder">Generating pathway plot...</div>';

        const result = await python.generateGseaPlot({
            topN: topN || 15,
            format: 'svg'
        });

        if (result.image) {
            container.innerHTML = result.image;
        } else if (result.error) {
            container.innerHTML = `<div class="mpl-placeholder">Error: ${result.error}</div>`;
        }
    } catch (error) {
        console.error('Error generating matplotlib GSEA:', error);
        container.innerHTML = `<div class="mpl-placeholder">Error: ${error.message}</div>`;
    }
}

/**
 * Generate matplotlib GSEA plot for visualization tab
 */
async function generateMplGseaViz(topN) {
    const container = document.getElementById('gsea-viz-plot-mpl');
    if (!container) return;

    try {
        container.innerHTML = '<div class="mpl-placeholder">Generating pathway plot...</div>';

        const result = await python.generateGseaPlot({
            topN: topN || 15,
            format: 'svg'
        });

        if (result.image) {
            container.innerHTML = result.image;
        } else if (result.error) {
            container.innerHTML = `<div class="mpl-placeholder">Error: ${result.error}</div>`;
        }
    } catch (error) {
        console.error('Error generating matplotlib GSEA:', error);
        container.innerHTML = `<div class="mpl-placeholder">Error: ${error.message}</div>`;
    }
}

/**
 * Generate Plotly PCA plot (interactive with hover tooltips)
 */
let _lastPcaSampleData = null;
async function generatePlotlyPca() {
    const container = document.getElementById('pca-plot-plotly');
    if (!container) return;

    try {
        container.innerHTML = '<div class="mpl-placeholder">Generating interactive PCA...</div>';

        // Get PCA data from Python (reuses same handler - returns sample_data)
        const colorBy = document.getElementById('pca-color-by')?.value || 'condition';
        const normalization = document.getElementById('pca-normalization')?.value || 'log2';
        const center = document.getElementById('pca-center')?.checked ?? true;
        const scale = document.getElementById('pca-scale')?.checked ?? true;
        const showLabels = document.getElementById('pca-show-labels')?.checked ?? true;
        const showLegend = document.getElementById('pca-show-legend')?.checked ?? true;
        const useSelectedGenes = document.getElementById('pca-use-selected-genes')?.checked ?? false;
        const genes = useSelectedGenes ? store.get('selectedGenes') : null;

        if (useSelectedGenes && (!genes || genes.length < 2)) {
            container.innerHTML = '<div class="mpl-placeholder">Select at least 2 genes for PCA on selected genes.</div>';
            return;
        }

        const result = await python.generatePcaPlot({
            colorBy, normalization, center, scale, showLabels, showLegend, genes, format: 'svg'
        });

        if (!result.sample_data || result.sample_data.length === 0) {
            container.innerHTML = '<div class="mpl-placeholder">No PCA data available</div>';
            return;
        }

        _lastPcaSampleData = result.sample_data;
        const varExplained = result.variance_explained || [0, 0];

        // Update variance info
        const infoContainer = document.getElementById('pca-variance-info');
        if (infoContainer && varExplained.length >= 2) {
            infoContainer.innerHTML = `
                <span class="variance-info">PC1: ${varExplained[0].toFixed(1)}% variance</span>
                <span class="variance-info">PC2: ${varExplained[1].toFixed(1)}% variance</span>
            `;
        }

        // Update dropdown
        if (result.metadata_columns) {
            const dropdown = document.getElementById('pca-color-by');
            if (dropdown && dropdown.options.length <= 1) {
                const current = dropdown.value;
                dropdown.innerHTML = result.metadata_columns.map(col =>
                    `<option value="${col}" ${col === current ? 'selected' : ''}>${col}</option>`
                ).join('');
            }
        }

        // Group samples by color_by attribute
        const groups = {};
        result.sample_data.forEach(s => {
            const key = s[colorBy] || 'unknown';
            if (!groups[key]) groups[key] = [];
            groups[key].push(s);
        });

        const colorMap = {
            'control': '#7BA3C9', 'treatment': '#C97B7B',
            'batch1': '#8FB996', 'batch2': '#B39DDB', 'batch3': '#E8C87D',
        };
        const defaultColors = ['#C97B7B', '#7BA3C9', '#8FB996', '#B39DDB', '#E8C87D', '#D4A5A5', '#A8D8EA'];

        const traces = [];
        let colorIdx = 0;
        Object.entries(groups).forEach(([groupName, samples]) => {
            const color = colorMap[groupName] || defaultColors[colorIdx % defaultColors.length];
            colorIdx++;

            // Build hover text with all metadata
            const hoverTexts = samples.map(s => {
                const lines = [`<b>${s.name}</b>`];
                Object.keys(s).forEach(k => {
                    if (k !== 'name' && k !== 'pc1' && k !== 'pc2') {
                        lines.push(`${k}: ${s[k]}`);
                    }
                });
                return lines.join('<br>');
            });

            traces.push({
                type: 'scatter',
                mode: showLabels ? 'markers+text' : 'markers',
                name: groupName,
                x: samples.map(s => s.pc1),
                y: samples.map(s => s.pc2),
                text: samples.map(s => s.name),
                textposition: 'top center',
                textfont: { size: 9 },
                hovertext: hoverTexts,
                hoverinfo: 'text',
                marker: {
                    color: color,
                    size: 12,
                    line: { color: '#333333', width: 1.5 },
                    opacity: 0.85
                }
            });
        });

        const layout = {
            title: { text: 'Principal Component Analysis', font: { size: 16, weight: 'bold' } },
            xaxis: { title: `PC1 (${varExplained[0].toFixed(1)}% variance)`, zeroline: true, zerolinecolor: '#ccc', gridcolor: '#f0f0f0' },
            yaxis: { title: `PC2 (${varExplained[1].toFixed(1)}% variance)`, zeroline: true, zerolinecolor: '#ccc', gridcolor: '#f0f0f0' },
            hovermode: 'closest',
            showlegend: showLegend,
            legend: { x: 1, xanchor: 'right', y: 1, bgcolor: 'rgba(255,255,255,0.9)', bordercolor: '#ccc', borderwidth: 1 },
            margin: { l: 60, r: 40, t: 60, b: 60 },
            paper_bgcolor: 'white',
            plot_bgcolor: 'white'
        };

        container.innerHTML = '';
        await Plotly.newPlot(container, traces, layout, { responsive: true, displayModeBar: true });

    } catch (error) {
        console.error('Error generating Plotly PCA:', error);
        container.innerHTML = `<div class="mpl-placeholder">Error: ${error.message}</div>`;
    }
}

/**
 * Generate Plotly heatmap (interactive)
 */
async function generatePlotlyHeatmap() {
    const container = document.getElementById('heatmap-plot-plotly');
    if (!container) return;

    try {
        container.innerHTML = '<div class="mpl-placeholder">Generating interactive heatmap...</div>';

        const geneSelection = document.getElementById('heatmap-gene-selection')?.value || 'variance';
        const topN = parseInt(document.getElementById('heatmap-topn')?.value || 50);
        const scale = document.getElementById('heatmap-scale')?.value || 'row';

        let genes = null;
        if (geneSelection === 'selected') {
            genes = store.get('selectedGenes');
            if (!genes || genes.length === 0) {
                container.innerHTML = '<div class="mpl-placeholder">Select genes from the volcano plot first.</div>';
                return;
            }
        }

        const data = await python.run('get_heatmap_data', {
            top_n: topN,
            genes: genes,
            scale: scale,
            annotation_columns: ['condition']
        });

        if (!data.genes || data.genes.length === 0) {
            container.innerHTML = '<div class="mpl-placeholder">No genes to display</div>';
            return;
        }

        const trace = {
            type: 'heatmap',
            z: data.z,
            x: data.samples,
            y: data.genes,
            colorscale: 'RdBu',
            reversescale: true,
            zmid: 0,
            hovertemplate: 'Gene: %{y}<br>Sample: %{x}<br>Value: %{z:.2f}<extra></extra>',
            colorbar: { title: scale === 'row' ? 'Z-score' : 'Log2(counts+1)', thickness: 15 }
        };

        const layout = {
            title: { text: 'Heatmap', font: { size: 16, weight: 'bold' } },
            xaxis: { title: 'Samples', tickangle: -45, tickfont: { size: 10 } },
            yaxis: { title: 'Genes', tickfont: { size: data.genes.length > 60 ? 6 : 9 }, autorange: 'reversed' },
            margin: { l: 120, r: 60, t: 60, b: 100 },
            paper_bgcolor: 'white',
            plot_bgcolor: 'white',
            height: Math.max(500, data.genes.length * 12 + 200)
        };

        container.innerHTML = '';
        await Plotly.newPlot(container, [trace], layout, { responsive: true, displayModeBar: true });

    } catch (error) {
        console.error('Error generating Plotly heatmap:', error);
        container.innerHTML = `<div class="mpl-placeholder">Error: ${error.message}</div>`;
    }
}

/**
 * Generate Plotly GSEA dot plot (interactive)
 */
async function generatePlotlyGsea(gseaResults, topN) {
    const container = document.getElementById('gsea-plot-plotly');
    if (!container) return;

    try {
        container.innerHTML = '<div class="mpl-placeholder">Generating interactive pathway plot...</div>';

        const results = gseaResults || store.get('gseaResults');
        if (!results || results.length === 0) {
            container.innerHTML = '<div class="mpl-placeholder">Run pathway analysis first</div>';
            return;
        }

        const n = topN || 15;
        const df = results.slice(0, n);

        const termCol = 'Term' in df[0] ? 'Term' : 'term';
        const nesCol = 'NES' in df[0] ? 'NES' : 'nes';
        const fdrCol = 'FDR q-val' in df[0] ? 'FDR q-val' : ('fdr' in df[0] ? 'fdr' : 'Adjusted P-value');

        const terms = df.map(r => {
            const t = r[termCol] || '';
            return t.length > 45 ? t.substring(0, 45) + '...' : t;
        });
        const nes = df.map(r => r[nesCol] || 0);
        const fdr = df.map(r => Math.max(r[fdrCol] || 0.1, 1e-10));
        const sizes = fdr.map(f => Math.max(8, Math.min(-Math.log10(f) * 4, 30)));
        const colors = nes.map(n => n > 0 ? '#C97B7B' : '#7BA3C9');

        const trace = {
            type: 'scatter',
            mode: 'markers',
            x: nes,
            y: terms.reverse(),
            marker: {
                size: sizes.reverse(),
                color: colors.reverse(),
                line: { color: '#333', width: 1 },
                opacity: 0.8
            },
            hovertemplate: df.reverse().map((r, i) =>
                `<b>${r[termCol]}</b><br>NES: ${(r[nesCol] || 0).toFixed(3)}<br>FDR: ${(r[fdrCol] || 0).toExponential(2)}<extra></extra>`
            ),
        };

        const layout = {
            title: { text: 'Pathway Enrichment', font: { size: 16, weight: 'bold' } },
            xaxis: { title: 'Normalized Enrichment Score (NES)', zeroline: true, zerolinecolor: '#999', gridcolor: '#f0f0f0' },
            yaxis: { tickfont: { size: 10 } },
            margin: { l: 280, r: 40, t: 60, b: 60 },
            paper_bgcolor: 'white',
            plot_bgcolor: 'white',
            shapes: [{
                type: 'line', x0: 0, x1: 0, y0: -0.5, y1: terms.length - 0.5,
                line: { color: '#999', width: 1, dash: 'dash' }
            }],
            height: Math.max(400, terms.length * 28 + 120)
        };

        container.innerHTML = '';
        await Plotly.newPlot(container, [trace], layout, { responsive: true, displayModeBar: true });

    } catch (error) {
        console.error('Error generating Plotly GSEA:', error);
        container.innerHTML = `<div class="mpl-placeholder">Error: ${error.message}</div>`;
    }
}

/**
 * Setup environment export panel (requirements.txt / conda yaml)
 */
function setupEnvironmentExport() {
    let envData = null;
    let currentFormat = 'requirements';

    const envTabs = document.querySelectorAll('.env-tab');
    const envOutput = document.getElementById('env-output');
    const btnLoad = document.getElementById('btn-load-env');
    const btnCopy = document.getElementById('btn-copy-env');

    // Format tabs
    envTabs.forEach(tab => {
        tab.addEventListener('click', () => {
            envTabs.forEach(t => t.classList.remove('active'));
            tab.classList.add('active');
            currentFormat = tab.dataset.env;
            if (envData) {
                envOutput.textContent = currentFormat === 'requirements'
                    ? envData.requirements_txt
                    : envData.conda_yaml;
            }
        });
    });

    // Load button
    btnLoad?.addEventListener('click', async () => {
        try {
            envOutput.textContent = '# Loading environment info...';
            envData = await python.getEnvironmentInfo();
            envOutput.textContent = currentFormat === 'requirements'
                ? envData.requirements_txt
                : envData.conda_yaml;
        } catch (error) {
            console.error('Error loading environment info:', error);
            envOutput.textContent = `# Error: ${error.message}`;
        }
    });

    // Copy button
    btnCopy?.addEventListener('click', async () => {
        const text = envOutput.textContent;
        try {
            await navigator.clipboard.writeText(text);
            const orig = btnCopy.textContent;
            btnCopy.textContent = 'Copied!';
            setTimeout(() => { btnCopy.textContent = orig; }, 2000);
        } catch (e) {
            console.error('Copy failed:', e);
        }
    });
}

/**
 * Render GSEA results table in visualization tab using Tabulator
 */
function renderGseaVizTable(results) {
    if (!results || results.length === 0) return;

    // Destroy existing table
    if (gseaVizTable) {
        gseaVizTable.destroy();
    }

    // Determine column names (varies by GSEA method)
    const firstRow = results[0];
    const termCol = 'Term' in firstRow ? 'Term' : 'term';
    const nesCol = 'NES' in firstRow ? 'NES' : ('nes' in firstRow ? 'nes' : null);
    const fdrCol = 'FDR q-val' in firstRow ? 'FDR q-val' : ('fdr' in firstRow ? 'fdr' : 'Adjusted P-value');
    const pvalCol = 'NOM p-val' in firstRow ? 'NOM p-val' : ('pval' in firstRow ? 'pval' : 'P-value');

    const columns = [
        {
            title: 'Pathway',
            field: termCol,
            headerFilter: 'input',
            headerFilterPlaceholder: 'Search...',
            width: 300,
            formatter: function(cell) {
                const val = cell.getValue();
                return val.length > 50 ? val.substring(0, 50) + '...' : val;
            },
            tooltip: true
        }
    ];

    if (nesCol && nesCol in firstRow) {
        columns.push({
            title: 'NES',
            field: nesCol,
            sorter: 'number',
            formatter: function(cell) {
                const val = cell.getValue();
                if (val == null) return '-';
                const color = val > 0 ? '#C97B7B' : '#7BA3C9';
                return `<span style="color: ${color}; font-weight: bold">${val.toFixed(3)}</span>`;
            }
        });
    }

    columns.push({
        title: 'P-value',
        field: pvalCol,
        sorter: 'number',
        formatter: function(cell) {
            const val = cell.getValue();
            return val != null ? val.toExponential(2) : '-';
        }
    });

    columns.push({
        title: 'FDR',
        field: fdrCol,
        sorter: 'number',
        headerFilter: 'number',
        headerFilterPlaceholder: 'max',
        headerFilterFunc: '<=',
        formatter: function(cell) {
            const val = cell.getValue();
            if (val == null) return '-';
            const style = val < 0.25 ? 'font-weight: bold; color: #8FB996' : '';
            return `<span style="${style}">${val.toExponential(2)}</span>`;
        }
    });

    // Add gene count column if available
    if ('Gene_set' in firstRow || 'geneset_size' in firstRow) {
        columns.push({
            title: 'Size',
            field: 'Gene_set' in firstRow ? 'Gene_set' : 'geneset_size',
            sorter: 'number'
        });
    }

    gseaVizTable = new Tabulator('#gsea-viz-table', {
        data: results,
        height: '400px',
        layout: 'fitColumns',
        pagination: true,
        paginationSize: 25,
        movableColumns: true,
        placeholder: 'No pathway results',
        columns: columns,
        initialSort: nesCol ? [{ column: nesCol, dir: 'desc' }] : []
    });
}


/**
 * Setup selected genes display
 */
function setupSelectedGenesDisplay() {
    const chipsContainer = document.getElementById('gene-chips');
    const clearBtn = document.getElementById('btn-clear-genes');

    // Update display when selection changes
    store.on('selectedGenes', async (genes) => {
        updateSelectedGenesDisplay(genes);
        updateSelectedCount(genes.length);

        // Update boxplot if visible
        if (document.getElementById('viz-boxplot')?.classList.contains('hidden') === false) {
            boxPlot.render();
            // Also refresh matplotlib boxplot
            if (genes.length > 0) {
                await generateMplBoxplot();
            }
        }
    });

    // Clear button
    clearBtn?.addEventListener('click', () => {
        store.clearSelectedGenes();
    });
}

/**
 * Update selected genes chips display
 */
function updateSelectedGenesDisplay(genes) {
    const container = document.getElementById('gene-chips');
    if (!container) return;

    if (genes.length === 0) {
        container.innerHTML = '<span class="placeholder">Click genes on volcano plot to select</span>';
        return;
    }

    container.innerHTML = genes.map(gene => `
        <span class="gene-chip">
            ${gene}
            <span class="remove" data-gene="${gene}">×</span>
        </span>
    `).join('');

    // Add remove handlers
    container.querySelectorAll('.remove').forEach(btn => {
        btn.addEventListener('click', () => {
            store.removeSelectedGene(btn.dataset.gene);
        });
    });
}

/**
 * Update status message in status bar
 */
function updateStatusMessage(message, isError = false) {
    const statusElement = document.getElementById('status-message');
    if (statusElement) {
        statusElement.textContent = message;
        statusElement.style.color = isError ? '#ffcccc' : 'white';
    }
}

/**
 * Show an inline status message under a tool section.
 * @param {string} elementId - ID of the .tool-status-msg element
 * @param {'loading'|'error'|'success'} type
 * @param {string} message
 * @param {number} [autoDismissMs] - auto-hide after ms (0 = never, default: success=5000, error=0)
 */
function showToolStatus(elementId, type, message, autoDismissMs) {
    const el = document.getElementById(elementId);
    if (!el) return;

    el.className = 'tool-status-msg visible ' + type;
    el.querySelector('.status-text').textContent = message;

    // Wire up dismiss button
    const dismissBtn = el.querySelector('.status-dismiss');
    if (dismissBtn) {
        dismissBtn.onclick = () => { el.className = 'tool-status-msg'; };
    }

    // Auto-dismiss for success messages
    const timeout = autoDismissMs !== undefined ? autoDismissMs : (type === 'success' ? 5000 : 0);
    if (timeout > 0) {
        setTimeout(() => {
            if (el.classList.contains(type)) {
                el.className = 'tool-status-msg';
            }
        }, timeout);
    }
}

function hideToolStatus(elementId) {
    const el = document.getElementById(elementId);
    if (el) el.className = 'tool-status-msg';
}

/**
 * Update selected genes count in status bar
 */
function updateSelectedCount(count) {
    const countElement = document.getElementById('selected-count');
    if (countElement) {
        countElement.textContent = `Selected: ${count} genes`;
    }
}

/**
 * Mark a workflow step as complete
 */
function markStepComplete(stepId) {
    const tab = document.querySelector(`.nav-tab[data-panel="${stepId}"]`);
    if (tab) {
        tab.classList.add('completed');
    }
}

// ---------------------------------------------------------------------------
// Citations Panel (Step 5)
// ---------------------------------------------------------------------------

function setupCitationsPanel() {
    const bibEntries = {
        kallisto: `@article{bray2016near,
  title={Near-optimal probabilistic RNA-seq quantification},
  author={Bray, Nicolas L and Pimentel, Harold and Melsted, P{\\'a}ll and Pachter, Lior},
  journal={Nature Biotechnology},
  volume={34},
  number={5},
  pages={525--527},
  year={2016},
  publisher={Nature Publishing Group},
  doi={10.1038/nbt.3519}
}`,
        pydeseq2: `@article{muzellec2023pydeseq2,
  title={PyDESeq2: a python package for bulk RNA-seq differential expression analysis},
  author={Muzellec, Boris and Telenczuk, Maria and Cabeli, Vincent and Andreux, Mathieu},
  journal={Bioinformatics},
  volume={39},
  number={9},
  pages={btad547},
  year={2023},
  publisher={Oxford University Press},
  doi={10.1093/bioinformatics/btad547}
}`,
        deseq2: `@article{love2014moderated,
  title={Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2},
  author={Love, Michael I and Huber, Wolfgang and Anders, Simon},
  journal={Genome Biology},
  volume={15},
  number={12},
  pages={550},
  year={2014},
  publisher={BioMed Central},
  doi={10.1186/s13059-014-0550-8}
}`,
        edgestar: `@article{pachter2026edgestar,
  title={Differential analysis of genomics count data with edge*},
  author={Pachter, Lior},
  journal={bioRxiv},
  year={2026},
  doi={10.64898/2026.02.16.706223}
}`,
        edger: `@article{robinson2010edger,
  title={edgeR: a Bioconductor package for differential expression analysis of digital gene expression data},
  author={Robinson, Mark D and McCarthy, Davis J and Smyth, Gordon K},
  journal={Bioinformatics},
  volume={26},
  number={1},
  pages={139--140},
  year={2010},
  publisher={Oxford University Press},
  doi={10.1093/bioinformatics/btp616}
}`,
        recount3: `@article{wilks2021recount3,
  title={recount3: summaries and queries for large-scale RNA-seq expression and splicing},
  author={Wilks, Christopher and Zheng, Shijie C and Chen, Feng Yong and Charles, Rone and Solomon, Brad and Ling, Jonathan P and Imada, Eddie Luidy and Zhang, David and Joseph, Lance and Leek, Jeffrey T and Jaffe, Andrew E and Nellore, Abhinav and Collado-Torres, Leonardo and Hansen, Kasper D and Langmead, Ben},
  journal={Genome Biology},
  volume={22},
  number={1},
  pages={323},
  year={2021},
  publisher={BioMed Central},
  doi={10.1186/s13059-021-02533-6}
}`,
        geoquery: `@article{davis2007geoquery,
  title={GEOquery: a bridge between the Gene Expression Omnibus (GEO) and BioConductor},
  author={Davis, Sean and Meltzer, Paul S},
  journal={Bioinformatics},
  volume={23},
  number={14},
  pages={1846--1847},
  year={2007},
  publisher={Oxford University Press},
  doi={10.1093/bioinformatics/btm254}
}`
    };

    function getAllBib() {
        return Object.values(bibEntries).join('\n\n');
    }

    function getAllNatureCitations() {
        const cards = document.querySelectorAll('.citation-card .citation-text');
        return Array.from(cards).map(el => el.textContent.trim()).join('\n\n');
    }

    document.getElementById('btn-download-bib')?.addEventListener('click', () => {
        const bib = getAllBib();
        const blob = new Blob([bib], { type: 'text/plain' });
        const url = URL.createObjectURL(blob);
        const a = document.createElement('a');
        a.href = url;
        a.download = 'rnaseq_analysis_citations.bib';
        document.body.appendChild(a);
        a.click();
        document.body.removeChild(a);
        URL.revokeObjectURL(url);
    });

    document.getElementById('btn-copy-citations')?.addEventListener('click', async () => {
        const text = getAllNatureCitations();
        try {
            await navigator.clipboard.writeText(text);
            const btn = document.getElementById('btn-copy-citations');
            const orig = btn.innerHTML;
            btn.innerHTML = '<span class="btn-icon">✓</span> Copied!';
            setTimeout(() => { btn.innerHTML = orig; }, 2000);
        } catch (e) {
            console.error('Copy failed:', e);
        }
    });
}

// ---------------------------------------------------------------------------
// FASTQ Processing Panel (Step 0)
// ---------------------------------------------------------------------------

function setupFastqPanel() {
    // Check platform — disable FASTQ on Windows (no kallisto binary)
    if (window.api && window.api.getPlatformInfo) {
        window.api.getPlatformInfo().then(info => {
            if (info.platform === 'win32') {
                const panel = document.getElementById('panel-fastq');
                if (panel) {
                    const warning = document.createElement('div');
                    warning.className = 'platform-warning';
                    warning.innerHTML = `
                        <div class="warning-banner">
                            <strong>FASTQ Processing Unavailable on Windows</strong>
                            <p>Kallisto pseudo-alignment requires macOS or Linux.
                            You can still load a pre-computed counts matrix in the Data tab.</p>
                        </div>
                    `;
                    const content = panel.querySelector('.panel-content');
                    if (content) content.prepend(warning);
                    // Disable all interactive elements
                    panel.querySelectorAll('button, input, select').forEach(el => {
                        el.disabled = true;
                    });
                }
                // Also mark FASTQ tab visually
                const fastqTab = document.querySelector('[data-panel="fastq"] .tab-label');
                if (fastqTab) fastqTab.textContent = 'FASTQ (N/A)';
                return; // Don't set up FASTQ event listeners
            }
        });
    }

    // State for the FASTQ panel
    const fastqState = {
        indexPath: null,
        fastaPath: null,
        metadataPath: null,
        samples: [],          // [{name, fastqFiles: [], status: 'pending'}]
        singleEnd: true,
    };

    // --- Index source radio buttons ---
    const indexRadios = document.querySelectorAll('input[name="index-source"]');
    const customFastaControls = document.getElementById('custom-fasta-controls');
    const customIndexControls = document.getElementById('custom-index-controls');

    indexRadios.forEach(radio => {
        radio.addEventListener('change', () => {
            customFastaControls?.classList.add('hidden');
            customIndexControls?.classList.add('hidden');
            if (radio.value === 'custom-fasta') {
                customFastaControls?.classList.remove('hidden');
            } else if (radio.value === 'custom-index') {
                customIndexControls?.classList.remove('hidden');
            }
        });
    });

    // --- Browse existing index file ---
    document.getElementById('btn-browse-index')?.addEventListener('click', async () => {
        try {
            const filePath = await python.openFileDialog({
                filters: [{ name: 'Kallisto Index', extensions: ['idx', 'kidx'] }]
            });
            if (filePath) {
                fastqState.customIndexPath = filePath;
                document.getElementById('index-path-display').textContent = filePath.split('/').pop();
            }
        } catch (e) {
            console.error('Browse index error:', e);
        }
    });

    // --- Check for prebuilt index on load ---
    python.onReady(async () => {
        try {
            const result = await python.kallistoCheckIndex();
            const statusEl = document.getElementById('prebuilt-index-status');
            const downloadBtn = document.getElementById('btn-download-index');

            if (result.has_index) {
                statusEl.textContent = 'Available';
                statusEl.classList.add('status-ok');
                fastqState.indexPath = result.index_path;
                downloadBtn?.classList.add('hidden');
            } else {
                statusEl.textContent = 'Not downloaded';
                statusEl.classList.add('status-missing');
                downloadBtn?.classList.remove('hidden');
            }
        } catch (e) {
            console.error('Failed to check index:', e);
        }
    });

    // --- Download index ---
    document.getElementById('btn-download-index')?.addEventListener('click', async () => {
        const statusEl = document.getElementById('prebuilt-index-status');
        const downloadBtn = document.getElementById('btn-download-index');
        statusEl.textContent = 'Downloading (138 MB)...';
        statusEl.classList.remove('status-missing');
        downloadBtn.disabled = true;

        try {
            const result = await python.kallistoDownloadIndex();
            fastqState.indexPath = result.index_path;
            statusEl.textContent = 'Available';
            statusEl.classList.add('status-ok');
            downloadBtn?.classList.add('hidden');
        } catch (e) {
            statusEl.textContent = 'Download failed: ' + e.message;
            statusEl.classList.add('status-missing');
            downloadBtn.disabled = false;
        }
    });

    // --- FASTQ mode toggle ---
    document.querySelectorAll('input[name="fastq-mode"]').forEach(radio => {
        radio.addEventListener('change', () => {
            fastqState.singleEnd = radio.value === 'single';
        });
    });

    // --- Load example FASTQ data ---
    document.getElementById('btn-load-example-fastq-bottom')?.addEventListener('click', async (e) => {
        e.preventDefault();
        try {
            updateStatusMessage('Loading example FASTQ data...');
            const result = await python.kallistoLoadExample();

            fastqState.fastaPath = result.fasta_path;
            fastqState.metadataPath = result.metadata_path;

            // Select custom-fasta radio and show controls
            const customRadio = document.querySelector('input[name="index-source"][value="custom-fasta"]');
            if (customRadio) customRadio.checked = true;
            customFastaControls?.classList.remove('hidden');

            // Display FASTA path
            document.getElementById('fasta-path-display').textContent = result.fasta_path;

            // Display metadata path
            document.getElementById('fastq-metadata-path-display').textContent = result.metadata_path;

            // Populate sample table
            fastqState.samples = result.samples.map(s => ({
                name: s.name,
                fastqFiles: [s.fastq_path],
                condition: s.condition,
                status: 'pending',
            }));
            renderSampleTable();

            updateStatusMessage('Example FASTQ data loaded');
        } catch (e) {
            console.error('Failed to load example FASTQ:', e);
            updateStatusMessage('Error: ' + e.message);
        }
    });

    // --- Browse FASTA ---
    document.getElementById('btn-browse-fasta')?.addEventListener('click', async () => {
        try {
            const filePath = await python.openFileDialog({
                filters: [{ name: 'FASTA', extensions: ['fa', 'fasta', 'fna', 'fa.gz', 'fasta.gz'] }]
            });
            if (filePath) {
                fastqState.fastaPath = filePath;
                document.getElementById('fasta-path-display').textContent = filePath;
            }
        } catch (e) {
            console.error('File dialog failed:', e);
        }
    });

    // --- Build index ---
    document.getElementById('btn-build-index')?.addEventListener('click', async () => {
        if (!fastqState.fastaPath) {
            updateStatusMessage('Select a FASTA file first');
            return;
        }

        const statusEl = document.getElementById('build-index-status');
        statusEl.textContent = 'Building index...';

        try {
            const result = await python.kallistoBuildIndex(fastqState.fastaPath);
            fastqState.indexPath = result.index_path;
            statusEl.textContent = 'Index built successfully';
            statusEl.classList.add('status-ok');
        } catch (e) {
            statusEl.textContent = 'Build failed: ' + e.message;
            statusEl.classList.add('status-missing');
        }
    });

    // --- Add sample ---
    document.getElementById('btn-add-sample')?.addEventListener('click', () => {
        const name = `Sample_${fastqState.samples.length + 1}`;
        fastqState.samples.push({ name, fastqFiles: [], status: 'pending' });
        renderSampleTable();
    });

    // --- Browse metadata ---
    document.getElementById('btn-browse-metadata-fastq')?.addEventListener('click', async () => {
        try {
            const filePath = await python.openFileDialog({
                filters: [{ name: 'CSV', extensions: ['csv', 'tsv', 'txt'] }]
            });
            if (filePath) {
                fastqState.metadataPath = filePath;
                document.getElementById('fastq-metadata-path-display').textContent = filePath;
                // Parse metadata and assign groups to matching samples
                try {
                    const meta = await python.loadMetadata(filePath);
                    if (meta && meta.metadata_preview && meta.metadata_columns && meta.metadata_columns.length > 0) {
                        const groupCol = meta.metadata_columns[0];
                        const groupData = meta.metadata_preview[groupCol] || {};
                        fastqState.samples.forEach(sample => {
                            if (groupData[sample.name] !== undefined) {
                                sample.condition = groupData[sample.name];
                            }
                        });
                    }
                } catch (parseErr) {
                    console.error('Could not parse metadata for groups:', parseErr);
                }
                renderSampleTable();
            }
        } catch (e) {
            console.error('File dialog failed:', e);
        }
    });

    // --- Run kallisto ---
    document.getElementById('btn-run-kallisto')?.addEventListener('click', async () => {
        // Resolve index path
        const selectedSource = document.querySelector('input[name="index-source"]:checked')?.value;
        let indexPath = null;

        if (selectedSource === 'custom-index') {
            // Use a user-provided index file
            indexPath = fastqState.customIndexPath;
            if (!indexPath) {
                updateStatusMessage('Select an index file first');
                return;
            }
        } else if (selectedSource === 'custom-fasta') {
            // Build from custom FASTA
            if (!fastqState.fastaPath) {
                updateStatusMessage('Select a FASTA file first');
                return;
            }
            updateStatusMessage('Building index from custom FASTA...');
            try {
                const result = await python.kallistoBuildIndex(fastqState.fastaPath);
                indexPath = result.index_path;
                fastqState.customIndexPath = indexPath;
                document.getElementById('build-index-status').textContent = 'Index built successfully';
            } catch (e) {
                updateStatusMessage('Index build failed: ' + e.message);
                return;
            }
        } else {
            // Use prebuilt index
            indexPath = fastqState.indexPath;
        }

        if (!indexPath) {
            updateStatusMessage('No index available. Download or build one first.');
            return;
        }

        if (fastqState.samples.length === 0) {
            updateStatusMessage('Add at least one sample');
            return;
        }

        // Show progress
        const progressContainer = document.getElementById('kallisto-progress');
        const progressBar = document.getElementById('kallisto-progress-bar');
        const progressText = document.getElementById('kallisto-progress-text');
        progressContainer?.classList.remove('hidden');

        const total = fastqState.samples.length;
        let completed = 0;

        for (const sample of fastqState.samples) {
            if (sample.fastqFiles.length === 0) {
                sample.status = 'error';
                renderSampleTable();
                continue;
            }

            sample.status = 'running';
            renderSampleTable();
            progressText.textContent = `Processing ${sample.name} (${completed + 1}/${total})...`;
            progressBar.style.width = `${(completed / total) * 100}%`;

            try {
                const result = await python.kallistoQuantSample(
                    sample.name,
                    sample.fastqFiles,
                    { singleEnd: fastqState.singleEnd, indexPath }
                );
                sample.status = 'done';
                sample.totalCounts = result.total_counts;
                sample.nTargets = result.n_targets;
            } catch (e) {
                sample.status = 'error';
                sample.error = e.message;
                console.error(`Kallisto quant failed for ${sample.name}:`, e);
            }

            completed++;
            renderSampleTable();
        }

        progressBar.style.width = '100%';
        progressText.textContent = `Done: ${completed}/${total} samples processed`;

        // Show results section
        const resultsSection = document.getElementById('fastq-results-section');
        resultsSection?.classList.remove('hidden');

        const summaryEl = document.getElementById('kallisto-results-summary');
        const doneSamples = fastqState.samples.filter(s => s.status === 'done');
        summaryEl.innerHTML = `
            <p>${doneSamples.length} of ${total} samples quantified successfully.</p>
            ${doneSamples.map(s => `<div>${s.name}: ${Math.round(s.totalCounts)} estimated counts across ${s.nTargets} targets</div>`).join('')}
        `;

        // Generate QC stats and plot
        try {
            const [qcStats, qcPlot] = await Promise.all([
                python.kallistoGetQcStats(),
                python.generateKallistoQcPlot(),
            ]);

            // Render QC plot
            const plotContainer = document.getElementById('kallisto-qc-plot');
            if (plotContainer && qcPlot.plot) {
                if (qcPlot.format === 'svg') {
                    plotContainer.innerHTML = qcPlot.plot;
                } else {
                    plotContainer.innerHTML = `<img src="data:image/png;base64,${qcPlot.plot}" style="max-width:100%">`;
                }
            }

            // Render QC stats table
            const tableContainer = document.getElementById('kallisto-qc-table');
            if (tableContainer && qcStats.samples && qcStats.samples.length > 0) {
                tableContainer.innerHTML = `
                    <table class="qc-table">
                        <thead>
                            <tr>
                                <th>Sample</th>
                                <th>Reads Processed</th>
                                <th>Pseudoaligned</th>
                                <th>Mapping Rate</th>
                                <th>Unique Rate</th>
                                <th>Est. Counts</th>
                                <th>Targets</th>
                            </tr>
                        </thead>
                        <tbody>
                            ${qcStats.samples.map(s => `
                                <tr>
                                    <td>${s.name}</td>
                                    <td>${s.n_processed.toLocaleString()}</td>
                                    <td>${s.n_pseudoaligned.toLocaleString()}</td>
                                    <td class="${s.p_pseudoaligned >= 50 ? 'qc-good' : s.p_pseudoaligned >= 20 ? 'qc-warn' : 'qc-bad'}">${s.p_pseudoaligned.toFixed(1)}%</td>
                                    <td>${s.p_unique.toFixed(1)}%</td>
                                    <td>${Math.round(s.total_counts).toLocaleString()}</td>
                                    <td>${s.n_targets.toLocaleString()}</td>
                                </tr>
                            `).join('')}
                        </tbody>
                    </table>
                `;
            }
        } catch (e) {
            console.error('QC stats generation failed:', e);
        }

        updateStatusMessage('Kallisto quantification complete');
    });

    // --- Download QC plot as PNG ---
    document.getElementById('btn-download-qc-plot')?.addEventListener('click', () => {
        const plotContainer = document.getElementById('kallisto-qc-plot');
        if (!plotContainer) return;
        const svg = plotContainer.querySelector('svg');
        const img = plotContainer.querySelector('img');
        if (svg) {
            // SVG → PNG via canvas
            const svgData = new XMLSerializer().serializeToString(svg);
            const canvas = document.createElement('canvas');
            const bbox = svg.getBoundingClientRect();
            canvas.width = bbox.width * 2;
            canvas.height = bbox.height * 2;
            const ctx = canvas.getContext('2d');
            ctx.scale(2, 2);
            const imgEl = new Image();
            imgEl.onload = () => {
                ctx.drawImage(imgEl, 0, 0);
                const a = document.createElement('a');
                a.download = 'kallisto_qc_summary.png';
                a.href = canvas.toDataURL('image/png');
                a.click();
            };
            imgEl.src = 'data:image/svg+xml;base64,' + btoa(unescape(encodeURIComponent(svgData)));
        } else if (img) {
            const a = document.createElement('a');
            a.download = 'kallisto_qc_summary.png';
            a.href = img.src;
            a.click();
        }
    });

    // --- Download QC table as CSV ---
    document.getElementById('btn-download-qc-table')?.addEventListener('click', () => {
        const table = document.querySelector('#kallisto-qc-table table');
        if (!table) return;
        const rows = Array.from(table.querySelectorAll('tr'));
        const csv = rows.map(row =>
            Array.from(row.querySelectorAll('th, td'))
                .map(cell => `"${cell.textContent.trim()}"`)
                .join(',')
        ).join('\n');
        const blob = new Blob([csv], { type: 'text/csv' });
        const a = document.createElement('a');
        a.download = 'kallisto_qc_stats.csv';
        a.href = URL.createObjectURL(blob);
        a.click();
        URL.revokeObjectURL(a.href);
    });

    // --- Combine & Continue ---
    document.getElementById('btn-combine-continue')?.addEventListener('click', async () => {
        try {
            updateStatusMessage('Combining counts...');
            const result = await python.kallistoCombineCounts({
                metadataPath: fastqState.metadataPath,
            });

            // Show Kallisto indicator in tool bar
            const kallistoBar = document.getElementById('tool-bar-kallisto');
            if (kallistoBar) kallistoBar.classList.remove('hidden');

            // Navigate to data panel
            navigateToPanel('data');
            showDataPreview();
            markStepComplete('fastq');
            updateStatusMessage(`Counts matrix created: ${result.counts_shape[0]} genes x ${result.counts_shape[1]} samples`);
        } catch (e) {
            console.error('Combine failed:', e);
            updateStatusMessage('Error combining counts: ' + e.message);
        }
    });

    // --- Render sample table ---
    function renderSampleTable() {
        const tbody = document.getElementById('fastq-sample-tbody');
        if (!tbody) return;

        tbody.innerHTML = fastqState.samples.map((sample, idx) => {
            const statusClass = sample.status === 'done' ? 'status-ok' :
                                sample.status === 'error' ? 'status-missing' :
                                sample.status === 'running' ? 'status-running' : '';
            const statusText = sample.status === 'done' ? 'Done' :
                               sample.status === 'error' ? 'Error' :
                               sample.status === 'running' ? 'Running...' : 'Pending';
            const filesText = sample.fastqFiles.length > 0
                ? sample.fastqFiles.map(f => f.split('/').pop()).join(', ')
                : '<em>No files</em>';
            const groupText = sample.condition || '<em>--</em>';

            return `
                <tr>
                    <td><input type="text" class="sample-name-input" value="${sample.name}" data-idx="${idx}"></td>
                    <td class="fastq-files-cell">
                        <span class="file-list">${filesText}</span>
                        <button class="btn btn-tiny btn-secondary btn-browse-sample-fastq" data-idx="${idx}">Browse</button>
                    </td>
                    <td><span class="sample-status ${statusClass}">${statusText}</span></td>
                    <td class="sample-group-cell">${groupText}</td>
                    <td><button class="btn btn-tiny btn-text btn-remove-sample" data-idx="${idx}">Remove</button></td>
                </tr>
            `;
        }).join('');

        // Re-attach event listeners
        tbody.querySelectorAll('.sample-name-input').forEach(input => {
            input.addEventListener('change', (e) => {
                const idx = parseInt(e.target.dataset.idx);
                fastqState.samples[idx].name = e.target.value;
            });
        });

        tbody.querySelectorAll('.btn-browse-sample-fastq').forEach(btn => {
            btn.addEventListener('click', async (e) => {
                const idx = parseInt(e.target.dataset.idx);
                try {
                    const filePath = await python.openFileDialog({
                        filters: [{ name: 'FASTQ', extensions: ['fastq', 'fq', 'fastq.gz', 'fq.gz'] }]
                    });
                    if (filePath) {
                        fastqState.samples[idx].fastqFiles = [filePath];
                        renderSampleTable();
                    }
                } catch (err) {
                    console.error('File dialog failed:', err);
                }
            });
        });

        tbody.querySelectorAll('.btn-remove-sample').forEach(btn => {
            btn.addEventListener('click', (e) => {
                const idx = parseInt(e.target.dataset.idx);
                fastqState.samples.splice(idx, 1);
                renderSampleTable();
            });
        });
    }
}

/**
 * Navigate to a specific panel programmatically
 */
function navigateToPanel(panelId) {
    // Update horizontal tab nav
    const navTabs = document.querySelectorAll('.nav-tab');
    navTabs.forEach(t => t.classList.remove('active'));
    const targetTab = document.querySelector(`.nav-tab[data-panel="${panelId}"]`);
    if (targetTab) targetTab.classList.add('active');

    // Show panel
    document.querySelectorAll('.panel').forEach(p => p.classList.add('hidden'));
    document.getElementById(`panel-${panelId}`)?.classList.remove('hidden');

    store.set('currentPanel', panelId);
}

// ---------------------------------------------------------------------------
// First-run setup flow
// ---------------------------------------------------------------------------

async function checkAndRunSetup() {
    const setupScreen = document.getElementById('setup-screen');
    const appDiv = document.getElementById('app');

    // If running in a browser (dev server), skip setup entirely
    if (!window.api || !window.api.checkPythonSetup) {
        initApp();
        return;
    }

    try {
        const status = await window.api.checkPythonSetup();

        if (status.ready) {
            // Python is already set up — start backend FIRST (non-blocking)
            window.api.startPythonBackend().catch(e => console.error('Failed to start Python backend:', e));
            setupScreen.style.display = 'none';
            appDiv.style.display = '';
            initApp();
            return;
        }

        // Need first-run setup
        setupScreen.style.display = 'flex';
        appDiv.style.display = 'none';

        const detail = document.getElementById('setup-detail');
        const progressBar = document.getElementById('setup-progress-bar');
        const note = document.getElementById('setup-note');

        // Listen for progress
        window.api.onPythonSetupProgress((data) => {
            detail.textContent = data.message;
            progressBar.style.width = `${Math.round(data.progress * 100)}%`;
        });

        detail.textContent = 'Downloading and installing Python environment...';
        note.textContent = 'This only needs to happen once. Please stay connected to the internet.';

        const result = await window.api.runPythonSetup();

        if (result.success) {
            detail.textContent = 'Setup complete! Starting application...';
            progressBar.style.width = '100%';

            // Brief pause so user sees "complete"
            await new Promise(r => setTimeout(r, 800));

            // Start backend FIRST (non-blocking)
            window.api.startPythonBackend().catch(e => console.error('Failed to start Python backend:', e));
            setupScreen.style.display = 'none';
            appDiv.style.display = '';
            initApp();
        } else {
            detail.textContent = 'Setup failed: ' + result.error;
            note.textContent = 'Please check your internet connection and try restarting the app.';
            progressBar.style.width = '0%';
        }

    } catch (err) {
        console.error('Setup check failed:', err);
        // Fallback: try to start Python backend first, then UI
        try {
            window.api.startPythonBackend().catch(e => console.error('Failed to start Python backend:', e));
        } catch (e) {
            console.error('Failed to start Python backend:', e);
        }
        setupScreen.style.display = 'none';
        appDiv.style.display = '';
        initApp();
    }
}

// Initialize when DOM is ready
document.addEventListener('DOMContentLoaded', checkAndRunSetup);
