"""
KallistoRunner — wraps the precompiled kallisto binary for index building,
quantification, and count-matrix assembly.
"""

import os
import sys
import platform
import subprocess
import tempfile
import tarfile
import lzma
import urllib.request
import shutil
from pathlib import Path
from typing import Dict, List, Optional, Callable

import pandas as pd


# ---------------------------------------------------------------------------
# App-level config (populated by main_handler on init_config from Electron)
# ---------------------------------------------------------------------------
APP_CONFIG: Dict = {
    'resources_path': None,
    'user_data_path': None,
    'platform': sys.platform,
    'is_packaged': False,
}


def _project_root() -> Path:
    """Return the project / resources root.

    When running inside a packaged Electron app the resources_path is set via
    init_config.  In development it falls back to two levels up from this file.
    """
    if APP_CONFIG.get('resources_path'):
        return Path(APP_CONFIG['resources_path'])
    return Path(__file__).resolve().parent.parent.parent


class KallistoRunner:
    def __init__(self):
        self.binary = self._get_binary_path()
        self.index_path: Optional[str] = None
        self.t2g: Optional[pd.DataFrame] = None

    # ------------------------------------------------------------------
    # Platform availability
    # ------------------------------------------------------------------
    @property
    def is_available(self) -> bool:
        return self.binary is not None

    # ------------------------------------------------------------------
    # Binary resolution
    # ------------------------------------------------------------------
    def _get_binary_path(self) -> Optional[str]:
        root = _project_root()
        system = sys.platform  # 'darwin', 'linux', 'win32'

        if system == 'darwin':
            binary = root / 'bin' / 'kallisto' / 'mac' / 'kallisto'
        elif system.startswith('linux'):
            binary = root / 'bin' / 'kallisto' / 'linux' / 'kallisto'
        else:
            # Windows or unsupported — kallisto not available
            return None

        if not binary.exists():
            return None

        # Ensure executable
        try:
            binary.chmod(0o755)
        except OSError:
            pass
        return str(binary)

    def _require_available(self):
        """Raise a friendly error if kallisto is not available."""
        if not self.is_available:
            raise RuntimeError(
                'Kallisto is not available on this platform. '
                'FASTQ processing requires macOS or Linux. '
                'You can still load a pre-computed counts matrix in the Data tab.'
            )

    # ------------------------------------------------------------------
    # Prebuilt index helpers
    # ------------------------------------------------------------------
    INDEX_URL = (
        'https://github.com/pachterlab/kallisto-transcriptome-indices'
        '/releases/download/v1/human_index_standard.tar.xz'
    )

    @staticmethod
    def _indices_dir() -> Path:
        """Writable directory for downloaded indices.

        In a packaged app, the app bundle is read-only, so we store indices
        in the user data directory instead.
        """
        if APP_CONFIG.get('user_data_path'):
            return Path(APP_CONFIG['user_data_path']) / 'kallisto_indices' / 'human_standard'
        return _project_root() / 'bin' / 'kallisto' / 'indices' / 'human_standard'

    def get_prebuilt_index_path(self) -> Optional[str]:
        idx = self._indices_dir() / 'index.idx'
        return str(idx) if idx.exists() else None

    def get_t2g_path(self) -> Optional[str]:
        t2g = self._indices_dir() / 't2g.txt'
        return str(t2g) if t2g.exists() else None

    def has_prebuilt_index(self) -> bool:
        return self.get_prebuilt_index_path() is not None

    def download_human_index(
        self,
        progress_callback: Optional[Callable[[float, str], None]] = None,
    ) -> str:
        """Download the prebuilt human transcriptome index (~138 MB).

        Args:
            progress_callback: optional fn(fraction, message) called during download.

        Returns:
            Path to the extracted index.idx file.
        """
        dest_dir = self._indices_dir()
        dest_dir.mkdir(parents=True, exist_ok=True)

        idx_path = dest_dir / 'index.idx'
        if idx_path.exists():
            return str(idx_path)

        if progress_callback:
            progress_callback(0.0, 'Downloading human transcriptome index (138 MB)...')

        tmp_xz = os.path.join(tempfile.gettempdir(), 'human_index_standard.tar.xz')

        # Stream download with progress
        req = urllib.request.urlopen(self.INDEX_URL)
        total = int(req.headers.get('Content-Length', 0))
        downloaded = 0
        chunk_size = 1024 * 256  # 256 KB

        with open(tmp_xz, 'wb') as f:
            while True:
                chunk = req.read(chunk_size)
                if not chunk:
                    break
                f.write(chunk)
                downloaded += len(chunk)
                if progress_callback and total > 0:
                    frac = min(downloaded / total, 0.9)
                    mb = downloaded / (1024 * 1024)
                    progress_callback(frac, f'Downloading... {mb:.0f} / {total / (1024*1024):.0f} MB')

        if progress_callback:
            progress_callback(0.9, 'Extracting index...')

        # Decompress .tar.xz  →  extract index.idx and t2g.txt
        with lzma.open(tmp_xz) as xz:
            with tarfile.open(fileobj=xz) as tar:
                for member in tar.getmembers():
                    basename = os.path.basename(member.name)
                    if basename in ('index.idx', 't2g.txt'):
                        member.name = basename
                        tar.extract(member, path=str(dest_dir))

        # Cleanup
        try:
            os.remove(tmp_xz)
        except OSError:
            pass

        if progress_callback:
            progress_callback(1.0, 'Human transcriptome index ready')

        return str(idx_path)

    def load_t2g(self, path: Optional[str] = None) -> pd.DataFrame:
        path = path or self.get_t2g_path()
        if path is None:
            return pd.DataFrame()
        df = pd.read_csv(path, sep='\t', header=None,
                         names=['transcript_id', 'gene_id', 'gene_name'])
        self.t2g = df
        return df

    def transcript_to_gene_map(self) -> Dict[str, str]:
        if self.t2g is None:
            self.load_t2g()
        if self.t2g is None or self.t2g.empty:
            return {}
        return dict(zip(self.t2g['transcript_id'], self.t2g['gene_name']))

    # ------------------------------------------------------------------
    # Index building
    # ------------------------------------------------------------------
    def build_index(self, fasta_path: str, output_path: Optional[str] = None) -> str:
        self._require_available()

        if output_path is None:
            output_path = os.path.join(
                tempfile.gettempdir(), 'kallisto_index.idx'
            )

        cmd = [self.binary, 'index', '-i', output_path, fasta_path]
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=600)

        if result.returncode != 0:
            raise RuntimeError(
                f'kallisto index failed:\n{result.stderr}'
            )

        self.index_path = output_path
        return output_path

    # ------------------------------------------------------------------
    # Quantification
    # ------------------------------------------------------------------
    def quant_sample(
        self,
        index_path: str,
        output_dir: str,
        fastq_files: List[str],
        single_end: bool = False,
        fragment_length: float = 200,
        sd: float = 20,
    ) -> str:
        self._require_available()

        os.makedirs(output_dir, exist_ok=True)

        cmd = [self.binary, 'quant', '-i', index_path, '-o', output_dir]

        if single_end:
            cmd += ['--single', '-l', str(fragment_length), '-s', str(sd)]

        cmd += fastq_files

        result = subprocess.run(cmd, capture_output=True, text=True, timeout=3600)

        if result.returncode != 0:
            raise RuntimeError(
                f'kallisto quant failed (exit {result.returncode}):\n{result.stderr}'
            )

        # Check that abundance.tsv was actually created
        abundance_file = os.path.join(output_dir, 'abundance.tsv')
        if not os.path.exists(abundance_file):
            raise RuntimeError(
                f'kallisto quant completed but no abundance.tsv produced:\n{result.stderr}'
            )

        return output_dir

    # ------------------------------------------------------------------
    # Abundance parsing
    # ------------------------------------------------------------------
    def parse_abundance(self, output_dir: str) -> pd.DataFrame:
        tsv = os.path.join(output_dir, 'abundance.tsv')
        if not os.path.exists(tsv):
            raise FileNotFoundError(f'No abundance.tsv in {output_dir}')

        df = pd.read_csv(tsv, sep='\t')
        df = df.set_index('target_id')
        return df

    # ------------------------------------------------------------------
    # Combine into counts matrix
    # ------------------------------------------------------------------
    def combine_counts(
        self,
        sample_dirs: Dict[str, str],
        use_gene_names: bool = True,
    ) -> pd.DataFrame:
        """Merge per-sample abundance.tsv files into a genes x samples int count matrix."""
        frames = {}
        for sample_name, output_dir in sample_dirs.items():
            ab = self.parse_abundance(output_dir)
            frames[sample_name] = ab['est_counts']

        counts = pd.DataFrame(frames)

        # Map transcript IDs to gene names if t2g is available
        if use_gene_names:
            t2g_map = self.transcript_to_gene_map()
            if t2g_map:
                counts.index = counts.index.map(
                    lambda tid: t2g_map.get(tid, tid)
                )
                # Sum counts for transcripts mapping to the same gene
                counts = counts.groupby(counts.index).sum()

        # Convert to integer counts
        counts = counts.round(0).astype(int)
        counts.index.name = 'gene'
        return counts
