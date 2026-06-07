"""Database backend abstraction and SQLite implementation for FragalyseQt."""
from __future__ import annotations
from hashlib import md5, sha1, sha256, sha3_256
from os.path import exists, getsize
from zlib import compress as zlib_compress, decompress as zlib_decompress
from numpy import (array as np_array, frombuffer as np_frombuffer, float32,
                   float64)
from abc import ABC, abstractmethod
from contextlib import contextmanager
from dataclasses import dataclass
from datetime import datetime, timezone
from typing import Optional


# ---------------------------------------------------------------------------
# Shared defaults
# ---------------------------------------------------------------------------

# Default `wlen` for scipy.signal.find_peaks (odd value — find_peaks rounds
# even values down internally, but keeping the canonical default odd keeps
# UI display and SQL DEFAULT consistent).  Used by:
#   * fragalyseqt.py — initial value of the "Detection window width" spinbox
#   * AnalysisRunRecord.peak_window — dataclass default
#   * analysis_run.peak_window column — SQL DEFAULT for both ALTER and CREATE
DEFAULT_PEAK_WINDOW = 51


# ---------------------------------------------------------------------------
# File hashing
# ---------------------------------------------------------------------------

def compute_hashes(path: str) -> dict:
    """Compute MD5/SHA-1/SHA-256/SHA-3-256 in a single pass. Returns dict
    with keys hash_md5, hash_sha1, hash_sha256, hash_sha3_256, size."""
    h = {
        'hash_md5': md5(),
        'hash_sha1': sha1(),
        'hash_sha256': sha256(),
        'hash_sha3_256': sha3_256(),
    }
    size = 0
    with open(path, 'rb') as f:
        for chunk in iter(lambda: f.read(1 << 20), b''):
            size += len(chunk)
            for d in h.values():
                d.update(chunk)
    return {k: v.hexdigest() for k, v in h.items()} | {'size': size}


# ---------------------------------------------------------------------------
# Channel signal compression (zlib + float32, stdlib only)
# ---------------------------------------------------------------------------

def compress_signal(signal) -> bytes:
    """Compress a signal channel list to a zlib BLOB (float32 precision)."""
    return zlib_compress(np_array(signal, dtype=float32).tobytes(), level=6)


def decompress_signal(data: bytes):
    """Decompress a channel signal BLOB back to a float64 ndarray."""
    return np_frombuffer(zlib_decompress(data), dtype=float32).astype(float64)


# ---------------------------------------------------------------------------
# File verification
# ---------------------------------------------------------------------------

@dataclass
class FileStatus:
    file_id: int
    file_name: str
    path: str
    # 'ok' | 'missing' | 'modified'
    status: str
    # empty string when status is 'ok'
    detail: str


def verify_file(stored: dict) -> FileStatus:
    """Verify one instrument_file row against the filesystem."""
    fid = stored['id']
    name = stored['file_name']
    path = stored['file_path']

    if not exists(path):
        return FileStatus(fid, name, path, 'missing', f'Not found at {path}')

    if getsize(path) != stored['file_size']:
        return FileStatus(fid, name, path, 'modified', 'File has been modified')

    h = {
        'hash_md5': md5(),
        'hash_sha1': sha1(),
        'hash_sha256': sha256(),
        'hash_sha3_256': sha3_256(),
    }
    with open(path, 'rb') as f:
        for chunk in iter(lambda: f.read(1 << 20), b''):
            for d in h.values():
                d.update(chunk)

    failed = [
        algo.replace('hash_', '')
        for algo, digest in h.items()
        if digest.hexdigest() != stored[algo]
    ]
    if failed:
        return FileStatus(fid, name, path, 'modified',
                          'Hash mismatch: ' + ', '.join(failed))

    return FileStatus(fid, name, path, 'ok', '')


def verify_session(backend: DatabaseBackend,
                   session_id: int) -> list[FileStatus]:
    """Return verification status for every unique file in a session."""
    seen: set[int] = set()
    results: list[FileStatus] = []
    for tab in backend.get_session_tabs(session_id):
        run = backend.get_run_info(tab['run_id'])
        fid = run['file_id']
        if fid in seen:
            continue
        seen.add(fid)
        results.append(verify_file(backend.get_file_info(fid)))
    return results


# ---------------------------------------------------------------------------
# Record dataclasses
# ---------------------------------------------------------------------------

@dataclass
class InstrumentFileRecord:
    created_by: str
    file_name: str
    file_path: str
    file_size: int
    hash_md5: str
    hash_sha1: str
    hash_sha256: str
    hash_sha3_256: str
    instrument: str = ''
    run_name: str = ''


@dataclass
class DyeChannelRecord:
    file_id: int
    channel: int
    dye_name: str
    wavelength: Optional[int] = None


@dataclass
class AnalysisRunRecord:
    file_id: int
    created_by: str
    min_height: float
    min_prominence: float
    min_width: float
    halfwindow_width: int
    baseline_correction: bool
    sizing_method: str
    size_standard: str
    panel: str
    peak_window: int = DEFAULT_PEAK_WINDOW
    supersedes_id: Optional[int] = None


@dataclass
class PeakCallRecord:
    run_id: int
    created_by: str
    channel: int
    dye_name: str
    position_dp: float
    position_bp: Optional[float]
    height: float
    area: Optional[float]
    fwhm: Optional[float]
    is_ladder: bool
    supersedes_id: Optional[int] = None
    # True when the peak clipped at the ADC ceiling and its height/area/fwhm
    # are Gaussian-flank estimates rather than measured values.
    saturated: bool = False


@dataclass
class AlleleCallRecord:
    peak_id: int
    created_by: str
    allele: str
    marker: str = ''
    bin_distance: Optional[float] = None
    note: Optional[str] = None
    supersedes_id: Optional[int] = None


@dataclass
class ChannelSignalRecord:
    file_id: int
    # 1-based channel index
    channel: int
    # zlib-compressed float32 array
    signal: bytes


@dataclass
class SavedSessionRecord:
    created_by: str
    name: str
    supersedes_id: Optional[int] = None


@dataclass
class SessionTabRecord:
    session_id: int
    tab_order: int
    run_id: int


# ---------------------------------------------------------------------------
# Abstract backend
# ---------------------------------------------------------------------------

class DatabaseBackend(ABC):

    @abstractmethod
    def store_file(self, record: InstrumentFileRecord) -> int:
        """Insert or reuse instrument_file row. Returns file id."""

    @abstractmethod
    def store_dye_channel(self, record: DyeChannelRecord) -> int:
        """Insert a dye_channel row. Returns new id."""

    @abstractmethod
    def store_analysis_run(self, record: AnalysisRunRecord) -> int:
        """Insert an analysis_run row. Returns new id."""

    @abstractmethod
    def store_peak_call(self, record: PeakCallRecord) -> int:
        """Insert a peak_call row. Returns new id."""

    @abstractmethod
    def store_allele_call(self, record: AlleleCallRecord) -> int:
        """Insert an allele_call row. Returns new id."""

    @abstractmethod
    def store_session(self, record: SavedSessionRecord) -> int:
        """Insert a saved_session row. Returns new id."""

    @abstractmethod
    def store_session_tab(self, record: SessionTabRecord) -> int:
        """Insert a session_tab row. Returns new id."""

    @abstractmethod
    def get_session_list(self) -> list[dict]:
        """Return current (non-superseded) saved sessions, newest first."""

    @abstractmethod
    def get_session_tabs(self, session_id: int) -> list[dict]:
        """Return session_tab rows ordered by tab_order."""

    @abstractmethod
    def get_file_info(self, file_id: int) -> dict:
        """Return instrument_file row for a given id."""

    @abstractmethod
    def get_dye_channels(self, file_id: int) -> list[dict]:
        """Return dye_channel rows ordered by channel."""

    @abstractmethod
    def get_run_info(self, run_id: int) -> dict:
        """Return analysis_run row for a given id."""

    @abstractmethod
    def get_peak_calls_for_run(self, run_id: int) -> list[dict]:
        """Return current peak_call rows for a run."""

    @abstractmethod
    def get_allele_calls_for_run(self, run_id: int) -> list[dict]:
        """Return allele_call rows joined with their peaks for a run."""

    @abstractmethod
    def find_file_by_path_and_hash(self, path: str,
                                   sha256: str) -> Optional[dict]:
        """Return existing instrument_file if path and sha256 both match."""

    @abstractmethod
    def store_channel_signal(self, record: ChannelSignalRecord) -> int:
        """Insert a channel_signal row. Returns new id."""

    @abstractmethod
    def get_channel_signals(self, file_id: int) -> list[dict]:
        """Return channel_signal rows for a file, ordered by channel."""

    @abstractmethod
    def find_session_by_name(self, name: str) -> Optional[dict]:
        """Return the current (non-superseded) session with this name, or None."""

    @abstractmethod
    def hide_session(self, session_id: int) -> None:
        """Mark a session as deleted by inserting a tombstone record.
        Append-only: no UPDATE or DELETE involved."""

    @abstractmethod
    def close(self) -> None:
        """Release connection/resources."""


# ---------------------------------------------------------------------------
# Shared SQLite infrastructure
# ---------------------------------------------------------------------------

class _SQLiteBase:
    """Common SQLite plumbing shared by all append-only backends."""

    def __init__(self, db_path: str, schema_sql: str) -> None:
        import sqlite3
        self._conn = sqlite3.connect(db_path, check_same_thread=False)
        self._conn.row_factory = sqlite3.Row
        self._conn.execute('PRAGMA foreign_keys = ON')
        self._auto_commit = True
        self._conn.executescript(schema_sql)
        # Lightweight schema migrations: add columns introduced after the
        # initial schema.  ALTER TABLE ADD COLUMN errors out if the column
        # already exists, so each migration is wrapped in its own try.
        for ddl in (
            'ALTER TABLE analysis_run '
            f'ADD COLUMN peak_window INTEGER NOT NULL '
            f'DEFAULT {DEFAULT_PEAK_WINDOW}',
            # Add halfwindow_width as a NEW column (no rename/drop — the
            # database is append-only).  Pre-existing rows keep their
            # legacy window_width values; new rows write only
            # halfwindow_width.  On read (`get_run_info`) the legacy
            # value is translated to the new half-width semantics so old
            # sessions still restore their saved smoothing setting.
            'ALTER TABLE analysis_run '
            'ADD COLUMN halfwindow_width INTEGER',
            # Saturated/clipped-peak flag (added after initial schema).
            # Append-only: ADD COLUMN never rewrites existing rows; pre-existing
            # peaks default to 0 (not saturated).
            'ALTER TABLE peak_call '
            'ADD COLUMN saturated INTEGER NOT NULL DEFAULT 0',
        ):
            try:
                self._conn.execute(ddl)
            except sqlite3.OperationalError:
                pass
        # Re-run the schema so the DROP/CREATE views are rebuilt AFTER the
        # ALTER migrations above: a `SELECT *` view created before a column was
        # added would otherwise not expose it on pre-existing databases.
        # (CREATE TABLE/INDEX IF NOT EXISTS are no-ops on the second pass.)
        self._conn.executescript(schema_sql)
        self._conn.commit()

    @contextmanager
    def transaction(self):
        self._conn.execute('BEGIN IMMEDIATE')
        self._auto_commit = False
        try:
            yield
            self._conn.commit()
        except Exception:
            self._conn.rollback()
            raise
        finally:
            self._auto_commit = True

    @staticmethod
    def _now() -> str:
        return datetime.now(timezone.utc).isoformat()

    def close(self) -> None:
        self._conn.close()


# ---------------------------------------------------------------------------
# SQLite backend (casework)
# ---------------------------------------------------------------------------

class SQLiteBackend(_SQLiteBase, DatabaseBackend):

    def __init__(self, db_path: str) -> None:
        _SQLiteBase.__init__(self, db_path, _SCHEMA_SQL)

    # --- write ---

    def store_file(self, record: InstrumentFileRecord) -> int:
        existing = self.find_file_by_path_and_hash(record.file_path,
                                                   record.hash_sha256)
        if existing:
            return existing['id']
        cur = self._conn.execute(
            'INSERT INTO instrument_file '
            '(created_at,created_by,file_name,file_path,file_size,'
            ' hash_md5,hash_sha1,hash_sha256,hash_sha3_256,instrument,run_name)'
            ' VALUES (?,?,?,?,?,?,?,?,?,?,?)',
            (self._now(), record.created_by, record.file_name,
             record.file_path, record.file_size,
             record.hash_md5, record.hash_sha1,
             record.hash_sha256, record.hash_sha3_256,
             record.instrument, record.run_name))
        if self._auto_commit:
            self._conn.commit()
        return cur.lastrowid

    def store_dye_channel(self, record: DyeChannelRecord) -> int:
        cur = self._conn.execute(
            'INSERT INTO dye_channel (file_id,channel,dye_name,wavelength)'
            ' VALUES (?,?,?,?)',
            (record.file_id, record.channel, record.dye_name,
             record.wavelength))
        if self._auto_commit:
            self._conn.commit()
        return cur.lastrowid

    def store_analysis_run(self, record: AnalysisRunRecord) -> int:
        cur = self._conn.execute(
            'INSERT INTO analysis_run '
            '(created_at,created_by,supersedes_id,file_id,'
            ' min_height,min_prominence,min_width,halfwindow_width,'
            ' baseline_correction,sizing_method,size_standard,panel,'
            ' peak_window)'
            ' VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?)',
            (self._now(), record.created_by, record.supersedes_id,
             record.file_id, record.min_height, record.min_prominence,
             record.min_width, record.halfwindow_width,
             int(record.baseline_correction),
             record.sizing_method, record.size_standard, record.panel,
             record.peak_window))
        if self._auto_commit:
            self._conn.commit()
        return cur.lastrowid

    def store_peak_call(self, record: PeakCallRecord) -> int:
        cur = self._conn.execute(
            'INSERT INTO peak_call '
            '(created_at,created_by,supersedes_id,run_id,channel,'
            ' dye_name,position_dp,position_bp,height,area,fwhm,is_ladder,'
            ' saturated)'
            ' VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?)',
            (self._now(), record.created_by, record.supersedes_id,
             record.run_id, record.channel, record.dye_name,
             record.position_dp, record.position_bp,
             record.height, record.area, record.fwhm,
             int(record.is_ladder), int(record.saturated)))
        if self._auto_commit:
            self._conn.commit()
        return cur.lastrowid

    def store_allele_call(self, record: AlleleCallRecord) -> int:
        cur = self._conn.execute(
            'INSERT INTO allele_call '
            '(created_at,created_by,supersedes_id,peak_id,'
            ' allele,marker,bin_distance,note)'
            ' VALUES (?,?,?,?,?,?,?,?)',
            (self._now(), record.created_by, record.supersedes_id,
             record.peak_id, record.allele, record.marker,
             record.bin_distance, record.note))
        if self._auto_commit:
            self._conn.commit()
        return cur.lastrowid

    def store_session(self, record: SavedSessionRecord) -> int:
        cur = self._conn.execute(
            'INSERT INTO saved_session (created_at,created_by,supersedes_id,name)'
            ' VALUES (?,?,?,?)',
            (self._now(), record.created_by,
             record.supersedes_id, record.name))
        if self._auto_commit:
            self._conn.commit()
        return cur.lastrowid

    def store_session_tab(self, record: SessionTabRecord) -> int:
        cur = self._conn.execute(
            'INSERT INTO session_tab (session_id,tab_order,run_id)'
            ' VALUES (?,?,?)',
            (record.session_id, record.tab_order, record.run_id))
        if self._auto_commit:
            self._conn.commit()
        return cur.lastrowid

    # --- read ---

    def get_session_list(self) -> list[dict]:
        cur = self._conn.execute(
            'SELECT * FROM current_saved_session ORDER BY created_at DESC')
        return [dict(r) for r in cur.fetchall()]

    def get_session_tabs(self, session_id: int) -> list[dict]:
        cur = self._conn.execute(
            'SELECT * FROM session_tab WHERE session_id=? ORDER BY tab_order',
            (session_id,))
        return [dict(r) for r in cur.fetchall()]

    def get_file_info(self, file_id: int) -> dict:
        cur = self._conn.execute(
            'SELECT * FROM instrument_file WHERE id=?', (file_id,))
        return dict(cur.fetchone())

    def get_dye_channels(self, file_id: int) -> list[dict]:
        cur = self._conn.execute(
            'SELECT * FROM dye_channel WHERE file_id=? ORDER BY channel',
            (file_id,))
        return [dict(r) for r in cur.fetchall()]

    def get_run_info(self, run_id: int) -> dict:
        cur = self._conn.execute(
            'SELECT * FROM analysis_run WHERE id=?', (run_id,))
        row = dict(cur.fetchone())
        # Backwards-compat read path.  New tables and new rows only carry
        # halfwindow_width; the schema is append-only, though, so old
        # databases still contain the legacy full-width window_width on
        # pre-rename rows.  Translate it to the new semantics on read so
        # that loading a historic session restores the user's stored
        # smoothing setting.  Old data is never modified — only the
        # in-memory dict surfaced to the caller is enriched.
        if (row.get('halfwindow_width') is None
                and row.get('window_width') is not None):
            row['halfwindow_width'] = (int(row['window_width']) - 1) // 2
        return row

    def get_peak_calls_for_run(self, run_id: int) -> list[dict]:
        cur = self._conn.execute(
            'SELECT * FROM current_peak_call WHERE run_id=? ORDER BY id',
            (run_id,))
        return [dict(r) for r in cur.fetchall()]

    def get_allele_calls_for_run(self, run_id: int) -> list[dict]:
        cur = self._conn.execute(
            'SELECT ac.* FROM current_allele_call ac'
            ' JOIN current_peak_call pc ON ac.peak_id=pc.id'
            ' WHERE pc.run_id=?', (run_id,))
        return [dict(r) for r in cur.fetchall()]

    def find_file_by_path_and_hash(self, path: str,
                                   sha256: str) -> Optional[dict]:
        cur = self._conn.execute(
            'SELECT * FROM instrument_file'
            ' WHERE file_path=? AND hash_sha256=? LIMIT 1',
            (path, sha256))
        row = cur.fetchone()
        return dict(row) if row else None

    def store_channel_signal(self, record: ChannelSignalRecord) -> int:
        cur = self._conn.execute(
            'INSERT INTO channel_signal (file_id, channel, signal) VALUES (?,?,?)',
            (record.file_id, record.channel, record.signal))
        if self._auto_commit:
            self._conn.commit()
        return cur.lastrowid

    def get_channel_signals(self, file_id: int) -> list[dict]:
        cur = self._conn.execute(
            'SELECT * FROM channel_signal WHERE file_id=? ORDER BY channel',
            (file_id,))
        return [dict(r) for r in cur.fetchall()]

    def find_session_by_name(self, name: str) -> Optional[dict]:
        cur = self._conn.execute(
            'SELECT * FROM current_saved_session WHERE name=? LIMIT 1', (name,))
        row = cur.fetchone()
        return dict(row) if row else None

    def hide_session(self, session_id: int) -> None:
        self._conn.execute(
            'INSERT INTO session_deletion (created_at, session_id) VALUES (?,?)',
            (self._now(), session_id))
        if self._auto_commit:
            self._conn.commit()


# ---------------------------------------------------------------------------
# SQLite backend — reference profiles (refprofiles.db)
# ---------------------------------------------------------------------------

_REFPROFILE_SCHEMA_SQL = """
CREATE TABLE IF NOT EXISTS reference_profile (
    id INTEGER PRIMARY KEY,
    created_at TEXT NOT NULL,
    supersedes_id INTEGER REFERENCES reference_profile(id),
    name TEXT NOT NULL,
    role TEXT,
    notes TEXT,
    session_id INTEGER
);
CREATE TABLE IF NOT EXISTS reference_allele (
    id INTEGER PRIMARY KEY,
    profile_id INTEGER NOT NULL REFERENCES reference_profile(id),
    marker TEXT NOT NULL,
    allele1 TEXT NOT NULL,
    allele2 TEXT
);
CREATE TABLE IF NOT EXISTS reference_profile_deletion (
    id INTEGER PRIMARY KEY,
    created_at TEXT NOT NULL,
    profile_id INTEGER NOT NULL REFERENCES reference_profile(id)
);
CREATE INDEX IF NOT EXISTS idx_ref_allele_profile
    ON reference_allele(profile_id);
DROP VIEW IF EXISTS current_reference_profile;
CREATE VIEW current_reference_profile AS
    SELECT * FROM reference_profile
    WHERE id NOT IN (
        SELECT supersedes_id FROM reference_profile
        WHERE supersedes_id IS NOT NULL)
      AND id NOT IN (
        SELECT profile_id FROM reference_profile_deletion);
"""


class RefProfileBackend(_SQLiteBase):

    def __init__(self, db_path: str) -> None:
        _SQLiteBase.__init__(self, db_path, _REFPROFILE_SCHEMA_SQL)

    def store_reference_profile(self, name: str, role, notes,
                                session_id, supersedes_id=None) -> int:
        cur = self._conn.execute(
            'INSERT INTO reference_profile'
            ' (created_at, supersedes_id, name, role, notes, session_id)'
            ' VALUES (?,?,?,?,?,?)',
            (self._now(), supersedes_id, name, role, notes, session_id))
        if self._auto_commit:
            self._conn.commit()
        return cur.lastrowid

    def store_reference_alleles(self, profile_id: int, alleles: list) -> None:
        self._conn.executemany(
            'INSERT INTO reference_allele (profile_id, marker, allele1, allele2)'
            ' VALUES (?,?,?,?)',
            [(profile_id, a['marker'], a['allele1'], a['allele2'])
             for a in alleles])
        if self._auto_commit:
            self._conn.commit()

    def delete_reference_profile(self, profile_id: int) -> None:
        self._conn.execute(
            'INSERT INTO reference_profile_deletion (created_at, profile_id)'
            ' VALUES (?,?)',
            (self._now(), profile_id))
        if self._auto_commit:
            self._conn.commit()

    def get_reference_profile(self, profile_id: int) -> dict:
        cur = self._conn.execute(
            'SELECT * FROM reference_profile WHERE id=?', (profile_id,))
        return dict(cur.fetchone())

    def get_reference_alleles(self, profile_id: int) -> list:
        cur = self._conn.execute(
            'SELECT marker, allele1, allele2 FROM reference_allele'
            ' WHERE profile_id=? ORDER BY rowid',
            (profile_id,))
        return [dict(r) for r in cur.fetchall()]

    def list_reference_profiles(self) -> list:
        cur = self._conn.execute(
            'SELECT * FROM current_reference_profile ORDER BY created_at DESC')
        return [dict(r) for r in cur.fetchall()]


# ---------------------------------------------------------------------------
# Backend factory
# ---------------------------------------------------------------------------

def open_backend(config: dict) -> DatabaseBackend:
    backend = config.get('backend', 'sqlite')
    if backend == 'sqlite':
        return SQLiteBackend(config['path'])
    if backend == 'postgresql':
        raise NotImplementedError('PostgreSQL backend not yet implemented')
    if backend == 'immudb':
        raise NotImplementedError('ImmuDB backend not yet implemented')
    raise ValueError(f'Unknown backend: {backend}')


# ---------------------------------------------------------------------------
# Schema DDL
# ---------------------------------------------------------------------------

_SCHEMA_SQL = f"""
CREATE TABLE IF NOT EXISTS instrument_file (
    id             INTEGER PRIMARY KEY,
    created_at     TEXT    NOT NULL,
    created_by     TEXT    NOT NULL,
    file_name      TEXT    NOT NULL,
    file_path      TEXT    NOT NULL,
    file_size      INTEGER NOT NULL,
    hash_md5       TEXT    NOT NULL,
    hash_sha1      TEXT    NOT NULL,
    hash_sha256    TEXT    NOT NULL,
    hash_sha3_256  TEXT    NOT NULL,
    instrument     TEXT,
    run_name       TEXT
);
CREATE TABLE IF NOT EXISTS dye_channel (
    id         INTEGER PRIMARY KEY,
    file_id    INTEGER NOT NULL REFERENCES instrument_file(id),
    channel    INTEGER NOT NULL,
    dye_name   TEXT    NOT NULL,
    wavelength INTEGER
);
CREATE TABLE IF NOT EXISTS channel_signal (
    id         INTEGER PRIMARY KEY,
    file_id    INTEGER NOT NULL REFERENCES instrument_file(id),
    channel    INTEGER NOT NULL,
    signal     BLOB    NOT NULL
);
CREATE TABLE IF NOT EXISTS analysis_run (
    id                  INTEGER PRIMARY KEY,
    created_at          TEXT    NOT NULL,
    created_by          TEXT    NOT NULL,
    supersedes_id       INTEGER REFERENCES analysis_run(id),
    file_id             INTEGER NOT NULL REFERENCES instrument_file(id),
    min_height          REAL,
    min_prominence      REAL,
    min_width           REAL,
    halfwindow_width        INTEGER,
    baseline_correction INTEGER NOT NULL DEFAULT 0,
    sizing_method       TEXT,
    size_standard       TEXT,
    panel               TEXT,
    peak_window         INTEGER NOT NULL DEFAULT {DEFAULT_PEAK_WINDOW}
);
CREATE TABLE IF NOT EXISTS peak_call (
    id              INTEGER PRIMARY KEY,
    created_at      TEXT    NOT NULL,
    created_by      TEXT    NOT NULL,
    supersedes_id   INTEGER REFERENCES peak_call(id),
    run_id          INTEGER NOT NULL REFERENCES analysis_run(id),
    channel         INTEGER NOT NULL,
    dye_name        TEXT,
    position_dp     REAL    NOT NULL,
    position_bp     REAL,
    height          REAL    NOT NULL,
    area            REAL,
    fwhm            REAL,
    is_ladder       INTEGER NOT NULL DEFAULT 0,
    saturated       INTEGER NOT NULL DEFAULT 0
);
CREATE TABLE IF NOT EXISTS allele_call (
    id              INTEGER PRIMARY KEY,
    created_at      TEXT    NOT NULL,
    created_by      TEXT    NOT NULL,
    supersedes_id   INTEGER REFERENCES allele_call(id),
    peak_id         INTEGER NOT NULL REFERENCES peak_call(id),
    allele          TEXT    NOT NULL,
    marker          TEXT,
    bin_distance    REAL,
    is_stutter      INTEGER NOT NULL DEFAULT 0,
    note            TEXT
);
CREATE TABLE IF NOT EXISTS saved_session (
    id              INTEGER PRIMARY KEY,
    created_at      TEXT    NOT NULL,
    created_by      TEXT    NOT NULL,
    supersedes_id   INTEGER REFERENCES saved_session(id),
    name            TEXT    NOT NULL
);
CREATE TABLE IF NOT EXISTS session_deletion (
    id          INTEGER PRIMARY KEY,
    created_at  TEXT    NOT NULL,
    session_id  INTEGER NOT NULL REFERENCES saved_session(id)
);
CREATE TABLE IF NOT EXISTS session_tab (
    id          INTEGER PRIMARY KEY,
    session_id  INTEGER NOT NULL REFERENCES saved_session(id),
    tab_order   INTEGER NOT NULL,
    run_id      INTEGER NOT NULL REFERENCES analysis_run(id)
);
DROP VIEW IF EXISTS current_allele_call;
DROP VIEW IF EXISTS current_peak_call;
DROP VIEW IF EXISTS current_analysis_run;
DROP VIEW IF EXISTS current_saved_session;
CREATE VIEW current_analysis_run AS
    SELECT * FROM analysis_run WHERE id NOT IN (
        SELECT supersedes_id FROM analysis_run WHERE supersedes_id IS NOT NULL);
CREATE VIEW current_peak_call AS
    SELECT * FROM peak_call
    WHERE id NOT IN (
        SELECT supersedes_id FROM peak_call WHERE supersedes_id IS NOT NULL)
      AND run_id IN (SELECT id FROM current_analysis_run);
CREATE VIEW current_allele_call AS
    SELECT * FROM allele_call
    WHERE id NOT IN (
        SELECT supersedes_id FROM allele_call WHERE supersedes_id IS NOT NULL)
      AND peak_id IN (SELECT id FROM current_peak_call);
CREATE VIEW current_saved_session AS
    SELECT * FROM saved_session
    WHERE id NOT IN (
        SELECT supersedes_id FROM saved_session WHERE supersedes_id IS NOT NULL)
      AND id NOT IN (
        SELECT session_id FROM session_deletion);
"""
