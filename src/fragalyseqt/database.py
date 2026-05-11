"""Database backend abstraction and SQLite implementation for FragalyseQt."""
from __future__ import annotations

import hashlib
import os
from abc import ABC, abstractmethod
from dataclasses import dataclass
from datetime import datetime, timezone
from typing import Optional


# ---------------------------------------------------------------------------
# File hashing
# ---------------------------------------------------------------------------

def compute_hashes(path: str) -> dict:
    """Compute MD5/SHA-1/SHA-256/SHA-3-256 in a single pass. Returns dict
    with keys hash_md5, hash_sha1, hash_sha256, hash_sha3_256, size."""
    h = {
        'hash_md5':      hashlib.md5(),
        'hash_sha1':     hashlib.sha1(),
        'hash_sha256':   hashlib.sha256(),
        'hash_sha3_256': hashlib.sha3_256(),
    }
    size = 0
    with open(path, 'rb') as f:
        for chunk in iter(lambda: f.read(1 << 20), b''):
            size += len(chunk)
            for d in h.values():
                d.update(chunk)
    return {k: v.hexdigest() for k, v in h.items()} | {'size': size}


# ---------------------------------------------------------------------------
# File verification
# ---------------------------------------------------------------------------

@dataclass
class FileStatus:
    file_id:   int
    file_name: str
    path:      str
    status:    str   # 'ok' | 'missing' | 'modified'
    detail:    str   # human-readable explanation, empty when ok


def verify_file(stored: dict) -> FileStatus:
    """Verify one instrument_file row against the filesystem."""
    fid  = stored['id']
    name = stored['file_name']
    path = stored['file_path']

    if not os.path.exists(path):
        return FileStatus(fid, name, path, 'missing', f'Not found at {path}')

    if os.path.getsize(path) != stored['file_size']:
        return FileStatus(fid, name, path, 'modified', 'File has been modified')

    h = {
        'hash_md5':      hashlib.md5(),
        'hash_sha1':     hashlib.sha1(),
        'hash_sha256':   hashlib.sha256(),
        'hash_sha3_256': hashlib.sha3_256(),
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


def verify_session(backend: DatabaseBackend, session_id: int) -> list[FileStatus]:
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
    created_by:    str
    file_name:     str
    file_path:     str
    file_size:     int
    hash_md5:      str
    hash_sha1:     str
    hash_sha256:   str
    hash_sha3_256: str
    instrument:    str = ''
    run_name:      str = ''


@dataclass
class DyeChannelRecord:
    file_id:    int
    channel:    int
    dye_name:   str
    wavelength: Optional[int] = None


@dataclass
class AnalysisRunRecord:
    file_id:             int
    created_by:          str
    min_height:          float
    min_prominence:      float
    min_width:           float
    window_width:        int
    baseline_correction: bool
    sizing_method:       str
    size_standard:       str
    panel:               str
    supersedes_id:       Optional[int] = None


@dataclass
class PeakCallRecord:
    run_id:        int
    created_by:    str
    channel:       int
    dye_name:      str
    position_dp:   float
    position_bp:   Optional[float]
    height:        float
    area:          Optional[float]
    fwhm:          Optional[float]
    is_ladder:     bool
    supersedes_id: Optional[int] = None


@dataclass
class AlleleCallRecord:
    peak_id:       int
    created_by:    str
    allele:        str
    marker:        str = ''
    bin_distance:  Optional[float] = None
    is_stutter:    bool = False
    note:          Optional[str] = None
    supersedes_id: Optional[int] = None


@dataclass
class SavedSessionRecord:
    created_by:    str
    name:          str
    supersedes_id: Optional[int] = None


@dataclass
class SessionTabRecord:
    session_id: int
    tab_order:  int
    run_id:     int


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
    def find_file_by_path_and_hash(self, path: str, sha256: str) -> Optional[dict]:
        """Return existing instrument_file if path and sha256 both match."""

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
# SQLite backend
# ---------------------------------------------------------------------------

class SQLiteBackend(DatabaseBackend):

    def __init__(self, db_path: str) -> None:
        import sqlite3
        self._conn = sqlite3.connect(db_path, check_same_thread=False)
        self._conn.row_factory = sqlite3.Row
        self._conn.execute('PRAGMA foreign_keys = ON')
        self._apply_schema()

    def _apply_schema(self) -> None:
        self._conn.executescript(_SCHEMA_SQL)
        self._conn.commit()

    @staticmethod
    def _now() -> str:
        return datetime.now(timezone.utc).isoformat()

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
        self._conn.commit()
        return cur.lastrowid

    def store_dye_channel(self, record: DyeChannelRecord) -> int:
        cur = self._conn.execute(
            'INSERT INTO dye_channel (file_id,channel,dye_name,wavelength)'
            ' VALUES (?,?,?,?)',
            (record.file_id, record.channel, record.dye_name, record.wavelength))
        self._conn.commit()
        return cur.lastrowid

    def store_analysis_run(self, record: AnalysisRunRecord) -> int:
        cur = self._conn.execute(
            'INSERT INTO analysis_run '
            '(created_at,created_by,supersedes_id,file_id,'
            ' min_height,min_prominence,min_width,window_width,'
            ' baseline_correction,sizing_method,size_standard,panel)'
            ' VALUES (?,?,?,?,?,?,?,?,?,?,?,?)',
            (self._now(), record.created_by, record.supersedes_id,
             record.file_id, record.min_height, record.min_prominence,
             record.min_width, record.window_width,
             int(record.baseline_correction),
             record.sizing_method, record.size_standard, record.panel))
        self._conn.commit()
        return cur.lastrowid

    def store_peak_call(self, record: PeakCallRecord) -> int:
        cur = self._conn.execute(
            'INSERT INTO peak_call '
            '(created_at,created_by,supersedes_id,run_id,channel,'
            ' dye_name,position_dp,position_bp,height,area,fwhm,is_ladder)'
            ' VALUES (?,?,?,?,?,?,?,?,?,?,?,?)',
            (self._now(), record.created_by, record.supersedes_id,
             record.run_id, record.channel, record.dye_name,
             record.position_dp, record.position_bp,
             record.height, record.area, record.fwhm,
             int(record.is_ladder)))
        self._conn.commit()
        return cur.lastrowid

    def store_allele_call(self, record: AlleleCallRecord) -> int:
        cur = self._conn.execute(
            'INSERT INTO allele_call '
            '(created_at,created_by,supersedes_id,peak_id,'
            ' allele,marker,bin_distance,is_stutter,note)'
            ' VALUES (?,?,?,?,?,?,?,?,?)',
            (self._now(), record.created_by, record.supersedes_id,
             record.peak_id, record.allele, record.marker,
             record.bin_distance, int(record.is_stutter), record.note))
        self._conn.commit()
        return cur.lastrowid

    def store_session(self, record: SavedSessionRecord) -> int:
        cur = self._conn.execute(
            'INSERT INTO saved_session (created_at,created_by,supersedes_id,name)'
            ' VALUES (?,?,?,?)',
            (self._now(), record.created_by,
             record.supersedes_id, record.name))
        self._conn.commit()
        return cur.lastrowid

    def store_session_tab(self, record: SessionTabRecord) -> int:
        cur = self._conn.execute(
            'INSERT INTO session_tab (session_id,tab_order,run_id)'
            ' VALUES (?,?,?)',
            (record.session_id, record.tab_order, record.run_id))
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
        return dict(cur.fetchone())

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

    def find_file_by_path_and_hash(self, path: str, sha256: str) -> Optional[dict]:
        cur = self._conn.execute(
            'SELECT * FROM instrument_file'
            ' WHERE file_path=? AND hash_sha256=? LIMIT 1',
            (path, sha256))
        row = cur.fetchone()
        return dict(row) if row else None

    def find_session_by_name(self, name: str) -> Optional[dict]:
        cur = self._conn.execute(
            'SELECT * FROM current_saved_session WHERE name=? LIMIT 1', (name,))
        row = cur.fetchone()
        return dict(row) if row else None

    def hide_session(self, session_id: int) -> None:
        self._conn.execute(
            'INSERT INTO session_deletion (created_at, session_id) VALUES (?,?)',
            (self._now(), session_id))
        self._conn.commit()

    def close(self) -> None:
        self._conn.close()


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

_SCHEMA_SQL = """
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
CREATE TABLE IF NOT EXISTS analysis_run (
    id                  INTEGER PRIMARY KEY,
    created_at          TEXT    NOT NULL,
    created_by          TEXT    NOT NULL,
    supersedes_id       INTEGER REFERENCES analysis_run(id),
    file_id             INTEGER NOT NULL REFERENCES instrument_file(id),
    min_height          REAL,
    min_prominence      REAL,
    min_width           REAL,
    window_width        INTEGER,
    baseline_correction INTEGER NOT NULL DEFAULT 0,
    sizing_method       TEXT,
    size_standard       TEXT,
    panel               TEXT
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
    is_ladder       INTEGER NOT NULL DEFAULT 0
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
