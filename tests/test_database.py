import os
import tempfile
import numpy as np
import pytest
from fragalyseqt.database import (
    compress_signal, decompress_signal, SQLiteBackend,
    SavedSessionRecord, InstrumentFileRecord, AnalysisRunRecord,
    PeakCallRecord, DyeChannelRecord, SessionTabRecord,
)


@pytest.fixture
def db():
    with tempfile.NamedTemporaryFile(suffix='.db', delete=False) as f:
        path = f.name
    backend = SQLiteBackend(path)
    yield backend
    backend.close()
    os.unlink(path)


# ---------------------------------------------------------------------------
# compress / decompress
# ---------------------------------------------------------------------------

def test_decompress_returns_ndarray():
    blob = compress_signal([1.0, 2.0, 3.0])
    result = decompress_signal(blob)
    assert isinstance(result, np.ndarray), "decompress_signal must return ndarray"


def test_compress_decompress_roundtrip():
    signal = np.array([0.0, 100.5, 1234.7, 65535.0], dtype=float)
    result = decompress_signal(compress_signal(signal))
    np.testing.assert_allclose(result, signal, rtol=1e-5)


def test_compress_accepts_ndarray():
    signal = np.array([1.0, 2.0, 3.0], dtype=np.float64)
    blob = compress_signal(signal)
    assert isinstance(blob, bytes)
    assert len(blob) > 0


def test_compress_accepts_list():
    blob = compress_signal([1.0, 2.0, 3.0])
    assert isinstance(blob, bytes)


# ---------------------------------------------------------------------------
# SQLiteBackend — basic operations
# ---------------------------------------------------------------------------

def test_store_and_retrieve_session(db):
    sid = db.store_session(SavedSessionRecord(created_by='test', name='my_session'))
    sessions = db.get_session_list()
    assert any(s['name'] == 'my_session' for s in sessions)


def test_find_session_by_name(db):
    db.store_session(SavedSessionRecord(created_by='test', name='find_me'))
    row = db.find_session_by_name('find_me')
    assert row is not None
    assert row['name'] == 'find_me'


def test_find_session_by_name_missing(db):
    assert db.find_session_by_name('nonexistent') is None


def test_hide_session(db):
    sid = db.store_session(SavedSessionRecord(created_by='test', name='to_hide'))
    db.hide_session(sid)
    sessions = db.get_session_list()
    assert not any(s['name'] == 'to_hide' for s in sessions)


def test_store_file_deduplicates(db):
    rec = InstrumentFileRecord(
        created_by='test', file_name='x.fsa', file_path='/tmp/x.fsa',
        file_size=100, hash_md5='a', hash_sha1='b',
        hash_sha256='c', hash_sha3_256='d',
    )
    id1 = db.store_file(rec)
    id2 = db.store_file(rec)
    assert id1 == id2


def _make_run(db):
    fid = db.store_file(InstrumentFileRecord(
        created_by='t', file_name='s.fsa', file_path='/tmp/sat.fsa',
        file_size=1, hash_md5='m', hash_sha1='s1', hash_sha256='s2',
        hash_sha3_256='s3'))
    return db.store_analysis_run(AnalysisRunRecord(
        file_id=fid, created_by='t', min_height=175, min_prominence=175,
        min_width=4, halfwindow_width=5, baseline_correction=False,
        sizing_method='Local Southern', size_standard='GS600LIZ', panel=''))


def _peak(run_id, saturated):
    return PeakCallRecord(
        run_id=run_id, created_by='t', channel=1, dye_name='FAM',
        position_dp=1234.5, position_bp=100.0, height=2500.0, area=5000.0,
        fwhm=8.0, is_ladder=False, saturated=saturated)


def test_saturated_flag_roundtrip(db):
    run_id = _make_run(db)
    db.store_peak_call(_peak(run_id, saturated=True))
    db.store_peak_call(_peak(run_id, saturated=False))
    rows = db.get_peak_calls_for_run(run_id)
    assert len(rows) == 2
    assert {bool(r['saturated']) for r in rows} == {True, False}


def test_saturated_migration_appends_column(db):
    # Simulate a pre-existing DB created before the 'saturated' column: drop it
    # via a rebuilt table without the column, then re-open and confirm the
    # append-only ALTER migration adds it (defaulting to 0) and the view
    # exposes it.
    path = db._conn.execute('PRAGMA database_list').fetchall()[0][2]
    run_id = _make_run(db)
    db.store_peak_call(_peak(run_id, saturated=True))
    # Emulate the old schema: drop the dependent views first, then the column.
    db._conn.executescript(
        'DROP VIEW IF EXISTS current_allele_call;'
        'DROP VIEW IF EXISTS current_peak_call;'
        'ALTER TABLE peak_call DROP COLUMN saturated;')
    db._conn.commit()
    db.close()
    reopened = SQLiteBackend(path)        # triggers the append-only migration
    try:
        rows = reopened.get_peak_calls_for_run(run_id)
        assert len(rows) == 1
        assert 'saturated' in rows[0]            # view exposes the new column
        assert bool(rows[0]['saturated']) is False   # old row defaults to 0
        # a fresh insert still records saturation correctly post-migration
        reopened.store_peak_call(_peak(run_id, saturated=True))
        sat = [bool(r['saturated'])
               for r in reopened.get_peak_calls_for_run(run_id)]
        assert sat.count(True) == 1 and sat.count(False) == 1
    finally:
        reopened.close()


def test_store_dye_channels(db):
    rec = InstrumentFileRecord(
        created_by='test', file_name='x.fsa', file_path='/tmp/x2.fsa',
        file_size=100, hash_md5='e', hash_sha1='f',
        hash_sha256='g', hash_sha3_256='h',
    )
    fid = db.store_file(rec)
    db.store_dye_channel(DyeChannelRecord(file_id=fid, channel=1, dye_name='FAM'))
    db.store_dye_channel(DyeChannelRecord(file_id=fid, channel=2, dye_name='VIC'))
    dyes = db.get_dye_channels(fid)
    assert len(dyes) == 2
    assert dyes[0]['dye_name'] == 'FAM'
    assert dyes[1]['dye_name'] == 'VIC'


def test_store_channel_signal_roundtrip(db):
    rec = InstrumentFileRecord(
        created_by='test', file_name='x.fsa', file_path='/tmp/x3.fsa',
        file_size=100, hash_md5='i', hash_sha1='j',
        hash_sha256='k', hash_sha3_256='l',
    )
    fid = db.store_file(rec)
    signal = np.array([1.0, 2.0, 3.0, 4.0])
    from fragalyseqt.database import ChannelSignalRecord
    db.store_channel_signal(ChannelSignalRecord(
        file_id=fid, channel=1, signal=compress_signal(signal)))
    rows = db.get_channel_signals(fid)
    assert len(rows) == 1
    result = decompress_signal(rows[0]['signal'])
    assert isinstance(result, np.ndarray)
    np.testing.assert_allclose(result, signal, rtol=1e-5)


# ---------------------------------------------------------------------------
# Transaction context manager
# ---------------------------------------------------------------------------

def test_transaction_commits_all(db):
    with db.transaction():
        sid = db.store_session(SavedSessionRecord(created_by='test', name='tx_session'))
    assert db.find_session_by_name('tx_session') is not None


def test_transaction_rolls_back_on_error(db):
    with pytest.raises(ValueError):
        with db.transaction():
            db.store_session(SavedSessionRecord(created_by='test', name='will_rollback'))
            raise ValueError("simulated error")
    assert db.find_session_by_name('will_rollback') is None


def test_auto_commit_restored_after_transaction(db):
    with db.transaction():
        pass
    assert db._auto_commit is True


def test_nested_independent_stores_without_transaction(db):
    # Individual store calls outside a transaction should still commit
    db.store_session(SavedSessionRecord(created_by='test', name='solo'))
    assert db.find_session_by_name('solo') is not None
