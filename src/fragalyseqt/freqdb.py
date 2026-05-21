# This file is part of FragalyseQt.
#
# FragalyseQt is free software: you can redistribute it and/or modify it under
# the terms of the GNU Affero General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option)
# any later version.
#
# FragalyseQt is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE. See the GNU Affero General Public License
# for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with FragalyseQt. If not, see <https://www.gnu.org/licenses/>.

# Allele frequency table storage and import for forensic STR statistics.
#
# Tables are stored as JSON in ~/.local/share/fragalyseqt/freqtables/.
# Supported import format:
#   CSV/TSV — long format: marker, allele, frequency (delimiter auto-detected)
#
# All allele keys are normalised: whole numbers without .0 ("12"), microvariants
# with one decimal place ("9.3"). Lookup applies the same normalisation so
# "12.0" and "12" resolve to the same entry.

from __future__ import annotations
from csv import DictReader, Sniffer
from dataclasses import dataclass, asdict
from json import load as json_load, dump as json_dump
from os import makedirs, listdir
from os.path import isdir, join, dirname

_TABLE_VERSION = 1
_DEFAULT_THETA = 0.01
_DEFAULT_MIN_FREQ = 0.01


@dataclass
class FrequencyTable:
    name: str
    panel: str
    population: str
    markers: dict[str, dict[str, float]]
    theta: float = _DEFAULT_THETA
    min_freq: float = _DEFAULT_MIN_FREQ


def _normalize_allele(v) -> str:
    try:
        f = float(v)
    except (ValueError, TypeError):
        return str(v)
    return str(int(f)) if f == int(f) else str(round(f, 1))


def import_freq_csv(path: str, name: str, panel: str, population: str,
                    theta: float = _DEFAULT_THETA,
                    min_freq: float = _DEFAULT_MIN_FREQ) -> FrequencyTable:
    with open(path, newline='', encoding='utf-8-sig') as fh:
        sample = fh.read(4096)
        fh.seek(0)
        dialect = Sniffer().sniff(sample, delimiters=',\t;')
        reader = DictReader(fh, dialect=dialect)
        first = next(reader)
        col_map = {k.strip().lower(): k for k in first.keys()}
        for required in ('marker', 'allele', 'frequency'):
            if required not in col_map:
                raise ValueError(f"CSV missing required column: '{required}'")

        markers: dict[str, dict[str, float]] = {}

        def _add(row):
            m = row[col_map['marker']].strip()
            a = _normalize_allele(row[col_map['allele']].strip())
            f = float(row[col_map['frequency']].strip())
            if not (0.0 <= f <= 1.0):
                raise ValueError(f"Frequency {f} out of [0,1] for {m}:{a}")
            markers.setdefault(m, {})[a] = f

        _add(first)
        for row in reader:
            _add(row)

    return FrequencyTable(name=name, panel=panel, population=population,
                          markers=markers, theta=theta, min_freq=min_freq)


def save_freq_table(table: FrequencyTable, path: str) -> None:
    makedirs(dirname(path) or '.', exist_ok=True)
    payload = {'_version': _TABLE_VERSION}
    payload.update(asdict(table))
    with open(path, 'w', encoding='utf-8') as fh:
        json_dump(payload, fh, indent=2, ensure_ascii=False)


def load_freq_table(path: str) -> FrequencyTable:
    with open(path, encoding='utf-8') as fh:
        data = json_load(fh)
    if data.get('_version') != _TABLE_VERSION:
        raise ValueError(f"Unsupported table version in {path}")
    data.pop('_version')
    return FrequencyTable(**data)


def list_freq_tables(data_dir: str) -> list[str]:
    d = join(data_dir, 'freqtables')
    if not isdir(d):
        return []
    return sorted(join(d, f) for f in listdir(d) if f.endswith('.json'))


def get_allele_freq(table: FrequencyTable, marker: str, allele: str) -> float:
    m = table.markers.get(marker)
    if m is None:
        return table.min_freq
    return m.get(_normalize_allele(allele), table.min_freq)
