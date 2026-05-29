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


def import_freq_fam(path: str, name: str, panel: str = '',
                    population: str = '',
                    theta: float = _DEFAULT_THETA,
                    min_freq: float = _DEFAULT_MIN_FREQ) -> FrequencyTable:
    with open(path, encoding='utf-8', errors='replace') as fh:
        lines = [ln.rstrip('\n\r') for ln in fh]

    def _unquote(s: str) -> str:
        s = s.strip()
        if len(s) >= 2 and s[0] == '"' and s[-1] == '"':
            return s[1:-1]
        return s

    pop_name = population
    i = 0
    n = len(lines)

    while i < n:
        if lines[i].strip().startswith('#TRUE#'):
            if i + 1 < n:
                candidate = _unquote(lines[i + 1])
                if candidate and not candidate.startswith('#'):
                    pop_name = candidate
                    i += 2
                else:
                    i += 1
            else:
                i += 1
            break
        i += 1

    if not population:
        population = pop_name

    markers: dict[str, dict[str, float]] = {}

    while i < n:
        line = lines[i].strip()
        if not (line.startswith('"') and line.endswith('"') and len(line) > 2):
            i += 1
            continue

        marker_name = _unquote(line)
        i += 12

        if i >= n:
            break

        count_line = lines[i].strip()
        i += 1
        try:
            allele_count = int(count_line.split()[0])
        except (ValueError, IndexError):
            continue

        alleles: dict[str, float] = {}
        for _ in range(allele_count):
            if i + 1 >= n:
                break
            allele_key = _normalize_allele(_unquote(lines[i].strip()))
            try:
                freq = float(lines[i + 1].strip())
                alleles[allele_key] = freq
            except ValueError:
                pass
            i += 2

        if alleles:
            markers[marker_name] = alleles

    if not markers:
        raise ValueError(f"No frequency data found in {path}")

    return FrequencyTable(name=name, panel=panel, population=population,
                          markers=markers, theta=theta, min_freq=min_freq)


def import_freq_gm(path: str, name: str, panel: str = '',
                   population: str = '',
                   theta: float = _DEFAULT_THETA,
                   min_freq: float = _DEFAULT_MIN_FREQ) -> FrequencyTable:
    # Parse GeneMarker / GeneMarkerHID allele frequency table format.
    # Header section (tab-separated key/value pairs):
    #   POPULATION  <population name>
    #   PANELNAME   <panel name>
    #   VERSION     <version>
    # Followed by one or more blank/ignored lines, then marker blocks:
    #   MARKER      <marker name>
    #   <allele>    <frequency>
    #   ...
    # Blank lines or MARKER lines signal the start of the next block.
    with open(path, encoding='utf-8-sig', errors='replace') as fh:
        lines = [ln.rstrip('\r\n') for ln in fh]

    pop_name = population
    panel_name = panel

    markers: dict[str, dict[str, float]] = {}
    current_marker: str | None = None

    for line in lines:
        if not line.strip():
            continue
        parts = line.split('\t')
        key = parts[0].strip().upper()

        if key == 'POPULATION' and len(parts) >= 2:
            if not population:
                pop_name = parts[1].strip()
            continue
        if key == 'PANELNAME' and len(parts) >= 2:
            if not panel:
                panel_name = parts[1].strip()
            continue
        if key == 'VERSION':
            continue
        if key == 'MARKER' and len(parts) >= 2:
            current_marker = parts[1].strip()
            markers.setdefault(current_marker, {})
            continue
        if current_marker is None:
            continue
        if len(parts) >= 2:
            try:
                allele = _normalize_allele(parts[0].strip())
                freq = float(parts[1].strip())
            except (ValueError, TypeError):
                continue
            # GeneMarker files include a trailing "N\t<count>" line per
            # marker that records the sample size — skip any row whose
            # value is outside [0, 1] rather than treating it as an error.
            if not (0.0 <= freq <= 1.0):
                continue
            markers[current_marker][allele] = freq

    if not markers:
        raise ValueError(f"No frequency data found in {path}")

    return FrequencyTable(name=name,
                          panel=panel_name,
                          population=pop_name,
                          markers=markers,
                          theta=theta,
                          min_freq=min_freq)


def _is_gm_format(path: str) -> bool:
    # Quick probe: GeneMarker files start with a POPULATION header line.
    try:
        with open(path, encoding='utf-8-sig', errors='replace') as fh:
            first = fh.readline()
        return first.upper().startswith('POPULATION')
    except OSError:
        return False


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
