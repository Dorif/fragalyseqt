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

# Reference profile storage and retrieval.
#
# Profiles are stored in an append-only SQLite database via DatabaseBackend.
# Logical updates use the supersedes_id pattern; logical deletes insert a
# tombstone into reference_profile_deletion. The current_reference_profile
# view hides superseded and deleted rows.

from __future__ import annotations
from dataclasses import dataclass, field
from .forensicstats import AlleleCall


@dataclass
class ReferenceProfile:
    name: str
    calls: list[AlleleCall]
    role: str | None = None
    notes: str | None = None
    session_id: int | None = None
    id: int | None = None
    created_at: str | None = None


def store_profile(db, profile: ReferenceProfile) -> int:
    with db.transaction():
        pid = db.store_reference_profile(profile.name, profile.role,
                                         profile.notes, profile.session_id)
        db.store_reference_alleles(pid, _calls_to_dicts(profile.calls))
    return pid


def update_profile(db, profile_id: int, **kwargs) -> int:
    current = db.get_reference_profile(profile_id)
    new_name = kwargs.get('name', current['name'])
    new_role = kwargs.get('role', current['role'])
    new_notes = kwargs.get('notes', current['notes'])
    calls = kwargs.get('calls')

    with db.transaction():
        new_id = db.store_reference_profile(
            new_name, new_role, new_notes, current['session_id'],
            supersedes_id=profile_id)
        alleles = (_calls_to_dicts(calls) if calls is not None
                   else db.get_reference_alleles(profile_id))
        db.store_reference_alleles(new_id, alleles)
    return new_id


def delete_profile(db, profile_id: int) -> None:
    db.delete_reference_profile(profile_id)


def get_profile(db, profile_id: int) -> ReferenceProfile:
    row = db.get_reference_profile(profile_id)
    alleles = db.get_reference_alleles(profile_id)
    return ReferenceProfile(id=row['id'], name=row['name'], role=row['role'],
                            notes=row['notes'], session_id=row['session_id'],
                            created_at=row['created_at'],
                            calls=[AlleleCall(marker=a['marker'],
                                              allele1=a['allele1'],
                                              allele2=a['allele2'])
                            for a in alleles],)


def list_profiles(db) -> list[dict]:
    return db.list_reference_profiles()


def profile_from_state(state, name: str,
                       role: str | None = None,
                       session_id: int | None = None) -> ReferenceProfile:
    from .comparison import allele_calls_from_state
    return ReferenceProfile(name=name, role=role, session_id=session_id,
                            calls=allele_calls_from_state(state),)


def profiles_from_codis_xml(path: str) -> list[ReferenceProfile]:
    from xml.etree.ElementTree import parse as _parse
    _NS = 'urn:CODISImportFile-schema'
    _N = f'{{{_NS}}}'

    root = _parse(path).getroot()
    ns = _N if root.find(f'.//{_N}SPECIMEN') is not None else ''

    def _t(name):
        return f'{ns}{name}'

    results = []
    for spec in root.iter(_t('SPECIMEN')):
        sid_el = spec.find(_t('SPECIMENID'))
        cat_el = spec.find(_t('SPECIMENCATEGORY'))
        name = sid_el.text.strip() if sid_el is not None and sid_el.text else 'Unknown'
        role = cat_el.text.strip() if cat_el is not None and cat_el.text else None

        calls = []
        for locus in spec.iter(_t('LOCUS')):
            ln_el = locus.find(_t('LOCUSNAME'))
            if ln_el is None or not ln_el.text:
                continue
            marker = ln_el.text.strip()
            alleles = []
            for al_el in locus.iter(_t('ALLELE')):
                val_el = al_el.find(_t('ALLELEVALUE'))
                if val_el is not None and val_el.text:
                    alleles.append(val_el.text.strip())
            if alleles:
                calls.append(AlleleCall(
                    marker=marker,
                    allele1=alleles[0],
                    allele2=alleles[1] if len(alleles) > 1 else None,
                ))

        if calls:
            results.append(ReferenceProfile(name=name, role=role, calls=calls))

    return results


def _calls_to_dicts(calls: list[AlleleCall]) -> list[dict]:
    return [{'marker': c.marker, 'allele1': c.allele1, 'allele2': c.allele2}
            for c in calls]
