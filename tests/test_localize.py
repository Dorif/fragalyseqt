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

import pytest

from fragalyseqt.localize import localizefq

# A missing or minimal locale is the normal state inside containers, cron jobs
# and minimal installs. The application must still start there, in English.
NEUTRAL_LOCALES = [None, '', 'C', 'c', 'POSIX', 'posix', 'C.UTF-8', 'C.utf8',
                   'C.UTF8', 'POSIX.UTF-8', 'en_US.UTF-8', 'en_GB.UTF-8']

# Locales the project has no translation for. They must fall back to
# English rather than leaving the interface dictionary empty.
UNSUPPORTED_LOCALES = ['es_ES.UTF-8', 'cs_CZ.UTF-8', 'it_IT.UTF-8',
                       'pl_PL.UTF-8', 'ja_JP.UTF-8', 'tr_TR.UTF-8',
                       'nl_NL.UTF-8', 'pt_BR.UTF-8']

TRANSLATED_LOCALES = ['ru_RU.UTF-8', 'uk_UA.UTF-8', 'ro_RO.UTF-8',
                      'fr_FR.UTF-8', 'bg_BG.UTF-8', 'de_DE.UTF-8',
                      'zh_CN.UTF-8']


def _localize(monkeypatch, lang):
    monkeypatch.setattr('os.name', 'posix')
    if lang is None:
        monkeypatch.delenv('LANG', raising=False)
    else:
        monkeypatch.setenv('LANG', lang)
    msg = {}
    localizefq(msg)
    return msg


@pytest.mark.parametrize('lang', NEUTRAL_LOCALES)
def test_neutral_locale_starts_up_in_english(lang, monkeypatch):
    # Regression: an unset $LANG used to crash on startup with
    # AttributeError: 'NoneType' object has no attribute 'lower', and a bare
    # LANG=C loaded no strings at all, so every lookup raised KeyError.
    msg = _localize(monkeypatch, lang)
    assert msg, f'no interface strings loaded for LANG={lang!r}'
    assert msg['aboutbtn'] == 'About'


@pytest.mark.parametrize('lang', TRANSLATED_LOCALES)
def test_translated_locale_is_not_swallowed_by_the_neutral_branch(lang,
                                                                 monkeypatch):
    # The neutral branch matches "c" and "posix" exactly, so real language
    # codes containing those letters must still reach their own block.
    msg = _localize(monkeypatch, lang)
    assert msg, f'no interface strings loaded for LANG={lang!r}'
    assert msg['aboutbtn'] != 'About', f'{lang} fell through to English'


@pytest.mark.parametrize('lang', NEUTRAL_LOCALES + TRANSLATED_LOCALES)
def test_keys_used_by_the_profile_dialogs_exist(lang, monkeypatch):
    # These are looked up unconditionally when the dialogs open, so a missing
    # key is a crash rather than an untranslated label.
    msg = _localize(monkeypatch, lang)
    for key in ('save_profile_db_casework', 'save_profile_db_refprofiles',
                'search_profile_dlg', 'search_profile_query',
                'search_profile_nothing_to_search',
                'search_profile_col_name', 'search_profile_col_role',
                'search_profile_col_match', 'search_profile_col_status'):
        assert key in msg, f'{key} missing for LANG={lang!r}'


@pytest.mark.parametrize('lang', UNSUPPORTED_LOCALES)
def test_unsupported_locale_falls_back_to_english(lang, monkeypatch):
    # Regression: locales without a translation block used to leave `iface`
    # empty, so the application crashed on the first interface lookup.
    msg = _localize(monkeypatch, lang)
    assert msg, f'no interface strings loaded for LANG={lang!r}'
    assert msg['aboutbtn'] == 'About'


@pytest.mark.parametrize('lang', UNSUPPORTED_LOCALES + NEUTRAL_LOCALES)
def test_fallback_provides_the_same_keys_as_english(lang, monkeypatch):
    # The fallback must be the complete English set, not a partial one.
    english = _localize(monkeypatch, 'en_US.UTF-8')
    assert set(_localize(monkeypatch, lang)) == set(english)
