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

# Guards for the XML entity expansion hardening.  The application reads XML
# it did not write -- panels, size standards, CODIS profiles, instrument
# files, SOAP bodies -- and the standard parser expands internal entities
# without any limit.  These tests pin both halves of the fix: hostile
# documents are refused, and every real document still parses.

import glob
import io
import os

import pytest

from fragalyseqt import safexml
from fragalyseqt.safexml import XMLSecurityError

_HERE = os.path.dirname(__file__)
_ROOT = os.path.join(_HERE, "..")

# A "billion laughs": each entity references the previous one ten times, so
# the document expands geometrically while staying a few hundred bytes on
# the wire.
BOMB = b"""<?xml version="1.0"?>
<!DOCTYPE lolz [
 <!ENTITY lol "lol">
 <!ENTITY lol1 "&lol;&lol;&lol;&lol;&lol;&lol;&lol;&lol;&lol;&lol;">
 <!ENTITY lol2 "&lol1;&lol1;&lol1;&lol1;&lol1;&lol1;&lol1;&lol1;&lol1;">
]>
<lolz>&lol2;</lolz>"""

# An external entity: the same declaration mechanism used to read a local
# file instead of to exhaust memory.
XXE = (b'<?xml version="1.0"?>\n'
       b'<!DOCTYPE r [ <!ENTITY x SYSTEM "file:///etc/passwd"> ]>'
       b'<r>&x;</r>')

GOOD = b'<?xml version="1.0"?><KitData><Kit><Name>GlobalFiler</Name></Kit>' \
       b'</KitData>'


@pytest.mark.parametrize("payload,label", [
    (BOMB, "internal entity bomb"),
    (XXE, "external entity reference"),
    (b"<!doctype r><r/>", "lowercase doctype"),
    (b"<!-- comment --><!DOCTYPE r><r/>", "doctype behind a comment"),
    (b'<?xml version="1.0"?>\n\n  <!DOCTYPE r>\n<r/>', "doctype after blanks"),
])
def test_hostile_documents_are_refused(payload, label):
    with pytest.raises(XMLSecurityError):
        safexml.fromstring(payload)


def test_refusal_is_a_valueerror():
    # Call sites already treat unreadable XML as a ValueError, so a refusal
    # must travel that same path instead of escaping as something unexpected.
    assert issubclass(XMLSecurityError, ValueError)
    with pytest.raises(ValueError):
        safexml.fromstring(BOMB)


@pytest.mark.parametrize("payload,label", [
    (GOOD, "plain document"),
    (b'<?xml version="1.0"?><!-- legal comment --><Kit><N>x</N></Kit>',
     "leading comment"),
    (b'<Kit><N>x</N></Kit>', "no prolog at all"),
    (b'\xef\xbb\xbf<?xml version="1.0"?><Kit><N>x</N></Kit>', "utf-8 BOM"),
])
def test_legitimate_documents_still_parse(payload, label):
    assert safexml.fromstring(payload) is not None


def test_str_input_is_accepted_like_the_stdlib():
    assert safexml.fromstring(GOOD.decode()).tag == "KitData"


def test_a_stream_is_left_where_it_was_found(tmp_path):
    # iterparse in fillfrf reads the file after the guard has inspected it,
    # so the guard must rewind: a stream left mid-prolog silently truncates
    # the document.
    buf = io.BytesIO(GOOD)
    safexml._check_stream(buf)
    assert buf.tell() == 0
    assert buf.read() == GOOD


def test_every_real_xml_asset_still_parses():
    # The blunt guard is only acceptable because nothing this program reads
    # uses a DTD.  If that ever stops being true, this test says so.
    from xml.etree.ElementTree import parse as stdparse

    files = sorted(glob.glob(os.path.join(
        _ROOT, "docs", "OSIRIS_PANELS", "*.xml")))
    files += sorted(glob.glob(os.path.join(
        _ROOT, "docs", "SPECS_AND_REFERENCES", "CODIS", "*.xml")))
    files += sorted(glob.glob(os.path.join(
        _ROOT, "src", "fragalyseqt", "**", "*.xml"), recursive=True))
    assert files, "no XML assets found -- the test would prove nothing"

    checked = 0
    for path in files:
        try:
            expected = stdparse(path).getroot().tag
        except Exception:
            continue        # already unparseable; not this guard's business
        assert safexml.parse(path).getroot().tag == expected, path
        checked += 1
    assert checked >= 50, f"only {checked} assets checked"
