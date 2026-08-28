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

"""Hardened XML entry points, built on the standard library alone.

Every XML document this application reads comes from outside it: analysis
panels exchanged between laboratories, size standards, CODIS profiles,
instrument files and SOAP request bodies.  The standard parser expands
internal entities without any limit, so a few hundred bytes of nested
declarations expand into gigabytes and take the process down with them --
the "billion laughs" denial of service.

The guard here is deliberately blunt: a document carrying a document type
declaration is refused before the parser ever sees it.  Entities can only
be declared inside a DTD, so removing the DTD removes internal expansion
and external entity references (XXE) in one step, with no dependency to
install and no parser internals to patch -- Python 3.14 no longer exposes
the expat handlers that the classic mitigation relied on.

The cost of being blunt is zero here: none of the formats this program
reads uses a DTD, which was verified against every XML asset shipped with
the project.

Function signatures mirror ``xml.etree.ElementTree`` exactly, so call sites
change only their import line.
"""

from xml.etree.ElementTree import (parse as _parse, fromstring as _fromstring,
                                   iterparse as _iterparse)

# How much of the document may be inspected while looking for a DTD.  A
# real prolog is a few dozen bytes; the allowance is generous purely so a
# large comment ahead of the root element is not mistaken for an attack.
_PROLOG_LIMIT = 1 << 20
_CHUNK = 1 << 16


class XMLSecurityError(ValueError):
    """Raised for a document refused on security grounds.

    Derived from ValueError because the call sites already treat malformed
    XML as a ValueError, so a refusal travels the same path as any other
    unreadable file instead of escaping as an unexpected error.
    """


def _as_text(head):
    """Decode enough of a document to scan its prolog.

    Errors are replaced rather than raised: this is a scan for one ASCII
    keyword, not a decode of the payload, and the real parser is the one
    entitled to complain about a broken encoding.
    """
    if isinstance(head, str):
        return head
    if head[:2] in (b'\xff\xfe', b'\xfe\xff'):
        return head.decode('utf-16', errors='replace')
    if head[:3] == b'\xef\xbb\xbf':
        return head[3:].decode('utf-8', errors='replace')
    return head.decode('utf-8', errors='replace')


def _verdict(head, eof):
    """Look for a DTD in the prolog.

    Returns True when a document type declaration is present, False when
    the root element was reached without one, and None when the text ran
    out mid-prolog and more of it is needed to decide.
    """
    text = _as_text(head)
    i, n = 0, len(text)
    while True:
        while i < n and text[i].isspace():
            i += 1
        if i >= n:
            return False if eof else None
        if text[i] != '<':
            # Character data before the root element: not valid XML at all.
            # Leave the diagnosis to the parser, which words it better.
            return False
        if text.startswith('<?', i):
            end = text.find('?>', i)
            if end < 0:
                return False if eof else None
            i = end + 2
        elif text.startswith('<!--', i):
            end = text.find('-->', i)
            if end < 0:
                return False if eof else None
            i = end + 3
        elif text[i:i + 9].upper() == '<!DOCTYPE':
            return True
        elif len(text) - i < 9 and not eof:
            return None
        else:
            # The root element starts here, so the prolog is over.
            return False


def _refuse():
    raise XMLSecurityError(
        'XML document type declarations are not accepted: they can define '
        'entities that expand without limit. Remove the <!DOCTYPE ...> '
        'declaration from the file.')


def _check_text(data):
    """Guard a document already held in memory."""
    if _verdict(data[:_PROLOG_LIMIT], eof=True):
        _refuse()


def _check_stream(fh):
    """Guard an open binary stream, leaving it ready to be parsed.

    The prolog is read in chunks and the stream is rewound afterwards, so
    the parser still sees the document from its first byte.  A stream that
    cannot seek is refused rather than silently parsed unguarded: that is
    a deliberate choice, since every source this program uses is either a
    real file or an in-memory buffer.
    """
    start = fh.tell()
    head = None
    answer = None
    while answer is None:
        chunk = fh.read(_CHUNK)
        if not chunk:
            answer = _verdict(head, eof=True) if head else False
            break
        head = chunk if head is None else head + chunk
        if len(head) >= _PROLOG_LIMIT:
            answer = _verdict(head, eof=True)
            break
        answer = _verdict(head, eof=False)
    fh.seek(start)
    if answer:
        _refuse()


def _check_source(source):
    """Guard whatever ElementTree accepts as a source."""
    if hasattr(source, 'read'):
        _check_stream(source)
        return
    with open(source, 'rb') as fh:
        _check_stream(fh)


def parse(source, parser=None):
    """As ``ElementTree.parse``, refusing documents that declare a DTD."""
    _check_source(source)
    return _parse(source, parser)


def fromstring(text, parser=None):
    """As ``ElementTree.fromstring``, refusing documents that declare a DTD."""
    _check_text(text)
    return _fromstring(text, parser)


def iterparse(source, events=None, parser=None):
    """As ``ElementTree.iterparse``, refusing documents that declare a DTD."""
    _check_source(source)
    if parser is None:
        return _iterparse(source, events)
    return _iterparse(source, events, parser)
