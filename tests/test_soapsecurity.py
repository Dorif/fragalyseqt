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

# The SOAP endpoint is the only way into this application from the network,
# and XML parsing is its most attackable step.  These tests pin the order of
# operations: a request is sized, then authenticated, and only then parsed.
# Counting parser entries is the point -- a status code alone cannot tell
# "refused before parsing" apart from "refused after parsing".

import threading
import urllib.error
import urllib.request
from http.server import HTTPServer

import pytest

from fragalyseqt import soapserver as ss

NS = 'http://fragalyseqt.dorif.dev/soap/v1'
TOKEN = 'S3CRET-token'

BOMB = b"""<?xml version="1.0"?>
<!DOCTYPE lolz [
 <!ENTITY lol "lol">
 <!ENTITY lol1 "&lol;&lol;&lol;&lol;&lol;&lol;&lol;&lol;&lol;&lol;">
]>
<soap:Envelope xmlns:soap="http://schemas.xmlsoap.org/soap/envelope/">
<soap:Body><fq:GetSessionList xmlns:fq="%s"/></soap:Body>
</soap:Envelope>""" % NS.encode()


def envelope(token=None, op='GetSessionList'):
    """A well-formed SOAP request, optionally carrying an auth token."""
    header = ''
    if token is not None:
        header = (f'<soap:Header><fq:Auth xmlns:fq="{NS}">'
                  f'<fq:token>{token}</fq:token></fq:Auth></soap:Header>')
    return (f'<?xml version="1.0"?><soap:Envelope '
            f'xmlns:soap="http://schemas.xmlsoap.org/soap/envelope/">{header}'
            f'<soap:Body><fq:{op} xmlns:fq="{NS}"/></soap:Body>'
            f'</soap:Envelope>').encode()


class _Bridge:
    """Enough of the real bridge for GetSessionList to succeed."""

    def get_session_list(self):
        return [{'session_id': 0, 'file_name': 'sample.hid'}]


class _Server:
    """A live endpoint plus a counter of how often the parser ran."""

    def __init__(self, token):
        self.parses = 0
        self._real = ss.parse_envelope

        def counting_parse(data):
            self.parses += 1
            return self._real(data)

        ss.parse_envelope = counting_parse
        ss._SOAPHandler.bridge = _Bridge()
        ss._SOAPHandler.token = token
        ss._SOAPHandler.wsdl = b'<wsdl/>'
        self._http = HTTPServer(('127.0.0.1', 0), ss._SOAPHandler)
        self.url = f'http://127.0.0.1:{self._http.server_address[1]}/svc'
        self._thread = threading.Thread(target=self._http.serve_forever,
                                        daemon=True)
        self._thread.start()

    def post(self, body, headers=None, content_length=None):
        head = {'Content-Type': 'text/xml'}
        head.update(headers or {})
        if content_length is not None:
            head['Content-Length'] = str(content_length)
        request = urllib.request.Request(self.url, data=body, headers=head,
                                         method='POST')
        try:
            with urllib.request.urlopen(request, timeout=10) as response:
                return response.status
        except urllib.error.HTTPError as exc:
            return exc.code

    def close(self):
        self._http.shutdown()
        self._http.server_close()
        ss.parse_envelope = self._real


@pytest.fixture
def guarded():
    """An endpoint that requires a token."""
    server = _Server(TOKEN)
    yield server
    server.close()


@pytest.fixture
def open_endpoint():
    """An endpoint with no token configured, as when the user sets none."""
    server = _Server(None)
    yield server
    server.close()


@pytest.mark.parametrize("body,label", [
    (BOMB, "entity bomb"),
    (envelope(), "well-formed, no token"),
    (envelope('WRONG'), "wrong token"),
    (envelope('pässwörd'), "non-ascii token"),
    (b'not xml at all', "not xml"),
])
def test_unauthenticated_requests_never_reach_the_parser(guarded, body, label):
    # The heart of the fix.  Whatever an anonymous caller sends -- a bomb, a
    # valid document or rubbish -- the parser must not run for it.
    assert guarded.post(body) == 401
    assert guarded.parses == 0, f"parser ran for an anonymous {label}"


@pytest.mark.parametrize("headers,label", [
    ({}, "token inside the SOAP header"),
    ({'X-Auth-Token': TOKEN}, "token in X-Auth-Token"),
    ({'Authorization': f'Bearer {TOKEN}'}, "token in Authorization"),
])
def test_authenticated_requests_are_served(guarded, headers, label):
    body = envelope(TOKEN) if not headers else envelope()
    assert guarded.post(body, headers) == 200, label
    assert guarded.parses == 1


def test_an_authenticated_bomb_is_still_refused(guarded):
    # Authentication decides who may be parsed, not what may be parsed: a
    # caller holding the token still cannot feed the parser a bomb.
    authenticated = BOMB.replace(
        b'<soap:Body>',
        f'<soap:Header><fq:Auth xmlns:fq="{NS}"><fq:token>{TOKEN}</fq:token>'
        f'</fq:Auth></soap:Header><soap:Body>'.encode())
    assert guarded.post(authenticated) == 400


def test_an_oversized_body_is_refused_before_it_is_read(guarded):
    # Announcing a body larger than the cap must be rejected outright, so an
    # anonymous caller cannot make the process receive an arbitrary payload.
    assert guarded.post(b'x', content_length=ss.MAX_BODY_BYTES + 1) == 413
    assert guarded.parses == 0


def test_a_bomb_is_refused_even_with_no_token_configured(open_endpoint):
    # Most installations run without a token on loopback.  The XML guard has
    # to protect them too -- authentication is a second line, not the only.
    assert open_endpoint.post(BOMB) == 400


def test_an_open_endpoint_still_serves_ordinary_requests(open_endpoint):
    assert open_endpoint.post(envelope()) == 200
    assert open_endpoint.parses == 1


@pytest.mark.parametrize("headers,expected", [
    ({'X-Auth-Token': TOKEN}, TOKEN),
    ({'Authorization': f'Bearer {TOKEN}'}, TOKEN),
    ({'Authorization': f'bearer {TOKEN}'}, TOKEN),
    ({'Authorization': f'Basic {TOKEN}'}, None),
    ({}, None),
])
def test_token_from_headers(headers, expected):
    assert ss.token_from_headers(headers) == expected


def test_token_from_body_scans_without_parsing():
    assert ss.token_from_body(envelope(TOKEN)) == TOKEN
    assert ss.token_from_body(envelope()) is None
    # The scan is bounded: a token buried past the limit is not honoured,
    # which keeps the scan itself from becoming a denial of service.
    padded = b'<pad>' + b'.' * ss._TOKEN_SCAN_BYTES + envelope(TOKEN)
    assert ss.token_from_body(padded) is None


def test_token_comparison_rejects_mismatches_without_raising():
    assert ss.token_matches(TOKEN, TOKEN) is True
    assert ss.token_matches('other', TOKEN) is False
    assert ss.token_matches(None, TOKEN) is False
    assert ss.token_matches('', TOKEN) is False
    # A non-ascii token must fail the comparison, not crash it: the standard
    # constant-time helper refuses non-ascii strings outright.
    assert ss.token_matches('pässwörd', TOKEN) is False
