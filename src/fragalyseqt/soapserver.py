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

from threading import Thread
from http.server import HTTPServer, BaseHTTPRequestHandler
from xml.etree.ElementTree import fromstring, tostring, Element, SubElement
from os.path import dirname, join

NS_SOAP = 'http://schemas.xmlsoap.org/soap/envelope/'
NS_FQ = 'http://fragalyseqt.dorif.dev/soap/v1'
_S = f'{{{NS_SOAP}}}'
_F = f'{{{NS_FQ}}}'


def _env():
    env = Element(f'{_S}Envelope',
                  {f'xmlns:soap': NS_SOAP, f'xmlns:fq': NS_FQ})
    body = SubElement(env, f'{_S}Body')
    return env, body


def _fault(code, message):
    env, body = _env()
    fault = SubElement(body, f'{_S}Fault')
    SubElement(fault, 'faultcode').text = code
    SubElement(fault, 'faultstring').text = str(message)
    return tostring(env, encoding='unicode', xml_declaration=False)


def parse_envelope(data):
    # Return (operation_name, params_dict, auth_token_or_None).
    root = fromstring(data)
    token = None
    header = root.find(f'{_S}Header')
    if header is not None:
        auth = header.find(f'{_F}Auth')
        if auth is not None:
            tok = auth.find(f'{_F}token')
            if tok is not None:
                token = tok.text
    body = root.find(f'{_S}Body')
    if body is None:
        raise ValueError('Missing SOAP Body')
    children = list(body)
    if not children:
        raise ValueError('Empty SOAP Body')
    op_el = children[0]
    local = op_el.tag.split('}')[1] if '}' in op_el.tag else op_el.tag
    params = {}
    for child in op_el:
        key = child.tag.split('}')[1] if '}' in child.tag else child.tag
        # Handle nested <params> element (SetAnalysisParams)
        if len(child):
            params[key] = {
                (c.tag.split('}')[1] if '}' in c.tag else c.tag): c.text
                for c in child
            }
        else:
            params[key] = child.text
    return local, params, token


class _SOAPHandler(BaseHTTPRequestHandler):
    # Class-level references set by SOAPServerThread before serve_forever()
    bridge = None
    token = None
    wsdl = None

    def log_message(self, fmt, *args):
        pass  # suppress per-request console output

    def do_GET(self):
        if '?wsdl' in self.path or '?WSDL' in self.path:
            self.send_response(200)
            self.send_header('Content-Type', 'text/xml; charset=utf-8')
            self.send_header('Content-Length', str(len(self.wsdl)))
            self.end_headers()
            self.wfile.write(self.wsdl)
        else:
            self.send_response(404)
            self.end_headers()

    def do_POST(self):
        length = int(self.headers.get('Content-Length', 0))
        raw = self.rfile.read(length)
        try:
            op, params, token = parse_envelope(raw)
        except Exception as exc:
            self._xml(500, _fault('soap:Server', f'Parse error: {exc}'))
            return
        if self.token and token != self.token:
            self.send_response(401)
            self.end_headers()
            return
        handler = getattr(self, f'_op_{op}', None)
        if handler is None:
            self._xml(500, _fault('soap:Client', f'Unknown operation: {op}'))
            return
        try:
            self._xml(200, handler(params))
        except IndexError:
            self._xml(500, _fault('soap:Client', 'session_id out of range'))
        except TimeoutError as exc:
            self._xml(500, _fault('soap:Server', str(exc)))
        except Exception as exc:
            self._xml(500, _fault('soap:Server', str(exc)))

    def _xml(self, code, body):
        data = body.encode('utf-8') if isinstance(body, str) else body
        self.send_response(code)
        self.send_header('Content-Type', 'text/xml; charset=utf-8')
        self.send_header('Content-Length', str(len(data)))
        self.end_headers()
        self.wfile.write(data)

    # ------------------------------------------------------------------
    # Helpers
    # ------------------------------------------------------------------

    @staticmethod
    def _resp(name):
        env, body = _env()
        resp = SubElement(body, f'{_F}{name}Response')
        return env, resp

    @staticmethod
    def _str(env):
        return tostring(env, encoding='unicode', xml_declaration=False)

    @staticmethod
    def _fields(parent, tag, d):
        el = SubElement(parent, f'{_F}{tag}')
        for k, v in d.items():
            SubElement(el, f'{_F}{k}').text = (
                'true' if v is True else 'false' if v is False else str(v))
        return el

    # ------------------------------------------------------------------
    # Read operations
    # ------------------------------------------------------------------

    def _op_GetSessionList(self, _params):
        env, resp = self._resp('GetSessionList')
        for info in self.bridge.get_session_list():
            self._fields(resp, 'session', info)
        return self._str(env)

    def _op_GetPeakTable(self, params):
        env, resp = self._resp('GetPeakTable')
        for peak in self.bridge.get_peak_table(int(params['session_id'])):
            self._fields(resp, 'peak', peak)
        return self._str(env)

    def _op_GetAnalysisParams(self, params):
        env, resp = self._resp('GetAnalysisParams')
        self._fields(resp, 'params',
                     self.bridge.get_analysis_params(int(params['session_id'])))
        return self._str(env)

    def _op_ExportCSV(self, params):
        env, resp = self._resp('ExportCSV')
        SubElement(resp, f'{_F}csv_content').text = (
            self.bridge.export_csv(int(params['session_id'])))
        return self._str(env)

    def _op_GetRawSignal(self, params):
        env, resp = self._resp('GetRawSignal')
        for v in self.bridge.get_raw_signal(int(params['session_id']),
                                            int(params['channel'])):
            SubElement(resp, f'{_F}sample').text = str(v)
        return self._str(env)

    # ------------------------------------------------------------------
    # Write operations
    # ------------------------------------------------------------------

    def _op_SubmitFile(self, params):
        sid = self.bridge.submit_file(params['file_name'],
                                      params['content_b64'])
        env, resp = self._resp('SubmitFile')
        SubElement(resp, f'{_F}session_id').text = str(sid)
        return self._str(env)

    def _op_SetAnalysisParams(self, params):
        sid = int(params['session_id'])
        nested = params.get('params', {})
        ok = self.bridge.set_analysis_params(sid, nested)
        env, resp = self._resp('SetAnalysisParams')
        SubElement(resp, f'{_F}success').text = 'true' if ok else 'false'
        return self._str(env)

    def _op_TriggerReanalysis(self, params):
        ok = self.bridge.trigger_reanalysis(int(params['session_id']))
        env, resp = self._resp('TriggerReanalysis')
        SubElement(resp, f'{_F}success').text = 'true' if ok else 'false'
        return self._str(env)

    def _op_TriggerBatchProcess(self, params):
        ok = self.bridge.trigger_batch_process(int(params['session_id']))
        env, resp = self._resp('TriggerBatchProcess')
        SubElement(resp, f'{_F}success').text = 'true' if ok else 'false'
        return self._str(env)

    def _op_CloseSession(self, params):
        ok = self.bridge.close_session(int(params['session_id']))
        env, resp = self._resp('CloseSession')
        SubElement(resp, f'{_F}success').text = 'true' if ok else 'false'
        return self._str(env)

    # ------------------------------------------------------------------
    # Panel operations
    # ------------------------------------------------------------------

    def _op_ListPanels(self, _params):
        env, resp = self._resp('ListPanels')
        for name in self.bridge.list_panels():
            SubElement(resp, f'{_F}panel_name').text = name
        return self._str(env)

    def _op_ImportPanel(self, params):
        names = self.bridge.import_panel(
            panels_name=params.get('panels_file_name', 'panel.txt'),
            panels_b64=params.get('panels_content_b64', ''),
            bins_name=params.get('bins_file_name', '') or '',
            bins_b64=params.get('bins_content_b64', '') or '',
            stutter_name=params.get('stutter_file_name', '') or '',
            stutter_b64=params.get('stutter_content_b64', '') or '',
        )
        env, resp = self._resp('ImportPanel')
        for name in names:
            SubElement(resp, f'{_F}panel_name').text = name
        return self._str(env)

    # ------------------------------------------------------------------
    # Database operations
    # ------------------------------------------------------------------

    def _op_SaveSession(self, params):
        sid = self.bridge.save_session(params.get('name', ''))
        env, resp = self._resp('SaveSession')
        SubElement(resp, f'{_F}session_id').text = str(sid)
        return self._str(env)

    def _op_ListSessions(self, _params):
        env, resp = self._resp('ListSessions')
        for s in self.bridge.list_sessions():
            self._fields(resp, 'session', {
                'session_id': s['id'],
                'name': s['name'],
                'created_at': s['created_at'],
                'created_by': s['created_by'],
            })
        return self._str(env)

    def _op_OpenSession(self, params):
        result = self.bridge.open_session(int(params['session_id']))
        env, resp = self._resp('OpenSession')
        SubElement(resp, f'{_F}success').text = 'true'
        SubElement(resp, f'{_F}readonly').text = 'true' if result['readonly'] else 'false'
        SubElement(resp, f'{_F}tabs').text = str(result['tabs'])
        return self._str(env)

    def _op_DeleteSession(self, params):
        ok = self.bridge.delete_session(int(params['session_id']))
        env, resp = self._resp('DeleteSession')
        SubElement(resp, f'{_F}success').text = 'true' if ok else 'false'
        return self._str(env)

    # ------------------------------------------------------------------
    # Frequency tables and reference profiles
    # ------------------------------------------------------------------

    def _op_ImportFreqTable(self, params):
        name = self.bridge.import_freq_table(
            params.get('file_name', 'table.csv'),
            params.get('content_b64', ''),
            params.get('table_name', ''),
        )
        env, resp = self._resp('ImportFreqTable')
        SubElement(resp, f'{_F}table_name').text = name
        return self._str(env)

    def _op_ListFreqTables(self, _params):
        env, resp = self._resp('ListFreqTables')
        for name in self.bridge.list_freq_tables():
            SubElement(resp, f'{_F}table_name').text = name
        return self._str(env)

    def _op_ListRefProfiles(self, _params):
        env, resp = self._resp('ListRefProfiles')
        for p in self.bridge.list_ref_profiles():
            self._fields(resp, 'profile', p)
        return self._str(env)

    def _op_ImportRefProfile(self, params):
        xml_b64 = params.get('xml_b64', '')
        role = params.get('role', '') or ''
        ids = self.bridge.import_ref_profile(xml_b64, role)
        env, resp = self._resp('ImportRefProfile')
        for pid in ids:
            SubElement(resp, f'{_F}profile_id').text = str(pid)
        return self._str(env)

    def _op_ExportIdentityPDF(self, params):
        pdf_b64 = self.bridge.export_identity_pdf(params)
        env, resp = self._resp('ExportIdentityPDF')
        SubElement(resp, f'{_F}pdf_b64').text = pdf_b64
        return self._str(env)

    def _op_ExportKinshipPDF(self, params):
        pdf_b64 = self.bridge.export_kinship_pdf(params)
        env, resp = self._resp('ExportKinshipPDF')
        SubElement(resp, f'{_F}pdf_b64').text = pdf_b64
        return self._str(env)

    def _op_CompareIdentity(self, params):
        result = self.bridge.compare_identity(params)
        env, resp = self._resp('CompareIdentity')
        self._comparison_result_to_xml(resp, result)
        return self._str(env)

    def _op_CompareKinship(self, params):
        result = self.bridge.compare_kinship(params)
        env, resp = self._resp('CompareKinship')
        self._comparison_result_to_xml(resp, result)
        return self._str(env)

    def _op_StoreRefProfile(self, params):
        pid = self.bridge.store_ref_profile(params)
        env, resp = self._resp('StoreRefProfile')
        SubElement(resp, f'{_F}profile_id').text = str(pid)
        return self._str(env)

    def _op_SearchProfiles(self, params):
        results = self.bridge.search_profiles_soap(params)
        env, resp = self._resp('SearchProfiles')
        for r in results:
            match_el = SubElement(resp, f'{_F}match')
            SubElement(match_el, f'{_F}profile_id').text = str(r['id'])
            SubElement(match_el, f'{_F}name').text = r['name']
            SubElement(match_el, f'{_F}role').text = r['role']
            SubElement(match_el, f'{_F}matched').text = str(r['matched'])
            SubElement(match_el, f'{_F}common').text = str(r['common'])
            status = 'exact' if r['matched'] == r['common'] else 'partial'
            SubElement(match_el, f'{_F}status').text = status
        return self._str(env)

    def _op_GetRefProfile(self, params):
        data = self.bridge.get_ref_profile(int(params['profile_id']))
        env, resp = self._resp('GetRefProfile')
        for k in ('id', 'name', 'role', 'notes', 'created_at'):
            SubElement(resp, f'{_F}{k}').text = str(data[k])
        for call in data['calls']:
            call_el = SubElement(resp, f'{_F}allele_call')
            SubElement(call_el, f'{_F}marker').text = call['marker']
            SubElement(call_el, f'{_F}allele1').text = call['allele1']
            SubElement(call_el, f'{_F}allele2').text = call['allele2']
        return self._str(env)

    def _op_UpdateRefProfile(self, params):
        profile_id = int(params.pop('profile_id'))
        new_id = self.bridge.update_ref_profile(profile_id, params)
        env, resp = self._resp('UpdateRefProfile')
        SubElement(resp, f'{_F}new_profile_id').text = str(new_id)
        return self._str(env)

    def _op_DeleteRefProfile(self, params):
        ok = self.bridge.delete_ref_profile(int(params['profile_id']))
        env, resp = self._resp('DeleteRefProfile')
        SubElement(resp, f'{_F}success').text = 'true' if ok else 'false'
        return self._str(env)

    @staticmethod
    def _comparison_result_to_xml(resp, result):
        SubElement(resp, f'{_F}combined_stat').text = str(result.combined_stat)
        SubElement(resp, f'{_F}log10_stat').text = str(result.log10_stat)
        SubElement(resp, f'{_F}verbal_scale').text = result.verbal_scale
        SubElement(resp, f'{_F}n_loci').text = str(result.n_loci)
        SubElement(resp, f'{_F}n_excluded').text = str(result.n_excluded)
        SubElement(resp, f'{_F}n_only_q').text = str(result.n_only_q)
        SubElement(resp, f'{_F}n_only_r').text = str(result.n_only_r)
        for locus in result.loci:
            loc = SubElement(resp, f'{_F}locus')
            SubElement(loc, f'{_F}marker').text = locus.marker
            SubElement(loc, f'{_F}alleles_q').text = '/'.join(
                str(a) for a in locus.alleles_q if a is not None)
            SubElement(loc, f'{_F}alleles_r').text = '/'.join(
                str(a) for a in locus.alleles_r if a is not None)
            SubElement(loc, f'{_F}freq_a1').text = f'{locus.freq_a1:.6f}'
            SubElement(loc, f'{_F}freq_a2').text = (
                f'{locus.freq_a2:.6f}' if locus.freq_a2 is not None else '')
            SubElement(loc, f'{_F}locus_stat').text = (
                f'{locus.locus_stat:.6f}' if locus.included else '')
            SubElement(loc, f'{_F}included').text = (
                'true' if locus.included else 'false')
            SubElement(loc, f'{_F}note').text = locus.note


class SOAPServerThread(Thread):
    def __init__(self, bridge, host='127.0.0.1', port=8742, token=None):
        super().__init__(daemon=True, name='SOAPServerThread')
        self._bridge = bridge
        self._host = host
        self._port = port
        self._token = token or None
        self._server = None
        wsdl_path = join(dirname(__file__), 'soap_service.wsdl')
        with open(wsdl_path, 'rb') as f:
            self._wsdl = f.read()

    def run(self):
        _SOAPHandler.bridge = self._bridge
        _SOAPHandler.token = self._token
        _SOAPHandler.wsdl = self._wsdl
        self._server = HTTPServer((self._host, self._port), _SOAPHandler)
        self._server.serve_forever()

    def stop(self):
        if self._server:
            self._server.shutdown()
            self._server = None
