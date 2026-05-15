from threading import Thread
from http.server import HTTPServer, BaseHTTPRequestHandler
from xml.etree.ElementTree import fromstring, tostring, Element, SubElement
from os.path import dirname, join

NS_SOAP = 'http://schemas.xmlsoap.org/soap/envelope/'
NS_FQ   = 'http://fragalyseqt.dorif.dev/soap/v1'
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
    """Return (operation_name, params_dict, auth_token_or_None)."""
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
    token  = None
    wsdl   = None

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
                'name':       s['name'],
                'created_at': s['created_at'],
                'created_by': s['created_by'],
            })
        return self._str(env)

    def _op_OpenSession(self, params):
        result = self.bridge.open_session(int(params['session_id']))
        env, resp = self._resp('OpenSession')
        SubElement(resp, f'{_F}success').text  = 'true'
        SubElement(resp, f'{_F}readonly').text = 'true' if result['readonly'] else 'false'
        SubElement(resp, f'{_F}tabs').text     = str(result['tabs'])
        return self._str(env)

    def _op_DeleteSession(self, params):
        ok = self.bridge.delete_session(int(params['session_id']))
        env, resp = self._resp('DeleteSession')
        SubElement(resp, f'{_F}success').text = 'true' if ok else 'false'
        return self._str(env)


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
        _SOAPHandler.token  = self._token
        _SOAPHandler.wsdl   = self._wsdl
        self._server = HTTPServer((self._host, self._port), _SOAPHandler)
        self._server.serve_forever()

    def stop(self):
        if self._server:
            self._server.shutdown()
            self._server = None
