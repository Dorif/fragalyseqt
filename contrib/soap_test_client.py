#!/usr/bin/env python3
"""FragalyseQt SOAP API test client.

Requires zeep:
  Debian/Ubuntu:  sudo apt install python3-zeep
  Fedora:         sudo dnf install python3-zeep
  Arch Linux:     sudo pacman -S python-zeep
  openSUSE:       sudo zypper install python-zeep
  other:          pip install zeep

Usage:
  1. Start FragalyseQt and open one or more files.
  2. Go to Settings -> SOAP API..., enable the API, click OK.
  3. Run this script:
       python3 contrib/soap_test_client.py
"""

import sys

BASE    = 'http://127.0.0.1:8742/FragalyseQtService'
WSDL    = f'{BASE}?wsdl'
BINDING = '{http://fragalyseqt.dorif.dev/soap/v1}FragalyseQtBinding'


def main():
    try:
        from zeep import Client
    except ImportError:
        print('zeep is not installed. See usage at the top of this file.')
        sys.exit(1)

    print(f'Connecting to {WSDL} ...')
    try:
        client = Client(WSDL)
        # Override endpoint so zeep uses the actual running port,
        # not the default 8742 hardcoded in the WSDL.
        svc = client.create_service(BINDING, BASE)
    except Exception as exc:
        print(f'Connection failed: {exc}')
        print('Is FragalyseQt running with SOAP API enabled?')
        sys.exit(1)

    print('\n=== GetSessionList ===')
    sessions = svc.GetSessionList()
    if not sessions:
        print('No files open in FragalyseQt.')
        return

    for s in sessions:
        print(f'  [{s.session_id}] {s.file_name}  |  {s.instrument}')
        print(f'       dyes={s.dye_count}  sizing={s.has_sizing}'
              f'  panel={s.panel_name if s.has_panel else "none"}')

    first_id = sessions[0].session_id

    print(f'\n=== GetAnalysisParams (session {first_id}) ===')
    params = svc.GetAnalysisParams(session_id=first_id)
    print(f'  height={params.min_height}  width={params.min_width}'
          f'  prominence={params.min_prominence}')
    print(f'  sizing={params.sizing_method}  standard={params.size_standard}')

    print(f'\n=== GetPeakTable (session {first_id}) ===')
    peaks = svc.GetPeakTable(session_id=first_id)
    print(f'  {len(peaks)} peaks found')
    for p in peaks[:5]:
        allele = f'  allele={p.allele}' if p.allele else ''
        size = f'  size={p.size_bp:.2f} bp' if p.size_bp == p.size_bp else ''
        print(f'  ch{p.channel} {p.dye_name:8s}'
              f'  pos={p.position_dp:8.2f}'
              f'  h={p.height:7.1f}{size}{allele}')
    if len(peaks) > 5:
        print(f'  ... and {len(peaks) - 5} more')

    print(f'\n=== ExportCSV (session {first_id}) ===')
    csv_content = svc.ExportCSV(session_id=first_id)
    lines = csv_content.strip().splitlines()
    print(f'  {len(lines)} lines (header + {len(lines) - 1} peaks)')
    print(f'  Header: {lines[0]}')

    print(f'\n=== GetRawSignal (session {first_id}, channel 1) ===')
    samples = svc.GetRawSignal(session_id=first_id, channel=1)
    if samples:
        print(f'  {len(samples)} samples, '
              f'min={min(samples):.1f}  max={max(samples):.1f}')
    else:
        print('  No signal data (analysis not yet run?)')

    print('\nAll read operations completed successfully.')


if __name__ == '__main__':
    main()
