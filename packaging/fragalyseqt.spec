Name:           fragalyseqt
Version:        0.5.4
Release:        1%{?dist}
Summary:        DNA fragment analysis tool
License:        AGPL-3.0
URL:            https://github.com/Dorif/fragalyseqt
Source0:        %{name}-%{version}.tar.gz

%define _buildhost reproducible
%define _buildtime 1778112000
%if "%{_vendor}" == "alt"
%global dist .altlinux
%else
%global dist %{nil}
%endif

BuildArch:      noarch
Group:          Sciences/Biology

Requires:       python3 >= 3.8
%if "%{_vendor}" == "alt"
%global __find_provides /bin/true
%global __find_requires /bin/true
Requires:       python3-module-pyqtgraph >= 0.11.0
Requires:       python3-module-scipy >= 1.5.0
Requires:       python3-module-charset-normalizer
Requires:       python3-module-platformdirs >= 2.0
Requires:       python3-module-numpy-testing
Requires:       python3-module-PyQt5
Requires:       gcc
Requires:       python3-module-pip
Requires:       python3-dev
%else
Requires:       (python3-pyqtgraph >= 0.11.0 or python-pyqtgraph >= 0.11.0 or python3-module-pyqtgraph >= 0.11.0)
Requires:       (python3-biopython >= 1.58 or python-biopython >= 1.58 or python3-module-biopython >= 1.58)
Requires:       (python3-scipy >= 1.5.0 or python-scipy >= 1.5.0 or python3-module-scipy >= 1.5.0)
Requires:       (python3-charset-normalizer or python-charset-normalizer or python3-module-charset-normalizer)
Requires:       (python3-platformdirs >= 2.0 or python-platformdirs >= 2.0 or python3-module-platformdirs >= 2.0)
Requires:       (python3-qt5 or python3-qt6 or python3-pyside6 or python3-module-PyQt5 or python3-module-PyQt6 or python3-module-pyside6)
Recommends:     (python3-zeep or python-zeep)
%endif

%description
CE-machine-independent tool for fragment analysis data processing.
Supports baseline correction and denoising, selective channel
hiding, peak detection and sizing, allele binning using GeneMapper,
GeneMarker and NCBI OSIRIS panels, stutter allelesfiltering, and
CSV or CODIS 3.2 CMF XML export.

%prep
%setup -q

%build

%install
install -dm755 %{buildroot}/usr/lib/fragalyseqt
cp -r src/fragalyseqt/* %{buildroot}/usr/lib/fragalyseqt/
find %{buildroot}/usr/lib/fragalyseqt -name "__pycache__" -type d -exec rm -rf {} + 2>/dev/null || true
python3 -m compileall -d /usr/lib/fragalyseqt -q %{buildroot}/usr/lib/fragalyseqt
python3 -O -m compileall -d /usr/lib/fragalyseqt -q %{buildroot}/usr/lib/fragalyseqt

cat > %{buildroot}/usr/lib/fragalyseqt/__main__.py << 'EOF'
from fragalyseqt.main import main
main()
EOF

install -dm755 %{buildroot}/usr/bin
cat > %{buildroot}/usr/bin/fragalyseqt << 'EOF'
#!/usr/bin/python3
import sys
sys.path.insert(0, '/usr/lib')
from fragalyseqt.main import main
main()
EOF
chmod 755 %{buildroot}/usr/bin/fragalyseqt

install -Dm644 packaging/fragalyseqt.desktop \
    %{buildroot}%{_datadir}/applications/fragalyseqt.desktop
install -Dm644 packaging/fragalyseqt.png \
    %{buildroot}%{_datadir}/icons/hicolor/512x512/apps/fragalyseqt.png
install -Dm644 README.md \
    %{buildroot}/usr/share/doc/fragalyseqt/README.md
install -Dm644 COPYING \
    %{buildroot}/usr/share/licenses/fragalyseqt/COPYING

%post
%if "%{_vendor}" == "alt"
python3 -m pip install --quiet biopython 2>/dev/null || true
%endif

%files
/usr/lib/fragalyseqt/
/usr/bin/fragalyseqt
/usr/share/doc/fragalyseqt/README.md
/usr/share/licenses/fragalyseqt/COPYING
%{_datadir}/applications/fragalyseqt.desktop
%{_datadir}/icons/hicolor/512x512/apps/fragalyseqt.png

%changelog
* Sat Jun 06 2026 Alexandr Dorif <dorif11@gmail.com> - 0.5.4-1

- Forensic STR statistics: identity LR, kinship index (Fung & Hu 2008).
- Allele frequency tables: NIST 1036 and STRidER 2025 bundled out of the box.
- Reference profile database with append-only storage.
- Profile comparison dialog with CODIS XML import and saved profiles.
- Search profiles in database, PDF and CSV export of results.
- Familias .fam frequency table parser.
- SOAP API extended with 11 new operations; WSDL updated.
- Frequency table and reference profile manager dialogs.
- CODIS CMF export and import extended to standard 3.3 and the Rapid
  Import CMF format (3.2 retained); CODIS XML import auto-detects the
  format.
- Smarter ILS ladder alignment: relative-spacing/ratio-of-ratios
  consensus matcher fixes mis-alignment of sub-ladders.
- Log-parabolic sub-datapoint peak refinement; recovery of true
  height/area for saturated (clipped) peaks from their flanks, flagged
  in the results table, CSV export and database.- Smart ILS ladder
  alignment via relative spacing pattern matching.
* Fri May 08 2026 Alexandr Dorif <dorif11@gmail.com> - 0.5.2-1
- Batch processing, session CSV export, SOAP API, database design.
* Wed May 06 2026 Alexandr Dorif <dorif11@gmail.com> - 0.5.1-2
- Ensuring compatibility with SUSE and AltLinux based distros, moving towards reproducible builds.
* Wed Apr 29 2026 Alexandr Dorif <dorif11@gmail.com> - 0.5.1-1
- Initial RPM packaging
