Name:           fragalyseqt
Version:        0.5.1
Release:        2
Summary:        DNA fragment analysis tool
License:        AGPL-3.0
URL:            https://github.com/Dorif/fragalyseqt
Source0:        %{name}-%{version}.tar.gz

%define _buildhost reproducible
%define _buildtime 1778112000

BuildArch:      noarch

Requires:       python3 >= 3.8
Requires:       (python3-pyqtgraph >= 0.11.0 or python-pyqtgraph >= 0.11.0)
Requires:       (python3-biopython >= 1.58 or python-biopython >= 1.58)
Requires:       (python3-scipy >= 1.5.0 or python-scipy >= 1.5.0)
Requires:       (python3-charset-normalizer or python-charset-normalizer)
Requires:       (python3-platformdirs >= 2.0 or python-platformdirs >= 2.0)
Requires:       (python3-qt5 or python3-qt6 or python3-pyside6)

%description
CE-machine-independent tool for fragment analysis data processing.
Supports baseline correction and denoising, selective channel
hiding, peak detection and sizing, allele binning using GeneMapper,
GeneMarker and NCBI OSIRIS panels, stutter allelesfiltering, and
CSV or CODIS 3.2 CMF XML export.

%prep
%autosetup

%build

%install
install -dm755 %{buildroot}/usr/lib/fragalyseqt
cp -r src/fragalyseqt/* %{buildroot}/usr/lib/fragalyseqt/
find %{buildroot}/usr/lib/fragalyseqt -name "__pycache__" -type d -exec rm -rf {} + 2>/dev/null || true

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
    %{buildroot}%{_datadir}/pixmaps/fragalyseqt.png
install -Dm644 README.md \
    %{buildroot}/usr/share/doc/fragalyseqt/README.md
install -Dm644 COPYING \
    %{buildroot}/usr/share/licenses/fragalyseqt/COPYING

%files
/usr/lib/fragalyseqt/
/usr/bin/fragalyseqt
/usr/share/doc/fragalyseqt/README.md
/usr/share/licenses/fragalyseqt/COPYING
%{_datadir}/applications/fragalyseqt.desktop
%{_datadir}/pixmaps/fragalyseqt.png

%changelog
* Wed May 06 2026 Alexandr Dorif <dorif11@gmail.com> - 0.5.1-2
- Ensuring compatibility with SUSE based distros, moving towards reproducible builds.
* Wed Apr 29 2026 Alexandr Dorif <dorif11@gmail.com> - 0.5.1-1
- Initial RPM packaging
