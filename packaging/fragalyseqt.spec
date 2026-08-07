Name:           fragalyseqt
Version:        0.5.4
Release:        1%{?dist}
Summary:        DNA fragment analysis tool
License:        AGPL-3.0-or-later
URL:            https://github.com/Dorif/fragalyseqt
Source0:        %{name}-%{version}.tar.gz

%define _buildhost reproducible
%define _buildtime 1778112000

# Two builds are produced, deliberately: one on Fedora that installs on every
# RPM distribution (openSUSE included), and one on ALT Linux, whose package
# naming and dependency generators differ irreconcilably.
#
# Distribution-specific macros (%%_docdir, %%_licensedir, %%__python3, ...) are
# NOT used anywhere in this spec on purpose: they resolve to different paths on
# different distributions, which would make the build non-reproducible and the
# resulting package non-portable. Paths are hardcoded instead.
%if "%{_vendor}" == "alt"
%global is_alt 1
%else
%global is_alt 0
%endif

BuildArch:      noarch
Group:          Sciences/Biology
BuildRequires:  python3
%if 0%{is_alt}
# ALT's rpm-build hard-depends on rpm-macros-python (Python 2 era), which
# registers a /usr/lib/rpm/python.prov generator whose python.prov.py helper
# is not shipped. Any package containing .py files therefore fails at
# "Finding Provides". The generators are switched off here; every dependency
# of this package is declared explicitly below anyway. Scoped to ALT so the
# universal RPM keeps its automatic detection.
AutoReq:        no
AutoProv:       no
Requires:       python3 >= 3.8
Requires:       python3-module-pyqtgraph >= 0.11.0
# ALT splits numpy.testing into a separate package, and the application does
# not work correctly without it. It hard-depends on python3-module-numpy of
# the same version, so numpy itself is pulled in by this line.
Requires:       python3-module-numpy-testing
Requires:       python3-module-scipy >= 1.5.0
Requires:       python3-module-charset-normalizer
Requires:       python3-module-platformdirs >= 2.0
Requires:       python3-module-PyQt5
# ALT splits parts of the Python standard library into separate packages.
Requires:       python3-modules-sqlite3
# ALT has no Biopython package; it is installed from PyPI in %%post, which
# keeps this package noarch instead of one build per architecture.
Requires:       python3-module-pip
# Needed only if PyPI has no matching Biopython wheel and pip has to build it.
Requires:       gcc
Requires:       python3-dev
%else
# Portable dependencies: one Fedora-built package must resolve on Fedora, RHEL
# and openSUSE alike, so every name is an alternation of the spellings used by
# those distributions. On openSUSE the versioned python313-* packages provide
# the unversioned python3-* names, so the python3-* spelling resolves there.
Requires:       python3 >= 3.8
Requires:       (python3-pyqtgraph >= 0.11.0 or python-pyqtgraph >= 0.11.0)
Requires:       (python3-biopython >= 1.58 or python-biopython >= 1.58)
Requires:       (python3-numpy or python-numpy)
Requires:       (python3-scipy >= 1.5.0 or python-scipy >= 1.5.0)
Requires:       (python3-charset-normalizer or python-charset-normalizer)
Requires:       (python3-platformdirs >= 2.0 or python-platformdirs >= 2.0)
Requires:       (python3-qt5 or python3-pyqt5 or python3-PyQt5 or python3-qt6 or python3-pyside6)
Recommends:     (python3-zeep or python-zeep)
%endif

%description
CE-machine-independent tool for fragment analysis data processing.
Supports baseline correction and denoising, selective channel hiding,
peak detection and sizing, allele binning using GeneMapper, GeneMarker
and NCBI OSIRIS panels, stutter filtering, CSV export, and CODIS CMF
export and import.

%prep
%setup -q

%build

%install
install -dm755 %{buildroot}/usr/lib/fragalyseqt
cp -a src/fragalyseqt/. %{buildroot}/usr/lib/fragalyseqt/
find %{buildroot}/usr/lib/fragalyseqt -type d -name '__pycache__' -prune -exec rm -rf {} +
find %{buildroot}/usr/lib/fragalyseqt -type f -name '*.pyc' -delete
python3 -m compileall -d /usr/lib/fragalyseqt -q %{buildroot}/usr/lib/fragalyseqt
python3 -O -m compileall -d /usr/lib/fragalyseqt -q %{buildroot}/usr/lib/fragalyseqt

install -dm755 %{buildroot}/usr/bin
cat > %{buildroot}/usr/bin/fragalyseqt << EOF
#!/usr/bin/python3
import sys
sys.path.insert(0, '/usr/lib/fragalyseqt')
sys.path.insert(0, '/usr/lib')
from fragalyseqt.main import main
main()
EOF
chmod 755 %{buildroot}/usr/bin/fragalyseqt

install -Dm644 packaging/fragalyseqt.desktop \
    %{buildroot}/usr/share/applications/fragalyseqt.desktop
install -Dm644 packaging/fragalyseqt.png \
    %{buildroot}/usr/share/icons/hicolor/512x512/apps/fragalyseqt.png
install -Dm644 README.md \
    %{buildroot}/usr/share/doc/fragalyseqt/README.md
install -Dm644 COPYING \
    %{buildroot}/usr/share/licenses/fragalyseqt/COPYING

%post
%if 0%{is_alt}
# ALT Linux does not package Biopython. Installing it here keeps this
# package architecture-independent instead of shipping one build per
# architecture with Biopython's compiled extensions bundled in.
if ! python3 -c 'import Bio' >/dev/null 2>&1; then
    python3 -m pip install --quiet biopython 2>/dev/null || true
    if ! python3 -c 'import Bio' >/dev/null 2>&1; then
        echo "WARNING: fragalyseqt could not install Biopython from PyPI." >&2
        echo "         The application will not start until it is available." >&2
        echo "         Install it manually: python3 -m pip install biopython" >&2
    fi
fi
%endif

%files
/usr/lib/fragalyseqt/
/usr/bin/fragalyseqt
/usr/share/doc/fragalyseqt/README.md
/usr/share/licenses/fragalyseqt/COPYING
/usr/share/applications/fragalyseqt.desktop
/usr/share/icons/hicolor/512x512/apps/fragalyseqt.png

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
- Containerised build matrix for Ubuntu, Fedora, openSUSE and ALT Linux.
- Distribution-specific Python dependency names.
* Tue May 19 2026 Alexandr Dorif <dorif11@gmail.com> - 0.5.3-1
- Smart ILS ladder alignment via relative spacing pattern matching.
- Performance optimisations and code quality improvements.
- Bug fixes: session restore, read-only tab rendering, stutter filter.
- UI improvements for HID users.
* Fri May 08 2026 Alexandr Dorif <dorif11@gmail.com> - 0.5.2-1
- Batch processing, session CSV export, SOAP API, database design.
* Wed May 06 2026 Alexandr Dorif <dorif11@gmail.com> - 0.5.1-2
- Ensuring compatibility with SUSE and AltLinux based distros, moving towards reproducible builds.
* Wed Apr 29 2026 Alexandr Dorif <dorif11@gmail.com> - 0.5.1-1
- Initial RPM packaging
