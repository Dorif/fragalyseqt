Name:           fragalyseqt
Version:        0.5.1
Release:        2
Summary:        DNA fragment analysis tool
License:        AGPL-3.0
URL:            https://github.com/Dorif/fragalyseqt
Source0:        %{name}-%{version}.tar.gz

%{?suse_version:%define pythons python313}
%global _bindir /usr/bin

BuildArch:      noarch
BuildRequires:  python3-devel >= 3.8
BuildRequires:  (python3-setuptools or python-setuptools)
BuildRequires:  (pyproject-rpm-macros or python-rpm-macros)

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
%pyproject_wheel

%install
%pyproject_install
install -Dm644 packaging/fragalyseqt.desktop \
    %{buildroot}%{_datadir}/applications/fragalyseqt.desktop
install -Dm644 packaging/fragalyseqt.png \
    %{buildroot}%{_datadir}/pixmaps/fragalyseqt.png

%files
%license COPYING
%doc README.md
%{python3_sitelib}/fragalyseqt/
%{python3_sitelib}/fragalyseqt-*.dist-info/
%{_bindir}/fragalyseqt
%{_datadir}/applications/fragalyseqt.desktop
%{_datadir}/pixmaps/fragalyseqt.png

%changelog
* Wed May 06 2026 Alexandr Dorif <dorif11@gmail.com> - 0.5.1-2
- Ensuring compatibility with SUSE based distros.
* Wed Apr 29 2026 Alexandr Dorif <dorif11@gmail.com> - 0.5.1-1
- Initial RPM packaging
