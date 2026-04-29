Name:           fragalyseqt
Version:        0.5.1
Release:        1
Summary:        DNA fragment analysis tool
License:        AGPL-3.0
URL:            https://github.com/Dorif/fragalyseqt
Source0:        %{name}-%{version}.tar.gz

BuildArch:      noarch
BuildRequires:  python3-devel >= 3.8
BuildRequires:  python3-setuptools
BuildRequires:  pyproject-rpm-macros

Requires:       python3 >= 3.8
Requires:       python3-pyqtgraph >= 0.11.0
Requires:       python3-biopython >= 1.58
Requires:       python3-scipy >= 1.5.0
Requires:       python3-charset-normalizer
Requires:       python3-platformdirs >= 2.0
Requires:       python3-qt5

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
* Wed Apr 29 2026 Alexandr Dorif <dorif11@gmail.com> - 0.5.1-1
- Initial RPM packaging
