# Copyright 1999-2024 Gentoo Authors
# Distributed under the terms of the GNU General Public License v2

EAPI=8

DISTUTILS_USE_PEP517=setuptools
PYTHON_COMPAT=( python3_{10..13} )

inherit distutils-r1 desktop

if [[ ${PV} == *9999 ]]; then
	inherit git-r3
	EGIT_REPO_URI="https://github.com/Dorif/${PN}.git"
else
	SRC_URI="https://codeload.github.com/Dorif/fragalyseqt/zip/refs/tags/jeffreys_bugfix -> ${P}.gh.zip"
	KEYWORDS="~amd64 ~x86"
fi

DESCRIPTION="Tool for DNA fragment analysis data handling."
HOMEPAGE="https://github.com/Dorif/fragalyseqt"

RESTRICT="mirror"

LICENSE="AGPL-3"
SLOT="0"
IUSE=""

PATCHES="${FILESDIR}/import.patch"

S="${WORKDIR}/${PN}-jeffreys_initial"

DOCS="${S}/README.md ${S}/ABIF_specs/ABIF_File_Format-2006.pdf ${S}/ABIF_specs/ABIF_File_Format-2009.pdf"

DEPEND="
	dev-python/pyside6[${PYTHON_USEDEP}]
	>=dev-python/pyqtgraph-0.11.0[${PYTHON_USEDEP}]
	>=sci-biology/biopython-1.58[${PYTHON_USEDEP}]
	>=dev-python/scipy-1.5.0[${PYTHON_USEDEP}]
	dev-python/charset-normalizer[${PYTHON_USEDEP}]
	>=dev-python/pybaselines-1.1.0[${PYTHON_USEDEP}]
"

RDEPEND="${DEPEND}"

src_prepare() {
	default

	# Put files to PEP517-compatible folder structure
	mkdir -p src/FragalyseApp
	mv *.py src/FragalyseApp
	touch src/FragalyseApp/__init__.py

	# Copy PEP517 setup data
	cp "${FILESDIR}"/pyproject.toml .
}

python_install_all() {
	distutils-r1_python_install_all

	# Install examples
	insinto /usr/share/${PN}
	doins -r TEST_FILES
	
	# Create desktop file
	doicon FragalyseQt.png
	make_desktop_entry ${PN} "FragalyseQt" FragalyseQt
}
