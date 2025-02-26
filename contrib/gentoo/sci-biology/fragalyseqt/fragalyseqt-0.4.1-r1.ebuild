# Copyright 1999-2025 Gentoo Authors
# Distributed under the terms of the GNU General Public License v2

EAPI=8

DISTUTILS_USE_PEP517=setuptools
PYTHON_COMPAT=( python3_{10..13} )

inherit distutils-r1 desktop

DESCRIPTION="Software for DNA fragment analysis (MLPA, QF-PCR etc.) data processing."
HOMEPAGE="https://github.com/Dorif/fragalyseqt"

MY_RELEASE="jeffreys_bugfix"

if [[ ${PV} == 9999 ]]; then
	EGIT_REPO_URI="https://github.com/Dorif/${PN}.git"
	inherit git-r3
else
	SRC_URI="https://codeload.github.com/Dorif/${PN}/zip/refs/tags/${MY_RELEASE} -> ${P}.gh.zip"
	KEYWORDS="~amd64 ~x86"
fi

LICENSE="AGPL-3"
SLOT="0"

PATCHES="${FILESDIR}/${PN}-fix_imports.patch"

S="${WORKDIR}/${PN}-${MY_RELEASE}"

DOCS="${S}/README.md ${S}/ABIF_specs/ABIF_File_Format-2006.pdf ${S}/ABIF_specs/ABIF_File_Format-2009.pdf"

DEPEND="
	dev-python/pyside:6[${PYTHON_USEDEP}]
	>=dev-python/pyqtgraph-0.11.0[qt6,${PYTHON_USEDEP}]
	>=sci-biology/biopython-1.58[${PYTHON_USEDEP}]
	>=dev-python/scipy-1.5.0[${PYTHON_USEDEP}]
	dev-python/charset-normalizer[${PYTHON_USEDEP}]
	>=dev-python/pybaselines-1.1.0[${PYTHON_USEDEP}]
"

RDEPEND="${DEPEND}"

MY_SRC="src/FragalyseApp"
MY_APP="FragalyseQt"

src_prepare() {
	default

	# Put files to PEP517-compatible folder structure
	mkdir -p "${MY_SRC}"
	mv *.py "${MY_SRC}"
	touch "${MY_SRC}"/__init__.py

	# Copy PEP517 setup data
	cp "${FILESDIR}"/pyproject.toml .
}

python_install_all() {
	distutils-r1_python_install_all

	# Install examples
	insinto /usr/share/${PN}
	doins -r TEST_FILES
	
	# Create desktop file
	doicon "${MY_APP}".png
	make_desktop_entry ${PN} "${MY_APP}" ${MY_APP}
}
