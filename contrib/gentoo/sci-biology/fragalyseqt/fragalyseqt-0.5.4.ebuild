# Copyright 1999-2026 Gentoo Authors
# Distributed under the terms of the GNU General Public License v2

EAPI=8

DISTUTILS_USE_PEP517=setuptools
PYTHON_COMPAT=( python3_{10..14} )

inherit distutils-r1 desktop

DESCRIPTION="Software for DNA fragment analysis (MLPA, QF-PCR etc.) data processing."
HOMEPAGE="https://github.com/Dorif/fragalyseqt"

if [[ ${PV} == 9999 ]]; then
	EGIT_REPO_URI="https://github.com/Dorif/${PN}.git"
	inherit git-r3
else
	SRC_URI="https://github.com/Dorif/${PN}/archive/refs/tags/v${PV}.tar.gz -> ${P}.tar.gz"
	KEYWORDS="~amd64 ~x86 ~arm64 ~riscv"
fi

LICENSE="AGPL-3+"
SLOT="0"

DOCS=( README.md docs/SPECS_AND_REFERENCES/ABIF_File_Format-2006.pdf docs/SPECS_AND_REFERENCES/ABIF_File_Format-2009.pdf )

DEPEND="
	|| ( dev-python/pyside:6[${PYTHON_USEDEP}]
	     dev-python/pyqt6[${PYTHON_USEDEP}]
	     dev-python/PyQt5[${PYTHON_USEDEP}] )
	>=dev-python/pyqtgraph-0.11.0[${PYTHON_USEDEP}]
	>=sci-biology/biopython-1.58[${PYTHON_USEDEP}]
	>=dev-python/scipy-1.5.0[${PYTHON_USEDEP}]
	dev-python/charset-normalizer[${PYTHON_USEDEP}]
	>=dev-python/platformdirs-2.0[${PYTHON_USEDEP}]
	dev-lang/python:*[sqlite]
"

RDEPEND="${DEPEND}"

# Optional: dev-python/zeep for SOAP API test client (contrib/soap_test_client.py)

MY_APP="FragalyseQt"

src_prepare() {
	default
}

python_install_all() {
	distutils-r1_python_install_all

	# Install examples and panels
	insinto /usr/share/${PN}
	doins -r docs/TEST_FILES
	doins -r docs/OSIRIS_PANELS

	# Create desktop file
	newicon -s 512 "${MY_APP}".png ${PN}.png
	make_desktop_entry ${PN} "${MY_APP}" ${MY_APP}
}
