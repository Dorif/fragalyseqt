# Copyright 1999-2024 Gentoo Authors
# Distributed under the terms of the GNU General Public License v2

EAPI=8

DISTUTILS_USE_PEP517=setuptools
PYTHON_COMPAT=( python3_{10..13} )

inherit distutils-r1 #pypi

DESCRIPTION="A library of algorithms for the baseline correction of experimental data."
HOMEPAGE="https://${PN}.readthedocs.io https://github.com/derb12/${PN} https://pypi.org/project/${PN}"

if [[ ${PV} == 9999 ]]; then
	EGIT_REPO_URI="https://github.com/derb12/${PN}.git"
	inherit git-r3
else
	#SRC_URI="$(pypi_sdist_url "${PN^}" "${PV}")"

	SRC_URI="https://github.com/derb12/pybaselines/archive/refs/tags/v${PV}.tar.gz -> ${P}.gh.tar.gz"
	KEYWORDS="~amd64 ~x86"
fi


LICENSE="BSD"
SLOT="0"

RDEPEND="
	>=dev-python/scipy-1.5.0[${PYTHON_USEDEP}]
	>=dev-python/numpy-1.20.0[${PYTHON_USEDEP}]
"
distutils_enable_tests pytest

S=${WORKDIR}/${P}
