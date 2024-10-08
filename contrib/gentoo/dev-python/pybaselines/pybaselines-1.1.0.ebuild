# Copyright 1999-2024 Gentoo Authors
# Distributed under the terms of the GNU General Public License v2

EAPI=8

DISTUTILS_USE_PEP517=setuptools
PYTHON_COMPAT=( python3_{10..13} )

inherit distutils-r1 pypi

if [[ ${PV} == *9999 ]]; then
	inherit git-r3
	EGIT_REPO_URI="https://github.com/derb12/${PN}.git"
else
	SRC_URI="$(pypi_sdist_url "${PN^}" "${PV}")"
	KEYWORDS="~amd64 ~x86"
fi

DESCRIPTION="A library of algorithms for the baseline correction of experimental data."
HOMEPAGE="https://${PN}.readthedocs.io https://github.com/derb12/${PN} https://pypi.org/project/${PN}"

LICENSE="BSD"
SLOT="0"
IUSE=""

RESTRICT="mirror"

RDEPEND="
	>=dev-python/scipy-1.5.0[${PYTHON_USEDEP}]
	>=dev-python/numpy-1.20.0[${PYTHON_USEDEP}]
"
distutils_enable_tests pytest

S=${WORKDIR}/${P}

# USE full - needs below - but seems to be only for faster processing, no new func
#>=dev-python/pentapy-1.1
#>=dev-python/numba-0.49


