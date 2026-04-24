<<<<<<< HEAD
# {# pkglts, pysetup.kwds
# format setup arguments

from setuptools import setup, find_packages


short_descr = "Simbal is a FSPM model for agroecological simulation"
readme = open('README.rst').read()
history = open('HISTORY.rst').read()

# version 
version = {}
with open("src/simbal/version.py") as fp:
    exec(fp.read(), version)
_version = version["__version__"]

# find packages
pkgs = find_packages('src')

authors = "Laurene perthame, Solen Farra, ..."
author_email="laurene.perthame at agroparistech.fr, solen.farra at inrae.fr, ..."

setup_kwds = dict(
    name='simbal',
    version=_version,
    description=short_descr,
    long_description=readme + '\n\n' + history,
    author=authors,
    author_email="alice.penanhoat@live.fr",
    url='https://github.com/openalea-incubator/caracol',
    license='cecill-c',
    zip_safe=False,

    packages=pkgs,
    
    package_dir={'': 'src'},
    install_requires=[],
    entry_points={},
    keywords='agroecology, FSPM, openalea',
    )
# #}
# change setup_kwds below before the next pkglts tag

# do not change things below
# {# pkglts, pysetup.call
setup(**setup_kwds)
# #}
=======
# -*- coding: latin-1 -*-
import sys
from setuptools import setup

"""

    setup
    ~~~~~

    Setup script for installation.

    See README.md for installing procedure.

    :copyright: Copyright 2023-2024 INRA-ECOSYS, see AUTHORS.
    :license: CeCILL-C, see LICENSE for details.

    **Acknowledgments**: The research leading these results has received funding through the 
    Investment for the Future programme managed by the Research National Agency 
    (BreedWheat project ANR-10-BTBR-03).

    .. seealso:: 1st article et al.
"""

"""
    Information about this versioned file:
        $LastChangedBy$
        $LastChangedDate$
        $LastChangedRevision$
        $URL$
        $Id$
"""

if sys.version_info < (3, 7):
    print('ERROR: Requires at least Python 3.7 to run.')
    sys.exit(1)

setup(
    name="simbal",
    version="0.0.1",
    packages=["simbal"],
    include_package_data=True,
    author="F.Rees, C.Pradal",
    author_email="frederic.rees@inrae.fr, christophe.pradal@cirad.fr",
    description="",
    long_description="""TODO""",
    license="CeCILL-C",
    keywords="",
    url="https://github.com/openalea-incubator/simbal.git",
    download_url="https://github.com/openalea-incubator/simbal.git"
)
>>>>>>> frederic
