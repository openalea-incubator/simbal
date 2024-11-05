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