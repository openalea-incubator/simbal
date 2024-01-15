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
