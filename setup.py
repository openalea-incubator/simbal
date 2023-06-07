# {# pkglts, pysetup.kwds
# format setup arguments

from setuptools import setup, find_packages


short_descr = "Caracol is a FSPM model for agroecological simulation"
readme = open('README.md').read()
history = open('HISTORY.rst').read()

# find packages
pkgs = find_packages('src')

setup_kwds = dict(
    name='caracol',
    version="0.0.2",
    description=short_descr,
    long_description=readme + '\n\n' + history,
    author="Solen, Laurene et al.",
    author_email="",
    url='https://github.com/openalea-incubator/caracol',
    license='cecill-c',
    zip_safe=False,

    packages=pkgs,
    
    package_dir={'': 'src'},
    install_requires=[
        ],
    tests_require=[
        "pytest",
        ],
    entry_points={},
    keywords='',
    )
# #}
# change setup_kwds below before the next pkglts tag

# do not change things below
# {# pkglts, pysetup.call
setup(**setup_kwds)
# #}
