#!/usr/bin/env python

from distutils.core import setup
from platform import python_version_tuple
import re

LICENSE = open("LICENSE").read()

VERSION = open("VERSION").read()

# strip links from the descripton on the PyPI
LONG_DESCRIPTION = open("README.rst").read().replace("`_", "`")

# strip Build Status from the PyPI package
if python_version_tuple()[:2] >= ('2', '7'):
    LONG_DESCRIPTION = re.sub("^Build status\n(.*\n){7}", "", LONG_DESCRIPTION, flags=re.M)

setup(name='The IHP-PING Package',
   version=VERSION,
   description=69*'*'+'\n* An integrated human protein-protein interaction network generator *\n'+69*'*',
   long_description=LONG_DESCRIPTION,
   author='Gaston K. Mazandu et al.',
   author_email='HPRCHR001@myuct.ac.za, babuken@gmail.com, funmite@aims.ac.za, vnembaware@gmail.com, nicholas.thomford@uct.ac.za, emile.chimusa@uct.ac.za, nicola.mulder@uct.ac.za, ambroise.wonkam@uct.ac.za, gaston.mazandu@uct.ac.za',
   maintainer = 'Gaston K. Mazandu',
   maintainer_email = 'gaston.mazandu@uct.ac.za, gmazandu@gmail.com, kuzamunu@aims.ac.za',
   url='http://web.cbio.uct.ac.za/ITGOM/post-analysis-tools/ihp-ping-dev/\nhttps://github.com/gkm-software-dev/post-analysis-tools',
   license=LICENSE,
   classifiers= [ "Development Status :: 4 - Beta",
                  "License :: OSI Approved :: GNU General Public License",
                  "Operating System :: OS Independent, but tested only on Linux",
                  "Programming Language :: Python :: Not tested on the version less than 2.7",
                  "Programming Language :: Python :: 2.7 or greater",
                  "Topic :: Software Development :: Libraries",
                  "Following libraries need to be installed prior to the installation and the\nuse of A-DaFO-Fun:",
                   "\t::ncbi-blast+: For protein sequence data\n\t::python-selenium: For DIP datasets\n\t::chromium-chromedriver: For DIP datasets\n\t::chromium-browser: For DIP datasets"],
   platforms = 'x86_64-pc-linux-gnu',
   requires = ['\nncbi_blast', 'python_selenium', 'chromium_chromedrive', 'chromium_browser\n'],
   provides = ['ihppinbuilder.py', 'PyPING.networkbuilder.py', 'PyPING.sequenceprocessing.py'],
   )
