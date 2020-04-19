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

setup(name='PySML Interface',
   version=VERSION,
   description='\n'+68*'*'+'\n* An open Python library implementing Semantic Similarity Measures *'+'\n* and common related applications'.ljust(68)+'*\n'+68*'*'+'\n',
   long_description=LONG_DESCRIPTION,
   author='Gaston K. Mazandu et al.',
   author_email='gmazandu@cbio.uct.ac.za, emile@cbio.uct.ac.za, mamana@aims.ac.za, nicola.mulder@uct.ac.za',
   maintainer = 'Gaston K. Mazandu',
   maintainer_email = 'gaston.mazandu@uct.ac.za, gmazandu@gmail.com, kuzamunu@aims.ac.za',
   url='http://web.cbio.uct.ac.za/ITGOM/post-analysis-tools/mysml-dev/\nhttps://github.com/gkm-software-dev/post-analysis-tools',
   license=LICENSE,
   classifiers= [ "Development Status :: 4 - Beta",
                  "License :: OSI Approved :: GNU General Public License",
                  "Operating System :: OS Independent, but tested only on Linux",
                  "Programming Language :: Python :: Not tested on the version less than 2.7",
                  "Programming Language :: Python :: 2.7 or greater",
                  "Topic :: Software Development :: Libraries",
                  "Following libraries need to be installed prior to the installation and the\nuse of PYSML:",
                   "\t::networkx\n\t::scipy\n\t::matplotlib"],
   platforms = 'x86_64-pc-linux-gnu',
   requires = ['\nscipy', 'scipy', 'matplotlib\n'],
   provides = ['procsemsim.py', 'PySML.informationcontent.py', 'PySML.conceptsimilarity.py', 'PySML.entitysimilarity.py', 'PySML.smlapps.conceptenrichment.py', 'PySML.smlapps.entityidentification.py', 'PySML.smlapps.entityclassification.py'],
   package_data={'PySML': ['tests/go-basic.obo'],})
