# coding: utf-8
"""A Python front-end to generating an integrated human 
   protein-protein interaction network.
"""

from __future__ import absolute_import, print_function, division

__version__ = 'v2.4.1'
__release__ = True
__status__ = "Production"
__author__ = """Gaston K. Mazandu & Christopher Hooper\n\t(c) 2020 All rights reserved."""
__author_email__ = """{christopher.hooper, gaston.mazandu}@uct.ac.za,\n\t gmazandu@gmail.com, kuzamunu@aims.ac.za"""
__license__ = "GPL (https://www.gnu.org/licenses/gpl-3.0.en.html)"

__all__ = ["sequenceprocessing", "networkgenerator"]

from .networkgenerator import *
from .sequenceprocessing import *

# Dynamically get the version of the installed module
try:
	import pkg_resources
	__version__ = pkg_resources.get_distribution(__name__).version
except Exception:
	pkg_resources = None
finally:
	del pkg_resources


