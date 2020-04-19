# coding: utf-8
"""a Python front-end to semantic similarity score retrieval and
   fundamental applications.
"""

from __future__ import absolute_import, print_function, division

__version__ = '2.5.1'
__release__ = True
__author__ = """Gaston K. Mazandu\n\t(c) 2020 All rights reserved."""
__author_email__ = """{gaston.mazandu@uct, kuzamunu@aims}.ac.za, gmazandu@gmail.com"""
__license__ = "GPL (https://www.gnu.org/licenses/gpl-3.0.en.html)"

__all__ = ["InformationContent", "ConceptSimilarity", "EntitySimilarity"]

from .informationcontent import InformationContent
from .conceptsimilarity import ConceptSimilarity
from .entitysimilarity import EntitySimilarity
from .smlapps.conceptenrichment import ConceptEnrichment
from .smlapps.entityidentification import EntityIdentification

# Dynamically get the version of the installed module
try:
	import pkg_resources
	__version__ = pkg_resources.get_distribution(__name__).version
except Exception:
	pkg_resources = None
finally:
	del pkg_resources

#meta['auto-generated-by'] = ['PySML v{}'.format(__version__)]
#meta['date'] = [datetime.datetime.now().strftime('%d:%m:%Y %H:%M')]
#https://github.com/gkm-software-dev/post-analysis-tools

