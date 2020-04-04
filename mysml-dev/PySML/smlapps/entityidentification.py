# coding: utf-8

from __future__ import print_function, division
from math import exp, log, tanh, atan, pi as PI

import sys
import networkx as nx

from .. import __version__
from ..imports.readontology import Ontology, output_str
from ..imports.tabulate import tabulate as tabs
from ..conceptsimilarity import ConceptSimilarity


class EntityIdentification(ConceptSimilarity):
	__cslot__ = []
	def __init__(self, ontofile = '', namespace='', is_a = None, part_of = None):
		ConceptSimilarity.__init__(self, ontofile, namespace, is_a, part_of)
		self.fouts = {}
		self.comments = {}
	
	def search(self, annot, ctargets, model=('nunivers','universal'), **kwargs):
		"""annot: is the annotaton file-dict(entity-annotation
		   ctargets: list, tuple or list of concepts for which
		   entities from annot should map to based on the model.
		   model is the concept similarity model used.
		"""
		self.fouts = dict([(t, {}) for t in ctargets])
		for t in ctargets:
			for ent in annot:
				Pairs = [(t, s) for s in annot[ent]]
				self.fouts[t][ent]=self.conceptInterface(Pairs, model, **kwargs)

	def retrieveEntity(self, annot, ctargets, model= None, **kwargs):
		kwargs = kwargs
		cutoff = kwargs['cutoff'] if cutoff in kwargs else 0.3



