# coding: utf-8

from __future__ import print_function, division
from math import exp, log, tanh, atan, pi as PI
from functools import reduce

import sys
import networkx as nx

from .. import __version__
from ..imports.readontology import Ontology, output_str
from ..imports.tabulate import tabulate as tabs
from ..conceptsimilarity import ConceptSimilarity


class ConceptEnrichment(ConceptSimilarity):
	__cslot__ = []
	def __init__(self, ontofile = '', namespace='', is_a = None, part_of = None):
		ConceptSimilarity.__init__(self, ontofile, namespace, is_a, part_of)
		self.fouts = {}
		self.comments = {}
		self.Targets = set()
		self.TargetMissing = set()
		self.Background = {}
	
	def search(self, model=('nunivers','universal'), **kwargs):
		"""annot: is the annotaton file-dict(entity-annotation
		   ctargets: list, tuple or list of concepts for which
		   entities from annot should map to based on the model.
		   model is the concept similarity model used.
		   
		Inputs: takes 5 parameters: the set of terms in the target 'targetterms', set of protein target 'targetproteins', set of proteins in the background proteins 'background', Approach under consideration 'app' and Agreement level 'agree'.
		 Outputs: returns fuzzy frequency of occurrences of each term in the target and reference gene sets given term semantic similarity approach under consideration.
		"""
		cutoff = kwargs['cutoff'] if 'cutoff' in kwargs else 0.3
		self.fouts = dict()
		Targets = reduce(lambda x, y: x | y, [self.Background[ent] for ent in self.Targets]); 
		for t in Targets:
			m = n = 0
			for ent in self.Background:
				Pairs = [(t, s) for s in self.Background[ent]]
				ds = self.conceptInterface(Pairs, model, **kwargs)
				if max(ds.values())==1: # The term occurs in ent
					n += 1
					if ent in self.Targets: m += 1
				else:
					for k in ds:
						if t in nx.ancestors(self.DagStr,k[1]): # t ancestor of s
							if ds[(t,k[1])] > cutoff: # Concept fuzzy-occurs in ent
								n += 1
								if ent in self.Targets: m += 1
								break
						elif k[1] in nx.ancestors(self.DagStr, t): # s descent from t
							if ds[(t,k[1])] > 0.7: # Concept fuzzy-occurs in ent
								n += 1
								if ent in self.Targets: m += 1
								break
			self.fouts[t] = (n,m) 

	def enrichedConcept(self, Annots, ctargets, model= None, **kwargs):
		"""This function computes concept similarity scores given:
           - List of concepts or concept pairs : TermPairs
           - List of tuples similarity model-IC approach: models
           or list of tuples similarity model-None if sim model does
           not require IC approach, e.g., edge-based models
           - Other measure parameter if requred
		"""
		if ctargets is None or not ctargets:
			print('\n\tYou should provide concept pairs for \n\twhich similarity scores should be\n\tcomputed are required\n')
			raise ValueError("\n\tIssue with the concept target set provided.\n\tPlease, provide the list or tuple or set of target concepts.")

		if model is None or not model: model = ('nunivers','universal')
		elif isinstance(model, tuple):
			if len(model)==2 and isinstance(model[0],str) and isinstance(model[1],str):
				if model[0].lower() in self.Models and model[1].lower() in self.AppMods:
					model = (model[0].lower(), model[1].lower())
				else: model = None
			elif len(model)==1 and isinstance(model[0],str):
				if model[0].lower() in self.Models: model = (model[0].lower(),)
				else: model = None
		elif isinstance(model, str) and model.lower() in self.Models:
			model = (model.lower(),)
		else:# Type Error! To still look at!
			raise TypeError('Please, include only implemented Information Content and Concept Similarity models as arguments. Check IC and/or CS provided and try again.')

		if model is None and not model: # proposed approach, not known 
			raise ValueError('Please, include only implemented Information Content and Concept Similarity models as arguments. Check IC and CS provided and try again.')

		self.Background = {}; self.EntityMissing = {}
		if isinstance(Annots, dict):
			for ent in Annots:
				self.Background[ent] = set(); self.EntityMissing[ent] = set()
				for tt in Annots[ent]:
					if tt in self.alt_id and self.alt_id[tt] in self.DagStr:
						self.Background[ent].add(self.alt_id[tt])
					elif tt in self.Dag and self.Dag.index(tt) in self.DagStr:
						self.Background[ent].add(self.Dag.index(tt))
					else: # Term is either obsolete or does not exist in the onto!
						self.EntityMissing[ent].add(tt)
				self.Background[ent].discard(self.oroot)
				if not self.Background[ent]: self.Background.pop(ent, False)
		else:
			raise TypeError('Please, the background set should be a dictionary with key-entities and value-ontology concept. Check the background set provided is a dictionary and try again.')
		if not self.Background:
			raise ValueError('Sorry the background map is empty. Please check the type of the set of concept targets.')

		self.Targets = set(); self.TargetMissing = set()
		if isinstance(ctargets, (list, set, tuple)):
			for p in ctargets:
				if p in self.Background: self.Targets.add(p)
				else: self.TargetMissing.add(p)
		else:
			raise TypeError('Please, the target set should be a set, tuple or list of concepts in the ontology. Check the type of the set of concept targets.')
		if not self.Targets:
			raise TypeError('Please, the target set should be a set, tuple or list of concepts in the ontology. Check the type of the set of concept targets.')

		if not self.ShortPath and model[0] in set(self.CatPaths): 
			# Get all the shortest path lengths for a to b in DAG
			self.ShortPath = nx.shortest_path_length(self.DagStr)
		if not self.DicLevels and model[0] in set(self.CatLevels):
			DicLevels = nx.bellman_ford(self.DagStr, self.oroot) 
			self.DicLevels = DicLevels[-1]
			del DicLevels
			self.deep = -min(list(set(self.DicLevels.values())))

		# Set IC and CS specific arguments
		kwargs = kwargs
		# Checking than setting!
		if 'sigma' in kwargs and not isinstance(kwargs['sigma'], float):
			if kwargs['sigma'] < 0 or kwargs['sigma'] > 1: 
				raise ValueError('Please, Zho IC sigma value should be a float between 0 and 1. Check this parameter and try again.')
			raise TypeError('Please, Zho IC sigma value should be a float between 0 and 1. Check this parameter and try again.')
		if 'TermStats' in kwargs and not isinstance(kwargs['TermStats'], dict):
			raise TypeError('Please, the potential term count values should be a dictionary in which keys are ontology terms and values are the count values of these terms in the dataset under consideration. Check this parameter and try again.')
		if 'TermIC' in kwargs and not isinstance(kwargs['TermIC'], dict):
			raise TypeError('Please, the potential term IC scores should be a dictionary in which keys are ontology terms and values are the IC values of these terms in the dataset under consideration. Check this parameter and try again.')
#		#Setting potential paraneters for CS 
		if 'cf' in kwargs and not kwargs['cf'] in [0,1,2,3]:
			raise TypeError('Please, the potential correction factor should be an integer between 0 and 3:\n0. No correction is applied\n1. Graph-based (XGraSM or EICA) \n2. Relevance\n3. Information coefficient (simIC)\nCheck this parameter and try again.')
		if 'gr' in kwargs and not kwargs['gr'] in [0,1]:
			raise TypeError('Please, the choice for the Graph-based (XGraSM or EICA) correction factor should be 0 or 1:\n0. For XGraSM model and \n1. For EICA model.\nCheck this parameter and try again.')
		if 'ak' in kwargs and not (isinstance(kwargs['ka'], float) or kwargs['ka'] < 0):
			#Al-Mubaid
			raise TypeError('Please, ...')
		if 'aa' in kwargs and not (isinstance(kwargs['aa'], float) or kwargs['aa'] < 0):
			raise TypeError('Please, ...')
		if 'ab' in kwargs and not (isinstance(kwargs['ab'], float) or kwargs['ab'] < 0):
			raise TypeError('Please, ...')
		if 'kz' in kwargs and not (isinstance(kwargs['kz'], int) or kwargs['kz'] < 2):
			# Zhong Model
			 raise TypeError('Please, ...')
		if 'alpha' in kwargs and not (isinstance(kwargs['alpha'], float) or kwargs['alpha'] < 0): # Edge-based Li alpha parameter
			raise TypeError('Please, edge-based Li alpha value should be a float greater or equal to 0. Check this parameter and try again.')
		if 'beta' in kwargs and not (isinstance(kwargs['beta'], float) or kwargs['beta'] <= 0): # Edge-based Li alpha parameter
			raise TypeError('Please, edge-based Li beta value should be a float strictly greater than 0. Check this parameter and try again.')
		if 'cutoff' in kwargs and not isinstance(kwargs['cutoff'], float):
			if kwargs['cutoff'] < 0 or kwargs['cutoff'] > 1: 
				raise ValueError('Please, concept similarity cutoff value should be a float between 0 and 1. Check this parameter and try again.')
			raise TypeError('Please, concept similarity cutoff value should be a float between 0 and 1. Check this parameter and try again.')

		self.search(model, **kwargs)

if __name__=='__main__':
	from PySML.smlapps import ConceptEnrichment
	cenriched = ConceptEnrichment()
	fp = open('PySML/tests/ReferenceSetTest.txt')
	Background = {}
	for line in fp:
		tline = line.split()
		Background[tline[0].strip()] = [t.strip() for t in tline[1].split(',')]
	Background['P30613']
	fp.close()

	fp = open('PySML/tests/TargetSetTest.txt')
	Targets = set()
	for line in fp:
		Targets.add(line.strip())
	list(Targets)[:4]
	cenriched.enrichedConcept(Background, Targets)
