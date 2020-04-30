#!/usr/bin/python
# coding: utf-8
"""
This python file is part of the PySML library, which is a tool for Gene 
Ontology-based functional analysis using term information content 
measures.
This python code implements fuzzy protein search or identification based
on semantic similarity concepts. Furthemore, this is context independent 
and can be used for any ontology annotated dataset as population background 
or reference and is not a context-based search.

The main website for the PySML library is:
 
http://web.cbio.uct.ac.za/ITGOM/post-analysis-tools/mysml-dev/
https://github.com/gkm-software-dev/post-analysis-tools

where users can find essential information about obtaining PySML. It is 
freely downloadable under GNU General Public License (GPL), pre-compiled 
for Linux version and protected by copyright laws. Users are free to copy, 
modify, merge, publish, distribute and display information contained in 
the package, provided that it is done with appropriate citation of the 
package and by including the permission notice in all copies or substantial 
portions of the module contained in this package.

PySML is distributed in the hope that it will be useful, but WITHOUT 
ANY WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED 
TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE 
AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS 
BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY WHETHER IN AN 
ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN 
CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
See <http://www.gnu.org/licenses/>.

This code was written by 
    Gaston K. Mazandu <gaston.mazandu@uct.ac.za, gmazandu@gmail.com, 
                       kuzamunu@aims.ac.za>
    (c) 2020 under free software (GPL), all rights reserved.
"""

from __future__ import print_function, division
from math import exp, log, tanh, atan, pi as PI

import sys, time, random


from .. import __version__
from ..error import InputError
from ..imports.readontology import Ontology, output_str
from ..imports.tabulate import tabulate as tabs
from ..conceptsimilarity import ConceptSimilarity

try: # It is sufficient to import "scipy" here once for some useful mathematical objects, functions and algorithm needed.
	from scipy import log, exp, sort, arange, array, mean, ceil
	import scipy.stats as dst
except ImportError:
	print(InputError('networkx library not installed', 'scipy library ImportError on line 53 in \nentityclassification.py source code\n\nInstall the networkx library and try again ...\n'))
	sys.exit(0)
	
try:
	import networkx as nx
except ImportError:
	print(InputError('networkx library not installed', 'networkx library ImportError on line 60 in \nentityclassification.py source code\n\nInstall the networkx library and try again ...\n'))
	sys.exit(1)

class EntityIdentification(ConceptSimilarity):
	__cslot__ = []
	def __init__(self, ontofile = '', namespace='', is_a = None, part_of = None):
		"""Instantiate a new object for a given ontology to retrieve term simila-
		   rity scores with a DAG inherited from ConceptSimilarity class and
		   from there, inheriting InformationContent class.

        Arguments:
            ontology (str): the ontology file
            namespace (str): the name of the subontology as named in the ontology
            is_a (float between 0 and 1, if provided): a semantic value of the
                topological relation is_a.
            part_of (float between 0 and 1, if provided): a semantic value of the
                topological relation part_of.

        Example:
            >>> fct = EntityIentification()
        """
		ConceptSimilarity.__init__(self, ontofile, namespace, is_a, part_of)
		self.Backgrpund = {}
		self.fouts = {}
		self.comments = {}
	
	def search(self, model, **kwargs):
		"""Parameter as defined in the retrieveEntity model
		"""
		dtype = [('cpts', int),('ss', float)]; self.fouts = {}
		for t in self.Targets:
			self.fouts[t] = []
			for ent in self.Background:
				Pairs = [(t, s) for s in self.Background[ent]]
				Sim = self.conceptInterface(Pairs, model, **kwargs)
				try: 
					csim = list(sort(array([(s[1],Sim[s]) for s in Sim], dtype=dtype), order= 'ss'))
					#Prot, high SS, AvgSS, termI
					if csim[-1][1] > 0.0: self.fouts[t].append((ent, csim[-1][1], mean([c[1] for c in csim]), csim[-1][0]))
				except:
					pass
		VideKey = set([t for t in self.fouts if not self.fouts[t]])
		for t in VideKey: self.fouts.pop(t, False)
		
		n = len(self.fouts); P03 = {}; Bonf = len(self.Targets); NewKeep = {}
		for t in self.fouts:
			l = 0; ss = 0.0; pvalue = 1
			if t in self.fouts:
				Temp = [p[1] for p in self.fouts[t] if p[1] >= kwargs['score']]
				l = len(Temp); ss = mean(Temp) if Temp else 0.0
				if l > 0: 
					pr = 1.0*l/n
					pvalue = 1-dst.binom.cdf(l-1, n, pr)
				else: pvalue = 1.0
			if pvalue <= 0.0: pvalue = 0.0
			if pvalue >= 1.0: pvalue = 1.0
			PvBf = Bonf*pvalue # Bonferroni correction
			if PvBf >= 1.0: PvBf = 1.0
			P03[t] = (pvalue, PvBf, l, ss)
			NewKeep[t] = sorted([p for p in self.fouts[t] if p[1] >= kwargs['score']], key = lambda x: x[1], reverse=True)
		return P03, NewKeep

	def retrieveEntity(self, AnnotationData, Tconcepts, model = ('nunivers','universal'), **kwargs):
		"""
		The retrieveEntity method for protein functional classification tool. This 
		method allows the identification of an entity (e.g. gene or protein) fuzzy 
		contributing to a given process using their semantic similarity scores 
		computed based on a selected concept semantic similarity model. 
		
		Arguments:
		
			AnnotationData (dict): A dictionary with entity as key and set of concepts
			as value.
			
			Tconcets: Set of concepts for which associated entities needs to be identified
			
			model (tuple): The entity semantic similarity measure to be used. 
			Refer to the Supplementary for more details on symbols used for different 
			measures.
			 
			**kwargs can be used to set different parameters needed for the model chosen
			as well as other parameters associated to entity retrieval model, including
			
			score (float > 0.0): The threshold score providing the semantic similarity 
			degree at which entities are considered to be semantically close or similar in 
			the ontology structure and it is set to 0.3 by default.
			
			stream (int 0 or 1): An Enum parameter taking values of 0 or 1. It is set to 1
			to output results on the screen, to 0 to output results in a file.			
		"""
		
		if not 'score' in kwargs: kwargs['score'] = 0.3
		stream = kwargs['stream'] if 'stream' in kwargs else 1
		
		# Check other Concept Semantic Similarity parameters here
		self.parameterChecks(**kwargs)
		
		self.Background = {}; self.EntityMissing = {}
		if isinstance(AnnotationData, dict):
			for ent in AnnotationData:
				self.Background[ent] = set(); self.EntityMissing[ent] = set()
				for tt in AnnotationData[ent]:
					if tt in self.alt_id and self.alt_id[tt] in self.DagStr:
						self.Background[ent].add(self.alt_id[tt])
					elif tt in self.Dag and self.Dag.index(tt) in self.DagStr:
						self.Background[ent].add(self.Dag.index(tt))
					else: # Term is either obsolete or does not exist in the onto!
						self.EntityMissing[ent].add(tt)
				self.Background[ent].discard(self.oroot)
				if not self.Background[ent]: self.Background.pop(ent, False)
		else:
			print(InputError('AnnototationData - Value Error', 'AnnotationData should be either a dictionary: key (entity)/value (ontology annotation) mapping\nor a string representing a full path to an annotation file.\n\nPlease refer to the tool documentation, fix this issue and try again ...\n'))
			sys.exit(2)
		if not self.Background:
			print(InputError('Annots - Type Error', 'Sorry the background map is empty. Please check the type of the set of concept targets.'))
			sys.exit(3)
		
		self.Targets = set()
		if isinstance(Tconcepts, (list, set, tuple)):
			for tt in Tconcepts:
				if tt in self.alt_id and self.alt_id[tt] in self.DagStr: self.Targets.add(self.alt_id[tt])
				elif tt in self.Dag and self.Dag.index(tt) in self.DagStr: self.Targets.add(self.Dag.index(tt))
		else:
			print(InputError('Tconcepts - Type Error', 'Tconcepts shoul be a list, set or a tuple of concepts.\n\nPlease refer to the tool documentation, fix this issue and try again ...\n'))
			sys.exit(4)
		if not self.Targets:
			print(InputError('Tconcepts - Value Error','Please, the target set should not be empty! Check the type of the set of concept targets.'))
			sys.exit(5)

		
		DicLevels = nx.bellman_ford(self.DagStr, self.oroot) 
		self.DicLevels = DicLevels[-1]
		del DicLevels
		self.deep = -min(list(set(self.DicLevels.values())))
		
		now = time.time()			
		P03, New = self.search(model, **kwargs)
		
		VideKey = set([t for t in New if not New[t]])
		for t in VideKey: Mew.pop(t, False)
		
		sts = [model[0].capitalize(), model[1].capitalize()] if len(model)==2 else [model[0].capitalize(), '']
				
		print("\nFuzzy Identification of proteins/genes using    %s:"%("[Running entityidentification.py]",))
		print("Total number of possible target GO IDs in the list: %d"%(len(self.Targets),))
		print("Semantic Similarity approaches used is  : {}-{}, ".format(*sts))

		outputfile = 'EntityIdentificationResults%d.txt'%(random.randint(0,100000),)
				
		outs = []
		if not New:#Prot, high SS, AvgSS, termI
			print("\n\nUnfortunately, no entity reached the threshold or cutoff of %.5f indicated. However,\nAll set of proteins and associated SS scores are in the file: %s"%(kwargs['score'], outputfile))
			for t in self.fouts: #self.fouts[t].append((ent, csim[-1][1], mean([c[1] for c in csim]), csim[-1][0]))
				outs.append([self.Dag[t], -self.DicLevels[t]] + list(self.fouts[t][0][:-1]))
				for i in range(1, len(self.fouts)):
					outs.append(['', ''] + list(self.fouts[t][i][:-1]))
			fw = open(outputfile, 'w')
			headers = ['Concept-ID', 'Level', 'Entity', 'high SS', 'Avg SS']
			fw.write(tabs(outs, headers, tablefmt = 'plain', floatfmt=".5f", stralign="center"))
			fw.close()
			print("\n***************************************************************\n")
			sys.exit()			
		
		for t in New: # GO ID, GO term, Term Level, Number of proteins, Average SS, p-value, Bonferroni correction
			if not New[t]: continue
			outs.append([self.Dag[t], -self.DicLevels[t], P03[t][-2], New[t][0][0], round(New[t][0][1],5), round(P03[t][-1],5), P03[t][0], P03[t][1]])
			for i in range(1, len(New[t])):
				outs.append(['', '', '', New[t][i][0], round(New[t][i][1],5), '', '', ''])
		

		headers = ['Concept-ID', 'Level', '# of Entities', 'Entity', 'ESSS', 'Avg SS', 'p-value', 'Cp-values']
		if stream: print(tabs(outs, headers, tablefmt = 'grid', floatfmt=".5f", stralign="center"))
		else:
			print("\nGenerale Statistics for each target GO ID can be found in the file: [%s]"%(outputfile,))
			fp = open(outputfile, 'w')
			fp.write(tabs(outs, headers, tablefmt = 'plain', floatfmt=".5f", stralign="center"))
			fw.close()
		print("\nLegend of Entity Identification:")
		print("----------------------------------------")
		print("Cp-values stand for corrected p-values done using Bonferroni multiple correction model.")
			
		print("\nProcessing accomplished on %s"%str(time.asctime(time.localtime())))
		print("Total time elapsed is approximately %.2f %s"%(time.time()-now, 'seconds'))
		print("\n***************************************************************\n")

if __name__=='__main__':
	pass


