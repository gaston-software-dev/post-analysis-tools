#!/usr/bin/python
# -*- coding: utf8 -*-

"""Important Note:
This python file is part of the PySML tool, which is a tool for semantic
similarity between objects annotated by ontology terms or concepts, e,g.,
Gene Ontology-based functional analysis using term information content 
measures.
This python code implements fuzzy entity search or identification based
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
from functools import reduce

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
	print(InputError('networkx library not installed', 'networkx library ImportError on line 62 in \nentityclassification.py source code\n\nInstall the networkx library and try again ...\n'))
	sys.exit(1)

class ConceptEnrichment(ConceptSimilarity):
	__cslot__ = []
	def __init__(self, ontofile = '', namespace='', is_a = None, part_of = None):
		ConceptSimilarity.__init__(self, ontofile, namespace, is_a, part_of)
		self.fouts = {}
		self.Targets = set()
		self.TargetConcepts = set()
		self.TargetMissing = set()
		self.Background = {}
	
	def search(self, model=('nunivers','universal'), **kwargs):
		"""Parameters as defined in the enrichedConcepts method.
		"""
		cutoff = kwargs['score'] if 'score' in kwargs else 0.3
		self.fouts = dict()
		self.TargetConcepts = reduce(lambda x, y: x | y, [self.Background[ent] for ent in self.Targets]);
		
		if cutoff == 0: # Using traditional approach 
			for t in self.TargetConcepts:
				initset = set(lcc[t][1]+[t])
				l = 0; m = 0
				for p in background:
					if t in self.Background[p]:   # The term occurs in p
						m += 1
						if p in self.Targets: l += 1
					else:       # Get here only when considering true-path rule when agree==0
						for s in self.Background[p]:
							if t in nx.ancestors(self.DagStr, s): # Indicating that the term t occurs as an ancestor of s
								m += 1
								if p in self.Targets: l += 1
								break
				self.fouts[t] = (n,m)
		else:
			for t in self.TargetConcepts:
				m = n = 0
				for ent in self.Background:
					if t in self.Background[ent]: # The term occurs in ent
						n += 1
						if ent in self.Targets: m += 1
					else:
						Pairs = [(t, s) for s in self.Background[ent]]
						ds = self.conceptInterface(Pairs, model, **kwargs)
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

	def enrichedConcepts(self, AnnotationData, Tentities, model= None, **kwargs):
		"""	
		The enrichedConcepts method for fuzzy concept enrichment analysis tool. This 
		method allows fuzzy enriched annotations contributing to a given process using 
		their semantic similarity scores computed based on a selected concept semantic 
		similarity model. 
		
		Arguments:
		
			AnnotationData (dict): A dictionary with entity as key and set of concepts
			as value, representing the reference or background dataset.
			
			Tentities: Set of targeted entities (targets) representing the system under
			consideration.
			
			model (tuple): The entity semantic similarity measure to be used. Refer to 
			the Supplementary for more details on symbols used for different measures.
			 
			**kwargs can be used to set different parameters needed for the model chosen
			as well as other parameters associated to entity retrieval model, including
			
			score (float >= 0.0): The threshold score providing the semantic similarity 
			degree at which entities are considered to be semantically close or similar in 
			the ontology structure and it is set to 0.3 by default.
			
			stream (int 0 or 1): An Enum parameter taking values of 0 or 1. It is set to 1
			to output results on the screen, to 0 to output results in a file.			
		"""
		# Check other Concept Semantic Similarity parameters here
		self.parameterChecks(**kwargs)
		stream = kwargs['stream'] if 'stream' in kwargs else 1
		if not Tentities or not isinstance(Tentities, (list, set, tuple)):
			print(InputError('Tentities - Type or Value Error', 'Tentities should be a no empty list, set or a tuple of concepts.\n\nPlease refer to the tool documentation, fix this issue and try again ...\n'))
			sys.exit(3)
		if not AnnotationData or not isinstance(AnnotationData, dict):
			print(InputError('AnnotationData - Type or Value Error', 'Tentities should be a no empty dict, mapping concepts to entities.\n\nPlease refer to the tool documentation, fix this issue and try again ...\n'))
			sys.exit(4)
		if not 'score' in kwargs: kwargs['score'] = 0.3
		elif not isinstance(kwargs['score'], (float, int)) or kwargs['score'] > 1.0 or kwargs['score'] < 0.0:
			print(InputError('score parameter - Value or Type Error', 'The value of the parameter score, if provided\nshould be a positive float <= 1, by default set to 0.3.\n\nPlease refer to the tool documentation, fix this issue and try again ...\n'))
			sys.exit(5)
		else:
			print(InputError('score parameter - Type Error', 'The value of the parameter score, if provided\nshould be a positive float <= 1, by default set to 0.3.\n\nPlease refer to the tool documentation, fix this issue and try again ...\n'))
			sys.exit(6)
		if not 'pvalue' in kwargs: kwargs['pvalue'] = 0.05
		elif not isinstance(kwargs['pvalue'], float) or kwargs['pvalue'] > 1.0 or kwargs['pvalue'] < 0.0:
			print(InputError('pvalue parameter - Value or Type Error', 'The value of the parameter pvalue, if provided\nshould be a positive float<=1, by default set to 0.05.\n\nPlease refer to the tool documentation, fix this issue and try again...\n'))
			sys.exit(7)
		else:
			print(InputError('score parameter - Type Error', 'The value of the parameter pvalue, if provided\nshould be a positive float <= 1, by default set to 0.05.\n\nPlease refer to the tool documentation, fix this issue and try again ...\n'))
			sys.exit(8)

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
			print(InputError('AnnototationData - Type Error', 'AnnotationData should be either a dictionary: key (entity)/value (ontology annotation) mapping\nor a string representing a full path to an annotation file.\n\nPlease refer to the tool documentation, fix this issue and try again ...\n'))
			sys.exit(9)
		if not self.Background:
			print(InputError('AnnotationData - Value Error', 'Sorry the background map is empty because concepts provided\ncould not be mapped to the ontology. Please check the type of the set of concept targets.'))
			sys.exit(10)

		self.Targets = set(); self.TargetMissing = set()
		for p in Tentities:
			if p in self.Background: self.Targets.add(p)
			else: self.TargetMissing.add(p)

		if not self.Targets:
			print(InputError('Tentities - Value Error','Please, the target entities provided not mapped! Check the type of the set of concept targets.'))
			sys.exit(9)

		now = time.time()
		print("\nComputing different parameters and concept frequency now, this may take time...")
		DicLevels = nx.bellman_ford(self.DagStr, self.oroot) 
		self.DicLevels = DicLevels[-1]
		del DicLevels
		self.deep = -min(list(set(self.DicLevels.values())))

		self.search(model, **kwargs)
		
		# computing p-values using hypergeometric distribution and corrected p-values by Bonferroni
		tn = time.time()
		print("\nComputing p-value now, this may take time...", end=',')
		Bonf = len(self.TargetConcepts); N = len(self.Background); n = len(self.Targets); P03 = {}
		for t in self.fouts:
			Pv = 1.0
			Pv = 1-dst.hypergeom.cdf(self.fouts[t][0]-1, N, self.fouts[t][1], n) if self.fouts[t][0] > 0 else 1.0
			if Pv <= 0.0: Pv = 0.0 # The Pv can become negative because of floating point
			if Pv >= 1.0: Pv = 1.0 # The Pv can go greater than 1.0 because of floating point
			PvBf = Bonf*Pv 		   # Bonferroni multiple correction
			if PvBf >= 1.0: PvBf = 1.0
			if Pv < kwargs['pvalue']: P03[t] = (Pv, PvBf)
		print("\rComputing p-value done, time elapsed: %d %s"%(int(round(time.time()-tn)), 'seconds...'))
		
		if not P03:
			print("\n\nUnfortunately, no enriched concept has been identified for\nthreshold or cutoff of %.5f indicated."%(kwargs['pvalue'],))
			print("\n***************************************************************\n")
			sys.exit(10)	

		print("\nFiltering identified enriched concepts now, this may take time...")
		# Filtering set of significant terms
		Fl03 = set(P03.keys())
		for t in Fl03.copy():
			if not t in Fl03: continue
			Fl03 -= nx.ancestors(self.DagStr, t) # Remove ancestors of the term t
		
		Fl03 = set([t for t in Fl03 if P03[t][1] < kwargs['pvalue']])
		if not Fl03: # No term passed the Bonferroni correction	
			pass
		
		# building output
		outs = []
		for t in Fl03: # Concept-ID, Term Level, p-value, Bonferroni correction
			outs.append([self.Dag[t], -self.DicLevels[t], P03[t][0], P03[t][1]])
		outs = sorted(outs, key = lambda x: (x[-1], -x[1]))
		sts = [model[0].capitalize(), model[1].capitalize()] if len(model)==2 else [model[0].capitalize(), '']
		
		print("\nIdentifying fuzzy enriched concepts using : %s"%("[conceptenrichment.py module]",))
		print("Total number of concets in the target set : %d"%(Bonf,))
		print("Number of enriched GO terms detected      : %d"%(len(outs),))
		print("Semantic Similarity approaches used is    : {}-{}\n".format(*sts))
		
		headers = ['Concept-ID', 'Level', 'p-value', 'Adj-p-values']
		if stream: print(tabs(outs, headers, tablefmt = 'grid', floatfmt="1.5e", stralign="center"))
		else:
			outputfile = 'EntityIdentificationResults%d.txt'%(random.randint(0,100000),)
			print("\nGenerale Statistics for each target concepts can be found in the file: [%s]"%(outputfile,))
			fp = open(outputfile, 'w')
			fp.write(tabs(outs, headers, tablefmt = 'plain', floatfmt="1.2e", stralign="center"))
			fw.close()
		
		print("\nProcessing accomplished on %s"%str(time.asctime(time.localtime())))
		print("Total time elapsed is approximately: %d %s"%(int(round((time.time()-now)/60)), 'minutes.'))
		print("\n************************************************************************************\n")

if __name__=='__main__':
	pass
