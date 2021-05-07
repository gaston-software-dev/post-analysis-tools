#!/usr/bin/python
# -*- coding: utf8 -*-

"""Important Note:
This python file is part of the PySML tool, which is a tool for semantic
similarity between objects annotated by ontology terms or concepts, e,g.,
Gene Ontology-based functional analysis using term information content 
measures.
This python module implements entity (or set of terms) semantic similarity 
or functional similarity measures. 
Furthemore, this is context independent and can be used for any ontology.
It is not only limited to the gene ontology or biomedical ontologies.

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
from math import exp, log, tanh, atan, sqrt, pi as PI

from functools import reduce

import sys, re


from . import __version__
from .error import InputError
from .imports.readontology import Ontology, output_str
from .imports.tabulate import tabulate as tabs
from .conceptsimilarity import ConceptSimilarity

try:
	import networkx as nx
except ImportError:
	print(InputError('networkx library not installed', 'networkx library ImportError on line 48 in \nconceptsimilarity.py source code\n\nInstall the networkx library and try again ...\n'))
	sys.exit(0)

class EntitySimilarity(ConceptSimilarity):
	"""Definition of the `EntitySimilarity` class for Entity semantic 
	   similarity score retrieval.

    The EntitySimilarity class actually behaves as a factory, creating
    a EntitySimilarity object via the default Python syntax only::
          >>> from PySML import EntitySimilarity
          >>> SScores = EntitySimilarity()
          >>> # This is very simple and interesting
          
	"""
	
	icf = {'avg':"Average", 'bma':"Best Match Average", 'abm':"Average Best Match", 'bmm':"Best Match Maximum", 'hdf':"Hausdorff (HDF) distance-based", 'vhdf':"HDF variant by Lerman et al.", 'max':"Maximum"}
	edf = {'aln':"Al-Mubaid & Nagar edge-based", 'intel':"IntelliGO by Benabderrahmane et al. edge-based ", 'spgk':"Shortest Path Graph Kernel edge-based", 'lp':'SimLP measure by Gentleman', 'ye': 'Normalized version of SimLP by Ye et al.'}
	icg_1 = {'simgic':"IC-%s-based Jaccard measure", 'simdic':"IC-%s-based Dice (Czekanowski or Lin like measure)", 'simuic': "IC-%s-based universal index measure", 'simcou':"IC-%s-based usual normalization Cosine", 'simcot': "IC-%s-based Tanimoto coefficient Cosine"}
	icg_2 = {'simui': "Ontology-based Jaccard measure", 'simub':"Ontology-based universal index", 'simdb':"Ontology-based Dice (Czekanowski or Lin like measure)", 'simnto': "Ontology-based normalized term overlap (NTO)", 'simcub':"Ontology-based usual normalization Cosine", 'simctb':"Ontology-based Tanimoto coefficient Cosine"}
	nnt = {'cho':"No-ontology based Chi et al.", 'ald':"No-ontology based Ali and Diane", 'kstats':"No-ontology based Kappa-statistics", 'nto':"No-ontology based normalized term overlap (NTO)", 'ub':"No-ontology based universal index", 'db':"No-ontology based Dice (Czekanowski or Lin like measure)", 'ui':"No-ontology-based Jaccard"}
	
	__fslot__ = ['fouts', 'EntityPairs', 'EntityMissing']
	
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
            >>> funsim = EntitySimilarity()
        """
		ConceptSimilarity.__init__(self, ontofile, namespace, is_a, part_of)
		self.fouts = {}
		self.EntityMissing = {}
		self.measures = []
		self.comments = {}

	def groupwise(self, annot, pairids, fct, **kwargs):
		"""Implements ontology group-wise based entity semantic similarity measures. 
		"""
		icapps = set([ft[1] for ft in fct if len(ft) > 1])
		self.getIC(list(icapps), **kwargs)
		Scores = self.AppScores.copy()
		if 'wang' in Scores: Scores['wang'] = dict([(t, sum(Scores['wang'][t].values())) for t in Scores['wang']])

		self.fouts.update(dict([(t, {}) for t in fct if not t in self.fouts]))
		for p, q in pairids:
			Pairs = [(s,t) for s in annot[p] for t in annot[q]]
			pt = reduce(lambda x, y: x | y, [nx.ancestors(self.DagStr,t) for t in annot[p]]); pt |= annot[p]
			qt = reduce(lambda x, y: x | y, [nx.ancestors(self.DagStr,t) for t in annot[q]]); qt |= annot[q]
			upq = list(pt|qt)
			for ft in fct:
				if ft[0]=='simgic':
					if p==q: self.fouts[ft][(p,q)] = 1.0
					else:
						sanc1 = sum([Scores[ft[1]][t] for t in pt & qt])
						sanc2 = sum([Scores[ft[1]][t] for t in pt | qt])
						self.fouts[ft][(p,q)] = round(sanc1/sanc2,5)
				elif ft[0]=='simdic':
					if p==q: self.fouts[ft][(p,q)] = 1.0
					else:
						sanc1 = sum([Scores[ft[1]][t] for t in pt & qt])
						sanc2 = sum([Scores[ft[1]][t] for t in pt])
						sanc3 = sum([Scores[ft[1]][t] for t in qt])
						self.fouts[ft][(p,q)] = round(2.0*sanc1/(sanc2+sanc3),5)
				elif ft[0]=='simuic':
					if p==q: self.fouts[ft][(p,q)] = 1.0
					else:
						sanc1 = sum([Scores[ft[1]][t] for t in pt & qt])
						sanc2 = sum([Scores[ft[1]][t] for t in pt])
						sanc3 = sum([Scores[ft[1]][t] for t in qt])
						self.fouts[ft][(p,q)] = round(sanc1/max(sanc2, sanc3),5)
				elif ft[0]=='simcou':
					if p==q: self.fouts[ft][(p,q)] = 1.0
					else:
						sanc1 = [Scores[ft[1]][t] if t in pt else 0.0 for t in upq]
						sanc2 = [Scores[ft[1]][t] if t in qt else 0.0 for t in upq]
						a = sum(map(lambda x, y: x*y, sanc1, sanc2))
						b = sqrt(sum(map(lambda x, y: x*y, sanc1, sanc1)))
						c = sqrt(sum(map(lambda x, y: x*y, sanc2, sanc2)))
						self.fouts[ft][(p,q)] = round(a/(b*c),5)
				elif ft[0]=='simcot':
					if p==q: self.fouts[ft][(p,q)] = 1.0
					else:
						sanc1 = [Scores[ft[1]][t] if t in pt else 0.0 for t in upq]
						sanc2 = [Scores[ft[1]][t] if t in qt else 0.0 for t in upq]
						a = sum(map(lambda x, y: x*y, sanc1, sanc2))
						b = sum(map(lambda x, y: x*y, sanc1, sanc1))
						c = sum(map(lambda x, y: x*y, sanc2, sanc2))
						self.fouts[ft][(p,q)] = round(a/(b+c-a),5)
				elif ft[0]=='simui':
					if p==q: self.fouts[ft][(p,q)] = 1.0
					else:
						self.fouts[ft][(p,q)] = round(len(pt&qt)/len(pt|qt),5)
				elif ft[0]=='simub':
					if p==q: self.fouts[ft][(p,q)] = 1.0
					else:
						self.fouts[ft][(p,q)] = round(len(pt&qt)/max(len(pt),len(qt)),5)
				elif ft[0]=='simdb':
					if p==q: self.fouts[ft][(p,q)] = 1.0
					else:
						self.fouts[ft][(p,q)] = round(2*len(pt&qt)/(len(pt)+len(qt)),5)
				elif ft[0]=='simnto':
					if p==q: self.fouts[ft][(p,q)] = 1.0
					else:
						self.fouts[ft][(p,q)] = round(len(pt&qt)/min(len(pt),len(qt)),5)
				elif ft[0]=='simcub':
					if p==q: self.fouts[ft][(p,q)] = 1.0
					else:
						sanc1 = [1.0 if t in pt else 0.0 for t in upq]
						sanc2 = [1.0 if t in qt else 0.0 for t in upq]
						a = sum(map(lambda x, y: x*y, sanc1, sanc2))
						b = sqrt(sum(map(lambda x, y: x*y, sanc1, sanc1)))
						c = sqrt(sum(map(lambda x, y: x*y, sanc2, sanc2)))
						self.fouts[ft][(p,q)] = round(a/(b*c),5)
				elif ft[0]=='simctb':
					if p==q: self.fouts[ft][(p,q)] = 1.0
					else:
						sanc1 = [1.0 if t in pt else 0.0 for t in upq]
						sanc2 = [1.0 if t in qt else 0.0 for t in upq]
						a = sum(map(lambda x, y: x*y, sanc1, sanc2))
						b = sqrt(sum(map(lambda x, y: x*y, sanc1, sanc1)))
						c = sqrt(sum(map(lambda x, y: x*y, sanc2, sanc2)))
						self.fouts[ft][(p,q)] = round(a/(b+c-a),5)

	def pairwise(self, annot, pairids, fct, **kwargs):
		"""Implements ontology pair-wise based entity semantic similarity measures. 
		"""
		self.fouts.update(dict([(t, {}) for t in fct if not t in self.fouts]))
		for p, q in pairids:
			Pairs = []
			for k in [(s,t) for s in annot[p] for t in annot[q]]:
				if not (k in Pairs or (k[1],k[0]) in Pairs): Pairs.append(k)
			for ft in fct: # Those using concept similarity score
				try: met = (ft[1], ft[2]) if len(ft)==3 else (ft[1],)
				except: pass
				au = {}
				if ft[0] in self.icf:
					tmp = self.conceptInterface(Pairs, models = met, **kwargs)
					for (t, s) in tmp: au[(t,s)] = au[(s,t)] = tmp[(t,s)]
				elif ft[0] in self.edf:
					if not self.ShortPath:
						self.ShortPath = nx.shortest_path_length(self.DagStr)
					if not self.DicLevels:
						DicLevels = nx.bellman_ford(self.DagStr, self.oroot) 
						self.DicLevels = DicLevels[-1]
						self.deep = -min(list(set(self.DicLevels.values())))
						del DicLevels

				if ft[0]=='avg':
					if p==q: self.fouts[ft][(p,q)] = 1.0
					else: self.fouts[ft][(p,q)] = round(sum(au.values())/len(au),5)
				elif ft[0]=='bma':
					if p==q: self.fouts[ft][(p,q)] = 1.0
					else:
						au1 = [max([au[(s,t)] for t in annot[q]]) for s in annot[p]]
						au2 = [max([au[(s,t)] for t in annot[p]]) for s in annot[q]]
						ss = sum(au1)/len(au1)+sum(au2)/len(au2)
						self.fouts[ft][(p,q)] = round(ss/2,5)
				elif ft[0]=='abm':
					if p==q: self.fouts[ft][(p,q)] = 1.0
					else:
						au1 = [max([au[(s,t)] for t in annot[q]]) for s in annot[p]]
						au2 = [max([au[(s,t)] for t in annot[p]]) for s in annot[q]]
						aut = au1 + au2
						self.fouts[ft][(p,q)] = round(sum(aut)/len(aut),5)
				elif ft[0]=='bmm': # Best Match Maximum corresponding to modified HDF.
					if p==q: self.fouts[ft][(p,q)] = 1.0
					else:
						au1 = [max([au[(s,t)] for t in annot[q]]) for s in annot[p]]
						au2 = [max([au[(s,t)] for t in annot[p]]) for s in annot[q]]
						aut = [sum(au1)/len(au1), sum(au2)/len(au2)]
						self.fouts[ft][(p,q)] = round(max(aut),5)
				elif ft[0]=='hdf':
					if p==q: self.fouts[ft][(p,q)] = 1.0
					else:
						au1 = [max([au[(s,t)] for t in annot[q]]) for s in annot[p]]
						au2 = [max([au[(s,t)] for t in annot[p]]) for s in annot[q]]
						ss = 1.0-max(1.0-min(au1), 1.0-min(au2))
						self.fouts[ft][(p,q)] = round(ss,5)
				elif ft[0]=='vhdf':
					if p==q: self.fouts[ft][(p,q)] = 1.0
					else:
						au1=[(1-max([au[(s,t)] for t in annot[q]]))**2 for s in annot[p]]
						au2=[(1-max([au[(s,t)] for t in annot[p]]))**2 for s in annot[q]]
						ss = 1-(sqrt(sum(au1)/len(au1))+sqrt(sum(au1)/len(au1)))/2
						self.fouts[ft][(p,q)] = round(ss,5)
				elif ft[0]=='max':
					if p==q: self.fouts[ft][(p,q)] = 1.0
					else: self.fouts[ft][(p,q)] = round(max(au.values()),5)
							
				elif ft[0]=='aln':# For Al-Mubaid&Nagar
					if p==q: self.fouts[ft][(p,q)] = 1.0
					else:
						try: aaln = kwargs['aaln']
						except KeyError: aaln = 1
						for (s,t) in Pairs:
							panc = nx.ancestors(self.DagStr, s); panc.add(s)
							qanc = nx.ancestors(self.DagStr, t); qanc.add(t)
							canc = panc & qanc
							dpl = min([self.DicLevels[a] for a in canc])
							slca = [a for a in canc if self.DicLevels[a]==dpl]
							vdist = sys.maxsize
							for c in slca:
								spl = self.ShortPath[c][s] + self.ShortPath[c][t]
								if spl < vdist: vdist = spl
							au[(s,t)] = au[(t,s)] = vdist
						self.fouts[ft][(p,q)] = exp(-aaln*sum(au.values())/len(au))
				elif ft[0]=='spgk': #  spgk
					if not self.deep: self.deep = -min(list(set(self.DicLevels.values())))
					if p==q: self.fouts[ft][(p,q)] = 2*self.deep*(self.deep+1)
					else:
						ss = 0.0
						for (s,t) in Pairs:
							panc = nx.ancestors(self.DagStr, s); panc.add(s)
							qanc = nx.ancestors(self.DagStr, t); qanc.add(t)
							canc = panc & qanc
							dpl = min([self.DicLevels[a] for a in canc])
							slca = [a for a in canc if self.DicLevels[a]==dpl]
							for c in slca: ss += (-self.DicLevels[c])*(-self.DicLevels[c]+1)
						self.fouts[ft][(p,q)] = ss
				elif ft[0]=='intel':# IntelliGO
					if p==q: self.fouts[ft][(p,q)] = 1.0
					else:
						sset1 = set(annot[p]); sset2 = set(annot[q])
						sset1 = list(sset1); sset2 = list(sset2)
						pairids = [(s, t) for s in sset1 for t in sset2]
						wu = self.swu(pairids); ss = sum([wu[(s,t)] for (s,t) in Pairs])
						pp = [(s,t) for s in sset1 for t in sset1]
						wu = self.swu(pp); ss1 = sum([wu[(s,t)] for (s,t) in pp])
						qq = [(s,t) for s in sset2 for t in sset2]
						wu = self.swu(qq); ss2 = sum([wu[(s,t)] for (s,t) in qq])
						self.fouts[ft][(p,q)] = ss/(sqrt(ss1)*sqrt(ss2))
				elif ft[0]=='lp': # SimLP not normalized!
					pt = reduce(lambda x, y: x | y, [nx.ancestors(self.DagStr,t) for t in annot[p]]); pt |= annot[p]
					qt = reduce(lambda x, y: x | y, [nx.ancestors(self.DagStr,t) for t in annot[q]]); qt |= annot[q]
					if not (pt & qt): self.fouts[ft][(p,q)] = 0.0
					else: self.fouts[ft][(p,q)] = -min([self.DicLevels[s] for s in pt & qt])
				elif ft[0]=='ye': # SimYE normalized version of SimLP
					if p==q: self.fouts[ft][(p,q)] = 1.0
					else:
						pt = reduce(lambda x, y: x | y, [nx.ancestors(self.DagStr,t) for t in annot[p]]); pt |= annot[p]
						qt = reduce(lambda x, y: x | y, [nx.ancestors(self.DagStr,t) for t in annot[q]]); qt |= annot[q]
						if not (pt & qt): self.fouts[ft][(p,q)] = 0.0
						else:
							dmin = -max([self.DicLevels[s] for s in pt & qt])
							dmax = -min([self.DicLevels[s] for s in pt & qt])
							if dmin==dmax: self.fouts[ft][(p,q)] = 0.0
							else: self.fouts[ft][(p,q)] = min([(-self.DicLevels[s]-dmin)/(dmax-dmin) for s in pt & qt])
						
	def ontology_indep(self, annot, pairids, fct, **kwargs):
		"""Implements measures that do not consider the ontology structure 
		"""
		self.fouts.update(dict([(t, {}) for t in fct if not t in self.fouts]))
		sdannot = reduce(lambda x, y: x | y, [set(annot[p]) for p in annot])
		
		for p, q in pairids:
			for ft in fct:
				if ft[0]=='cho': # SimCHO
					if p==q: self.fouts[ft][(p,q)] = 1.0
					else:
						interpq = set(annot[p]) & set(annot[q])
						if not interpq: self.fouts[ft][(p,q)] = 0.0
						else:
							cpq = min([sum([1 if s in annot[z] else 0 for z in annot]) for s in interpq])
							minmax = [sum([1 if s in annot[z] else 0 for z in annot]) for s in sdannot]
							cmin = min(minmax); cmax = max(minmax)
							self.fouts[ft][(p,q)] = log(cpq/cmax)/log(cmin/cmax)
				elif ft[0]=='ald': # SimALD
					if p==q: self.fouts[ft][(p,q)] = 1.0
					else:
						interpq = set(annot[p]) & set(annot[q])
						if not interpq: self.fouts[ft][(p,q)] = 0.0
						else:
							ss = dict([(s, sum([1 if s in annot[z] else 0 for z in annot])) for s in sdannot])
							sss = sum(ss.values())
							self.fouts[ft][(p,q)] = max([1 - ss[s]/sss for s in interpq])
				elif ft[0]=='kstats': # Kappa-stats
					if p==q: self.fouts[ft][(p,q)] = 1.0
					else:
						z = list(set(annot[p]) | set(annot[q])); n = len(z)
						pp = [1 if s in annot[p] else 0 for s in z]; qq = [1 if s in annot[q] else 0 for s in z]
						rho = sum([pp[i]==qq[i] for i in range(n)])/n
						alpha = sum([pp[i]==0 for i in range(n)])*sum([qq[i]==0 for i in range(n)])
						alpha += sum([pp[i]==1 for i in range(n)])*sum([qq[i]==1 for i in range(n)])
						alpha /= (n*n)
						self.fouts[ft][(p,q)] = (rho - alpha)/(1-alpha)
				elif ft[0]=='nto': # NTO-like
					if p==q: self.fouts[ft][(p,q)] = 1.0
					else:
						pp = set(annot[p]); qq = set(annot[q])
						self.fouts[ft][(p,q)] = len(pp & qq)/min(len(pp), len(qq))
				elif ft[0]=='ub': # UB-like
					if p==q: self.fouts[ft][(p,q)] = 1.0
					else:
						pp = set(annot[p]); qq = set(annot[q])
						self.fouts[ft][(p,q)] = len(pp & qq)/max(len(pp), len(qq))
				elif ft[0]=='db': # DB-like
					if p==q: self.fouts[ft][(p,q)] = 1.0
					else:
						pp = set(annot[p]); qq = set(annot[q])
						self.fouts[ft][(p,q)] = 2*len(pp & qq)/(len(pp) + len(qq))
				elif ft[0]=='ui': # UI-like
					if p==q: self.fouts[ft][(p,q)] = 1.0
					else:
						pp = set(annot[p]); qq = set(annot[q])
						self.fouts[ft][(p,q)] = len(pp & qq)/len(pp | qq)
	
	def processEntityAnnot(self, Annots):
		"""Processing entity-annotation map transforming physical dataset
		to logical dataset preparing entity semantic similarity scores
		"""
		CurrentAnnots = {}; self.EntityMissing = {}
		if isinstance(Annots, dict):
			for ent in Annots:
				CurrentAnnots[ent] = set(); self.EntityMissing[ent] = set()
				for tt in Annots[ent]:
					if tt in self.alt_id and self.alt_id[tt] in self.DagStr:
						CurrentAnnots[ent].add(self.alt_id[tt])
					elif tt in self.Dag and self.Dag.index(tt) in self.DagStr:
						CurrentAnnots[ent].add(self.Dag.index(tt))
					else: # Term is either obsolete or does not exist in the onto!
						self.EntityMissing[ent].add(tt)
				CurrentAnnots[ent].discard(self.oroot)
		elif isinstance(Annots, str): # Probably annots are in a file! Check if the path/to/file exists
			try:
				pathfile = os.path.abspath(Annots.strip())
				# Read into the file and read concepts
				fp = open(pathfile)			
				for line in fp:
					tline = line.strip()
					if not tline: continue
					tline = re.split("\s+|,|;|", tline); ent = tline[0]
					CurrentAnnots[ent] = set(); self.EntityMissing[ent] = set()
					for tt in [s.strip() for s in tline[1:]]:
						if tt in self.alt_id and self.alt_id[tt] in self.DagStr:
							CurrentAnnots[ent].add(self.alt_id[tt])
						elif tt in self.Dag and self.Dag.index(tt) in self.DagStr:
							CurrentAnnots[ent].add(self.Dag.index(tt))
						else: # Term is either obsolete or does not exist in the onto!
							self.EntityMissing[ent].add(tt)
					CurrentAnnots[ent].discard(self.oroot)
			except:
				try: # Possibly list of terms is provided as comma, space, semi column or colum separated string
					fp = [line.strinp() for line in Annots.split('\n')]
					for line in fp:
						tline = line.strip()
						if not tline: continue
						tline = re.split("\s+|,|;|", tline); ent = tline[0]
						CurrentAnnots[ent] = set(); self.EntityMissing[ent] = set()
						for tt in [s.strip() for s in tline[1:]]:
							if tt in self.alt_id and self.alt_id[tt] in self.DagStr:
								CurrentAnnots[ent].add(self.alt_id[tt])
							elif tt in self.Dag and self.Dag.index(tt) in self.DagStr:
								CurrentAnnots[ent].add(self.Dag.index(tt))
							else: # Term is either obsolete or does not exist in the onto!
								self.EntityMissing[ent].add(tt)
						CurrentAnnots[ent].discard(self.oroot)
				except:
					print(InputError('Path to concept file - not found', 'When ontology concepts are provided in a file, the full path to the file needs to be provided.\n\nPlease refer to the tool documentation, fix this issue and try again ...\n'))
					sys.exit(1)
		else:
			print(InputError('Annots - Type Error', 'Annots should be either a dictionary: key (entity)/value (ontology annotation) mapping\nor a string representing a full path to an annotation file.\n\nPlease refer to the tool documentation, fix this issue and try again ...\n'))
			sys.exit(2)
		if CurrentAnnots: return CurrentAnnots
		else:
			print(InputError('Annots - Value Error', 'Annots should be either a dictionary: key (entity)/value (ontology annotation) mapping\nor a string representing a full path to an annotation file.\n\nPlease refer to the tool documentation, fix this issue and try again ...\n'))
			sys.exit(3)
	
	def cleanEntityMeasures(self, measures):
		"""Transforming physical measures to be processed to logical measures for
		preparing the computation of entity semantic similarity scores.
		"""
		self.measures = []
		if measures is None or not measures:
			met = ('bma', 'nunivers','universal')
			self.measures.append(met)
			self.comments[met] = "runs %s with %s concept similarity model under %s IC-approach"%(self.icf[met[0]], met[1].capitalize(), met[2].capitalize())
		elif isinstance(measures, (tuple, list)):
			for meas in measures:
				if isinstance(meas, (list, tuple)) and len(meas)==3:
					if meas[0].lower() in self.icf and meas[1].lower() in self.Nodebased and meas[2].lower() in self.AppMods:
						met = (meas[0].lower(), meas[1].lower(), meas[2].lower())
						if not met in self.measures:
							self.measures.append(met)
							self.comments[met] = "runs %s with %s concept similarity model under %s IC-approach"%(self.icf[met[0]], met[1].capitalize(), met[2].capitalize())
					else:
						print(InputError('Model::%s - Structure Error'%(meas[0],), 'In the current setting, this model cannot be combined with %s and %s, as suggested.\n\nPlease refer to the tool documentation, fix this issue and try again ...\n'%(meas[1], meas[2])))
						sys.exit(4)
				elif isinstance(meas, (list, tuple)) and len(meas)==2:
					if meas[0].lower() in self.icf and meas[1].lower() in self.Nodebased:
						met = (meas[0].lower(), meas[1].lower(), 'universal')
						if not met in self.measures:
							self.measures.append(met)
							self.comments[met] = "runs %s with %s concept similarity model under %s IC-approach"%(self.icf[met[0]], met[1].capitalize(), met[2].capitalize())
					elif meas[0].lower() in self.icf and meas[1].lower() in self.AppMods:
						met = (meas[0].lower(), 'nunivers', meas[1].lower())
						if not met in self.measures:
							self.measures.append(met)
							self.comments[met] = "runs %s with %s concept similarity model under %s IC-approach"%(self.icf[met[0]], met[1].capitalize(), met[2].capitalize())
					elif meas[0].lower() in self.icf and meas[1].lower() in self.CatLevels|self.CatPaths: # For edge-based CS
						met = (meas[0].lower(), meas[1].lower())
						if not met in self.measures:
							self.measures.append(met)
							self.comments[met] = "runs %s with edge-based %s concept similarity model"%(self.icf[met[0]], met[1].capitalize())
					elif meas[0].lower() in self.icg_1 and meas[1].lower() in self.AppMods:
						#Case where you have groupwise and IC approach
						met = (meas[0].lower(), meas[1].lower())
						if not met in self.measures:
							self.measures.append(met)
							self.comments[met] = "runs %s"%(self.icg_1[met[0]]%(met[1].capitalize(),),)
					else:
						print(InputError('Model::%s - Structure Error'%(meas[0],), 'In the current setting, this model cannot be combined with %s, as suggested.\n\nPlease refer to the tool documentation, fix this issue and try again ...\n'%(meas[1],)))
						sys.exit(5)
				elif isinstance(meas, (list, tuple)) and len(meas)==1:
					if meas[0].lower() in self.icf:
						met = (meas[0].lower(), 'nunivers', 'universal')
						if not met in self.measures:
							self.measures.append(met)
							self.comments[met] = "runs %s with %s concept similarity model under %s IC-approach"%(self.icf[met[0]], met[1].capitalize(), met[2].capitalize())
					elif meas[0].lower() in self.AppMods:
						met = ('bma', 'nunivers', meas[0].lower())
						if not met in self.measures:
							self.measures.append(met)
							self.comments[met] = "runs %s with %s concept similarity model under %s IC-approach"%(self.icf[met[0]], met[1].capitalize(), met[2].capitalize())
					elif meas[0].lower() in self.CatLevels|self.CatPaths:
						met = ('bma', meas[0].lower())
						if not met in self.measures:
							self.measures.append(met)
							self.comments[met] = "runs %s with edge-based %s concept similarity model"%(self.icf[met[0]], met[1].capitalize())
					elif meas[0].lower() in self.Nodebased:
						met = ('bma', meas[0].lower(), 'universal')
						if not met in self.measures:
							self.measures.append(met)
							self.comments[met] = "runs %s with %s concept similarity model under %s IC-approach"%(self.icf[met[0]], met[1].capitalize(), met[2].capitalize())
					elif meas[0].lower() in self.edf:
						met = (meas[0].lower(), )
						if not met in self.measures:
							self.measures.append(met)
							self.comments[met] = "runs %s edge-based concept similarity model"%(self.edf[met[0]],)
					elif meas[0].lower() in self.icg_1:
						met = (meas[0].lower(), 'universal')
						if not met in self.measures:
							self.measures.append(met)
							self.comments[met] = "runs %s"%(self.icg_1[met[0]]%(met[1].capitalize(),),)
					elif meas[0].lower() in self.icg_2 or meas[0].lower() in self.nnt:
						met = (meas[0].lower(), )
						cc = self.icg_2[met[0]] if met[0] in self.icg_2 else self.nnt[met[0]]
						if not met in self.measures:
							self.measures.append(met)
							self.comments[met] = "runs %s concept semantic similarity model"%(cc,)
					else:
						print(InputError('Model::%s - Structure Error'%(meas[0],), 'In the current setting, this model does not exist or combination format issue.\n\nPlease refer to the tool documentation, fix this issue and try again ...\n'))
						sys.exit(6)
						
				elif isinstance(meas, str):
					if meas.lower() in self.icf:
						met = (meas.lower(), 'nunivers', 'universal')
						if not met in self.measures:
							self.measures.append(met)
							self.comments[met] = "runs %s with %s concept similarity model under %s IC-approach"%(self.icf[met[0]], met[1].capitalize(), met[2].capitalize())
					elif meas.lower() in self.AppMods:
						met = ('bma', 'nunivers', meas.lower())
						if not met in self.measures:
							self.measures.append(met)
							self.comments[met] = "runs %s with %s concept similarity model under %s IC-approach"%(self.icf[met[0]], met[1].capitalize(), met[2].capitalize())
					elif meas.lower() in self.CatLevels|self.CatPaths:
						met = ('bma', meas.lower())
						if not met in self.measures:
							self.measures.append(met)
							self.comments[met] = "runs %s with edge-based %s concept similarity model"%(self.icf[met[0]], met[1].capitalize())
					elif meas.lower() in self.Nodebased:
						met = ('bma', meas.lower(), 'universal')
						if not met in self.measures:
							self.measures.append(met)
							self.comments[met] = "runs %s with %s concept similarity model under %s IC-approach"%(self.icf[met[0]], met[1].capitalize(), met[2].capitalize())
					elif meas.lower() in self.edf:
						met = (meas.lower(), )
						if not met in self.measures:
							self.measures.append(met)
							self.comments[met] = "runs %s edge-based concept similarity model"%(self.edf[met[0]],)
					elif meas.lower() in self.icg_1:
						met = (meas.lower(), 'universal')
						if not met in self.measures:
							self.measures.append(met)
							self.comments[met] = "runs %s"%(self.icg_1[met[0]]%(met[1].capitalize(),),)
					elif meas.lower() in self.icg_2 or meas.lower() in self.nnt:
						met = (meas.lower(), )
						cc = self.icg_2[met[0]] if met[0] in self.icg_2 else self.nnt[met[0]]
						if not met in self.measures:
							self.measures.append(met)
							self.comments[met] = "runs %s concept semantic similarity model"%(cc,)
					else:
						print(InputError('Model::%s - Value Error'%(meas,), 'This model does not exist, therefore it is not implemented.\n\nPlease refer to the tool documentation, fix this issue and try again ...\n'))
						sys.exit(7)
				else:
					print(InputError('Model::%s - Type Error'%(str(meas),), 'The model presentation format has not been recognized.\n\nPlease refer to the tool documentation, fix this issue and try again ...\n'))
					sys.exit(8)
		elif isinstance(measures, str): # A single model presented as 
			meas = measures
			if meas.lower() in self.icf:
				met = (meas.lower(), 'nunivers', 'universal')
				if not met in self.measures:
					self.measures.append(met)
					self.comments[met] = "runs %s with %s concept similarity model under %s IC-approach"%(self.icf[met[0]], met[1].capitalize(), met[2].capitalize())
			elif meas.lower() in self.AppMods:
				met = ('bma', 'nunivers', meas.lower())
				if not met in self.measures:
					self.measures.append(met)
					self.comments[met] = "runs %s with %s concept similarity model under %s IC-approach"%(self.icf[met[0]], met[1].capitalize(), met[2].capitalize())
			elif meas.lower() in self.CatLevels|self.CatPaths:
				met = ('bma', meas.lower())
				if not met in self.measures:
					self.measures.append(met)
					self.comments[met] = "runs %s with edge-based %s concept similarity model"%(self.icf[met[0]], met[1].capitalize())
			elif meas.lower() in self.Nodebased:
				met = ('bma', meas.lower(), 'universal')
				if not met in self.measures:
					self.measures.append(met)
					self.comments[met] = "runs %s with %s concept similarity model under %s IC-approach"%(self.icf[met[0]], met[1].capitalize(), met[2].capitalize())
			elif meas.lower() in self.edf:
				met = (meas.lower(), )
				if not met in self.measures:
					self.measures.append(met)
					self.comments[met] = "runs %s edge-based concept similarity model"%(self.edf[met[0]],)
			elif meas.lower() in self.icg_1:
				met = (meas.lower(), 'universal')
				if not met in self.measures:
					self.measures.append(met)
					self.comments[met] = "runs %s"%(self.icg_1[met[0]]%(met[1].capitalize(),),)
			elif meas.lower() in self.icg_2 or meas.lower() in self.nnt:
				met = (meas.lower(), )
				cc = self.icg_2[met[0]] if met[0] in self.icg_2 else self.nnt[met[0]]
				if not met in self.measures:
					self.measures.append(met)
					self.comments[met] = "runs %s concept semantic similarity model"%(cc,)
			else:
				print(InputError('Model::%s - Value Error'%(meas,), 'This model does not exist, therefore it is not implemented.\n\nPlease refer to the tool documentation, fix this issue and try again ...\n'))
				sys.exit(9)		
		else:
			print(InputError('Model::%s - Type Error'%(str(meas),), 'The model presentation format has not been recognized.\n\nPlease refer to the tool documentation, fix this issue and try again ...\n'))
			sys.exit(10)
		
	def entitySim(self, Annots, Pairs = [], measures = None, **kwargs):
		"""This function computes entity similarity scores given:
		   - Dictionary of entities (key) and concpts (values)
           - List of entites or entity pairs : Pairs
           - List of tuples of entity similarity model-IC approach,
             for groupe-wise, entity similarity model-CS model-IC
             for IC term-based!
             or list of tuples entity similarity model- if sim model does
           not require IC approach, e.g., edge-based models
           - Other measure parameter if requred
		"""
		#Checking parameters:
		if 'aaln' in kwargs and not (isinstance(kwargs['aaln'], float) or kwargs['aaln'] <= 0): 
			print(InputError('alpha (aaln) parameter - Value or Type Error', 'The value of the parameter alpha of the Al-Mubaid&Nagar edge-based SS measure, if provided\nshould be a no zeros positive float, by default it is 1.0.\n\nPlease refer to the tool manual, fix this issue and try again ...\n'))
			sys.exit(1)
		if 'c' in kwargs and not (isinstance(kwargs['c'], int) or kwargs['c'] <= 1):
			print(InputError('alpha (aaln) parameter - Value or Type Error', 'The value of the parameter alpha of the Al-Mubaid&Nagar edge-based SS measure, if provided\nshould be an integer greater than 1, by default it is 2.\n\nPlease refer to the tool manual, fix this issue and try again ...\n'))
		
		# Check other Concept Semantic Similarity parameters here
		self.parameterChecks(**kwargs)
		
		self.cleanEntityMeasures(measures) # From physical measure presented to logical measure!
		if not self.measures: # proposed measure, not known 
			print(InputError('Value Error: Measures to be retrieved cannot be found','Please, include only implemented Entity Similarity measures, \nInformation Content and Concept Similarity models as arguments. \n\nCheck ESM and CS and IC provided and try again.'))
			sys.exit(11)
			
		CurrentAnnots = self.processEntityAnnot(Annots) # From physical annotations presented to logical measure!
		if not CurrentAnnots:
			print(InputError('Value Error: Entity-ontology annotation issue','Please provide an appropriate Entity-ontology annotation mapping.\n\nPlease refer to the tool documentation, fix this issue and try again ...\n'))
			sys.exit(12)

		EntityPairs = []
		if Pairs:
			for p, q in Pairs:
				if p in CurrentAnnots and q in CurrentAnnots: EntityPairs.append((p,q))
				else:
					pass # impossibility of retrieving scores, add to an appropriate set
		else:
			tmp = list(CurrentAnnots.keys()); n = len(tmp)
			EntityPairs = [(tmp[i],tmp[j]) for i in range(n) for j in range(i,n)]
		if not EntityPairs:
			print(InputError('Value Error: Entity-Entity pair issue','Please provide an appropriate set of entity pairs for which SS scores are needed.\n\nPlease refer to the tool documentation, fix this issue and try again ...\n'))
			sys.exit(13)

		ICg = [meas for meas in self.measures if meas[0] in self.icg_1 or meas[0] in self.icg_2] #Group-wise cat
		CSf = [meas for meas in self.measures if (meas[0] in self.icf or meas[0] in self.edf)] # Pairwise cat
		NOs = [meas for meas in self.measures if meas[0] in self.nnt]
		
		Fpairs = []
		for k in EntityPairs:
			if not (k[1], k[0]) in Fpairs: Fpairs.append(k)
		del EntityPairs
		
		if ICg: self.groupwise(CurrentAnnots, Fpairs, ICg) # Run for pairwise
		if CSf: self.pairwise(CurrentAnnots, Fpairs, CSf) # Run for groupwise
		if NOs: self.ontology_indep(CurrentAnnots, Fpairs, NOs) # Run no-ontology based

	@output_str
	def __str__(self):
		"""Return a string representation of the entity ot functional-semantic similarity scores.
		For printing on the screen
        """
		meth = self.fouts[self.measures[0]].keys()
		headers = ['Entity1', 'Entity2']+['Measure-{}'.format(i+1) for i in range(len(self.measures))] 
		results = [[k[0], k[1]] + [self.fouts[s][k] for s in self.measures] for k in meth]
		return '\n'+tabs(results, headers, tablefmt = 'grid', floatfmt=".5f", stralign="center")+'\nLegends:\n==========\n'+'\n'.join(["Measure-%d: %s"%(i+1, self.comments[meas]) for (i,meas) in enumerate(self.measures)])+'\n----------\n'

	@output_str
	def __repr__(self):
		""" Return repr(self) a string representation of the entity semantic similarity
		scores. For printing in a file
		"""
		meth = self.fouts[self.measures[0]].keys()
		headers = ['Entity1', 'Entity2']+['Measure-{}'.format(i+1) for i in range(len(self.measures))] 
		results = [[k[0], k[1]] + [self.fouts[s][k] for s in self.measures] for k in meth]
		return '\n'+tabs(results, headers, tablefmt = 'plain', floatfmt=".5f", stralign="center")+'\nLegends:\n==========\n'+'\n'.join(["Measure-%d: %s"%(i+1, self.comments[meas]) for (i,meas) in enumerate(self.measures)])+'\n----------\n'

