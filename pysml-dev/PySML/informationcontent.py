#!/usr/bin/python
# -*- coding: utf8 -*-

"""Important Note:
This python file is part of the PySML tool, which is a tool for semantic
similarity between objects annotated by ontology terms or concepts, e,g.,
Gene Ontology-based functional analysis using term information content 
measures.
This python module implements term or concept information content scores. 
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

from __future__ import absolute_import, print_function, division

import time, math, os, sys, re

from . import __version__
from .error import InputError
from .imports.readontology import Ontology, output_str
from .imports.tabulate import tabulate as tabs

try:
	import networkx as nx
except ImportError:
	print(InputError('networkx library not installed', 'networkx library ImportError on line 56 in \ninformationcontent.py source code\n\nInstall the networkx library and try again ...\n'))
	sys.exit(0)

LOG = math.log

class InformationContent(object):
	"""Definition of the `InformationContent` class for information 
	   content (IC) retrieval scores.

    The InformationContent class actually behaves as a factory, creating
    a InformationContent via the default Python syntax only::
          >>> from MySML import InformationContent
          >>> 
          >>> # This is very simple and interesting

	"""
	_ontofile = os.getcwd() + "/tests/go-basic.obo"
	_namespace = 'biological_process'
	AppMods = {'universal':'GOuniversal', 'wang':'processWang', 'zhang':'processZhang', 'seco':'processSeco', 'zhou':'processZhou', 'seddiqui':'processSeddiqui', 'zanchez': 'processZanchez', 'meng': 'processMeng', 'stats': 'processStats', 'ic':'virtualIC', 'assdd': 'processSSDD'}
	__slot__ = ['DicLevels','current','approach','AppScores','TargetMissing','alt_id', 'ontodata', 'Dag', 'DagStr', 'oroot']

	def __init__(self, ontofile = '', namespace='', is_a = None, part_of = None):
		"""Instantiate a new object for a given ontology to retrieve IC.

        Arguments:
            ontology (str): the ontology file
            namespace (str): the name of the subontology as named in the ontology
            is_a (float between 0 and 1, if provided): a semantic value of the
                topological relation is_a.
            part_of (float between 0 and 1, if provided): a semantic value of the
                topological relation part_of.

        Example:
        	>>> from PySML import InformationContent
            >>> ICScores = InformationContent()
        """
        
		if is_a is None: self.is_a = 0.8
		elif isinstance(is_a, (int, float)): self.is_a = is_a
		else: # Inconsistent case to be solved: break and send type error 
			print(InputError('is_a relation semantic Value Error', 'is_a semantic value on line 46, if provided, ranges between 0 and 1\n\nPlease refer to the tool documentation, fix this value and try again ...\n'))
			sys.exit(1)
		
		if part_of is None: self.part_of = 0.6
		elif isinstance(part_of, (int, float)): self.part_of = part_of
		else: # Inconsistent case to be solved: break and send type error
			print(InputError('part_of semantic relation Value Error', 'parf_of semantic value on line 46, if provided, ranges between 0 and 1\n\nPlease refer to the tool documentation, fix this value and try again ...\n'))
			sys.exit(2)

		self.ontofile = ontofile if ontofile else self._ontofile
		self.namespace = namespace if namespace else self._namespace

		self.alt_id = {}

		self.TargetTerms = {}
		self.TargetMissing = set()
		self.approach = []
		self.AppScores = {}
		self.DicLevels = {}
		self.ontodata = None
		self.Dag = None
		self.DagStr = None
		self.oroot = None

		self.getOntoFeatures(self.ontofile)

	def getOntoFeatures(self, ontofile):
		"""Building the logical DAG of the ontology from the file
		"""
		
		ontodata = Ontology(ontofile)

		if self.namespace:
			Dag = [str(k.id) for k in ontodata if 'namespace' in ontodata[k.id].other and self.namespace in ontodata[k.id].other['namespace']]
		else:
			Dag = [str(k.id) for k in ontodata]

		DagStr = nx.DiGraph(); toprel = {'is_a': self.is_a, 'part_of':self.part_of} 
		for j in range(len(Dag)):
			child = Dag[j]; croot = False
			relations = ontodata[child].relations
			for rel in relations:
				if rel.obo_name=='is_a' or rel.obo_name=='part_of':
					DagStr.add_edges_from([(Dag.index(p), j, {'cap': toprel[rel.obo_name]}) for p in relations[rel].id], weight=-1)
					croot = True
			if not croot and not 'is_obsolete' in ontodata[child].other: oroot = j

		for k in ontodata:
			try:
				for tt in ontodata[k.id].other['alt_id']: 
					self.alt_id[str(tt)] = Dag.index(str(k.id))
			except: pass

		self.ontodata = ontodata; self.Dag = Dag
		self.DagStr = DagStr; self.oroot = oroot

	def updateDec(self, pinfo):
		"""Adjust topological scores in the context of GO universal
		   approach to avoid log 0 when computing term IC.
		"""
		prod = 1.0e+00; k = 0.0
		for i in pinfo:
			cc = i[0]; k += i[1]
			while cc < 0.1: # Making position bigger in order to avoid log 0
				cc *= 10; k += 1
			prod *= cc
		return prod, k

	def GOuniversal(self, **kwargs):
		"""Using the GO universal approach to compute the IC scores of ontology 
		   terms as a dictionary: key (term)/value (score) mapping.
        """
        
		tinfo = {self.oroot:(1.0e+00, 0)}; infoc = {self.oroot: 0.0}; Log10 = LOG(10)
		rlevel = sorted(list(set(self.DicLevels.values())), reverse=True)
		for i in rlevel[1:]:
			for j in [c for c in self.DicLevels if self.DicLevels[c]==i]:
				parents = self.ontodata[self.Dag[j]].parents.id
				pinfo = self.updateDec([(tinfo[p][0]/len(self.ontodata[self.Dag[p]].children), tinfo[p][1]) for p in [self.Dag.index(s) for s in parents]])
				tinfo[j] = pinfo
				infoc[j] = -LOG(pinfo[0]) + pinfo[1]*Log10
		return infoc

	def processWang(self, **kwargs):
		"""Using the Wang et al. node-based approach to compute the IC scores
		   of ontology terms as a dictionary: key (term)/value (score) mapping.
        """
        
		Contribution = {}
		for ter in self.DagStr:
			Dict = dict(); Dict[ter] = 1.0
			ancest = [(self.DicLevels[a], a) for a in nx.ancestors(self.DagStr, ter)]
			ancest = [x[1] for x in sorted(ancest)]; tersub = set(ancest + [ter])
			for anc in ancest:
				Dict[anc] = max(Dict[c]*self.DagStr.get_edge_data(anc, c)['cap'] for c in set(self.DagStr[anc]) & tersub)
			Contribution[ter] = Dict.copy()
		return Contribution

	def processZhang(self, **kwargs):
		"""Using the Zhang et al. node-based approach to compute the IC scores
		   of ontology terms as a dictionary: key (term)/value (score) mapping.
        """
		tinfo = dict([(t, 1) for t in self.DagStr])
		rlevel = sorted(list(set(self.DicLevels.values())))
		for i in rlevel:
			for j in [c for c in self.DicLevels if self.DicLevels[c]==i]:
				if self.DagStr[j]: tinfo[j] = sum([tinfo[k] for k in self.DagStr[j]])
		troot = tinfo[self.oroot]
		for t in tinfo:
			 tinfo[t]= -LOG(1.0*tinfo[t]/troot)
		return tinfo

	def processSeco(self, **kwargs):
		"""Using the Seco et al. node-based approach to compute the IC scores
		   of ontology terms as a dictionary: key (term)/value (score) mapping.
        """
		tinfo = dict([(t, 1) for t in self.DagStr])
		rlevel = sorted(list(set(self.DicLevels.values())))
		for i in rlevel:
			for j in [c for c in self.DicLevels if self.DicLevels[c]==i]:
				if self.DagStr[j]: tinfo[j] += sum([tinfo[k] for k in self.DagStr[j]])
		troot = tinfo[self.oroot]-1
		for t in tinfo:
			 tinfo[t]= 1.0 - LOG(tinfo[t])/LOG(troot)
		tinfo[self.oroot] = 0.0
		return tinfo

	def processZhou(self, **kwargs):
		"""Using the Zho et al. based approach to compute the IC scores of 
		   ontology terms as a dictionary: key (term)/value (score) mapping.
        """
		sigma = kwargs['sigma'] if 'sigma' in kwargs else 0.5

		sinfo = dict([(t, 1) for t in self.DagStr])
		rlevel = sorted(list(set(self.DicLevels.values())))
		for i in rlevel:
			for j in [c for c in self.DicLevels if self.DicLevels[c]==i]:
				if self.DagStr[j]: sinfo[j] += sum([sinfo[k] for k in self.DagStr[j]])
		troot = sinfo[self.oroot]-1
		for t in sinfo:
			 sinfo[t]= 1.0 - LOG(sinfo[t])/LOG(troot)
		 
		tinfo = {}; tinfo[self.oroot] = 0.0;
		sinfo.pop(self.oroot)
		for t in sinfo:
			tinfo[t]=sigma*sinfo[t]+(1-sigma)*LOG(-self.DicLevels[t])/LOG(-rlevel[0])
		return tinfo

	def processSeddiqui(self, **kwargs):
		"""Using the Sediqqui et al. based approach to compute the IC scores of
		   ontology terms as a dictionary: key (term)/value (score) mapping.
        """
		sinfo = dict([(t, 1) for t in self.DagStr])
		rlevel = sorted(list(set(self.DicLevels.values())))
		for i in rlevel:
			for j in [c for c in self.DicLevels if self.DicLevels[c]==i]:
				if self.DagStr[j]: sinfo[j] += sum([sinfo[k] for k in self.DagStr[j]])
		troot = sinfo[self.oroot]-1
		for t in sinfo:
			 sinfo[t]= 1.0 - LOG(sinfo[t])/LOG(troot)
		 
		tinfo = {}; tinfo[self.oroot] = 0.0;
		sinfo.pop(self.oroot)
		nedge = 1.0*self.DagStr.size(); nnode = 1.0*self.DagStr.order()
		sigma = LOG(nedge +1)/(LOG(nedge) + LOG(nnode))
		for t in sinfo:
			tinfo[t]=(1-sigma)*sinfo[t]+sigma*LOG(len(self.DagStr[t])+1.0)/LOG(nedge)
		return tinfo

	def processZanchez(self, **kwargs):
		"""Using the Zanchez et al. based approach to compute the IC scores of 
		   ontology terms as a dictionary: key (term)/value (score) mapping.
        """
		tinfo = dict(); leaves = set([a for a in self.DagStr if not self.DagStr[a]])
		for t in self.DagStr:
			tleaves = nx.descendants(self.DagStr, t) & leaves
			tancest = nx.ancestors(self.DagStr, t); tancest.add(t)
			tinfo[t]= -LOG((1.0*len(tleaves)/len(tancest) + 1.0)/(len(leaves)+1.0))
		return tinfo

	def processMeng(self, **kwargs):
		"""Using the Meng et al. based approach to compute the IC scores of 
		   ontology terms as a dictionary: key (term)/value (score) mapping.
        """
		rlevel = min(self.DicLevels.values())
		tinfo = dict(); tinfo[self.oroot] = 0.0; nnode = 1.0*self.DagStr.order()
		rDag = set(self.DagStr.nodes()); rDag.remove(self.oroot)
		for t in rDag:
			tdesc = nx.descendants(self.DagStr, t)
			tinfo[t] = (LOG(-self.DicLevels[t])/LOG(-rlevel))*(1.0-LOG(sum([-1.0/self.DicLevels[d] for d in tdesc])+1.0)/LOG(nnode))
		return tinfo

	def processStats(self, **kwargs):
		"""Using the statistics based approach to compute the IC scores of 
		   ontology terms as a dictionary: key (term)/value (score) mapping.
		   Here it is assumed that the user provides the statistics or 
		   counts of each terms in the ontology as a dictionary. 
        """
		TermStats = kwargs['TermStats'] if 'TermStats' in kwargs else {}
		if not TermStats or not isinstance(TermStats, dict): # Break it as TermStats is empty or not a dictionary
			print(InputError('TermStats Value or Type Error', 'TermStats type and value on line 299, if provided, should be a dictionary\nkey (term)/value (statistical count) mapping\n\nPlease refer to the tool documentation, fix this issue and try again ...\n'))
			sys.exit(3)
		ctinfo = {}
		rlevel = sorted(list(set(self.DicLevels.values())))
		
		if(isinstance(list(TermStats.values())[0],int)):
			
			tinfo = {}
			for t in self.DagStr:
				tinfo[t] = TermStats[self.Dag[t]] if self.Dag[t] in TermStats else 0
			
			if not tinfo:
				print(InputError('Ontology terms error or not found', 'TermStats dictionary: key (term)/value (statistical count or entity list) mapping\nNo term provided in TermStats is in the current version of the ontology.\n\nPlease refer to the tool documentation, fix this issue and try again ...\n'))
				sys.exit(4)
			troot = tinfo[self.oroot]
			
			for t in tinfo:
				try: ctinfo[t]= -LOG(1.0*tinfo[t]/troot)
				except: pass
		else:
			
			tinfo = {}
			for t in self.DagStr:
				tinfo[t] = set(TermStats[self.Dag[t]]) if self.Dag[t] in TermStats else set()
			
			for i in rlevel:
				for j in [c for c in self.DicLevels if self.DicLevels[c]==i]:
					for k in self.DagStr[j]:
						tinfo[j] = tinfo[j].union(tinfo[k])
			if not tinfo:
					print(InputError('Ontology terms error or not found', 'TermStats dictionary: key (term)/value (statistical count or entity list) mapping\nNo term provided in TermStats is in the current version of the ontology.\n\nPlease refer to the tool documentation, fix this issue and try again ...\n'))
					sys.exit(4)
			troot = 1.0*len(tinfo[self.oroot])
			print(troot)
			for t in tinfo:
				 try: ctinfo[t]= -LOG(1.0*len(tinfo[t])/troot)
				 except: pass
		
		return ctinfo

	def virtualIC(self, **kwargs):
		"""Assuming that the IC scores were retrieved from external platform, 
		   provided here as a dictionary: key (term)/value (score) mapping.
		   For retrieving other scores, e.g, term similarity scores or
		   entity semantic similarity scores. 
        """
		TermIC = kwargs['TermIC'] if 'TermIC' in kwargs else {}
		if not TermIC or not isinstance(TermIC, dict):
			print(InputError('TermIC Value or Type Error', 'TermStats type and value on line 320, if provided, should be a dictionary\nkey (term)/value (IC scores) mapping\nPlease refer to the tool documentation, fix this issue and try again ...\n'))
			sys.exit(5)

		TermIndexIC = {}
		for t in TermIC:
			try: TermIndexIC[self.Dag.index(t)] = TermIC[t]
			except: pass
		if not TermIndexIC:
			print(InputError('Ontology term error or not found', 'TermIC dictionary: key (term)/value (IC scores) mapping\nNo term provided in TermIC is in the current version of the ontology.\n\nPlease refer to the tool documentation, fix this issue and try again ...\n'))
			sys.exit(6)
		return TermIndexIC

	def processSSDD(self, **kwargs): # Computing once for all T-values for SSDD approach
		"""Using the SSDD approach by Xu et al. to compute the T-values of 
		   ontology terms as a dictionary: key (term)/value (score) mapping.
        """
		data = {}; tinfo = {self.oroot:1.0e+00}
		rlevel = sorted(list(set(self.DicLevels.values())), reverse=True)
		for i in rlevel[1:]:
			for j in [c for c in self.DicLevels if self.DicLevels[c]==i]:
				par = [self.Dag.index(s) for s in self.ontodata[self.Dag[j]].parents.id]
				w = 1+len(nx.descendants(self.DagStr, j))
				ss = 0.0
				for p in par:
					wp = 1+len(nx.descendants(self.DagStr, p))
					ss += w*tinfo[p]/wp
				tinfo[j] = ss/len(par)
		return tinfo

	def getIC(self, approach = None, TermList = None, **kwargs):
		self.approach = []
		if approach is None or not approach: # No approach, GOuniversal as a default
			self.approach.append('universal')
		elif isinstance(approach, str):
			approach = approach.strip().lower()
			if not approach:
				self.approach.append('universal')
			elif approach in self.AppMods: 
				self.approach.append(approach)
			else: # proposed approach, not known
				print(InputError('IC approach - Value Error', 'IC approaches provided were not found among existing approaches\n\nPlease refer to the tool documentation, fix this issue and try again ...\n'))
				sys.exit(6)
		elif isinstance(approach, (tuple, list, set)):
			approach = list(approach)
			self.approach = [(app.lower(), approach.index(app)) for app in set(approach) if app.lower() in self.AppMods]
			self.approach = [app[0] for app in sorted(self.approach, key=lambda x: x[1])]
			if not self.approach: # proposed approach, not known 
				print(InputError('IC approach - Value Error', 'IC approaches provided were not found among existing approaches\n\nPlease refer to the tool documentation, fix this issue and try again ...\n'))
				sys.exit(7)
		else: # Type Error
			print(InputError('IC approach input - Type Error', 'IC approaches should be provided as a list or tuple,\nIf one approach, then it can be provided as a string.\n\nPlease refer to the tool documentation, fix this issue and try again ...\n'))
			sys.exit(8)

		tmp = set(self.approach); tmp.discard('zanchez')
		if not self.DicLevels and tmp: # Getting terms levels for approaches requiring them
			# Get longest paths [short neg] as levels
			DicLevels = nx.bellman_ford(self.DagStr, self.oroot) 
			self.DicLevels = DicLevels[-1]
			del DicLevels

		# Parametrizing IC score retrieval: We need to check instead
		kwargs = kwargs
		if not 'sigma' in kwargs: kwargs['sigma'] = 0.5
		elif isinstance(kwargs['sigma'], float):
			if kwargs['sigma'] < 0 or kwargs['sigma'] > 1:# Impossible case, take an action
				print(InputError('sigma parameter - Value Error', 'The value of sigma for the Zho et al. approach, if provided\nshould be between 0 and 1, by default it is 0.5.\n\nPlease refer to the tool documentation, fix this issue and try again ...\n'))
				sys.exit(8)
		else: # Value error, not appropriate format: take an action
			print(InputError('sigma parameter - Type Error', 'The sigma value for the Zho et al. approach should be a float\nranging between 0 and 1.\n\nPlease refer to the tool documentation, fix this issue and try again ...\n'))
			sys.exit(8)
		
		if not 'TermStats' in kwargs: kwargs['TermStats'] = {}
		elif not isinstance(kwargs['TermStats'], dict):# Value error, not appropriate format: take an action
			print(InputError('TermStats - Type Error', 'TermStats should be a dictionary: key (term)/value (statistical count) mapping.\n\nPlease refer to the tool documentation, fix this issue and try again ...\n'))
			sys.exit(9)
		if not 'TermIC' in kwargs: kwargs['TermIC'] = {}
		elif not isinstance(kwargs['TermIC'], dict):# Value error, not appropriate format: take an action
			print(InputError('TermIC - Type Error', 'TermIC should be a dictionary: key (term)/value (IC scores) mapping.\n\nPlease refer to the tool documentation, fix this issue and try again ...\n'))
			sys.exit(10)

		# Processing IC approaches provided go here
		for app in self.approach:
			if not app in self.AppScores:
				self.AppScores[app]=getattr(self, self.AppMods[app])(**kwargs)
				
		if isinstance(TermList, str):# Probably terms in a file! Check if the path/to/file exists
			try:
				pathfile = os.path.abspath(TermList)
				# Read into the file and read concepts
				fp = open(pathfile)
				TermList = re.split("\s+|,|;|", fp.read().strip())
			except:
				try: # Possibly list of terms is provided as comma, space, semi column or colum separated string
					TermList = re.split("\s+|,|;|", TermList)
				except:
					print(InputError('Path to concept file - not found', 'When ontology concepts are provided in a file, the full path to the file needs to be provided.\n\nPlease refer to the tool documentation, fix this issue and try again ...\n'))
					sys.exit(11)
		
		if TermList: # Selected IC scores for terms in the list
			print(os.getcwd(), TermList)
			if isinstance(TermList, (list, set, tuple)):
				TermList = [s.strip() for s in TermList if s.strip()] # Leave out empty string
				self.TargetTerms = {}; self.TargetMissing = {}
				for tt in TermList:
					if tt in self.alt_id and self.alt_id[tt] in self.DagStr:
						self.TargetTerms[tt] = self.alt_id[tt]
					elif tt in self.Dag and self.Dag.index(tt) in self.DagStr:
						self.TargetTerms[tt] = self.Dag.index(tt)
					else: # Term obsolete or does not exist in the current onto!
						self.TargetMissing.add(tt)
			else: # Type error! break
				print(InputError('TermList - Type Error', 'ontology concepts should be in a list, tuple or set\n or provided in a file or comma, space, semi column or colum separated string.\n\nPlease refer to the tool documentation, fix this issue and try again ...\n'))
				sys.exit(12)
		else:
			self.TargetTerms = dict([(self.Dag[j], j) for j in self.DagStr])
			

	@output_str
	def __str__(self):
		"""Return a string reprensentation of the concept-IC score mapping.
        """
		results = []; App = self.AppScores.keys()
		headers = ["Concept"] + [app.capitalize() for app in App]
		for t in self.TargetTerms:
			st = [t]
			for app in App:
				if not self.TargetTerms[t] in self.AppScores[app]:
					st.append('Undef')
				elif isinstance(self.AppScores[app][self.TargetTerms[t]], dict):
					st.append(sum(self.AppScores[app][self.TargetTerms[t]].values()))
				elif isinstance(self.AppScores[app][self.TargetTerms[t]], (float, str)): 
					st.append(self.AppScores[app][self.TargetTerms[t]])
			results.append(st)
		if results: 
			st = '\n'+tabs(results, headers, tablefmt = 'grid', floatfmt=".5f", stralign="center")
			return st if not self.TargetMissing else st + '\n\nFollowing terms are missing in the current version of the ontology: ' + '\t'.join(self.TargetMissing)
		else: return "Please perform first getIC before attempting to print IC scores"

	@output_str
	def __repr__(self):
		""" Return repr(self)
		"""
		results = []; App = self.AppScores.keys()
		headers = ["Concept"] + [app.capitalize() for app in App]
		for t in self.TargetTerms:
			st = [t]
			for app in App:
				if not self.TargetTerms[t] in self.AppScores[app]:
					st.append('Undef')
				elif isinstance(self.AppScores[app][self.TargetTerms[t]], dict):
					st.append(sum(self.AppScores[app][self.TargetTerms[t]].values()))
				elif isinstance(self.AppScores[app][self.TargetTerms[t]], (float, str)): 
					st.append(self.AppScores[app][self.TargetTerms[t]])
			results.append(st)
		if results: 
			st = '\n'+tabs(results, headers, tablefmt = 'plain', floatfmt=".5f", stralign="center")
			return st if not self.TargetMissing else st + '\n\nFollowing terms are missing in the current version of the ontology: ' + '\t'.join(self.TargetMissing)
		else: return "Please perform first getIC before attempting to print IC scores"
