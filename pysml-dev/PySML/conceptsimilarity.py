#!/usr/bin/python
# -*- coding: utf8 -*-

"""Important Note:
This python file is part of the PySML tool, which is a tool for semantic
similarity between objects annotated by ontology terms or concepts, e,g.,
Gene Ontology-based functional analysis using term information content 
measures.
This python code implements fuzzy entity search or identification based
on semantic similarity concepts. Furthemore, this is context independent 
and can be used for any GO annotated dataset as population background or
reference and is not a context-based search.

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

import sys, os

from . import __version__
from .error import InputError
from .imports.readontology import Ontology, output_str
from .imports.tabulate import tabulate as tabs
from .informationcontent import InformationContent

try:
	import networkx as nx
except ImportError:
	print(InputError('networkx library not installed', 'networkx library ImportError on line 48 in \nconceptsimilarity.py source code\n\nInstall the networkx library and try again ...\n'))
	sys.exit(0)

class ConceptSimilarity(InformationContent):
	"""Definition of the `ConceptSimilarity` class for concept semantic 
	   similarity score retrieval.

    The ConceptSimilarity class actually behaves as a factory, creating
    a ConceptSimilarity object via the default Python syntax only::
          >>> from PySML import ConceptSimilarity
          >>> 
          >>> # This is very simple and interesting
	"""
	
	Models = {'resnik': 'sresnik', 'lin': 'slin', 'nunivers': 'snunivers', 'wang': 'swang', 'jiang':'sjiang', 'faith':'sfaith', 'ps':'sps', 'aic':'saic', 'rada':'srada', 'resnik_edge':'sresnik_edge', 'leacock':'sleacock', 'wu':'swu', 'pekar':'spekar', 'li_edge':'sli_edge', 'slimani':'sslimani', 'shenoy':'sshenoy', 'wang_edge':'swang_edge', 'zhong':'szhong', 'almubaid': 'salmubaid', 'rss':'srss', 'ssdd':'sssdd', 'shen':'sshen', 'hrss':'shrss'}
	CatLevels = set(['wu', 'pekar', 'slimani', 'shenoy', 'wang_edge', 'zhong', 'almubaid', 'rss', 'ssdd'])
	CatPaths = set(['rada', 'resnik_edge', 'leacock', 'li_edge', 'shenoy', 'almubaid', 'rss'])
	Nodebased = set(['resnik', 'lin', 'nunivers', 'wang', 'jiang', 'faith', 'aic', 'hrss', 'ps', 'shen'])
	
	__cslot__ = ['ShortPath', 'deep', 'outputs', 'TargetPairs', 'models']
	
	def __init__(self, ontofile = '', namespace='', is_a = None, part_of = None):
		"""Instantiate a new object for a given ontology to retrieve term simila-
		   rity scores with a DAG inherited from InformationContent class.

        Arguments:
            ontology (str): the ontology file
            namespace (str): the name of the subontology as named in the ontology
            is_a (float between 0 and 1, if provided): a semantic value of the
                topological relation is_a.
            part_of (float between 0 and 1, if provided): a semantic value of the
                topological relation part_of.

        Example:
        	>>> from PySML import ConceptSimilarity
            >>> CScores = ConceptSimilarity()
        """
		InformationContent.__init__(self, ontofile, namespace, is_a, part_of)
		self.ShortPath = {}
		self.deep = None
		self.outputs = {}
		self.TargetPairs = {}
		self.keepterms = None
		self.models = []
		self.interface = {}

	def cfact(self, icanc, cf = 0): # Produce correction factor
		"""Given set of common ancestors and correction index to be applied,
		   the correction score is output.
		      cf = 0 when no correction is applied
		         = 1 when graph-based correction
		         = 2 for Relevance correction
		         = 3 for SimIC correction
		         
		"""
		if cf==0: return 1.0								
		if cf==1: return sum(icanc)/(len(icanc)*max(icanc))	
		elif cf==2: return 1.0 - exp(-max(icanc))			
		elif cf==3: return 1.0 - 1.0/(1.0 + max(icanc))		

	def sresnik(self, pairids, **kwargs):
		"""Using the Resnik-based approach to compute the ontology concept similarity scores
		
		Arguments:
		   app    : The IC approach to be used and if not provided, then 'universal' is used
		   cf     : The correction score and if not provided, then it set to 0, 
		   gr     : Used when graph-based correction is used, set to 0 by default if all
		            informative commun ancestors (ICA) model is applied and to 1 if 
		            exclusively inherited shared common ancestors (EICA) is applied.
		   sigma  : This is the Zho et al. parameter, set to 0.5 by default, otherwise it is
		            a float value ranging between 0 and 1.
		   pairids: Concept ID or concept ID pair list for which semantic similarity scores
		            are computed.
		   
		"""
		kwargs = kwargs
		app = kwargs['app'] if 'app' in kwargs else 'universal' 
		cf = kwargs['cf'] if 'cf' in kwargs else 0
		gr = kwargs['gr'] if 'gr' in kwargs else 0
		# Setting potential parameters for information content
		if not 'sigma' in kwargs: kwargs['sigma'] = 0.5
		if not 'TermStats' in kwargs: kwargs['TermStats'] = {}
		if not 'TermIC' in kwargs: kwargs['TermIC'] = {}
		
		if not app in self.AppScores: self.getIC(approach = app, **kwargs)
		if app != 'wang': ics = self.AppScores[app]
		else: ics = dict([(t, sum(self.AppScores[app][t].values())) for t in self.AppScores[app]])
		MaxValue = max(ics.values())
		data = {}
		for p, q in pairids:
			if p==q: data[(p, q)] = 1.0
			else: 
				panc = nx.ancestors(self.DagStr, p); panc.add(p)
				qanc = nx.ancestors(self.DagStr, q); qanc.add(q)
				canc = panc & qanc; canc.discard(self.oroot)
				if gr: canc = [a for a in canc if set(self.DagStr[a])&((panc|qanc)-canc)]
				icanc = [ics[t] for t in canc]
				vic = self.cfact(icanc, cf)*abs(max(icanc)) if icanc else 0.0
				data[(p,q)] = vic/MaxValue if vic else 0.0
		return data

	def slin(self, pairids, **kwargs):
		"""Using the Lin-based approach to compute the ontology concept similarity scores
		
		Arguments: As in the sresnik method above
		
		"""
		kwargs = kwargs
		app = kwargs['app'] if 'app' in kwargs else 'universal' 
		cf = kwargs['cf'] if 'cf' in kwargs else 0
		gr = kwargs['gr'] if 'gr' in kwargs else 0
		# Setting potential parameters for information content
		if not 'sigma' in kwargs: kwargs['sigma'] = 0.5
		if not 'TermStats' in kwargs: kwargs['TermStats'] = {}
		if not 'TermIC' in kwargs: kwargs['TermIC'] = {}
		if not app in self.AppScores: self.getIC(approach = app, **kwargs)
		if app != 'wang': ics = self.AppScores[app]
		else: ics = dict([(t, sum(self.AppScores[app][t].values())) for t in self.AppScores[app]])
		data = {}
		for p, q in pairids:
			if p==q: data[(p, q)] = 1.0
			else: 
				panc = nx.ancestors(self.DagStr, p); panc.add(p)
				qanc = nx.ancestors(self.DagStr, q); qanc.add(q)
				canc = panc & qanc; canc.discard(self.oroot)
				if gr: canc = [a for a in canc if set(self.DagStr[a])&((panc|qanc)-canc)]
				icanc = [ics[t] for t in canc]
				vic = self.cfact(icanc, cf)*abs(max(icanc)) if icanc else 0.0
				data[(p,q)] = 2*vic/(ics[p]+ics[q]) if vic else 0.0
		return data

	def sjiang(self, pairids, **kwargs): 
		"""Using any variant of J&C-based approach to compute the ontology concept 
		similarity scores
		
		Arguments: As in the sresnik method above and addition jv described above
		   jv = 0 for the Resnik-based normalization
		      = 1 for the Couto-based normalization
		      = 2 for the Leacock & Chodorow normalization
		      = 3 for the Garla & Brandt normalization
		      = 4 for the Rada normalization
		      = 5 for the canonical normalization
		
		"""
		# Setting potential parameters for concept similarity: approach and correct fact
		kwargs = kwargs
		app = kwargs['app'] if 'app' in kwargs else 'universal' 
		cf = kwargs['cf'] if 'cf' in kwargs else 0
		gr = kwargs['gr'] if 'gr' in kwargs else 0
		jv = kwargs['jv'] if 'jv' in kwargs else 0
		# Setting potential parameters for information content
		if not 'sigma' in kwargs: kwargs['sigma'] = 0.5
		if not 'TermStats' in kwargs: kwargs['TermStats'] = {}
		if not 'TermIC' in kwargs: kwargs['TermIC'] = {}
		if not app in self.AppScores: self.getIC(approach = app, **kwargs)
		if app != 'wang': ics = self.AppScores[app]
		else: ics = dict([(t, sum(self.AppScores[app][t].values())) for t in self.AppScores[app]])
		data = {}; icmax = max(ics.values())
		if jv==0: # Resnik-based normalization
			for p, q in pairids:
				if p==q: data[(p, q)] = 1.0
				else:
					panc = nx.ancestors(self.DagStr, p); panc.add(p)
					qanc = nx.ancestors(self.DagStr, q); qanc.add(q)
					canc = panc & qanc; canc.discard(self.oroot)
					if gr: canc = [a for a in canc if set(self.DagStr[a])&((panc|qanc)-canc)]
					icanc = [ics[t] for t in canc]
					vic = self.cfact(icanc, cf) if icanc else 0.0
					jcnn = ics[p]+ics[q]-2*abs(max(icanc)) if icanc else 1.0
					data[(p,q)] = vic*(1.0 - jcnn/(2*icmax)) if jcnn != 1.0 else 0.0
		elif jv==1: # Couto-based normalization
			for p, q in pairids:
				if p==q: data[(p, q)] = 1.0
				else:
					panc = nx.ancestors(self.DagStr, p); panc.add(p)
					qanc = nx.ancestors(self.DagStr, q); qanc.add(q)
					canc = panc & qanc; canc.discard(self.oroot)
					if gr: canc = [a for a in canc if set(self.DagStr[a])&((panc|qanc)-canc)]
					icanc = [ics[t] for t in canc]
					vic = self.cfact(icanc, cf) if icanc else 0.0
					jcnn = ics[p]+ics[q]-2*abs(max(icanc)) if icanc else 1.0
					data[(p,q)] = vic*(1.0 - min(1,jcnn/icmax)) if jcnn != 1.0 else 0.0
		elif jv==2: # Leacock & Chodorow normalization
			for p, q in pairids:
				if p==q: data[(p, q)] = 1.0
				else:
					panc = nx.ancestors(self.DagStr, p); panc.add(p)
					qanc = nx.ancestors(self.DagStr, q); qanc.add(q)
					canc = panc & qanc; canc.discard(self.oroot)
					if gr: canc = [a for a in canc if set(self.DagStr[a])&((panc|qanc)-canc)]
					icanc = [ics[t] for t in canc]
					vic = self.cfact(icanc, cf) if icanc else 0.0
					jcnn = ics[p]+ics[q]-2*abs(max(icanc)) if icanc else 1.0
					data[(p,q)] = vic*(1.0 - log(jcnn+1)/log(icmax))
		elif jv==3: # Garla & Brandt normalization
			for p, q in pairids:
				if p==q: data[(p, q)] = 1.0
				else:
					panc = nx.ancestors(self.DagStr, p); panc.add(p)
					qanc = nx.ancestors(self.DagStr, q); qanc.add(q)
					canc = panc & qanc; canc.discard(self.oroot)
					if gr: canc = [a for a in canc if set(self.DagStr[a])&((panc|qanc)-canc)]
					icanc = [ics[t] for t in canc]
					vic = self.cfact(icanc, cf) if icanc else 0.0
					jcnn = ics[p]+ics[q]-2*abs(max(icanc)) if icanc else 1.0
					data[(p,q)] = vic*(1.0 - log(jcnn+1)/log(icmax+1))
		elif jv==4: # Rada normalization
			for p, q in pairids:
				if p==q: data[(p, q)] = 1.0
				else:
					panc = nx.ancestors(self.DagStr, p); panc.add(p)
					qanc = nx.ancestors(self.DagStr, q); qanc.add(q)
					canc = panc & qanc; canc.discard(self.oroot)
					if gr: canc = [a for a in canc if set(self.DagStr[a])&((panc|qanc)-canc)]
					icanc = [ics[t] for t in canc]
					vic = self.cfact(icanc, cf) if icanc else 0.0
					jcnn = ics[p]+ics[q]-2*abs(max(icanc)) if icanc else 1.0
					data[(p,q)] = vic/(1.0 + jcnn)
		else: # canonical normalization
			for p, q in pairids:
				if p==q: data[(p, q)] = 1.0
				else:
					panc = nx.ancestors(self.DagStr, p); panc.add(p)
					qanc = nx.ancestors(self.DagStr, q); qanc.add(q)
					canc = panc & qanc; canc.discard(self.oroot)
					if gr: canc = [a for a in canc if set(self.DagStr[a])&((panc|qanc)-canc)]
					icanc = [ics[t] for t in canc]
					vic = self.cfact(icanc, cf) if icanc else 0.0
					jcnn = ics[p]+ics[q]-2*abs(max(icanc)) if icanc else 1.0
					data[(p,q)] = vic*(1.0 - log(jcnn+1)/log(ics[p]+ics[q]+1))
		return data

	def snunivers(self, pairids, **kwargs):
		"""Using the universal-based normalization approach to compute the ontology 
		concept similarity scores
		
		Arguments: As in the sresnik method above
		
		"""
		# Setting potential parameters for concept similarity: approach and correct fact
		kwargs = kwargs
		app = kwargs['app'] if 'app' in kwargs else 'universal' 
		cf = kwargs['cf'] if 'cf' in kwargs else 0
		gr = kwargs['gr'] if 'gr' in kwargs else 0
		# Setting potential parameters for information content
		if not 'sigma' in kwargs: kwargs['sigma'] = 0.5
		if not 'TermStats' in kwargs: kwargs['TermStats'] = {}
		if not 'TermIC' in kwargs: kwargs['TermIC'] = {}
		if not app in self.AppScores: self.getIC(approach = app, **kwargs)
		if app != 'wang': ics = self.AppScores[app]
		else: ics = dict([(t, sum(self.AppScores[app][t].values())) for t in self.AppScores[app]])
		data = {}
		for p, q in pairids:
			if p==q: data[(p, q)] = 1.0
			else: 
				panc = nx.ancestors(self.DagStr, p); panc.add(p)
				qanc = nx.ancestors(self.DagStr, q); qanc.add(q)
				canc = panc & qanc; canc.discard(self.oroot)
				if gr: canc = [a for a in canc if set(self.DagStr[a])&((panc|qanc)-canc)]
				icanc = [ics[t] for t in canc]
				vic = self.cfact(icanc, cf)*abs(max(icanc)) if icanc else 0.0
				data[(p,q)] = vic/max(ics[p], ics[q]) if vic else 0.0
		return data

	def swang(self, pairids, **kwargs):
		"""Using the Wang-based approach to compute the ontology concept similarity scores
		
		Arguments: As in the sresnik method above
		
		"""
		# Setting potential parameters for concept similarity: approach and correct fact
		kwargs = kwargs
		app = 'wang' 
		cf = kwargs['cf'] if 'cf' in kwargs else 0
		gr = kwargs['gr'] if 'gr' in kwargs else 0
		# Setting potential parameters for information content
		if not 'sigma' in kwargs: kwargs['sigma'] = 0.5
		if not 'TermStats' in kwargs: kwargs['TermStats'] = {}
		if not 'TermIC' in kwargs: kwargs['TermIC'] = {}
		if not app in self.AppScores: self.getIC(approach = app, **kwargs)
		ics = self.AppScores[app]
		data = {}
		for p, q in pairids:
			if p==q: data[(p, q)] = 1.0
			else:
				panc = nx.ancestors(self.DagStr, p); panc.add(p)
				qanc = nx.ancestors(self.DagStr, q); qanc.add(q)
				canc = panc & qanc
				if gr: canc = [a for a in canc if set(self.DagStr[a])&((panc|qanc)-canc)]
				icanc = [sum(ics[t].values()) for t in canc]
				cic = self.cfact(icanc, cf) if icanc else 0.0
				pqanc = set(ics[p]) & set(ics[q])
				data[(p,q)] = cic*abs(sum([ics[p][a]+ics[q][a] for a in pqanc])/(sum(ics[p].values())+sum(ics[q].values())))
		return data

	def sfaith(self, pairids, **kwargs):
		"""Using the FaiTH-based approach to compute the ontology concept similarity scores
		
		Arguments: As in the sresnik method above
		
		"""
		# Setting potential parameters for concept similarity: approach and correct fact
		kwargs = kwargs
		app = kwargs['app'] if 'app' in kwargs else 'universal' 
		cf = kwargs['cf'] if 'cf' in kwargs else 0
		gr = kwargs['gr'] if 'gr' in kwargs else 0
		# Setting potential parameters for information content
		if not 'sigma' in kwargs: kwargs['sigma'] = 0.5
		if not 'TermStats' in kwargs: kwargs['TermStats'] = {}
		if not 'TermIC' in kwargs: kwargs['TermIC'] = {}
		if not app in self.AppScores: self.getIC(approach = app, **kwargs)
		if app != 'wang': ics = self.AppScores[app]
		else: ics = dict([(t, sum(self.AppScores[app][t].values())) for t in self.AppScores[app]])
		data = {}
		for p, q in pairids:
			if p==q: data[(p, q)] = 1.0
			else: 
				panc = nx.ancestors(self.DagStr, p); panc.add(p)
				qanc = nx.ancestors(self.DagStr, q); qanc.add(q)
				canc = panc & qanc; canc.discard(self.oroot)
				icanc = [ics[t] for t in canc]
				if gr: canc = [a for a in canc if set(self.DagStr[a])&((panc|qanc)-canc)]
				vic = self.cfact(icanc, cf)*abs(max(icanc)) if icanc else 0.0
				data[(p,q)] = vic/(ics[p]+ics[q]-abs(max(icanc))) if vic else 0.0
		return data

	def sps(self, pairids, **kwargs):
		"""Using the P&S-based approach to compute the ontology concept similarity scores
		
		Arguments: As in the sresnik method above
		
		"""
		# Setting potential parameters for concept similarity: approach and correct fact
		kwargs = kwargs
		app = kwargs['app'] if 'app' in kwargs else 'universal' 
		cf = kwargs['cf'] if 'cf' in kwargs else 0
		gr = kwargs['gr'] if 'gr' in kwargs else 0
		# Setting potential parameters for information content
		if not 'sigma' in kwargs: kwargs['sigma'] = 0.5
		if not 'TermStats' in kwargs: kwargs['TermStats'] = {}
		if not 'TermIC' in kwargs: kwargs['TermIC'] = {}
		if not app in self.AppScores: self.getIC(approach = app, **kwargs)
		if app != 'wang': ics = self.AppScores[app]
		else: ics = dict([(t, sum(self.AppScores[app][t].values())) for t in self.AppScores[app]])
		data = {}
		for p, q in pairids:
			if p==q: data[(p, q)] = 1.0
			else: 
				panc = nx.ancestors(self.DagStr, p); panc.add(p)
				qanc = nx.ancestors(self.DagStr, q); qanc.add(q)
				canc = panc & qanc; canc.discard(self.oroot)
				if gr: canc = [a for a in canc if set(self.DagStr[a])&((panc|qanc)-canc)]
				icanc = [ics[t] for t in canc]
				vic = self.cfact(icanc, cf) if icanc else 0.0
				data[(p,q)] = vic*max(0,3*abs(max(icanc))-ics[p]-ics[q]) if vic else 0.0
		return data

	def saic(self, pairids, **kwargs):
		"""Using the Agregate Information content-based approach by Song et al.
		to compute the ontology concept similarity scores
		
		Arguments: As in the sresnik method above
		
		"""
		# Setting potential parameters for concept similarity: approach and correct fact
		kwargs = kwargs
		app = kwargs['app'] if 'app' in kwargs else 'universal' 
		cf = kwargs['cf'] if 'cf' in kwargs else 0
		gr = kwargs['gr'] if 'gr' in kwargs else 0
		# Setting potential parameters for information content
		if not 'sigma' in kwargs: kwargs['sigma'] = 0.5
		if not 'TermStats' in kwargs: kwargs['TermStats'] = {}
		if not 'TermIC' in kwargs: kwargs['TermIC'] = {}

		if not app in self.AppScores: self.getIC(approach = app, **kwargs)
		if app != 'wang': ics = self.AppScores[app]
		else: ics = dict([(t, sum(self.AppScores[app][t].values())) for t in self.AppScores[app]])
		tics = {}
		for t in ics:
			panc = nx.ancestors(self.DagStr, t); panc.add(t)
			tics[t] = sum([1.0/(1.0+exp(-1.0/ics[t])) if ics[t] else 1.0 for t in panc])
		data = {}
		for p, q in pairids:
			if p==q: data[(p, q)] = 1.0
			else: 
				panc = nx.ancestors(self.DagStr, p); panc.add(p)
				qanc = nx.ancestors(self.DagStr, q); qanc.add(q)
				canc = panc & qanc; canc.discard(self.oroot)
				if gr: canc = [a for a in canc if set(self.DagStr[a])&((panc|qanc)-canc)]
				icanc = [tics[t] for t in canc]
				vic = self.cfact(icanc, cf) if icanc else 0.0
				data[(p,q)] = vic*2.0*sum([1.0/(1.0+exp(-1.0/ics[t])) if ics[t] else 1.0 for t in canc])/(tics[p] + tics[q])
		return data

	def srada(self, pairids, **kwargs):
		"""Using the Rada-based approach to compute the ontology concept similarity scores
		
		Arguments: As in the sresnik method above
		
		"""
		data = {}
		for p, q in pairids:
			if p==q: data[(p, q)] = 1.0
			else: 
				panc = nx.ancestors(self.DagStr, p); panc.add(p)
				qanc = nx.ancestors(self.DagStr, q); qanc.add(q)
				canc = panc & qanc
				spl = min([self.ShortPath[a][p] for a in canc])
				spl += min([self.ShortPath[a][q] for a in canc])
				data[(p,q)] = 1.0/(1.0 + spl)
		return data
		
	def sresnik_edge(self, pairids, **kwargs):
		"""Using the Edge-based Resnik approach to compute the ontology concept similarity scores
		
		Arguments: As in the sresnik method above
		
		"""
		data = {}
		for p, q in pairids:
			if p==q: data[(p, q)] = 1.0
			else:
				panc = nx.ancestors(self.DagStr, p); panc.add(p)
				qanc = nx.ancestors(self.DagStr, q); qanc.add(q)
				canc = panc & qanc
				spl = min([self.ShortPath[a][p] for a in canc])
				spl += min([self.ShortPath[a][q] for a in canc])
				data[(p,q)] = 1.0 - spl/(2.0*self.deep)
		return data
	
	def sleacock(self, pairids, **kwargs):
		"""Using the Leacock & Chodorow approach to compute the ontology concept similarity scores
		
		Arguments: As in the sresnik method above
		
		"""
		data = {}
		for p, q in pairids:
			if p==q: data[(p, q)] = 1.0
			else:
				panc = nx.ancestors(self.DagStr, p); panc.add(p)
				qanc = nx.ancestors(self.DagStr, q); qanc.add(q)
				canc = panc & qanc
				spl = min([self.ShortPath[a][p] for a in canc])
				spl += min([self.ShortPath[a][q] for a in canc])
				data[(p,q)] = 1.0 - log(spl)/log(2.0*self.deep)
		return data

	def swu(self, pairids, **kwargs):
		"""Using the Wu and Palmer approach to compute the ontology concept similarity scores
		
		Arguments: As in the sresnik method above
		
		"""
		data = {}
		for p, q in pairids:
			if p==q: data[(p, q)] = 1.0
			else:
				if (q,p) in data:
					data[(p,q)] = data[(q,p)]
					continue
				panc = nx.ancestors(self.DagStr, p); panc.add(p)
				qanc = nx.ancestors(self.DagStr, q); qanc.add(q)
				canc = panc & qanc
				spl = min([self.DicLevels[a] for a in canc])
				data[(p,q)] = 2*spl/(self.DicLevels[p] + self.DicLevels[q])
		return data
	
	def spekar(self, pairids, **kwargs):
		"""Using the Pekar & Staab approach to compute the ontology concept similarity scores
		
		Arguments: As in the sresnik method above
		
		"""
		data = {}
		for p, q in pairids:
			if p==q: data[(p, q)] = 1.0
			else:
				panc = nx.ancestors(self.DagStr, p); panc.add(p)
				qanc = nx.ancestors(self.DagStr, q); qanc.add(q)
				canc = panc & qanc
				spl = min([self.DicLevels[a] for a in canc])
				data[(p,q)] = spl/(self.DicLevels[p] + self.DicLevels[q]-spl)
		return data
		
	def sli_edge(self, pairids, **kwargs):
		"""Using the Li et al. edge-based approach to compute the ontology concept similarity scores
		
		Arguments: As in the sresnik method above and following additional parameters:
		    alpha: A float, which is greater or equal to 0, set to 0.2 as default if not provided
		    beta : A float, which is stricty greater to 0, set to 0.6 as default if not provided
		
		"""
		alpha = kwargs['alpha'] if 'alpha' in kwargs else 0.2
		beta = kwargs['beta'] if 'beta' in kwargs else 0.6
		data = {}
		for p, q in pairids:
			if p==q: data[(p, q)] = 1.0
			else:
				panc = nx.ancestors(self.DagStr, p); panc.add(p)
				qanc = nx.ancestors(self.DagStr, q); qanc.add(q)
				canc = panc & qanc
				spl = min([self.ShortPath[a][p] for a in canc])
				spl += min([self.ShortPath[a][q] for a in canc])
				deepc = -min([self.DicLevels[a] for a in canc])
				data[(p,q)] = exp(-alpha*spl)*tanh(beta*deepc)
		return data

	def sslimani(self, pairids, **kwargs):
		"""Using the Slimani et al. edge-based approach to compute the ontology concept similarity scores	
		Arguments: As in the sresnik method above
		
		"""
		data = {}
		for p, q in pairids:
			if p==q: data[(p, q)] = 1.0
			else:
				panc = nx.ancestors(self.DagStr, p); panc.add(p)
				qanc = nx.ancestors(self.DagStr, q); qanc.add(q)
				canc = panc & qanc
				dpl = -min([self.DicLevels[a] for a in canc])
				pp = set([self.Dag.index(s) for s in self.ontodata[self.Dag[p]].parents.id]) 
				pq = set([self.Dag.index(s) for s in self.ontodata[self.Dag[p]].parents.id]) 
				lda = 1.0 if set(pp) & set(pq) else 0.0
				cf = (1-lda)/(min(-self.DicLevels[p], -self.DicLevels[q])-dpl+1)+lda/(abs(-self.DicLevels[p]+self.DicLevels[q])+1)
				data[(p,q)] = 2*cf*dpl/(-self.DicLevels[p] - self.DicLevels[q])
		return data

	def sshenoy(self, pairids, **kwargs):
		"""Using the Shenoy et al.-based approach to compute the ontology concept similarity scores	
		Arguments: As in the sresnik method above
		
		"""
		data = {}
		for p, q in pairids:
			if p==q: data[(p, q)] = 1.0
			else:
				panc = nx.ancestors(self.DagStr, p); panc.add(p)
				qanc = nx.ancestors(self.DagStr, q); qanc.add(q)
				canc = panc & qanc
				dpl = -min([self.DicLevels[a] for a in canc])
				pp = set([self.Dag.index(s) for s in self.ontodata[self.Dag[p]].parents.id]) 
				pq = set([self.Dag.index(s) for s in self.ontodata[self.Dag[p]].parents.id]) 
				lda = 1.0 if set(pp) & set(pq) else 0.0
				spl = min([self.ShortPath[a][p] for a in canc])
				spl += min([self.ShortPath[a][q] for a in canc])
				cf = exp(-lda*spl/self.deep)
				data[(p,q)] = 2*cf*dpl/(-self.DicLevels[p] - self.DicLevels[q])
		return data

	def swang_edge(self, pairids, **kwargs):
		"""Using the Wang et al. edge-based approach to compute the ontology concept similarity scores	
		Arguments: As in the sresnik method above
		
		"""
		data = {}
		for p, q in pairids:
			if p==q: data[(p, q)] = 1.0
			else:
				panc = nx.ancestors(self.DagStr, p); panc.add(p)
				qanc = nx.ancestors(self.DagStr, q); qanc.add(q)
				canc = panc & qanc
				dpl = min([self.DicLevels[a] for a in canc])
				slca = [a for a in canc if self.DicLevels[a]==dpl]
				ss = 0.0
				for c in slca:
					numer = [len(pth) for pth in nx.all_simple_paths(self.DagStr, self.oroot, c)]
					deno1 = [len(pth) for pth in nx.all_simple_paths(self.DagStr, self.oroot, p) if c in pth]
					deno2 = [len(pth) for pth in nx.all_simple_paths(self.DagStr, self.oroot, q) if c in pth]
					ss += (sum(numer)/len(numer))**2/((sum(deno1)/len(deno1))*(sum(deno2)/len(deno2)))
				data[(p,q)] = ss/len(slca)
		return data

	def szhong(self, pairids, **kwargs):
		"""Using the Zhong-based approach to compute the ontology concept similarity scores
		
		Arguments: As in the sresnik method above and following additional parameters:
		       zk: k parameter, a float greater than to 1, set to 2 by default
		
		"""
		kz = kwargs['zk'] if 'zk' in kwargs else 2
		data = {}
		for p, q in pairids:
			if p==q: data[(p, q)] = 1.0
			else:
				panc = nx.ancestors(self.DagStr, p); panc.add(p)
				qanc = nx.ancestors(self.DagStr, q); qanc.add(q)
				canc = panc & qanc
				dpl = -min([self.DicLevels[a] for a in canc])
				dp, dq = -self.DicLevels[p], -self.DicLevels[q]
				data[(p,q)] = 1.0 - (1.0/kz**dpl - 0.5/kz**dp - 0.5/kz**dq)
		return data

	def salmubaid(self, pairids, **kwargs):
		"""Using the Almubaid-based approach to compute the ontology concept similarity scores
		
		Arguments: As in the sresnik method above and following additional parameters:
		    ak: k parameter, a float great or equal to 1, set to 1 by default
		    aa: alpha parameter, a float which is strictly greater to 0, set to 1 by default
		    ab: beta parameter, a float which is stricty greater to 0, set to 1 by default
		
		"""
		ak = kwargs['ak'] if 'ak' in kwargs else 1
		aa = kwargs['aa'] if 'aa' in kwargs else 1
		ab = kwargs['ab'] if 'ab' in kwargs else 1
		data = {}
		for p, q in pairids:
			if p==q: data[(p, q)] = 1.0
			else:
				panc = nx.ancestors(self.DagStr, p); panc.add(p)
				qanc = nx.ancestors(self.DagStr, q); qanc.add(q)
				canc = panc & qanc
				dpl = -min([self.DicLevels[a] for a in canc])
				spl = min([self.ShortPath[a][p] for a in canc])
				spl += min([self.ShortPath[a][q] for a in canc])
				ds = log(ak + (spl - 1)**aa*(self.deep - dpl)**ab)
				data[(p,q)] = 1.0/(1.0 + ds) # Normalized using Rada model
		return data

	def srss(self, pairids, **kwargs):
		"""Using the relative specificity similarity-based approach by Wu et al.
		to compute the ontology concept similarity scores
		
		Arguments: As in the sresnik method above
		
		"""
		data = {}; leaves = set([a for a in self.DagStr if not self.DagStr[a]])
		for p, q in pairids:
			if p==q: data[(p, q)] = 1.0
			else:
				panc = nx.ancestors(self.DagStr, p); panc.add(p)
				qanc = nx.ancestors(self.DagStr, q); qanc.add(q)
				canc = panc & qanc
				dpl = -min([self.DicLevels[a] for a in canc])
				spl = min([self.ShortPath[a][p] for a in canc])
				spl += min([self.ShortPath[a][q] for a in canc])
				pleaves = nx.descendants(self.DagStr, p) & leaves
				qleaves = nx.descendants(self.DagStr, q) & leaves
				bp = min([self.ShortPath[p][a] for a in pleaves]) if pleaves else 0.01
				bq = min([self.ShortPath[q][a] for a in qleaves]) if qleaves else 0.01
				data[(p,q)] = (self.deep/(self.deep+spl))*dpl/(dpl+min(bp,bq))
		return data

	def sssdd(self, pairids, **kwargs):
		"""Using the shortest semantic differentiation distance (SSDD) approach by Wu et al.
		to compute the ontology concept similarity scores
		
		Arguments: As in the sresnik method above
		
		"""
		data = {}
		if not 'assdd' in self.AppScores: self.getIC(approach='assdd')
		tinfo = self.AppScores['assdd']
		for p, q in pairids:
			if p==q: data[(p, q)] = 1.0
			else:
				panc = nx.ancestors(self.DagStr, p); panc.add(p)
				qanc = nx.ancestors(self.DagStr, q); qanc.add(q)
				canc = panc & qanc
				dpl = min([self.DicLevels[a] for a in canc])
				slca = [a for a in canc if self.DicLevels[a]==dpl]
				pPath = []; qPath = []
				for c in slca:
					for pth in nx.all_shortest_paths(self.DagStr, c, p):
						pPath.append(pth)
					for pth in nx.all_shortest_paths(self.DagStr, c, q):
						qPath.append(pth)
				Path = [set(pp+qq) for pp in pPath for qq in qPath]
				data[p,q] = 1 - atan(min([sum([tinfo[c] for c in pth]) for pth in Path]))/(PI/2)
		return data 
				
	def sshen(self, pairids, **kwargs):
		"""Using the Shen et al.-based approach to compute the ontology concept similarity scores
		
		Arguments: As in the sresnik method above
		
		"""
		# Setting potential parameters for concept similarity: approach and correct fact
		kwargs = kwargs
		app = kwargs['app'] if 'app' in kwargs else 'universal' 
		cf = kwargs['cf'] if 'cf' in kwargs else 0
		gr = kwargs['gr'] if 'gr' in kwargs else 0
		# Setting potential parameters for information content
		if not 'sigma' in kwargs: kwargs['sigma'] = 0.5
		if not 'TermStats' in kwargs: kwargs['TermStats'] = {}
		if not 'TermIC' in kwargs: kwargs['TermIC'] = {}
		data = {}
		if not app in self.AppScores: self.getIC(approach = app, **kwargs)
		if app != 'wang': tinfo = self.AppScores[app]
		else: tinfo = dict([(t, sum(self.AppScores[app][t].values())) for t in self.AppScores[app]])
		for p, q in pairids:
			if p==q: data[(p, q)] = 1.0
			else:
				panc = nx.ancestors(self.DagStr, p); panc.add(p)
				qanc = nx.ancestors(self.DagStr, q); qanc.add(q)
				canc = panc & qanc; c = mica = -1
				for a in canc:
					if tinfo[a] > mica: mica, c = tinfo[a], a 
				pPath = list(nx.all_shortest_paths(self.DagStr, c, p))
				qPath = list(nx.all_shortest_paths(self.DagStr, c, q))
				Path = [set(pp+qq) for pp in pPath for qq in qPath]
				data[p,q] = 1 - atan(min([sum([1/tinfo[c] if tinfo[c] > 0 else 0.0 for c in pth]) for pth in Path]))/(PI/2)
		return data 

	def shrss(self, pairids, **kwargs):
		"""Using the hybrid relative specificity similarity (HRSS)-based approach by Wu et al.
		to compute the ontology concept similarity scores
		
		Arguments: As in the sresnik method above
		
		"""
		# Setting potential parameters for concept similarity: approach and correct fact
		kwargs = kwargs
		app = kwargs['app'] if 'app' in kwargs else 'universal' 
		cf = kwargs['cf'] if 'cf' in kwargs else 0
		gr = kwargs['gr'] if 'gr' in kwargs else 0
		# Setting potential parameters for information content
		if not 'sigma' in kwargs: kwargs['sigma'] = 0.5
		if not 'TermStats' in kwargs: kwargs['TermStats'] = {}
		if not 'TermIC' in kwargs: kwargs['TermIC'] = {}
		data = {}; leaves = set([a for a in self.DagStr if not self.DagStr[a]])
		if not app in self.AppScores: self.getIC(approach = app, **kwargs)
		if app != 'wang': tinfo = self.AppScores[app]
		else: tinfo = dict([(t, sum(self.AppScores[app][t].values())) for t in self.AppScores[app]])
		MaxValue = max(tinfo.values())
		tinfo = dict([(t, tinfo[t]/MaxValue) for t in tinfo])
		for p, q in pairids:
			if p==q: data[(p, q)] = 1.0
			else:
				panc = nx.ancestors(self.DagStr, p); panc.add(p)
				qanc = nx.ancestors(self.DagStr, q); qanc.add(q)
				canc = panc & qanc
				mica = max([tinfo[c] for c in canc])
				dic = abs(tinfo[p]-mica) + abs(tinfo[q]-mica)
				bp = max([tinfo[d] for d in nx.descendants(self.DagStr, p) & leaves])
				bq = max([tinfo[d] for d in nx.descendants(self.DagStr, q) & leaves])
				bic = (abs(bp - tinfo[p]) + abs(bq - tinfo[q]))/2
				data[(p,q)] = (1/(1+dic))*(mica/(mica+bic))
		return data
	
	def parameterChecks(**kwargs):
		""" Checking IC and CS specific arguments.
		Check whether different model parameters meet requirements, i.e.,
		all model parameters set in the required range and pass if true,
		otherwise, exit processing
		"""
		if not kwargs: pass
		else:		
			if 'sigma' in kwargs and not (isinstance(kwargs['sigma'], float) or 0 < kwargs['sigma'] < 1): 
				print(InputError('sigma parameter - Value or Type Error', 'The value of sigma for the Zho et al. approach, if provided\nshould be a float ranging between 0 and 1, by default it is 0.5.\n\nPlease refer to the tool documentation, fix this issue and try again ...\n'))
				sys.exit(3)
			if 'TermStats' in kwargs and not isinstance(kwargs['TermStats'], dict):
				print(InputError('TermStats - Type Error', 'TermStats should be a dictionary: key (term)/value (statistical count) mapping.\n\nPlease refer to the tool documentation, fix this issue and try again ...\n'))
				sys.exit(4)
			if 'TermIC' in kwargs and not isinstance(kwargs['TermIC'], dict):
				print(InputError('TermIC - Type Error', 'TermIC should be a dictionary: key (term)/value (IC scores) mapping.\n\nPlease refer to the tool documentation, fix this issue and try again ...\n'))
				sys.exit(5)
	#		#Setting potential parameters for CS 
			if 'cf' in kwargs and not kwargs['cf'] in [0,1,2,3]:
				print(InputError('Correction factor - Value Error', 'Please, the potential correction factor should be an integer between 0 and 3:\n0. No correction is applied\n1. Graph-based (XGraSM or EICA) \n2. Relevance\n3. Information coefficient (simIC)\nCheck this parameter and try again...\n'))
				sys.exit(6)
			if 'gr' in kwargs and not kwargs['gr'] in [0,1]:
				print(InputError('Please, the choice for the Graph-based (XGraSM or EICA) correction factor should be 0 or 1:\n0. For XGraSM model and \n1. For EICA model.\nCheck this parameter and try again.'))
				sys.exit(7)
			if 'jv' in kwargs and not kwargs['jv'] in [0,1,2,3,4,5]: # For Jiang and Conrath approach
				print(InputError('Correction factor - Value Error', 'Please, the potential correction factor should be an integer between 0 and 5: jv = \n\t0 for the Resnik-based normalization\n\t1 for the Couto-based normalization\n\t2 for the Leacock & Chodorow normalization\n\t3 for the Garla & Brandt normalization\n\t4 for the Rada normalization\n\t5 for the canonical normalization\nCheck this parameter and try again...\n'))
				sys.exit(8)
			if 'ak' in kwargs and not (isinstance(kwargs['ka'], int) or kwargs['ka'] < 1): # Al-Mubaid parameter
				print(InputError('k parameter - Value or Type Error', 'The value of the parameter, k, for the Al-Mubaid et al. approach, if provided\nshould be an interger greater or equal to 1, by default it is 1.\n\nPlease refer to the tool documentation, fix this issue and try again ...\n'))
				sys.exit(9)
			if 'aa' in kwargs and not (isinstance(kwargs['aa'], float) or kwargs['aa'] <= 0): # Al-Mubaid parameter
				print(InputError('alpha (aa) parameter - Value or Type Error', 'The value of the parameter alpha, aa, for the Al-Mubaid et al. approach, if provided\nshould be a no zeros positive float, by default it is 1.\n\nPlease refer to the tool documentation, fix this issue and try again ...\n'))
				sys.exit(10)
			if 'ab' in kwargs and not (isinstance(kwargs['ab'], float) or kwargs['ab'] < 0): # Al-Mubaid parameter
				print(InputError('beta (ab) parameter - Value or Type Error', 'The value of the parameter beta, ab, for the Al-Mubaid et al. approach, if provided\nshould be a no zeros positive float, by default it is 1.\n\nPlease refer to the tool documentation, fix this issue and try again ...\n'))
				sys.exit(11)
			if 'zk' in kwargs and not (isinstance(kwargs['zk'], int) or kwargs['zk'] < 2): # Zhong parameter
				print(InputError('k parameter - Value or Type Error', 'The value of the parameter, k, for the Zhong et al. approach, if provided\nshould be an integer strictly greater than 1, by default it is 2.\n\nPlease refer to the tool documentation, fix this issue and try again ...\n'))
				sys.exit(12)
			if 'alpha' in kwargs and not (isinstance(kwargs['alpha'], float) or kwargs['alpha'] < 0): # Edge-based Li parameter
				print(InputError('alpha parameter - Value or Type Error', 'The value of the parameter alpha for the edge-based Li et al. approach, if provided\nshould be a positive float, by default it is 0.2.\n\nPlease refer to the tool documentation, fix this issue and try again ...\n'))
				sys.exit(13)
			if 'beta' in kwargs and not (isinstance(kwargs['beta'], float) or kwargs['beta'] <= 0): # Edge-based Li parameter
				print(InputError('beta (ab) parameter - Value or Type Error', 'The value of the parameter beta for the edge-based Li et al. approach, if provided\nshould be a no zeros positive float, by default it is 0.6.\n\nPlease refer to the tool documentation, fix this issue and try again ...\n'))
				sys.exit(14)
	
	def conceptInterface(self, TermPairs, models = ('nunivers', 'universal'), **kwargs):
		try: kwargs['app'] = models[1]
		except: pass
		return getattr(self, self.Models[models[0]])(TermPairs, **kwargs)

	def computeSim(self, TermPairs, models = None, **kwargs): 
		"""This function computes concept similarity scores given:
           - List of concepts or concept pairs : TermPairs
           - List of tuples similarity model-IC approach: models
           or list of tuples similarity model-None if sim model does
           not require IC approach, e.g., edge-based models
           - Other measure parameter if requred
		"""
		#Checking parameters:
		self.parameterChecks(**kwargs)
		
		if not TermPairs:
			print(InputError('TermPairs - Value Error', 'Provide the list or tuple of concepts or pairs of concepts.\n\nPlease refer to the tool documentation, fix this issue and try again ...\n'))
			sys.exit(1)

		self.models = []
		if models is None or not models: 
			self.models.append(('nunivers','universal'))
		elif isinstance(models, (tuple, list, set)):
			for mod in models:
				if isinstance(mod, (list, tuple)) and len(mod)==2:
					if mod[0].lower() in self.Models and mod[1].lower() in self.AppMods:
						self.models.append((mod[0].lower(), mod[1].lower()))
				elif isinstance(mod, (list, tuple)) and len(mod)==1:
					if mod[0].lower() in self.Models:
						self.models.append((mod[0].lower(),))
				elif isinstance(mod, str) and mod.lower() in self.Models:
					self.models.append((mod.lower(),))
		elif isinstance(models, str) and models.lower() in self.Models:
			self.models.append((models.lower(),))

		if not self.models: # proposed approach, not known 
			print(InputError('Concept semantic similarity model- Value Error', 'Please provide the list or tuple of semantic similarity model and IC approach, if appliable, (ssm, ic), pairs.\n\nPlease refer to the tool documentation, fix this issue and try again ...\n'))
			sys.exit(2)

		self.TargetPairs = {}; self.TargetMissing = {}
		if isinstance(TermPairs, (list, set, tuple)):
			for tt in TermPairs:
				if tt in self.alt_id and self.alt_id[tt] in self.DagStr:
					self.TargetPairs[self.alt_id[tt]] = tt
				elif tt in self.Dag and self.Dag.index(tt) in self.DagStr:
					self.TargetPairs[self.Dag.index(tt)] = tt
				else: # Term is either obsolete or does not exist in the current onto!
					self.TargetMissing.add(tt)

		tmodels = set([s[0] for s in self.models])
		if not self.ShortPath and set(self.CatPaths) & tmodels: 
			# Get all the shortest path lengths for a to b in DAG
			self.ShortPath = nx.shortest_path_length(self.DagStr)
		if not self.DicLevels and set(self.CatLevels) & tmodels:
			DicLevels = nx.bellman_ford(self.DagStr, self.oroot) 
			self.DicLevels = DicLevels[-1]
			del DicLevels
			self.deep = -min(list(set(self.DicLevels.values())))
		del tmodels
		
		self.outputs = {}; pq = []; self.keepterms = sorted(list(self.TargetPairs.keys()))
		for k in [(p, q) for p in self.keepterms for q in self.keepterms]:
			if not (k in pq or (k[1],k[0]) in pq): pq.append(k)

		for fct in self.models:
			try: kwargs['app'] = fct[1]
			except: pass
			self.outputs[fct] = getattr(self,self.Models[fct[0]])(pq, **kwargs)

	@output_str
	def __str__(self):
		"""Return a string reprensentation of the concept-semantic similarity scores.
        """
		headers = ['Concept1', 'Concept2']; tmp = list(self.outputs.keys())
		for ms in tmp:
			if len(ms)==2: headers.append('{}-{}'.format(ms[0].capitalize(), ms[1][:4] if len(ms[1])>3 else ms[1]))
			else: headers.append('{}-{}'.format(ms[0].capitalize(), ''))
		nn = len(self.keepterms)
		results = [[self.Dag[k[0]], self.Dag[k[1]]] + [self.outputs[s][k] for s in tmp] for k in [(self.keepterms[i], self.keepterms[j]) for i in range(nn) for j in range(i,nn)]]
		return '\n'+tabs(results, headers, tablefmt = 'grid', floatfmt=".5f", stralign="center")

	def __repr__(self):
		""" Return repr(self)
		"""
		return self.outputs
