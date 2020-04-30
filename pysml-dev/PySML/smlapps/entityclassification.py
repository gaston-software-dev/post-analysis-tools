#!/usr/bin/python
# -*- coding: utf8 -*-

"""Important Note:
This python module is part of the PySML library, which is a tool for Gene 
Ontology-based functional analysis using term information content 
measures.
This python code implements hierarchical, Graph spectral (kmeans) and 
model or communuity-based clustering of genes and proteins based on their 
Gene Ontology annotations. Hierarchical and Graph spectral approaches are 
used as implemented under scipy and community-based approach is implemen-
ted using following modules as written by 
Thomas Aynaud <thomas.aynaud@lip6.fr>, 2009: "partition_at_level", 
"best_partition", "generate_dendogram", "induced_graph" .

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
from functools import reduce

import sys, time, random
import networkx as nx

from .. import __version__
from ..error import InputError
from ..imports.readontology import Ontology, output_str
from ..imports.tabulate import tabulate as tabs
from ..entitysimilarity import EntitySimilarity

try:
	import networkx as nx
except ImportError:
	print(InputError('networkx library not installed', 'networkx library ImportError on line 60 in \nentityclassification.py source code\n\nInstall the networkx library and try again ...\n'))
	sys.exit(0)

try: # It is sufficient to import "scipy" here once for some useful mathematical objects, functions and algorithm needed.
	from scipy import exp, array, zeros, diag, unique, mean, sqrt, dot, ceil
	from scipy.cluster import hierarchy
	from scipy.spatial import distance
	from scipy.cluster.vq import kmeans2 
	from scipy.linalg import svd
except ImportError:
	print(InputError('The scipy library not installed', 'scipy library ImportError on line 70 in \nentityclassification.py source code\n\nInstall the scipy library and try again ...\n'))
	sys.exit(1)

try: # It is sufficient to import "matplotlib" here once for plotting.
	import matplotlib.pylab as plt 
except ImportError:
	print(InputError('The matplotlib library not installed', 'matplotlib library ImportError on line 73 in \nentityclassification.py source code\n\nInstall the matplotlib library and try again ...\n'))
	sys.exit(2)

__PASS_MAX = -1
__MIN = 0.0000001

class Status:
	""" Handling several data in one structure:
	    A class as implemented by Thomas Aynaud for model-based
	    classification.
	"""
	node2com = {}
	total_weight = 0
	internals = {}
	degrees = {}
	gdegrees = {}
		
	def __init__(self):
		self.node2com = dict([])
		self.total_weight = 0
		self.degrees = dict([])
		self.gdegrees = dict([])
		self.internals = dict([])
		self.loops = dict([])
	
	def copy(self) :
		"""Perform a deep copy of status"""
		new_status = Status()
		new_status.node2com = self.node2com.copy()
		new_status.internals = self.internals.copy()
		new_status.degrees = self.degrees.copy()
		new_status.gdegrees = self.gdegrees.copy()
		new_status.total_weight = self.total_weight

	def init(self, graph, part = None) :
		"""Initialize the status of a graph with every node in one community"""
		count = 0
		self.node2com = dict([])
		self.total_weight = 0
		self.degrees = dict([])
		self.gdegrees = dict([])
		self.internals = dict([])
		self.total_weight = graph.size(weight = 'weight')
		if part == None :
			for node in graph.nodes() :
				self.node2com[node] = count
				deg = float(graph.degree(node, weight = 'weight'))
				self.degrees[count] = deg
				self.gdegrees[node] = deg
				self.loops[node] = float(graph.get_edge_data(node, node, {"weight":0}).get("weight", 1))
				self.internals[count] = self.loops[node]
				count = count + 1
		else :
			for node in graph.nodes() :
				com = part[node]
				self.node2com[node] = com
				deg = float(graph.degree(node, weigh = 'weight'))
				self.degrees[com] = self.degrees.get(com, 0) + deg
				self.gdegrees[node] = deg
				inc = 0.
				for neighbor, datas in graph[node].items() :
					weight = datas.get("weight", 1)
					if part[neighbor] == com :
						if neighbor == node : inc += float(weight)
						else : inc += float(weight) / 2.
				self.internals[com] = self.internals.get(com, 0) + inc

def __modularity(status) :
	""" Compute the modularity of the partition of the graph faslty using status precomputed """
	links = float(status.total_weight)
	result = 0.
	for community in set(status.node2com.values()) :
		in_degree = status.internals.get(community, 0.)
		degree = status.degrees.get(community, 0.)
		if links > 0 : result += in_degree / links - ((degree / (2.*links))**2)
	return result

def __renumber(dictionary) :
	""" Renumber the values of the dictionary from 0 to n """
	count = 0
	ret = dictionary.copy()
	new_values = dict([])

	for key in dictionary.keys() :
		value = dictionary[key]
		new_value = new_values.get(value, -1)
		if new_value == -1 :
			new_values[value] = count
			new_value = count
			count += 1
		ret[key] = new_value
	return ret

def induced_graph(partition, graph):
	"""Produce the graph where nodes are the communities """
	ret = nx.Graph()
	ret.add_nodes_from(partition.values())

	for node1, node2, datas in graph.edges_iter(data = True) :
		weight = datas.get("weight", 1)
		com1 = partition[node1]
		com2 = partition[node2]
		w_prec = ret.get_edge_data(com1, com2, {"weight":0}).get("weight", 1)
		ret.add_edge(com1, com2, weight = w_prec + weight)
	return ret

def __neighcom(node, graph, status) :
	""" Compute the communities in the neighborood of node in the graph given
    with the decomposition node2com
	"""
	weights = {}
	for neighbor, datas in graph[node].items() :
		if neighbor != node :
			weight = datas.get("weight", 1)
			neighborcom = status.node2com[neighbor]
			weights[neighborcom] = weights.get(neighborcom, 0) + weight

	return weights


def __remove(node, com, weight, status) :
	""" Remove node from community com and modify status"""
	status.degrees[com] = ( status.degrees.get(com, 0.) - status.gdegrees.get(node, 0.) )
	status.internals[com] = float( status.internals.get(com, 0.) - weight - status.loops.get(node, 0.) )
	status.node2com[node] = -1


def __insert(node, com, weight, status) :
	""" Insert node into community and modify status"""
	status.node2com[node] = com
	status.degrees[com] = ( status.degrees.get(com, 0.) + status.gdegrees.get(node, 0.) )
	status.internals[com] = float( status.internals.get(com, 0.) + weight + status.loops.get(node, 0.) )


def __one_level(graph, status) :
	""" Compute one level of communities """
	modif = True
	nb_pass_done = 0
	cur_mod = __modularity(status)
	new_mod = cur_mod

	while modif  and nb_pass_done != __PASS_MAX :
		cur_mod = new_mod
		modif = False
		nb_pass_done += 1

		for node in graph.nodes() :
			com_node = status.node2com[node]
			degc_totw = status.gdegrees.get(node, 0.) / (status.total_weight*2.)
			neigh_communities = __neighcom(node, graph, status)
			__remove(node, com_node, neigh_communities.get(com_node, 0.), status)
			best_com = com_node
			best_increase = 0
			for com, dnc in neigh_communities.items() :
				incr =  dnc  - status.degrees.get(com, 0.) * degc_totw
				if incr > best_increase :
					best_increase = incr
					best_com = com                    
			__insert(node, best_com, neigh_communities.get(best_com, 0.), status)
			if best_com != com_node : modif = True                
		new_mod = __modularity(status)
		if new_mod - cur_mod < __MIN : break

def generate_dendogram(graph, part_init = None):
	""" Finding the communities in the graph and return the associated dendogram """
	current_graph = graph.copy()
	status = Status()
	status.init(current_graph, part_init)
	mod = __modularity(status)
	status_list = list()
	__one_level(current_graph, status)
	new_mod = __modularity(status)
	partition = __renumber(status.node2com)
	status_list.append(partition)
	mod = new_mod
	current_graph = induced_graph(partition, current_graph)
	status.init(current_graph)

	while True :
		__one_level(current_graph, status)
		new_mod = __modularity(status)
		if new_mod - mod < __MIN : break
		partition = __renumber(status.node2com)
		status_list.append(partition)
		mod = new_mod
		current_graph = induced_graph(partition, current_graph)
		status.init(current_graph)
	return status_list[:]

def partition_at_level(dendogram, level) :
	""" Return the partition of the nodes at the given level """
	partition = dendogram[0].copy()
	for index in range(1, level + 1) :
		for node, community in partition.items() :
			partition[node] = dendogram[index][community]
	return partition

def best_partition(graph, partition = None):
	dendo = generate_dendogram(graph, partition)
	return partition_at_level(dendo, len(dendo) - 1)


class EntityClassification(EntitySimilarity):
	__cslot__ = []
	def __init__(self, ontofile = '', namespace='', is_a = None, part_of = None):
		"""Instantiate a new object for a given ontology to retrieve term simila-
		   rity scores with a DAG inherited from EntitySimilarity class and
		   from there, inheriting ConceptSimilarity class.

        Arguments:
            ontology (str): the ontology file
            namespace (str): the name of the subontology as named in the ontology
            is_a (float between 0 and 1, if provided): a semantic value of the
                topological relation is_a.
            part_of (float between 0 and 1, if provided): a semantic value of the
                topological relation part_of.

        Example:
            >>> fct = EntityClassification()
        """
		EntitySimilarity.__init__(self, ontofile, namespace, is_a, part_of)
		self.fouts = {}
		self.EntityMissing = {}
		self.measures = []
		self.comments = {}
		
	def entityfct(self, AnnotationData, **kwargs):
		"""
		The entityfct method for entity functional classification tool. This method 
		allows the partitioning of an entity (e.g. gene or protein) set into a set of 
		meaningful sub-class patterns using their semantic similarity scores computed
		based on entity-associated concepts and derived from a selected semantic sim-
		ilarity model. 
		
		Arguments:
		
			AnnotationData (dict): A dictionary with entity as key and set of concepts
			as value.
			 
			**kwargs can be used to set different parameters needed for processing the
			classification, including measure, mclust and nclust:
			
			measure (str or tuple): The entity semantic similarity measure to be used. 
			Refer to the Supplementary for more details on symbols used for different 
			measures.
			mclust (int 1, 2, 3): Classification model under consideration, and this 
			method implements three different models (mclust):
		  		- hierarchical clustering (mclust = 1)
		  		- Graph spectral clustering or kmeans (mclust = 2)
		  		- community detecting model by Thomas Aynaud, 2009 (mclust = 3)
		  	nclust (int): Number of clusters (nclust) applies only for the kmeans model 
		  	and it is set to 0 by default. In this case, if mclust is less than 2 then
			the comminity detecting model is applied instead of kmeans!
			
			score (float > 0.0): The threshold score providing the semantic similarity 
			degree at which entities are considered to be semantically close or similar in 
			the ontology structure and it is set to 0.3 by default.
			
			stream (int 0 or 1): An Enum parameter taking values of 0 or 1. It is set to 1
			to output results on the screen, to 0 to output results in a file.
			
			Other parameters which can be required depending on the entity semantic simila-
			rity measure used.

		Usage:
		------
			>>> entityfct(AnnotationData, measure=('bma', 'nunivers','unuversal'), score=0.3, mclust=1, nclust=0)

		Examples:
		---------
		    >>> background = {'Q5H9L2':['GO:0006355','GO:0006351'], 'P03891':['GO:0022904','GO:0044281','GO:0044237','GO:0006120'], 'Prot1':['GO:0006355', 'GO:0022904', 'GO:0044281'], 'Prot2':['GO:0044237', 'GO:0006120']}
			>>> entityfct(background, measure=('bma', 'nunivers','unuversal'), score=0.0)
		"""
		# Clustering process starts here
		measure = kwargs['measure'] if 'measure' in kwargs else None
		if not 'mclust' in kwargs: kwargs['mclust'] = 1 # hiererchical-based
		if not 'stream' in kwargs: kwargs['stream'] = 1
		if not 'nclust' in kwargs: kwargs['nclust'] = 0
		
		self.entitySim(AnnotationData, measures = measure, **kwargs)
		g = nx.Graph(); 
		agree = kwargs['score'] if 'score' in kwargs else 0.3
		data = self.fouts[self.measures[-1]]
		for p in data: # Constructing graph goes here!
			if agree==0.0 and data[p] > 0.0: 
				g.add_edge(p[0], p[1], weight=1.0-data[p])
			elif 0 < agree < 1.0 and data[p] >= agree: g.add_edge(p[0], p[1], weight = 1.0-data[p])
			elif agree==1.0 and data[p]==1.0: g.add_edge(p[0], p[1], weight = 1.0)
		models = ['Hierarchical', 'Graph spectral based (kmeans)', 'Model-based']
		now = time.time()
		if g:
			# Outputting different results
			print("Entity classification based on functional similarity")
			print("# The number of entities and entity pairs detected are %d and %d, respectively."%(len(g),g.size()))
			print("The distance is based on functional similarity measure  : %s"%self.comments[self.measures[0]])
			print("The clustering model used is                            : %s %s"%(models[kwargs['mclust']-1], 'approach'))
			if kwargs['stream']:
				print("\nDifferent clusters are displayed in the table below or in the following figure.\nIf possible, use full screen mode for more convenient visualization:")
			else:
				outputfile = 'ConceptSSFile%d'%(random.randint(0,100000),)
				if kwargs['mclust']: outputfile += '.png'	
				else: outputfile += '.txt'
				print("Different clusters can be found in the file       : [%s]"%(outputfile,))
			classes = {}
			if kwargs['mclust']==3:
				partition = best_partition(g); j = 0
				for i in set(partition.values()):
					classes[j] = [nodes for nodes in partition.keys() if partition[nodes] == i]
					j += 1
			elif kwargs['mclust']==1: # Hierarchical approach
				Index = g.nodes()
				n = len(Index)
				distances = zeros((n,n))
				path_length = nx.all_pairs_dijkstra_path_length(g)	
				for u,p in path_length.items():
					for v,d in p.items():
						distances[Index.index(u)][Index.index(v)] = d
						distances[Index.index(v)][Index.index(u)] = d
				sd = distance.squareform(distances)
				hier = hierarchy.average(sd)
				fig = plt.figure()
				plt.clf()
				hierarchy.dendrogram(hier, orientation='right', labels=Index[:])
				plt.grid()
				if kwargs['stream']: plt.show()
				else: plt.savefig(outputfile, format="png")
			elif kwargs['mclust']==2: # Graph spectral based (kmeans)
				d =  kwargs['nclust']; classes = {}# Number of presumed clusters
				if d < 2:
					partition = best_partition(g); j = 0
					for i in set(partition.values()):
						classes[j] = [nodes for nodes in partition.keys() if partition[nodes] == i]
						j += 1
				else:
					W = nx.adj_matrix(g)
					D = diag([sum(sum(array(w))) for w in W])
					L = D - W
					S, V, D = svd(L)
					N = g.nodes()
					test = True; d = 2
					while test:
						cidx = []
						res, idx = kmeans2(S[:,-d+1:], d, minit='random') #-d+1:-1
						cidx = unique(idx)
						classes = {}
						for i in xrange(len(N)):
							if idx[i] in classes: classes[idx[i]].append(N[i])
							else: classes[idx[i]] = [N[i]]

						k = 0 # Checking if nodes in the identified class are connected!
						for iclass in classes.itervalues():
							if not nx.is_connected(nx.subgraph(g,iclass)):
								k = 1
								break
						if not k: test = False
			if classes:
				outs = []
				for i in sorted(classes.keys()):
					st = str(classes[i])[1:-1]
					outs.append(('%d'%(i+1,), '%d'%(len(classes[i]),), st[:78]))
					for i in range(1, int(ceil(len(st)/78.0))):
						try: outs.append(('','',st[i*78:(i+1)*78+1]))
						except: pass
					outs.append(('','','\n'))
				headers = ['# Cluster', '# of Proteins', 'Protein identifiers']
				if kwargs['stream']: # Print on the screen
					print('%s'%tabs(outs[:-1], headers, tablefmt = 'grid', floatfmt=".5f", stralign="center"))
				else:
					try:
						fp = open(outputfile, 'w')
						fp.write("# Clustering proteins based on functional similarity using [%s] approach\n"%(models[kwargs['mclust']-1],))
						fp.write("# The number of possible entities and entity pairs detected are %d and %d, respectively.\n"%(len(g),g.size()))
						fp.write("# The distance is based on [%s] functional similarity measure\n\n"%(appnames[kwargs['measure']],))
						fp.write('%s'%tabs(outs[:-1], headers, tablefmt = 'plain', floatfmt=".5f", stralign="left"))
						fp.close()
					except IOError:
						print("File cannot be opened in writing. Check possibly the writing permission and try again ...")
						sys.exit(8)
		else:
			print('Trying %s clustering approach using the distance inferred from\n %s functional similarity measure has failed. Please, check your presumed\n list and options you have selected, and try again!'%(models[kwargs['mclust']-1], appnames[kwargs['measure']]))
			sys.exit(9)
		print("\nProcessing accomplished on %s"%str(time.asctime(time.localtime())))
		print("Total time elapsed is approximately %.2f %s"%(time.time()-now, 'seconds'))
		print("\n*****************************************************************\n")
		
if __name__=='__main__':
	pass
