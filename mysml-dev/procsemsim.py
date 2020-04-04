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
 
http://web.cbio.uct.ac.za/ITGOM/PySML
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

import sys, os, re, time, datetime, random

from PySML import __author__, __author_email__, __version__, __license__
from PySML import *

def main():
	from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
	parser = ArgumentParser(description = __doc__, formatter_class = ArgumentDefaultsHelpFormatter)
	parser.add_argument("-a", "--annotationfile", dest = "annot", type = str, help="Full path to the file appropriate Entity-term mapping ")
	parser.add_argument("-f", "--ontologyfile", dest = "ontology", type = str, help="Full path to the file appropriate Entity-term mapping ")
	parser.add_argument("-m", "--models", dest="models", nargs = '*', type = str, help="SS models to be considered ")
	parser.add_argument("-o", "--out", dest="outfile", type = str, help = "Naming the SS scores output file ", default = os.getcwd(), metavar="FILE")
	parser.add_argument("-d", "--data", dest = "data", type = str, nargs = '+', help="Full path to the file containing list of terms, term-term pairs, entity-entity ", metavar="str")
	parser.add_argument("-t", "--type", dest="mtype", type = str, help="Type of SS: ic (IC), cs (Term SS), es (Entity SS) ", required = True, metavar='str')
	parser.add_argument("-s", "--stream", dest="stream", type = int, help="Output on (1) screen or (0) file", default = 1, metavar='int')
	
	argss = parser.parse_args()

	# Providing system and platform requirements
	print('\n'+74*'*')
	print(('PySML v{}'.format(__version__)+', developed by {}'.format(__author__)).center(71))
	print(('E-mail: {}'.format(__author_email__)).center(71))
	print('\n'+74*'*'+'\n*'+'These are requirements for running Integrated Human PPI Generator tool'.center(71)+' *\n'+74*'*')
	print('\nPlatform: Python version >= 2.7.x\n\nOS      : Independent, but tested only on Linux (Ubuntu)')
	print("\nThe current PyPML library version retrieves semantic similarity scores\nof all known SS models:\n")
	print("For more information, please refer to the manual at:\n\t1. http://web.cbio.uct.ac.za/ITGOM/pysml/PySML_Documentation.pdf\n\t2. http://web.cbio.uct.ac.za/ITGOM/adagofun/PKG-INFO\n\t3. https://github.com/gkm-software-dev/post-analysis-tools\n\nOr go to the PySML directory and type the following command line:")
	print("\n\tpython setup.py --long-description")
	print("\nRequires: python networkx library")
	print('\nOUTPUT  : SS scores are displayed on the screen or written to a file.\n')
	print(('PySML is being run on {}, under'.format(datetime.datetime.now().strftime('%d-%m-%Y %H:%M'))).center(71))
	print(('{}'.format(__license__)).center(71))
	print(74*'*')
	print('\nEnter 1 to continue and 2 to exit\n')
	
	if sys.version_info.major < 3: raw_pass = raw_input
	else: raw_pass = input
	while True:
		a = raw_pass('> ')
		try: 
			a = int(a)
			if a in [1,2]: break
			print('Please enter 1 to continue and 2 to exit')
		except: print('Please enter 1 to continue and 2 to exit')
	if a==2:
		print("\nExiting PySML. Thanks for using this tool! You can try again at any time.\n")
		print(74*'*'+'\n')
		sys.exit(2)
		
	# Check whether python network library required for creating ontology DAG exists
	import imp
	spec = None
	try:
		imp.find_module('networkx')
		spec = True
	except ImportError: spec = None
	if spec is None:
		print('\n'+74*'*'+'\n'+"python-networkx library under Python has not been found.")
		print("install python-networkx and try again. Execution cannot be pursued, now exiting ...\n"+74*'*'+'\n')
		sys.exit(2)
	
	# Checking whether required arguments are provided
	allentsim = set(['avg', 'bma', 'abm', 'bmm', 'hdf', 'vhdf', 'max', 'aln', 'intel', 'spgk', 'lp', 'ye', 'simgic', 'simdic', 'simuic', 'simcou', 'simcot', 'simui', 'simub', 'simdb', 'simnto', 'simcub', 'simctb', 'cho', 'ald', 'kstats', 'nto', 'ub', 'db', 'ui'])
	allconsim = set(['resnik', 'lin', 'nunivers', 'wang', 'jiang', 'faith', 'aic', 'hrss', 'ps', 'shen', 'wu', 'pekar', 'slimani', 'shenoy', 'wang_edge', 'zhong', 'almubaid', 'rss', 'ssdd', 'rada', 'resnik_edge', 'leacock', 'li_edge', 'shenoy', 'almubaid', 'rss'])
	allconic = set(['universal', 'wang', 'zhang', 'seco', 'zho', 'seddiqui', 'zanchez', 'meng', 'stats', 'vic', 'assdd'])
	measures = {'ic':(allconic, [('universal',)]), 'cs': (allconsim, [('nunivers', 'universal')]), 'es': (allentsim, [('bma', 'nunivers', 'universal')])}
	
	nnt = set(['cho', 'ald', 'kstats', 'nto', 'ub', 'db', 'ui']); models = []
	if argss.models is None: models = measures[measures[argss.mtype][1]]
	elif isinstance(argss.models, str) and argss.models.strip():
		if argss.models.strip() in measures[argss.mtype][0]: models = [(argss.models.strip(),)]
		else:
			print()
			parser.error("\nModel Error -- Model or model type provided does not exist\nPlease read the library manual for more information\nargument -m/--models: expected an appropriate argument.\n")
	elif isinstance(argss.models, list):
		for mod in argss.models:
			tline = re.split(":|,|-", mod.strip())
			if len(tline) in [1, 2, 3] and tline[0] in measures[argss.mtype][0]: models.append(tuple(tline))
			else:
				print()
				parser.error("\nModel Error -- Model provided is not appropriate or do not match Model Type\nPlease read the library manual for more information\nargument -m/--models: expected an appropriate argument.\n")
	else:
		print()
		parser.error("\nModel structure format Error\nPlease read the library manual for more information\nargument -f/--ontologyfile: expected ontology file.\n")

	if not models:
		print()
		parser.error("\nModel Error -- Model format provided is not valid\nPlease read the library manual for more information\nargument -m/--models: expected an appropriate argument.\n")
		
	# This means that the ontology must be provided
	Pairs = []; Annots = {}
	if argss.mtype in ['ic', 'cs']: 
		if argss.ontology is None or (isinstance(argss.ontology, str) and not argss.ontology.strip()):
			print()
			parser.error("\nOntology file: Value Error\nFor the type of measure chosen, an ontology file should be provided\nargument -f/--ontologyfile: expected ontology file.\n")
		if argss.data is None and argss.mtype=='cs':
			print()
			parser.error("\nConcept or concept pair list : Value Error\nFor the type of measure chosen, a list of concepts or concept pairs should be provided\nargument -d/--data: expected list pof concepts or concept pairs.\n")
		elif isinstance(argss.data, str) and argss.data.strip():
			try:
				arg = os.path.abspath(argss.data)
				fp = open(arg)
				for line in fp:
					tline = line.strip()
					if not tline: continue
					tline = [s.strip() for s in re.split("\s+|,|;", tline)]
					if len(tline)==2: Pairs.append(tuple(tline))
					elif len(tline)==1: Pairs.append(tline[0])
				fp.close()
			except:
				for mod in argss.data:
					tline = re.split(",|-", mod.strip())
					if len(tline)==2: Pairs.append(tuple(tline))
					elif len(tline)==1: Pairs.append(tline[0])
		if not Pairs and argss.mtype=='cs':
			print()
			parser.error("\nConcept or concept pair list : Value Error\nA list of concepts or concept pairs provided is not valid\nargument -d/--data: expected list of concepts or concept pairs.\n")

	elif argss.mtype=='es':
		if set([s[0] for s in models])-nnt: # When the ontology is required, in case where ontology es is provided
			if argss.ontology is None or (isinstance(argss.ontology, str) and not argss.ontology.strip()):
				print()
				parser.error("\nOntology file: Value Error\nFor the type of measure chosen, an ontology file should be provided\nargument -f/--ontologyfile: expected ontology file.\n")
		if argss.annot is None or (isinstance(argss.annot, str) and not argss.annot.strip()):
			print()
			parser.error("\nAnnotation file: Value Error\nFor the type of measure chosen, an annotation file, entity-concept mapping should be provided\nargument -a/--annotationfile: expected annotation file.\n")
		else:
			try:
				arg = os.path.abspath(argss.annot)
				fp = open(arg)
				for line in fp:
					tline = line.strip()
					if not tline: continue
					tline = [s.strip() for s in re.split("\s+", tline)]
					if len(tline)==2:
						Annots[tline[0]] = set(re.split(";|:|,|-", tline))
				fp.close()
			except:
				try: Annots = eval(argss.annot)
				except:
					parser.error("\nAnnotation variable: Value Error\nFor the type of measure chosen, the annotation variable entity-concept\n\nmapping should be a dictionary-like string\nargument -a/--annotationfile: expected annotation file.\n")

		if not isinstance(Annots, dict) or not Annots:
			parser.error("\nAnnotation variable: Value or Type Error\nFor the type of measure chosen, the annotation variable entity-concept mapping should be a no empty dictionary\nargument -a/--annotationfile: expected annotation file.\n")
		if isinstance(argss.data, str) and argss.data.strip():
			try:
				arg = os.path.abspath(argss.data)
				fp = open(arg)
				for line in fp:
					tline = line.strip()
					if not tline: continue
					tline = [s.strip() for s in re.split("\s+|,|;", tline)]
					if len(tline)==2: Pairs.append(tuple(tline))
					elif len(tline)==1: Pairs.append(tline[0])
				fp.close()
			except:
				for mod in argss.data:
					tline = re.split(",|-", mod.strip())
					if len(tline)==2: Pairs.append(tuple(tline))
					elif len(tline)==1: Pairs.append(tline[0])
				print(argss.data, Pairs)
		print(Pairs)
		if not Pairs: Pairs = [(p, q) for p in Annots for q in Annots]

	# Now processing different scores
	now = time.time()
	if argss.mtype=='ic':
		simscore = InformationContent() # getIC(self, approach = None, TermList = None, **kwargs)
		simscore.getIC([s[0] for s in models])
		ScoreFile = 'InformationContentFile%d.txt'%(random.randint(0,100000),)
	elif argss.mtype=='cs':
		pass
		ScoreFile = 'EntitySSFile%d.txt'%(random.randint(0,100000),)
	elif argss.mtype=='es':
		simscore = EntitySimilarity(ontofile = argss.ontology, namespace = 'biological_process')
		simscore.entitySim(Annots, Pairs, models)
		ScoreFile = 'EntitySSFile%d.txt'%(random.randint(0,100000),)
	
	# Finally, outputting the score on screen or written into a file! 
	if argss.stream:
		print(simscore)
	else:# Print in a file
		arg = '/'.join([argss.outfile, ScoreFile])
		try:
			fw = open(arg, 'w')
			fw.write(repr(simscore)+'\n')
			print("Scores are reported in the file: %s"%(arg,))
		except:
			parser.error("\nWriting output error\nImpossible to write in %s\nargument -o/--out: Output error, scores cannot be written in the file.\n"%(arg,))
			
	print("\nThanks for choosing PySML. Start processing on %s"%str(time.asctime(time.localtime())))
	
	print("Processing accomplished on %s. Thanks!"%str(time.asctime(time.localtime())))
	tt = time.time()-now
	nh = tt//3600; rs = tt-nh*3600
	nm = rs//60; rs -= nm*60 
	print("Total time elapsed is approximately %02d:%02d:%02d"%(nh, nm, rs))
	print(74*'*'+'\n')

if __name__=='__main__':
	main()
