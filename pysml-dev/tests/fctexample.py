#!/usr/bin/python
# -*- coding: utf8 -*-

from __future__ import print_function
import sys, os, time

try:
	__IPYTHON__
except NameError:
	print("\nNote that this illustration should be run using Interactive Python (ipython)\nand not just python. Please run this under as follow:\n\tipython tests/fctexample.py\n\nOr alternatively, run ipython, by typing ipython from the shell and then type:\n\t run tests/idexample.py\n\nRemember that this is run under the pysml directory and \nipython can be ipython2 or ipython3 or just ipython, depending on the version installed.\n")
	sys.exit(0)

dcd = os.getcwd()
sys.path.append('/'.join(dcd.split('/')[:-1]))


from PySML.smlapps import EntityClassification

if __name__=='__main__':
	fct = EntityClassification()
	background = {'Q5H9L2':['GO:0006355','GO:0006351'], 'P03891':['GO:0022904','GO:0044281','GO:0044237','GO:0006120'], 'Prot1':['GO:0006355', 'GO:0022904', 'GO:0044281'], 'Prot2':['GO:0044237', 'GO:0006120']}
	data = dict(mclust=1)
	
	fct.entityfct(background, **data)

