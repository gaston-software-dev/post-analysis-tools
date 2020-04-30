#!/usr/bin/python
# -*- coding: utf8 -*-

from __future__ import print_function
import sys, os, time

try:
	__IPYTHON__
except NameError:
	print("\nNote that this illustration should be run using Interactive Python (ipython)\nand not just python. Please run this under as follow:\n\tipython tests/example2.py\n\nOr alternatively, run ipython, by typing ipython from the shell and then type:\n\t run tests/idexample.py\n\nRemember that this is run under the pysml directory and \nipython can be ipython2 or ipython3 or just ipython, depending on the version installed.\n")
	sys.exit(0)

dcd = os.getcwd()
sys.path.append('/'.join(dcd.split('/')[:-1]))

from PySML import *

if __name__=='__main__':
	now = time.time()
	semsim = ConceptSimilarity()
	now = time.time()
	semsim.computeSim(['GO:0000022', 'GO:0051231', 'GO:1903047', 'GO:0000278', 'GO:0007052', 'GO:0000023', 'GO:0005984'], ['nunivers', ('resnik', 'zhang'), 'wang', 'wang_edge', ('lin', 'zanchez'), 'aic', 'wu', 'hrss', 'jiang'])
	print(semsim)
	print("The time elapsed is:", time.time()-now)
