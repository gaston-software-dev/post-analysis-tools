#!/usr/bin/python
# -*- coding: utf8 -*-

from __future__ import print_function
import sys, os, time, random

try:
	__IPYTHON__
except NameError:
	print("\nNote that this illustration should be run using Interactive Python (ipython)\nand not just python. Please run this under as follow:\n\tipython tests/cenrichexample.py\n\nOr alternatively, run ipython, by typing ipython from the shell and then type:\n\t run tests/cenrichexample.py\n\nRemember that this is run under the pysml directory and \nipython can be ipython2 or ipython3 or just ipython, depending on the version installed.\n")
	sys.exit(0)

dcd = os.getcwd() 

sys.path.append('/'.join(dcd.split('/')[:-1]))


from PySML.smlapps import ConceptEnrichment

if __name__=='__main__':
	fp = open('tests/ExampleRefSet.txt') # Read reference annotation dataset
	background = {}
	for line in fp:
		tline = [s.strip() for s in line.strip().split('\t')]
		background[tline[0]] = set(tline[1].split(','))
	fp.close()
	
	targets = set(); d = list(background.keys())
	while len(targets) != 10: # Generating randomly 10 targets 
		targets.add(d[random.randint(0, len(d)-1)])
	
	idt = ConceptEnrichment()
	idt.enrichedConcepts(background, targets, ('nunivers', 'universal'))

