#!/usr/bin/python
# -*- coding: utf8 -*-

from __future__ import print_function
import sys, os, time

dcd = os.getcwd()
sys.path.append('/'.join(dcd.split('/')[:-1]))

from PySML import *
#from PySML import EntitySimilarity

if __name__=='__main__':
	now = time.time()
	funsim = EntitySimilarity()
	EntityAnnots = dict([("Prot1", ['GO:0000022', 'GO:0051231', 'GO:1903047', 'GO:0000278', 'GO:0007052', 'GO:0000023', 'GO:0005984']), ("Prot2", ['GO:0000022', 'GO:0051231', 'GO:1903047', 'GO:0000278', 'GO:0007052']), ("Prot3", ['GO:1903047', 'GO:0000278', 'GO:0007052', 'GO:0000023', 'GO:0005984'])])
	#funsim.entitySim(EntityAnnots, measures=[('bma', 'nunivers','universal'), ('simgic', 'zhang'), ('simdic', 'seco'), ('bmm', 'lin', 'zanchez'), 'ub'])
	funsim.entitySim(EntityAnnots, measures = ('simgic', 'simui', 'ssdd'))
	print(funsim)
	print("The time elapsed is:", time.time()-now)
	
	#"Prot1", ['GO:0000022', 'GO:0051231', 'GO:1903047', 'GO:0000278', 'GO:0007052', 'GO:0000023', 'GO:0005984']), ("Prot2", ['GO:0000022', 'GO:0051231', 'GO:1903047', 'GO:0000278', 'GO:0007052']), ("Prot3", ['GO:1903047', 'GO:0000278', 'GO:0007052', 'GO:0000023', 'GO:0005984']
	# "{'Prot1':['GO:0000022', 'GO:0051231', 'GO:1903047', 'GO:0000278', 'GO:0007052', 'GO:0000023', 'GO:0005984'], 'Prot2':['GO:0000022', 'GO:0051231', 'GO:1903047', 'GO:0000278', 'GO:0007052'], 'Prot3':['GO:1903047', 'GO:0000278', 'GO:0007052', 'GO:0000023', 'GO:0005984']}"
	# python procsemsim.py -t es -m bma kstats ui simgic intel -d mama -a "{'Prot1':['GO:0000022', 'GO:0051231', 'GO:1903047', 'GO:0000278', 'GO:0007052', 'GO:0000023', 'GO:0005984'], 'Prot2':['GO:0000022', 'GO:0051231', 'GO:1903047', 'GO:0000278', 'GO:0007052'], 'Prot3':['GO:1903047', 'GO:0000278', 'GO:0007052', 'GO:0000023', 'GO:0005984']}" -f /home/user/SADaCC-SPARCoProject/go-basic.obo
