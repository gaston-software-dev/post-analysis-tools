#!/usr/bin/python
# -*- coding: utf8 -*-

from __future__ import print_function
import sys, os, time

dcd = os.getcwd()
sys.path.append('/'.join(dcd.split('/')[:-1]))

from PySML import *
#from PySML import ConceptSimilarity

if __name__=='__main__':
	now = time.time()
#	http://purl.obolibrary.org/obo/go/go-basic.obo
#	http://purl.obolibrary.org/obo/go.obo
#	http://purl.obolibrary.org/obo/go.owl
	semsim = ConceptSimilarity()
	now = time.time()
	semsim.computeSim(['GO:0000022', 'GO:0051231', 'GO:1903047', 'GO:0000278', 'GO:0007052', 'GO:0000023', 'GO:0005984'], ['nunivers', ('resnik', 'zhang'), 'wang', 'wang_edge', ('lin', 'zanchez'), 'aic', 'wu', 'hrss', 'jiang'])
	print(semsim)
	print("The time elapsed is:", time.time()-now)

#		self.outputs['slimani'] = self.sslimani_shenoy(pq)
#		self.outputs['shenoy'] = self.sslimani_shenoy(pq, 0)
#		self.outputs['sli_edge'] = self.sli_edge(pq)

#		self.outputs['swang_edge'] = self.wang_edge(pq)
#		self.outputs['szhong'] = self.szhong(pq)
#		self.outputs['salmubaid'] = self.salmubaid(pq)
#		self.outputs['srss'] = self.srss(pq)
#		self.outputs['ssdd'] = self.sssdd(pq)
#		self.outputs['sshen'] = self.sshen(pq)
#		self.outputs['shrss'] = self.shrss(pq)

#		self.outputs['sresnik'] = self.saic(pq, 'zhang')
#		self.outputs['sresnik_xgr'] = self.saic(pq, 'zhang', cf = 1)
#		self.outputs['sresnik_xsi'] = self.saic(pq, 'zhang', cf = 1, gr = 1)
#		self.outputs['sresnik_rel'] = self.saic(pq, 'zhang', cf = 2)
#		self.outputs['sresnik_sic'] = self.saic(pq, 'zhang', cf = 3)

