#!/usr/bin/python
# -*- coding: utf8 -*-

from __future__ import print_function
import sys, os, time

dcd = os.getcwd()
sys.path.append('/'.join(dcd.split('/')[:-1]))

from PySML import *
#from PySML import InformationContent

if __name__=='__main__':
	ICScores = InformationContent()

	GOIds = ['GO:0000022', 'GO:0051231', 'GO:1903047', 'GO:0000278', 'GO:0007052', 'GO:0000023', 'GO:0005984']
	ICScores.getIC() # By default, it run the Universal approach
	ICScores.getIC('zho')
	ICScores.getIC('seddiqui')
	ICScores.getIC(['meng', 'zanchez', 'zhang', 'wang', 'seco'], GOIds)
	print(ICScores)

