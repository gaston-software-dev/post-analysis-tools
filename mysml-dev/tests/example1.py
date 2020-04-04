#!/usr/bin/python
# -*- coding: utf8 -*-

from __future__ import print_function
import sys, os, time

dcd = os.getcwd()
sys.path.append('/'.join(dcd.split('/')[:-1]))

from PySML import *
#from PySML import InformationContent

if __name__=='__main__':
	now = time.time()
	ICScores = InformationContent()

	GOIds = ['GO:0000022', 'GO:0051231', 'GO:1903047', 'GO:0000278', 'GO:0007052', 'GO:0000023', 'GO:0005984']
	ICScores.getIC() # By default, it run the Universal approach
	ICScores.getIC('zho')
	ICScores.getIC('seddiqui')
	ICScores.getIC(['meng', 'zanchez', 'zhang', 'wang', 'seco'], GOIds)
	print(ICScores)
	print("The time elapsed is:", time.time()-now)

#python procsemsim.py -t ic -m universal zanchez zhang wang seco -f /home/user/SADaCC-SPARCoProject/go-basic.obo -s 0 
#https://wtgrants.wellcome.ac.uk/Forms/timeout?returnUrl=%2Fforms%2Fen%2FApply%2FNewApplication

# QC metric
#https://docs.google.com/document/d/1EtsbhJA5K2-Fsoqyb2UueVnZQPkTAw6HH6ylmQM-1HY/edit?ts=5e538afc
# 4th SCDO meeting
# https://docs.google.com/document/d/15ZGv9FZoNR9jKn5RF4wqmA61D9sIi0iNNWiNxMwbfZI/edit#heading=h.z0syentme099
# https://github.com/jackmo375/Database
#https://github.com/gkm-software-dev/post-analysis-tools
