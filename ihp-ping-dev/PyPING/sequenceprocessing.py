#!/usr/local/bin/python
# -*- coding: utf8 -*-

# Copyright (C) 2020 Gaston K. Mazandu @HumanGenetics-UCT

from __future__ import print_function

import subprocess
from math import erf, sqrt

try:
	from math import log2
except ImportError:
	from math import log
	log2 = lambda x: log(x)/log(2)

def computeScoreBlast(Input):
	""" This program computes interaction reliability using similarity score produced from
	    BLAST and return all result in the file named yournamefileSeqSim.txt whose fields
	    are source_id, dest_id, interaction reliability obtained.
	    
	    Arguments:
	    	Input (FILE): The file from the BLAST ouputs containing protein and bit scores
	"""
	Control = {}
	fp = open(Input)
	count = 0
	stringp = "\rReading Blast Scores: %s %d%s"%(45*'.',0,'%')
	nl = subprocess.check_output(["/bin/sh", "-c", "wc -l %s ; exit 0"%(Input,)]).decode('utf-8')
	nl = int(nl.split()[0])
	for line in fp:
		print(stringp, end=',')
		count += 1
		ligne = line.strip() 
		if not ligne: continue  # skip empty
		ligne = ligne.split('\t')
		key = ligne[0].split('|')[1].strip(),ligne[1].split('|')[1].strip()
		if key in Control:
			if Control[key] < float(ligne[-1].strip()): Control[key] = float(ligne[-1].strip())
		else: Control[key] = float(ligne[-1].strip())
		perc = (count*45)//nl
		stringp = "\rReading Blast Scores: %s %d%s"%(perc*'='+(45-perc)*'.',count*100//nl,'%')
	print("\rReading Blast Scores: %s %d%s"%(45*'=',100,'%'))
	print("Reading successfully completed ...")
	fp.close()
	for key in Control.keys(): # Removing no symetrical existence
		if not (key[1],key[0]) in Control: Control.pop(key, None)

	KeepScores = {}
	count = 0
	stringp = "\rComputing Reliability Scores: %s %d%s"%(38*'.',0,'%')
	nl = len(Control)
	for key in Control:
		print(stringp, end=',')
		count += 1
		if key in KeepScores or (key[1], key[0]) in KeepScores: continue
		if key[0] != key[1]:
			try:
				Lrel1 = Control[key]+Control[key[1],key[0]]
				Lrel1 /= 2*max(Control[key[0],key[0]],Control[key[1],key[1]])
				KeepScores[tuple(sorted(key))] = round(Lrel1,5)
			except:
				pass
		perc = (count*38)//nl
		stringp = "\rComputing Reliability Scores: %s %d%s"%(perc*'='+(38-perc)*'.',count*100//nl,'%')
	print("\rComputing Reliability Scores: %s %d%s"%(38*'=',100,'%'))
	print("Process successfully completed ...") 
	return KeepScores
	
SignatureNumber = {}; ControlProt = {}

def quartile(x, n = 2):
	'''This function computes the quartile 1, 2 (median) and 3 for a given no empty array 
	   or list of data x and by default, it will be computing the median
	   
	   Arguments:
	      x (list, vector): A list or a vector of numbers
	      n (int): Type of quartile 1, 2, 3
	'''
	ax = sorted(x)
	if n==1: r = .25
	elif n==2: r = .5
	else: r = .75
	p = 1.0 + (len(ax)-1)*r
	ip = int(p)
	delta = (ax[ip]-ax[ip-1])*(p-ip)
	return ax[ip-1]+delta

def computeFamilyScore(SignatureNumber, ControlProt, alpha = 0.75, penal = 0.675):
	'''This function computes the scores of between proteins sharing common signature and 
	   output file containing these interactions, with calibration parameter alpha = 0.75 
	   and sigma penalized with a factor 0.675
	   
	   Arguments:
	      SignatureNumber (dict): The map InterPro signature-number of occurrence
	      ControlProt (dict): The map Protein-InterPro signature
	      aplha (float): The calibration factor between 0 and 1
	      penal (float): The standard deviation adjustment factor between 0 and 1
	'''
	
	x = SignatureNumber.values()
	Q1 = quartile(x,1); Q3 = quartile(x,3); IQR = Q3-Q1
	data = [b for b in x if (b >= Q1-1.5*IQR and b <= Q3+1.5*IQR)]
	mu = float(sum(data))/len(data)
	ss = sqrt(sum([(x-mu)**2 for x in data])/len(data))

	KeepScores = {}
	Protein = ControlProt.keys()
	count = 0
	stringp = "\rComputing InterPro Scores: %s %d%s"%(40*'.',0,'%')
	nl = len(ControlProt)
	for i in range(len(Protein)):
		print(stringp, end=',')
		count += 1
		origin = Protein[i]
		for j in range(i+1,len(Protein)):
			dest = Protein[j]
			CommonSig = list(set(ControlProt[origin]) & set(ControlProt[dest]))
			if len(CommonSig) != 0: # Meaning that the 2 proteins have commun signature!
				eta = 0.0
				for com in CommonSig:
					eta += min(ControlProt[origin].count(com), ControlProt[dest].count(com))
				f = 0.5*(1 + erf(pow(1.0*eta,alpha)/(penal*ss*sqrt(2))))
				uncert = -f*log2(f)-(1.0-f)*log2(1.0-f) if 0<f<1 else 0.0
				score = 1.0 - uncert
				if score: KeepScores[tuple(sorted((origin, dest)))] = score
				del eta, f, uncert, score
			del dest, CommonSig
		del origin
		perc = (count*40)//nl
		stringp = "\rComputing InterPro Scores: %s %d%s"%(perc*'='+(40-perc)*'.',count*100//nl,'%')
	print("\rComputing InterPro Scores: %s %d%s"%(40*'=',100,'%'))
	print("Process successfully completed ...")
	return KeepScores

if __name__=='__main__':
	pass
