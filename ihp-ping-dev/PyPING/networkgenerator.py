#!/usr/bin/python
# -*- coding: utf8 -*-

"""
*******************************************************************************
Important Note:
This python file is part of the IHP-PING tool, which is a tool for 
generating an integrated human protein-protein interaction network.

The main website for the IHP-PING package is:
 
http://web.cbio.uct.ac.za/ITGOM/post-analysis-tools/ihp-ping-dev
https://github.com/gkm-software-dev/post-analysis-tools

where users can find essential information about obtaining IHP-PING. It is 
freely downloadable under GNU General Public License (GPL), pre-compiled 
for Linux version and protected by copyright laws. Users are free to copy, 
modify, merge, publish, distribute and display information contained in 
the package, provided that it is done with appropriate citation of the 
package and by including the permission notice in all copies or substantial 
portions of the module contained in this package.

IHP-PING is distributed in the hope that it will be useful, but WITHOUT 
ANY WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED 
TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE 
AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS 
BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY WHETHER IN AN 
ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN 
CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
See <http://www.gnu.org/licenses/>.

This code was written by Gaston K. Mazandu and Christopher Hooper
    E-mail: <gaston.mazandu@uct.ac.za, gmazandu@gmail.com, 
                       kuzamunu@aims.ac.za>
    (c) 2020 under free software (GPL), all rights reserved.\n                  
*******************************************************************************   
"""

from __future__ import print_function
import xml.etree.ElementTree as et
import os, sys, subprocess, gzip, zipfile, tarfile, io, ssl, random, time, shutil, base64, json

try:
	import urllib.request as accessurl
except ImportError:
	import urllib2 as accessurl

try:
	import urllib.parse as urlparse
except ImportError:
	import urllib as urlparse

try:
	from urllib.request import urlretrieve
except ImportError:
	from urllib import urlretrieve

from functools import reduce
from . import __author__, __author_email__, __version__, __license__, __status__
from .sequenceprocessing import *

_STRINGMAP = {}; _GENEIDMAP = {}
_REVIEWEDID = set()
_GENEMAP = {}
_REFSEQ = {}
_REQ = rtvhttp = GetAllScores = None

def progressbar(sstr, dsize = 8192):
	"""This is controlling the downloading process provided that the size
       of the file being downloaded is known.
       
       Arguments:
           sstr (str): The title of the file being downloaded
           dsize (int): The read size
	"""
	global _REQ
	response = accessurl.urlopen(_REQ)
	try: total_size = response.info().getheader('Content-Length').strip()
	except: total_size = 0
	total_size = int(total_size)
	bytes_so_far = 0; data_bytes = ''; chunk = 1
	chunk_size = dsize
	if total_size > 0:
		stringp = "\rDownload %s progress:%s|%d%s"%(sstr,20*'.',0,'%')
		while chunk:
			print(stringp, end=',')
			chunk = response.read(chunk_size)
			bytes_so_far += len(chunk); data_bytes += chunk
			perc = (bytes_so_far*20)//total_size
			stringp = "\rDownload %s progress:%s|%d%s (%d/%d)"%(sstr,perc*'='+(20-perc)*'.',bytes_so_far*100//total_size,'%',bytes_so_far,total_size)
			if bytes_so_far >= total_size: 
				print("\rDownload %s progress:%s|%d%s (%d/%d)"%(sstr,20*'=',100,'%', total_size,total_size), end=',')
		print("\rDownload %s progress:%s|%d%s (%d/%d)"%(sstr,20*'=',100,'%', total_size,total_size))
	else:
		print("\rDownload %s file with unknown size... It may take time!"%(sstr,))
		data_bytes = response.read()
	return data_bytes

def reportbar(bn, bs, ts):
	"""This is controlling the downloading process provided that the size
       of the file being downloaded is unknown.
       
       Arguments:
           bn (int): The block size of the URL target
           bs (int): The read size
           ts (int): The total file size of the URL target. 
	"""
	global rtvhttp
	rsf = bn * bs
	if ts > 0:
		perc = rsf * 20 // ts
		s = "\rDownload %s:%s|%d%s (%d/%d)"%(str(rtvhttp), perc*'='+(20-perc)*'.',rsf*100//ts,'%',rsf,ts)
		sys.stderr.write(s)
		if rsf >= ts: # near the end
			sys.stderr.write("\n")
	else: # total size is unknown
		sys.stderr.write("\rDownload %s %d ..." % (str(rtvhttp), rsf))

def get_uniprotid():
	"""Downloading all the reviewed human UniProt identifiers (IDs) 
	"""
	global _REVIEWEDID, _REQ, _GENEMAP
	url_idrev = "https://www.uniprot.org/uniprot/?query=reviewed:yes+AND+organism:9606&format=tab"
	try:
		_REQ = accessurl.Request(url_idrev)
	except accessurl.HTTPError as error:
		print(error.code)
		print(error.read())
	except accessurl.URLError as error:
		if hasattr(error, 'reason'):
			print('We failed to reach a server.')
			print('Reason: ', error.reason)
		elif hasattr(error, 'code'):
			print('The server couldn\'t fulfill the request.')
			print('Error code: ', error.code)
	else:
		response = progressbar("UniProt reviewed IDs")
		response = response.decode().split('\n')
		for line in response[1:]:
			tline = line.strip()
			if not tline: continue
			tline = [s.strip() for s in tline.split('\t')]
			_REVIEWEDID.add(tline[0])
			try: _GENEMAP[tline[0]] = tline[4].split()[0].strip()
			except: pass 

def id_select(RequiredID):
	"""Selecting all the identifier (ID) mapping required from data provided to UniProt IDs
	
	   Arguments:
	      RequiredID (list): List of the dataset requiring ID mapping. 
	"""
	global _STRINGMAP, _GENEIDMAP, _REFSEQ, rtvhttp
	rtvhttp_ref = "Request %s - UniProt ID map"
	dictmap = dict(stringdb = "https://www.uniprot.org/uniprot/?query=reviewed:yes+AND+organism:9606+AND+database:string&format=tab&compress=yes&columns=id,database(string)", biogrid = "https://www.uniprot.org/uniprot/?query=reviewed:yes+AND+organism:9606+AND+database:geneid&format=tab&compress=yes&columns=id,database(geneid)", hprd = "https://www.uniprot.org/uniprot/?query=reviewed:yes+AND+organism:9606+AND+database:refseq&format=tab&compress=yes&columns=id,database(refseq)")
	Req = {}
	try:
		for r in RequiredID:
			Req[r] = accessurl.Request(dictmap[r])
	except accessurl.HTTPError as error:
		print(error.code)
		print(error.read())
	except accessurl.URLError as error:
		if hasattr(error, 'reason'):
			print('We failed to reach a server.')
			print('Reason: ', error.reason)
		elif hasattr(error, 'code'):
			print('The server couldn\'t fulfill the request.')
			print('Error code: ', error.code)
	else:
		FileNames = {}
		try:
			for r in RequiredID:
				FileNames[r] = "%s_idmapping%d.dat.gz"%(r, random.randint(0,100000))
				rtvhttp = rtvhttp_ref%(r.upper(),)
				urlretrieve(dictmap[r], FileNames[r], reportbar)
				print()
		except IOError as error:
			print(error)
		else:
			SetData = dict(stringdb = _STRINGMAP, biogrid = _GENEIDMAP, hprd = _REFSEQ)
			for idmap in RequiredID:
				with gzip.open(FileNames[idmap], "rb") as data:
					for line in data:
						tline = line.strip()
						if not tline: continue
						tline = [s.strip() for s in tline.split()]
						for maps in tline[1].split(';'):
							try:
								tmp = maps.split()[0].strip()
								if tmp: SetData[idmap][tmp] = tline[0].strip()
							except: pass
			for r in RequiredID:
				if os.path.isfile(FileNames[r]): os.remove(FileNames[r])
			del SetData
			return

def get_protseq():
	"""Downloading human protein sequence dataset and processing
	   Predicting interactions from protein sequence
	"""
	url_protseq = "https://www.uniprot.org/uniprot/?query=reviewed:yes+AND+organism:9606&format=fasta"
	global _REQ, GetAllScores
	_REQ = None
	try:
		_REQ = accessurl.Request(url_protseq)
	except accessurl.HTTPError as error:
		print(error.code)
		print(error.read())
	except accessurl.URLError as error:
		if hasattr(error, 'reason'):
			print('We failed to reach a server.')
			print('Reason: ', error.reason)
		elif hasattr(error, 'code'):
			print('The server couldn\'t fulfill the request.')
			print('Error code: ', error.code)
	else:
		DiffFiles = set()
		response = progressbar("Protein Sequences")
		FileNamePrefix = "Sequence%d"%(random.randint(0,100000),)
		fw = open(FileNamePrefix+'.fasta', 'w'); DiffFiles.add(FileNamePrefix+'.fasta')
		fw.write(response); fw.close()
		print("\rNow building local blast database and perform blastall ...")

		ccb = 0
		if sys.platform.startswith('linux') or sys.platform.startswith('darwin'):
			locate_app = 'which'
		elif sys.platform.startswith('win32'): locate_app = 'where'

		if not subprocess.call([locate_app, 'makeblastdb']): ccb = 1
		elif not subprocess.call([locate_app, 'formatdb']): ccb = 2

		if ccb==1:
			os.popen("makeblastdb -in %s.fasta -dbtype prot -out %s"%(FileNamePrefix, FileNamePrefix))
			os.popen('blastp -query %s.fasta -db %s -out %sBlast.txt -outfmt "6 qseqid sseqid bitscore" -evalue 0.01'%(FileNamePrefix, FileNamePrefix, FileNamePrefix))
			DiffFiles.add(FileNamePrefix+'.phr'); DiffFiles.add(FileNamePrefix+'.pin') 
			DiffFiles.add(FileNamePrefix+'.psq') #; DiffFiles.add("formatdb.log")
			os.system('clear')
			print('\n'+74*'*')
		elif ccb==2:
			os.popen("formatdb -i %s.fasta -p T -n %s"%(FileNamePrefix, FileNamePrefix))
			os.popen('blastall -p blastp -d %s -i %s.fasta -o %sBlast.txt -m 8 -e 0.01'%(FileNamePrefix, FileNamePrefix, FileNamePrefix))
			DiffFiles.add(FileNamePrefix+'.phr'); DiffFiles.add(FileNamePrefix+'.pin') 
			DiffFiles.add(FileNamePrefix+'.psq') #; DiffFiles.add("formatdb.log")
			os.system('clear')
			print('\n'+74*'*')
		else:
			print("Blast is not installed which is required when using sequence data. Please install blast and try again, or remove sequence from the resources parameter when running IHP-PING.")
			sys.exit() # Message that Blast is not installed and exit
		print("\rNow computing sequence blast scores, this may take time ...")

		GetAllScores = computeScoreBlast(FileNamePrefix+'Blast.txt')
		DiffFiles.add(FileNamePrefix+'Blast.txt')
		# Deleting all files generated to produce Sequence Blast Scores
		for fp in DiffFiles:
			if os.path.isfile(fp): os.remove(fp)
		del DiffFiles
		return

def get_protintepro():
	"""Downloading reviewed human InterPro signature- UniProt identifier mapping and
	   processing the InterPro data-Predicting interactions from InterPro dignature
	"""
	global _REQ, GetAllScores
	url_interpro = "https://www.uniprot.org/uniprot/?query=reviewed:yes+AND+organism:9606+AND+database:interpro&format=tab&compress=yes&columns=id,database(interpro)"
	try:
		_REQ = accessurl.Request(url_interpro)
	except accessurl.HTTPError as error:
		print(error.code)
		print(error.read())
	except accessurl.URLError as error:
		if hasattr(error, 'reason'):
			print('We failed to reach a server.')
			print('Reason: ', error.reason)
		elif hasattr(error, 'code'):
			print('The server couldn\'t fulfill the request.')
			print('Error code: ', error.code)
	else:
		response = progressbar("InterPro")
		response = gzip.GzipFile(fileobj = io.BytesIO(response))
		Sign = {}; ProtList = {}
		for line in response:
			sline = line.strip()
			if not sline: continue
			sline = sline.split('\t')
			ProteinId, InterProSig = sline[0].strip(), sline[1].strip()
			if InterProSig:
				ProtList[ProteinId] = [s.strip() for s in InterProSig.split(';')[:-1]]
				for sig in ProtList[ProteinId]:
					try: Sign[InterProSig] += 1
					except: Sign[InterProSig] = 1
			del ProteinId, InterProSig
	print("\nNow processing InterPro interaction scores. This may take time ...")

	InterProScores = computeFamilyScore(Sign, ProtList)
	for key in InterProScores:
		if key in GetAllScores: # There was a bug here! InterPro instead of All
			GetAllScores[key] = 1.0-(1.0-GetAllScores[key])*(1.0-InterProScores[key])
		else: GetAllScores[key] = InterProScores[key]
	del InterProScores
	return

def get_mint(ns, ind = 0):
	"""Retrieving interactions from the MINT online database
	   
	   Arguments:
	       ns (int): The current number of Columns
	       ind (ind): The column position of interaction score
	"""
	global _REQ, _REVIEWEDID, GetAllScores
	url_mint = "http://www.ebi.ac.uk/Tools/webservices/psicquic/mint/webservices/current/search/query/species:human"
	try:
		_REQ = accessurl.Request(url_mint)
	except accessurl.HTTPError as error:
		print(error.code)
		print(error.read())
	except accessurl.URLError as error:
		if hasattr(error, 'reason'):
			print('We failed to reach a server.')
			print('Reason: ', error.reason)
		elif hasattr(error, 'code'):
			print('The server couldn\'t fulfill the request.')
			print('Error code: ', error.code)
	else:
		response = progressbar("MINT")
		print("Now reading MINT PPI file of size %s bytes. This may take time ..."%(sys.getsizeof(response),))
		count = 0; response = io.BytesIO(response)
		for line in response:
			tline = line.strip()
			count += 1
			print("\rCurrently %d lines are already read ..."%(count,), end=',')
			if not tline: continue
			tline = [s.strip() for s in tline.split('\t')]
			if not tline[0].startswith('uniprot') or not tline[1].startswith('uniprot') or not 'intact' in tline[-1]: continue
			key = tline[0].split(':')[-1].split('-')[0].strip(), tline[1].split(':')[-1].split('-')[0].strip()
			if key[0] in _REVIEWEDID and key[1] in _REVIEWEDID:
				key = tuple(sorted(key))
				try: GetAllScores[key][ind] = float(tline[-1].split(':')[-1])
				except:
					GetAllScores[key] = ns*[0]
					GetAllScores[key][ind] = float(tline[-1].split(':')[-1])
		print("\nThe MINT extraction process successfully completed ...")
		return

def get_dip(ns, ind = 0):
	"""Retrieving interactions from the DIP online database. The user
	   is required to have a Credential Identifier from the DIP team 
	   for free access to non-Profit organization.
	   
	   It also requires to have python-selenium package and Chrome-
	   drivers to be installed.
	     
	   Arguments:
	       ns (int): The current number of Columns
	       ind (ind): The column position of interaction score
	"""
	try:
		from selenium import webdriver
		from selenium.webdriver.common.keys import Keys
		from selenium.common.exceptions import TimeoutException
	except ImportError:
		print("\n==========\nInstall selenium package appropriate to your python version,\nas well as chromedriver and chromium-browser.\n\nPlease install packages required and try again ...\nExiting  IHP-PING now ...\n==========\n")
		sys.exit()
		
	global _REVIEWEDID, GetAllScores
	page_load_timeout = 30
	query_url = "https://dip.doe-mbi.ucla.edu/dip/File.cgi?FN=2017/tab25/Hsapi20170205.txt.gz"
	
	options = webdriver.ChromeOptions()
	options.add_argument("--headless")
	
	driver = webdriver.Chrome(chrome_options = options)
	driver.set_page_load_timeout(page_load_timeout)
	
	login_url = 'https://dip.doe-mbi.ucla.edu/dip/Login.cgi'
	driver.get(login_url)
	
	catchuser = driver.find_element_by_name("login")
	catchaccs = driver.find_element_by_name("pass")

	catchuser.send_keys("WoodenTable")
	catchaccs.send_keys("Dip")
	driver.find_element_by_name("Login").click()
	
	driver.execute_script("""
         window.file_contents = null;
         var xhr = new XMLHttpRequest();
         xhr.responseType = 'blob';
         xhr.onload = function() {
             var reader  = new FileReader();
             reader.onloadend = function() {
                 window.file_contents = reader.result;
             };
             reader.readAsDataURL(xhr.response);
         };
         xhr.open('GET', %(download_url)s);
         xhr.send();
     """.replace('\r\n', ' ').replace('\r', ' ').replace('\n', ' ') % {
    'download_url': json.dumps(query_url),
    })

	print('Downloading DIP dataset... It may take time! Please')
	downloaded_file = None
	while downloaded_file is None: # Returns the file retrieved base64 encoded (perfect for downloading binary)
		downloaded_file = driver.execute_script('return (window.file_contents !== null ? window.file_contents.split(\',\')[1] : null);')

	FileName = "DIP%d.txt.gz"%(random.randint(0,100000),)
	fp = open(FileName, 'wb')
	fp.write(base64.b64decode(downloaded_file))
	fp.close(); driver.close() # close the downloaded and the web browser, or it'll persist after python exits.

	with gzip.open(FileName, "rb") as fp:
		count = 0
		for line in fp:
			tline = line.decode().strip()
			count += 1
			print("\rCurrently %d lines are already read ..."%(count,), end=',')
			if not tline: continue
			tline = [s.strip() for s in tline.split('\t')]
			if not 'uniprot' in tline[0] or not 'uniprot' in tline[1]: continue
			key = tline[0].split(':')[-1].strip(), tline[1].split(':')[-1].strip()
			if key[0] in _REVIEWEDID and key[1] in _REVIEWEDID:
				key = tuple(sorted(key))
				try: GetAllScores[key][ind] = 0.7
				except:
					GetAllScores[key] = ns*[0]
					GetAllScores[key][ind] = 0.7
	print("\nThe DIP extraction process successfully completed ...")
	if os.path.isfile(FileName): os.remove(FileName)
	return

def get_mips(ns, ind = 0): # This is xml file!!!
	"""Retrieving interactions from the MPPI-MIPS online database
	   
	   Arguments:
	       ns (int): The current number of Columns
	       ind (ind): The column position of interaction score
	"""
	global _REQ, _REVIEWEDID, GetAllScores
	url_mips = "http://mips.helmholtz-muenchen.de/proj/ppi/data/mppi.gz"
	try:
		_REQ = accessurl.Request(url_mips)
	except accessurl.HTTPError as error:
		print(error.code)
		print(error.read())
	except accessurl.URLError as error:
		if hasattr(error, 'reason'):
			print('We failed to reach a server.')
			print('Reason: ', error.reason)
		elif hasattr(error, 'code'):
			print('The server couldn\'t fulfill the request.')
			print('Error code: ', error.code)
	else:
		response = progressbar("MIPS")
		print("Reading MIPS gz file of size %s bytes. This may take time ..."%(sys.getsizeof(response),))
		response = io.BytesIO(response)
		response = gzip.GzipFile(fileobj = response)

		tree = et.parse(response)
		root = tree.getroot(); count = 0
		human = root.findall(".//*[@ncbiTaxId='9606']/../../..")

		for interaction in human:
			prot_ids = interaction.findall(".//{net:sf:psidev:mi}primaryRef")
			count += 1
			print("\rCurrently %d lines are already read ..."%(count,), end=',')
			key = [prot.attrib['id'] for prot in prot_ids if 'id' in prot.attrib]
			if len(key) != 2: continue
			key = key[0].strip(), key[1].strip()
			if key[0] in _REVIEWEDID and key[1] in _REVIEWEDID:
				key = tuple(sorted(key))
				try: GetAllScores[key][ind] = 0.6
				except:
					GetAllScores[key] = ns*[0]
					GetAllScores[key][ind] = 0.6
		print("\nThe MIPS extraction process successfully completed ...")
		return

def get_stringdb(ns, ind = 0):
	"""Retrieving interactions from the STRING online database
	   
	   Arguments:
	       ns (int): The current number of Columns
	       ind (ind): The column position of interaction score
	"""
	global _STRINGMAP, _REVIEWEDID, _REQ, GetAllScores 
	url_stringdb = "https://stringdb-static.org/download/protein.links.v11.0/9606.protein.links.v11.0.txt.gz"
	try:
		_REQ = accessurl.Request(url_stringdb)
		_REQ.add_header("User-Agent","Mozilla/5.0 (X11; U; Linux i686) Gecko/20071127 Firefox/2.0.0.11")
	except accessurl.HTTPError as error:
		print(error.code)
		print(error.read())
	except accessurl.URLError as error:
		if hasattr(error, 'reason'):
			print('We failed to reach a server.')
			print('Reason: ', error.reason)
		elif hasattr(error, 'code'):
			print('The server couldn\'t fulfill the request.')
			print('Error code: ', error.code)
	else:
		response = progressbar("STRING")
		print("Reading STRING zip file of size %s bytes. This may take time ..."%(sys.getsizeof(response),))
		response = gzip.GzipFile(fileobj = io.BytesIO(response))
		response.readline(); count = 0
		for line in response:
			sline = line.strip()
			count += 1
			print("\rCurrently %d lines are already read ..."%(count,), end=',')
			if not sline: continue
			sline = sline.split()
			Prot1, Prot2, Score = str(sline[0].decode()), str(sline[1].decode()), str(sline[-1].decode())
			if Prot1 in _STRINGMAP and Prot2 in _STRINGMAP:
				if _STRINGMAP[Prot1] in _REVIEWEDID and _STRINGMAP[Prot1] in _REVIEWEDID:
					key = tuple(sorted([_STRINGMAP[Prot1], _STRINGMAP[Prot2]]))
					try: GetAllScores[key][ind] = float(Score)/1000
					except:
						GetAllScores[key] = ns*[0]
						GetAllScores[key][ind] = float(Score)/1000
			del Prot1, Prot2, Score
		print("\nThe STRING extraction process successfully completed ...")
		return

def get_hprd(ns, ind = 0):
	"""Retrieving interactions from the HPRD online database
	   
	   Arguments:
	       ns (int): The current number of Columns
	       ind (ind): The column position of interaction score
	"""
	global _REFSEQ, _REVIEWEDID, _REQ, GetAllScores
	url_hprd = "http://www.hprd.org/RELEASE9/HPRD_Release9_041310.tar.gz"
	try:
		_REQ = accessurl.Request(url_hprd)
	except accessurl.HTTPError as error:
		print(error.code)
		print(error.read())
	except URLError as error:
		if hasattr(error, 'reason'):
			print('We failed to reach a server.')
			print('Reason: ', error.reason)
		elif hasattr(error, 'code'):
			print('The server couldn\'t fulfill the request.')
			print('Error code: ', error.code)
	else:
		response = progressbar("HPRD")
		print("Reading HPRD tar.gz file of size %s bytes. This may take time ..."%(sys.getsizeof(response),))
		fdir = io.BytesIO(response)
		tar = tarfile.open(mode="r:gz", fileobj = fdir)
		FolderName = "HprdDataset%d"%(random.randint(0,100000),)
		os.mkdir(FolderName)
		tar.extractall(FolderName)
		datasets = tar.getnames(); count = 0
		msize = [os.path.getsize('./%s/%s'%(FolderName, d)) for d in datasets]
		fp = open('./%s/%s'%(FolderName, datasets[msize.index(max(msize))]))
		for line in fp:
			sline = line.strip()
			count += 1
			print("\rCurrently %d lines are already read ..."%(count,), end=',')
			if not sline: continue
			sline = sline.split()
			Prot1, Prot2 = str(sline[2].decode()), str(sline[5].decode()) #0-3GeneSymbol
			Prot1, Prot2 = Prot1.strip(), Prot2.strip()
			if Prot1 in _REFSEQ and Prot2 in _REFSEQ:
				if _REFSEQ[Prot1] in _REVIEWEDID and _REFSEQ[Prot1] in _REVIEWEDID:
					dline = str(sline[-2].decode()), str(sline[-1].decode())
					dscore = len([s.strip() for s in dline[0].strip().split(';')]) + len(set([s.strip() for s in dline[1].strip().split(',')]))
					key = tuple(sorted([_REFSEQ[Prot1], _REFSEQ[Prot2]]))
					try: GetAllScores[key][ind] = 1.0 - 1.0/dscore
					except:
						GetAllScores[key] = ns*[0]
						GetAllScores[key][ind] = 1.0 - 1.0/dscore
			del Prot1, Prot2
		fp.close()
		shutil.rmtree(FolderName, ignore_errors=True)
		print("\nThe HPRD extraction process successfully completed ...")
		return

def get_intact(ns, ind = 0):
	"""Retrieving interactions from the IntAct online database
	   
	   Arguments:
	       ns (int): The current number of Columns
	       ind (ind): The column position of interaction score
	"""
	global _REVIEWEDID, _REQ, GetAllScores
	url_intact = "ftp://ftp.ebi.ac.uk/pub/databases/intact/current/psimitab/intact.zip"
	try:
		_REQ = accessurl.Request(url_intact)
	except accessurl.HTTPError as error:
		print(error.code)
		print(error.read())
	except accessurl.URLError as error:
		if hasattr(error, 'reason'):
			print('We failed to reach a server.')
			print('Reason: ', error.reason)
		elif hasattr(error, 'code'):
			print('The server couldn\'t fulfill the request.')
			print('Error code: ', error.code)
	else:
		response = progressbar("IntAct")
		print("Reading IntAct zip file of size %s bytes. This may take time ..."%(sys.getsizeof(response),))
		response = io.BytesIO(response)
		try:
			response = zipfile.ZipFile(response)
		except zipfile.BadZipFile:
			print('Error: The zip file is corrupted')
		else:
			for name in response.namelist():
				fpath = os.getcwd()+'/IntactDataset%d.txt'%(random.randint(0,100000),)
				with open(fpath, 'w') as fw:
					fw.write(response.read(name))
				fp = open(fpath, 'r'); count = 0; DataInter = {}
				for line in fp:
					ligne = line.strip()
					count += 1
					print("\rCurrently %d lines are already read ..."%(count,), end=',') 
					if not ligne: continue  
					ligne = ligne.split('\t')
					Prot1, Prot2 = ligne[0].strip(), ligne[1].strip()
					nat = ligne[35].strip() #positive (if false) or negative (if true)
					org =  ligne[28].strip()
					if not org: continue
					if not (Prot1.startswith("uniprotkb") and Prot2.startswith("uniprotkb") and nat=='false'): continue
					Prot1, Prot2 = Prot1.split(':')[-1].strip().split('-')[0].strip(), Prot2.split(':')[-1].strip().split('-')[0].strip()
					if Prot1 in _REVIEWEDID and Prot2 in _REVIEWEDID:
						key = tuple(sorted([Prot1, Prot2]))
						try: GetAllScores[key][ind] = 0.6
						except:
							GetAllScores[key] = ns*[0]
							GetAllScores[key][ind] = 0.6
				if os.path.isfile(fpath): os.remove(fpath)
				print("\nThe IntAct extraction process successfully completed ...",)
				return

def get_biogrid(ns, ind = 0):
	"""Retrieving interactions from the BioGRID online database
	   
	   Arguments:
	       ns (int): The current number of Columns
	       ind (ind): The column position of interaction score
	"""
	global _GENEIDMAP, _REVIEWEDID, _REQ, GetAllScores
	url_biogrid = "https://downloads.thebiogrid.org/Download/BioGRID/Release-Archive/BIOGRID-3.5.174/BIOGRID-ORGANISM-3.5.174.tab2.zip" # tab2 or mitab
	try:
		_REQ = accessurl.Request(url_biogrid)
	except accessurl.HTTPError as error:
		print(error.code)
		print(error.read())
	except accessurl.URLError as error:
		if hasattr(error, 'reason'):
			print('We failed to reach a server.')
			print('Reason: ', error.reason)
		elif hasattr(error, 'code'):
			print('The server couldn\'t fulfill the request.')
			print('Error code: ', error.code)
	else:
		response = progressbar("BioGRID")
		print("Reading BioGRID zip file of size %s bytes. This may take time ..."%(sys.getsizeof(response),))
		response = io.BytesIO(response)
		try:
			response = zipfile.ZipFile(response)
		except zipfile.BadZipFile:
			print('Error: The zip file is corrupted')
		else:
			for name in response.namelist():
				if 'Homo_sapiens' in name:
					fpath = os.getcwd()+'/BioGridDataset%d.txt'%(random.randint(0,100000),)
					with open(fpath, 'w') as fw:
						fw.write(response.read(name))
					fp = open(fpath, 'r'); count = 0; DataInter = {}
					for line in fp:
						ligne = line.strip()
						count += 1
						print("\rCurrently %d lines are already read ..."%(count,), end=',') 
						if not ligne: continue  
						ligne = [x.strip() for x in ligne.split('\t')]
						Prot1, Prot2 = ligne[1], ligne[2]
						if Prot1 in _GENEIDMAP and Prot2 in _GENEIDMAP:
							if _GENEIDMAP[Prot1] in _REVIEWEDID and _GENEIDMAP[Prot2] in _REVIEWEDID:
								key = tuple(sorted([_GENEIDMAP[Prot1], _GENEIDMAP[Prot2]]))
								try: GetAllScores[key][ind] = 0.6
								except:
									GetAllScores[key] = ns*[0]
									GetAllScores[key][ind] = 0.6
					if os.path.isfile(fpath): os.remove(fpath)
					print("\nThe BioGRID extraction process successfully completed ...",)
					break
			return
			
def is_valid_file(parser, arg):
	"""
	Check if a valid path to a file has been passed.

	Parameters
	parser, file path
	return the file path or fail with an error"""
	arg = os.path.abspath(arg)
	return arg if os.path.isdir(arg) else parser.error("The folder %s does not exist!" % arg)

def main():
	"""Integrating all human PPI datasets requested into a unified network
	"""
	global GetAllScores, _GENEMAP, _REVIEWEDID
	from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
	parser = ArgumentParser(description = __doc__, formatter_class = ArgumentDefaultsHelpFormatter)
	parser.add_argument("-r", "--resources", dest="databases", type = str, nargs = '+', help="Database to be integrated or considered ", default = 'all', metavar="list")
	parser.add_argument("-o", "--dir", dest = "folder", type = lambda x: is_valid_file(parser, x), default = os.getcwd(), help="Folder which will contain the PPI produced ", metavar="FILE")
	parser.add_argument("-i", "--identifiers", dest="identifiers", type = str, help="Identifier outputs: uniprot or genename", default = "uniprot", metavar = "str")
	parser.add_argument("-f", "--outformat", dest="outformat", type = str, help="Output format tsv, csv or csv2", default = "tsv", metavar = "str")
	
	args = parser.parse_args()
	
	# Providing system and platform requirements
	print('\n'+74*'*')
	print(('IHP-PING %s'%(__version__)+', developed by %s'%(__author__)).center(71))
	print(('E-mail: %s'%(__author_email__)).center(71))
	print('\n'+74*'*'+'\n*'+'These are requirements for running Integrated Human PPI Generator tool'.center(71)+' *\n'+74*'*')
	print('\nPlatform: Python version >= 2.7.x\n\nOS      : Independent, but only tested on Linux')
	print("\nThe current IHP-PING version retrieves data from eight following resources\nwith symbols used as shown:\n")
	print("  stringdb: For the STRING database\n"+"  sequence: For interactions inferred from sequence data from UniProt and\n"+12*" "+"InterPro databases\n"+"  mint    : For the MINT database\n"+"  dip     : For the DIP database\n"+"  biogrid : For the BioGrid database\n"+"  mips    : For the MIPS database\n"+"  intact  : For the INtact database\n"+"  hprd    : For the HPRD database\n"+"  all     : To produce integrated PPI network from all sources above.\n")
	print("Currently, the IHP-PING tool requires the math Python package to be in for\ninferring PPI from sequence data. For this, also ensure that local [BLAST]\nsoftware is installed.\nPlease refer to the manual for more information on these resources\n")
	print("Furthermore, IHP-PING also requires the selenium Python package for infer-\nring DIP dataset, as well as chromedriver and chromium-browser. So, ensure\nthe selenium package, chromedriver and chromium-browser are installed.\nAlso bear in mind that you would need a credential Identifier from the DIP\nTeam for free access for non-Profit organization\nPlease refer to the manual for more information on these resources\n")
	print('\nPlease ensure that different packages are installed and appropriate for\n your python version\n')
	print('\nOUTPUT FORMAT: IHP-PING provides a certain flexibility in the output \nformat. Choose between csv, csv2 or tsv\n')
	print(74*'*')
	print('\nEnter 1 to continue and 2 to exit\n')
	
	if sys.version_info.major < 3: raw_pass = raw_input
	else: raw_pass = input
	while True:
		a = raw_pass('> ')
		try: 
			a = int(a)
			if a in [1, 2]: break
			print('Please enter 1 to continue and 2 to exit')
		except: print('Please enter 1 to continue and 2 to exit')
	if a==2:
		print("\nExit IHP-PING. Thanks for using this tool! You can try again at any time.\n")
		print(74*'*'+'\n')
		sys.exit(2)
	if args.outformat=='tsv': frmt, sprt = 'tsv', '\t'
	elif args.outformat=='csv': frmt, sprt = 'csv', ','
	elif args.outformat=='csv2': frmt, sprt = 'csv2', ';'
	else: # Format not implemented
		print('\nFILE-Output File Format Error: Only three formats are currntly considered-\n\n\ttsv: For tab separated value file format, \n\tcsv: For comma separated value file format and \n\tcsv2: For semi-column separated value file format.\n\nPlease refer to the tool manual, fix this issue and try again ...\n')
		print(74*'*'+'\n')
		sys.exit(3)
	
	if isinstance(args.databases, str): args.databases = ['all']
		 
	allr = ['sequence', 'intact', 'stringdb', 'biogrid', 'dip', 'mint', 'hprd', 'mips']
	diffsource = dict([(s.strip().lower(), s) for s in args.databases])
	if 'all' in diffsource:
		existsource = [(s, allr.index(s)) for s in allr]
	else: 
		existsource = [(s, allr.index(s)) for s in diffsource if s in allr]
		existsource = sorted(existsource, key = lambda x: x[1])
	UserReq = [x[0] for x in existsource]

	nexist = set([s for s in diffsource])-set(UserReq + ['all'])

	# Check whether math package required for sequence exists
	locate_app = 'which' if sys.platform.startswith('linux') or sys.platform.startswith('darwin') else 'where'
	if 'sequence' in UserReq:
		ccb = 0
		if not subprocess.call([locate_app, 'makeblastdb']): ccb = 1
		elif not subprocess.call([locate_app, 'formatdb']): ccb = 2
		if not ccb:
			print('\n'+74*'*'+'\n'+"The local NCBI BLAST is required, when using SEQUENCE DATA, but has not \nbeen found. Please install local NCBI BLAST and try again.\n\nExecution cannot be pursued, now exiting ...\n"+74*'*'+'\n')
			sys.exit(2)
	if 'dip' in UserReq:
		import imp
		spec = None
		try:
			imp.find_module('selenium')
			spec = True
		except ImportError: pass
		if spec is None:
			print('\n'+74*'*'+'\n'+"selenium package under Python has not been found. Please ensure that the package is")
			print("installed and try again. Execution cannot be pursued, now exiting ...\n"+74*'*'+'\n')
			sys.exit(2)
		
		ccb = False if subprocess.call([locate_app, 'chromedriver']) or subprocess.call([locate_app, 'chromium-browser']) else True
		if not ccb:
			print('\n'+74*'*'+'\n'+"Retrieving the DIP dataset requires to have chromedriver and chromium-browser.\nPlease install chrome drivers and chromium-browser and try again ...\n\nExecution cannot be pursued, now exiting ...\n"+74*'*'+'\n')
			sys.exit(2)
		
	if nexist:
		print('\n'+74*'*'+'\n'+"WARNING: Note that the current IHP-PING version retrieves data from eight\nfollowing resources and symbols used\n")
		print("  stringdb: For the STRING database\n"+"  sequence: For interactions inferred from sequence data from UniProt and\n"+12*" "+"InterPro databases\n"+"  mint    : For the MINT database\n"+"  dip     : For the DIP database\n"+"  biogrid : For the BioGrid database\n"+"  mips    : For the MIPS database\n"+"  intact  : For the INtact database\n"+"  hprd    : For the HPRD database\n"+"  all     : To produce integrated PPI network from all sources above.\n")
		print("Please refer to the manual for more information on these resources\n")
		print("Thus, %s provided are not in the range and will not be considered in\nthe PPI file\n"%(', '.join([diffsource[a] for a in nexist]),)+74*'*'+'\n')
	
	if not existsource: # There is no user resource providedd --> tell user
		print("Since no valid resource has been provided, Execution cannot be pursued,\nnow exiting ...\n"+74*'*'+'\n') 
		sys.exit(3)

	print("\nThanks for choosing IHP-PING. Start processing on %s"%str(time.asctime(time.localtime())))
	now = time.time()

	treq = UserReq[:]
	if 'sequence' in treq: # This is the longest process
		treq.pop(treq.index('sequence'))
		print("\nWarning: Retrieving interactions from sequence may go over 2 hours.\nAdditionally, the IHP-PING tool requires the local [BLAST] software to be\ninstalled.\n\nEnter 1 to continue and 2 to exit ...\n")
		while True:
			a = raw_pass('> ')
			try: 
				a = int(a)
				if a in [1,2]: break
				print('Please enter 1 to continue and 2 to exit')
			except: print('Please enter 1 to continue and 2 to exit')
		if a==2:
			print("\nExit IHP-PING, please you can try again without sequence data option...\n")
			print(74*'*'+'\n')
			sys.exit(2)	

	if len(treq) >= 1:
		get_uniprotid()
		reqidmap = set(UserReq) & set(['stringdb', 'biogrid', 'hprd']) # need ID mapping
		if reqidmap: id_select(reqidmap)
	
	if args.identifiers == 'uniprot': _GENEMAP = dict([(p, p) for p in _REVIEWEDID])
	elif args.identifiers == 'genename': pass
	else: #ID not considered
		print('\nIdentifier (ID) System Error: Only two ID system are currntly considered-\n\n\tuniprot: For the UniProt ID system and, \n\tgenename: For the gene name ID system.\n\nPlease refer to the tool manual, fix this issue and try again ...\n')
		print(74*'*'+'\n')
		sys.exit(3)

	nowseq = time.time()
	print('\n*******')
	n = len(UserReq)
	if UserReq[0]=='sequence':
		get_protseq()
		get_protintepro()
		for key in GetAllScores:
			GetAllScores[key] = [GetAllScores[key]] + (n-1)*[0]
		print('\n'+74*'*')
		tt = time.time()-nowseq
		nh = tt//3600; rs = tt-nh*3600
		nm = rs//60; rs -= nm*60 
		print("Time elapsed in retrieving interactions from sequence data is %d:%02d:%02d"%(nh, nm, rs))
	else:
		GetAllScores = {}
		extractdata = eval("get_%s"%(UserReq[0],))
		extractdata(n)

	tind = 1
	for res in UserReq[1:]:
		extractdata = eval("get_%s"%(res,))
		print('\n*******')
		extractdata(n, tind)
		tind += 1

	print('\n*******')
	nn = len(GetAllScores)
	print("Downloading all %d interactions into a file. This may take time ..."%(nn,))
	if n==len(allr):
		InteractionFile = 'AllSourceInteractionFile%d.%s'%(random.randint(0,100000),frmt)
	else:
		InteractionFile = '%s%d.%s'%(''.join([s.capitalize() for s in UserReq]), random.randint(0,100000), frmt)

	argfile = os.path.abspath('/'.join([args.folder, InteractionFile]))
	fp = open(argfile,'w'); count = 0
	
	if n==1:
		fp.write('#Prot1%sProt2%s%s\n'%(sprt, sprt, UserReq[0]))
		for p in GetAllScores:
			try: 
				if _GENEMAP[p[0]]!=_GENEMAP[p[1]]:
					fp.write(_GENEMAP[p[0]]+'%s'%(sprt,)+_GENEMAP[p[1]]+'%s%.5f\n'%(sprt, GetAllScores[p][0]))
			except: pass
			count += 1
			perc = (count*20)//nn
			print("\rDownloading interactions progress:%s|%d%s (%d/%d)"%(perc*'='+(20-perc)*'.',count*100//nn,'%', count, nn), end=',')
	else:
		fp.write('#Prot1%sProt2%s'%(sprt, sprt)+sprt.join(UserReq[i] for i in range(n))+'%sScore\n'%(sprt,))
		for p in GetAllScores:
			try: 
				if _GENEMAP[p[0]]!=_GENEMAP[p[1]]:
					fp.write(_GENEMAP[p[0]]+'%s'%(sprt,)+_GENEMAP[p[1]]+'%s'%(sprt,)+sprt.join(['%.5f'%(s,) for s in GetAllScores[p]])+'%s%.5f\n'%(sprt, 1.0-reduce(lambda x, y: x*y,[1.0-d for d in GetAllScores[p]])))
			except: pass
			count += 1
			perc = (count*20)//nn
			print("\rDownloading interactions progress:%s|%d%s (%d/%d)"%(perc*'='+(20-perc)*'.',count*100//nn,'%', count, nn), end=',')
	print("\rDownloading interactions progress:%s|%d%s (%d/%d)"%(20*'=',100,'%',nn,nn))
	fp.close()
	print("Total number of interactions in the file is:", nn)
	print("Interactions are reported in the file: %s"%(argfile,))
	print("Processing accomplished on %s. Thanks!"%str(time.asctime(time.localtime())))
	tt = time.time()-now
	nh = tt//3600; rs = tt-nh*3600
	nm = rs//60; rs -= nm*60 
	print("Total time elapsed is approximately %d:%02d:%02d"%(nh, nm, rs))
	print(74*'*'+'\n')

