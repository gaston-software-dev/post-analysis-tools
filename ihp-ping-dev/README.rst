==============================
The IHP-PING Python Package
==============================

			THE BASIC README

(See the packge "PDF/RESOURCES" for more detailed information at the following 
links: 
                          - For the IHP-PING package at
http://web.cbio.uct.ac.za/ITGOM/post-analysis-tools/ihp-ping-dev/IHP-PING_Manual_2020.pdf
				                 
where tool administration and usage are also provided)				            

1. Description
--------------

.. See the PDF package documentation for more information on the use of the
   tool and different protein-protein interaction datasets.                  

IHP-PING is a repository of python modules for easing integration of existing 
human protein-protein interaction (PPI) datasets from multiple sources into a 
unified PPI network on-the-fly (in real time), which is stored locally for 
further user applications. This provides a platform that enables the retrieval 
of human PPI datasets stored in different online resources and produced a real 
time up-to-date integrated human PPI network with increased accuracy, confidence 
and coverage. IHP-PING retrieves human PPI datasets from eight different 
resources (Refer to the [PDF manual] 


2. Environment
--------------

The IHP-PING system is composed of one main high level folder: PyPING and one 
main python module, ihppinbuilder.py, which serves as an interface processing 
the human PPI interaction network and storing this network locally for further 
user applications. The PyPING folder includes two Python modules: 
	(1) networkgenerator.py, which downloads PPI datasets requested by a user 
	    and builds the integrated PPI network and 
	(2) sequenceprocessing.py, processing and scoring sequence dataset, which 
	    includes protein sequences and InterPro domains. 

As a package, IHP-PING can be imported in another Python module, however, a 
user can directly produce a unified human PPI network in three logical steps: 
PPI Extraction, Mapping Process and Integration Process, via user interface 
and input processing using a simple single command-line terminal on any 
computer or any operating system running Python. The only required argument is 
the list of human PPI datasets to be incorporated into the unified PPI network, 
with output format parameters according to user preferences. Each user-requested 
database is downloaded from specific uniform resource locations (URLs). 
Eight different resources used include STRING, IntAct, MINT, BioGrid, DIP, HPRD, 
and MPPI-MIPS, with additional interactions predicted from protein sequences 
receive a score according to sequence similarity and shared protein signatures 
computed using an information theory-based scheme.

3. Administration
-----------------

IHP-PING v2.4.1 requires Linux operating system and Python version >= 2.7, 
requiring the installation of the NCBI BLAST software locally when retrieving 
interactions predicted from sequence data. It also requires the selenium Python 
package, as well as chromedriver and chromium-browser, for retrieving the DIP 
dataset. These needs to be installed prior to the use of IHP-PING. It is worth
mentioning that IHP-PING can run in any operating system running Python, but it
has been tested only on Linux machine.

To use IHP-PING, the user needs to download the `tar.gz' file from the main URL
(see point x below) and extract all files as follows: 
::
     tar xzf post-analysis-tools.tar.gz

or alternatively, it can also be retrieved from the github public platform using 
git clone command line as follows.
::
	 git clone https://github.com/gkm-software-dev/post-analysis-tools.git

After downloading and/or uncompressing, move to the appropriate working directory  
where IHP-PING and related commands are executed using the following terminal 
command: 
::
	 cd post-analysis-tools/ihp-ping-dev/


It is worth noting that this package is free to use under GNU General Public Li-
cense. You are free to copy, distribute and display information contained herein, 
provided that it is done with appropriate citation of the tool. Thus, by using 
the IHP-PING package, it is assumed that you have read and accepted the agreement 
provided and that you agreed to be bound to all terms and conditions of this 
agreement. Please, use the following command line to see the package licence::

::	 python setup.py --licence    ::

It is worth mentioning that Python enables easy transitioning of codes between 
computers, rendering IHP-PING portable and an effective framework for human for
retrieving and integrating human PPI networks into a unifed PPI netwok. 

4. Build status
---------------

The main website for the IHP-PING package is 
http://web.cbio.uct.ac.za/ITGOM/post-analysis-tools/ihp-ping-dev where users can 
find essential information about IHP-PING. It is freely downloadable under GNU 
General Public License (GPL), pre-compiled for Linux version and protected by 
copyright laws. Users are free to copy, modify, merge, publish, distribute and 
display informamation contained in the package, provided that it is done with 
appropriate citation of the package and by including the permission notice in all 
copies or substantial portions of the module contained in this package.

It is currently maintained by one member of the core-development team, Gaston K. 
Mazandu <gmazandu@gmail.com, gmazandu@cbio.uct.ac.za, kuzamunu@aims.ac.za, who 
regularly updates the information available in this package and makes every effort 
to ensure the quality of this information.

5. Quick start guide
--------------------
IHP-PING can be processed through one main python module, ihppinbuilder.py, which 
serves as an interface. Get help on how to run IHP-PING through this interface 
module using the following command:
::
	 python ihppinbuilder.py -h
or
::
	 python ihppinbuilder.py --help

From the help option above, IHP-PING is run using the following one line command:
::
 python ihppinbuilder.py -r resources -o outputfolder -i outputProtID -f outputfileformat

As illustrations:
::
     python ihppinbuilder.py  -f csv
     python ihppinbuilder.py  -r stringdb biogrid dip -i genename -f csv2

The first command generates a unified human PPI network derived from all sources 
under consideration currently, and save under the working directory (default) in 
csv format. Similarly, for the second command, only STRING, BioGrid and DIP 
databases are used and the network is saved as a csv2 (semi-column separated 
value) file with the gene name ID system.

Finally, as any python library or package, IHP-PING can be imported and used in 
another Python models. For accessing and learning about different classes of the 
two main modules under PyPING, networkgenerator providing functions for downloding 
and integrating the human PPI network and sequenceprocessing for processing 
sequence datasets (protein sequences and InterPro datasets). Please access the 
python interpreter or the command shell for interactive computing (IPython) and run 
following commands:
::
	 >>> from PyPING import *
	 >>> help(networkgenerator)
	 >>> help(sequenceprocessing)

6. Version history
------------------

- 2.4.1: Initial IHP-PING release in April 2020.
::
	 python setup --version

7. Package URL
--------------

http://web.cbio.uct.ac.za/ITGOM/post-analysis-tools/ihp-ping-dev/
https://github.com/gkm-software-dev/post-analysis-tools
::
	 python setup --url

8. Maintainer
-------------

Gaston K. Mazandu
Email: gaston.mazandu@uct.ac.za, gmazandu@gmail.com, 
       kuzamunu@aims.ac.za

9. Contributors
---------------

Hooper C, Opap K, Makinde FL, Nembaware V, Thomford NE, Chimusa ER, Wonkam A, 
Mulder NJ, Mazandu GK
Emails: HPRCHR001@myuct.ac.za, babuken@gmail.com, funmite@aims.ac.za, 
vnembaware@gmail.com, nicholas.thomford@uct.ac.za, emile.chimusa@uct.ac.za, 
nicola.mulder@uct.ac.za, ambroise.wonkam@uct.ac.za, gaston.mazandu@uct.ac.za
::
	 python setup --classifiers
Classifier: License :: GPL (>= 2)
Classifier: Operating System :: OS Independent, but tested only on Linux
Classifier: Programming Language :: Python :: >= 2.7
Classifier: Topic :: Software Development :: Libraries

Sincerely,

Gaston K. Mazandu
