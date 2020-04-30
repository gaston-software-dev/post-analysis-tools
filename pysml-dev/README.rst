==============================
The Python SSM library
==============================

			THE BASIC README

(See the library "PDF/RESOURCES" for more detailed information at the following 
links: 
                          - For the PySML library at
http://web.cbio.uct.ac.za/ITGOM/post-analysis-tools/pysml-dev/PySML_Manual_2020.pdf
				                 
where tool administration and usage are also provided)				            

1. Description
--------------

.. See the PDF package documentation for more information on the use of the tool 
   and different semantic similarity measures

PySML is a portable and expandable Python library enabling the retrieval of 
Semantic Similarity (SS) scores for any ontology to overcome issues related to 
computation, reproducibility and reusability of SS scores for any ontology in 
any application. PySML implements a large collection of SS measures consisting 
of 9 information content (IC) approaches, yielding over 624 ontology concept 
and 4430 entity pair-wise SS measures, providing a platform that eases the 
manipulation of existing SS measures for interested users and a broad 
computational audience (Refer to the PDF manual for details).


2. Environment
--------------

The PySML system is composed of one main high level folder: PySML and one main 
python module, procsemsim.py, which serves as an interface running different SS 
measures and applications implemented. The tests folder contains an illustrative 
Python modules for testing the PySML interface. The PySML folder includes three 
Python modules: 
..  informationcontent.py: implementing a class for retrieving IC scores
..  conceptsimilarity.py: inherits InformationContent class to implement a class
                          for calculating concept SS scores.
..  entitysimilarity.py: inheriting ConceptSimilarity class to build a class for
                         computing entity semantic similarity scores. 

It also contains two sub-folders, namely smlapps and imports: 
.. smpapps sub-folder: containing source codes common applications related to 
    semantic similarity measures implemented under PySML. 
.. imports sub-folder: containing imported modules for reading an ontology and 
    outputting different results.

As a library, PySML can be imported in a Python module, however, a user can also 
directly retrieve SS scores or run embedded applications in two main steps: User 
interface and input processing via a simple single command-line terminal (please
refer to the PDF manual). SS scores produced are presented in a table format, 
displayed on the screen or directed into a file.

3. Administration
-----------------

PySML v2.5.1 requires Linux operating system and Python version >= 2.7 and one 
package, python-networkx for any application implemnted. This needs to be 
installed prior to the use of PySML and for running common related applications 
to SS measures: Entity Fuzzy classification, Entity Fuzzy Identification and 
Concept Fuzzy Enrichment, additional python-scipy and python-matplotlib should be 
installed. It is worth mentioning that PySML can run in any operating system 
running Python, but it has been tested only on Linux machine.

To use PySML, the user needs to download the tar.gz file from the main URL
(see point 7 below) and extract all files as follows: 

:: tar xzf post-analysis-tools.tar.gz                                        ::

or alternatively, it can also be retrieved from the github public platform using 
git clone command line as follows.

.. git clone https://github.com/gkm-software-dev/post-analysis-tools.git     ::

After downloading and/or uncompressing, move to the appropriate working directory  
where PySML and related commands are executed using the following terminal 
command: 

:: cd post-analysis-tools/ihp-ping-dev/ ::


It is worth noting that this package is free to use under GNU General Public Li-
cense. You are free to copy, distribute and display information contained herein, 
provided that it is done with appropriate citation of the tool. Thus, by using 
the PySML package, it is assumed that you have read and accepted the agreement 
provided and that you agreed to be bound to all terms and conditions of this 
agreement. Please, use the following command line to see the package licence:

:: python setup.py --licence ::

It is worth mentioning that Python enables easy transitioning of codes between 
computers, rendering PySML portable and an effective framework for human for
retrieving and integrating human PPI networks into a unifed PPI netwok. 

4. Build status
---------------

The main website for the PySML library is: 
.. http://web.cbio.uct.ac.za/ITGOM/post-analysis-tools/pysml-dev          :: 
where users can find essential information about PySML. It is freely downloadable 
under GNU General Public License (GPL), pre-compiled for Linux version and 
protected by copyright laws. Users are free to copy, modify, merge, publish, 
distribute and  display informamation contained in the package, provided that it 
is done with appropriate citation of the package and by including the permission 
notice in all copies or substantial portions of the module contained in this 
package.

It is currently maintained by one member of the core-development team, Gaston K. 
Mazandu <gmazandu@gmail.com, gmazandu@cbio.uct.ac.za, kuzamunu@aims.ac.za, who 
regularly updates the information available in this package and makes every effort 
to ensure the quality of this information.

5. Quick start guide
--------------------
IHP-PING can be processed through one main python module, ihppinbuilder.py, which 
serves as an interface. Get help on how to run IHP-PING through this interface 
module using the following command:

::   python procsemsim.py -h ::

or

::   python procsemsim.py --help ::

From the help option above, IHP-PING is run using the following one line command:

:: python procsemsim.py -t ss-model -m models -p parameters -d dataset -a annotationfile -f ontologyfile -n namespace -o outputfile -s value ::

As illustrations:
..  python procsemsim.py -t es -m bma kstats ui -a "{'Prot1':['GO:0000022', 'GO:0051231', 'GO:1903047', 'GO:0000278', 'GO:0007052', 'GO:0000023', 'GO:0005984'], 'Prot2':['GO:0000022', 'GO:0051231', 'GO:1903047', 'GO:0000278', 'GO:0007052'], 'Prot3':['GO:1903047', 'GO:0000278', 'GO:0007052', 'GO:0000023', 'GO:0005984']}"         ::

..  python procsemsim.py -t ic -m meng universal zanchez zhang wang seco -f tests/go-basic.obo -s 0 ::

..  python procsemsim.py -t ic -m meng universal zanchez zhang wang seco -d GO:1900309 GO:1900308 GO:1900303 GO:1900302 GO:0019990                                                                          ::

.. The first command produces a table of entity pairwise 'Prot1', 'Prot2' and 'Prot3'
SS scores for models (BMA, Nunivers, Universal), non-ontology Kappa-Statistics 
(SimKPS) and Jaccard-like (UI-like) using the ontology GO biological process by 
default and displaying the result on the screen. 
.. The second command will process IC scores for Meng et al., Universal, Zanchez et al. 
Zhang et al. Wang et al. and Seco et al. models using the ontology provided under the
 **tests** folder with biological_process as a default ontology namespace, writing 
all ontology concept IC scores in a file whose the name is printed on the screen and 
located in the current working directory (by default). 
.. The last command is similar to the the second, but it uses a default ontology, which 
is provided in the **tests** folder with biological_process as ontology namespace by 
default,  displaying on the screen (-s 1) by default, IC scores only for the five 
concepts provided.

As any python library or package, PySML can be imported and used in another Python 
models. For accessing and learning about different classes of the three main classes 
under PySML, InformationContent, ConceptSimilarity and EntitySimilarity, Please 
access the python interpreter or the command shell for interactive computing (IPython) 
and run following commands:
.. >>> PySML import *                             ::
.. >>> help(InformationContent)                   ::
.. >>> help(ConceptSimilarity)                    ::
.. >>> help(EntitySimilarity)                     ::

6. Version history
------------------

- 2.4.1: Initial IHP-PING release in April 2020.
::
	 python setup --version

7. Package URL
--------------

.. http://web.cbio.uct.ac.za/ITGOM/post-analysis-tools/pysml-dev/          ::
.. https://github.com/gkm-software-dev/post-analysis-tools                 ::

::	 python setup --url           ::

8. Maintainer
-------------

Gaston K. Mazandu
Email: gaston.mazandu@uct.ac.za, gmazandu@gmail.com, 
       kuzamunu@aims.ac.za

9. Contributors
---------------

Mazandu GK, Opap K, Makinde FL, Nembaware V, Agamah F, Bope C, Chimusa ER, Wonkam A, 
Mulder NJ
Emails: gaston.mazandu@uct.ac.za, babuken@gmail.com, funmite@aims.ac.za, 
vnembaware@gmail.com, francisagamahh@gmail.com, christian.bope@gmail.com, emile.chimusa@uct.ac.za, 
nicola.mulder@uct.ac.za, ambroise.wonkam@uct.ac.za

::  python setup --classifiers        ::
Classifier: License :: GPL (>= 2)
Classifier: Operating System :: OS Independent, but tested only on Linux
Classifier: Programming Language :: Python :: >= 2.7
Classifier: Topic :: Software Development :: Libraries

Sincerely,

Gaston K. Mazandu
