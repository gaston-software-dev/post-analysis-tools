A repository of python modules for retrieving semantic similarity (SS) scores for any ontology in any format (OWL, OBO, RDF), resolving issues related to computation, reproducibility and reusability of SS scores for any ontology in any application. This provides a platform that eases the manipulation of existing SS measures for interested users and a broad computational audience. **PySML** PySML implements a large collection of SS measures consisting of 9 information content (IC) approaches, yielding over 624 ontology concept and 4430 entity pair-wise SS measures. (Please refer to the [PDF manual](http://web.cbio.uct.ac.za/ITGOM/post-analysis-tools/pysml-dev/PySML_Manual_2020.pdf) for details).

### Tool administration and usage
**PySML** system is composed of one main high level folder: _PySML_ and one main python module, _procsemsim.py_, which serves as an interface running different SS measures and applications implemented. The **tests** folder contains an illustrative Python modules for testing the _PySML_ interface. 

#### 1. Requirements
**PySML** v2.5.1 requires Python &ge; 2.7) and one package, _python-networkx_ is required for any application implemented. This needs to be installed prior to the use of PySML and for running common related applications to SS measures: Entity Fuzzy classification, Entity Fuzzy Identification and Concept Fuzzy Enrichment, additional _python-scipy_ and _python-matplotlib_ should be installed. It is worth noting that this library is an adaptable, portable and expandable Python resource that can run on any computer and any operating system provided that the computer runs Python and satisfies the requirements, but it has been tested only on Linux (ubuntu).

#### 2. Quick start guide
Please refer to the [PDF reference manual](http://web.cbio.uct.ac.za/ITGOM/post-analysis-tools/pysml-dev/PySML_Manual_2020.pdf) on how to download or clone the **PySML** package. **PySML** is run using the following one line command:
<pre>
     python procsemsim.py -t ss-model -m models -p parameters -d dataset -a annotationfile -f ontologyfile -n namespace -o outputfile -s value
</pre>
As illustrations:
<pre>
    python procsemsim.py -t es -m bma kstats ui -a "{'Prot1':['GO:0000022', 'GO:0051231', 'GO:1903047', 'GO:0000278', 'GO:0007052', 'GO:0000023', 'GO:0005984'], 'Prot2':['GO:0000022', 'GO:0051231', 'GO:1903047', 'GO:0000278', 'GO:0007052'], 'Prot3':['GO:1903047', 'GO:0000278', 'GO:0007052', 'GO:0000023', 'GO:0005984']}"
    python procsemsim.py -t ic -m meng universal zanchez zhang wang seco -f tests/go-basic.obo -s 0
    python procsemsim.py -t ic -m meng universal zanchez zhang wang seco -d GO:1900309 GO:1900308 GO:1900303 GO:1900302 GO:0019990
</pre>
The first command produces a table of entity pairwise 'Prot1', 'Prot2' and 'Prot3' SS scores for models (BMA, Nunivers, Universal), non-ontology Kappa-Statistics (SimKPS) and Jaccard-like (UI-like) using the ontology GO biological process by default and displaying the result on the screen. The second command will process IC scores for Meng et al., Universal, Zanchez et al. Zhang et al. Wang et al. and Seco et al. models using the ontology provided under the **tests** folder with biological\_process as a default ontology namespace, writing all ontology concept IC scores in a file whose the name is printed on the screen and located in the current working directory (by default). The last command is similar to the the second, but it uses a default ontology, which is provided in the **tests** folder with biological\_process as ontology namespace by default,  displaying on the screen (-s 1) by default, IC scores only for the five concepts provided.

### Contact
Please use the [pysml-dev](http://web.cbio.uct.ac.za/ITGOM/post-analysis-tools/pysml-dev/) specific website link to contact the maintainer for related suggestions or for reporting potential bugs and errors or specific concerns related to the post-analysis-tools related codes. 

### Acknowledgements
Any work dependent on open-source software owes debt to those who developed these tools. The authors thank everyone involved with free software, from the core developers to those who contributed to the documentation. Many thanks to the authors of the freely available libraries for making this work possible. This study is supported by the National Institutes of Health (NIH), USA, under Common Fund under H3ABioNet (U24HG006941) and SADaCC (1U01HG007459-01).
