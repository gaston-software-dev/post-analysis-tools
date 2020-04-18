### The IHP-PING Package

A repository of python modules for easing integration of existing human PPI datasets from multiple sources into a unified PPI network on-the-fly (in real time), which is stored locally for further user applications. This provides a platform that enables the retrieval of human PPI datasets stored in different online resources and produced a real time up-to-date integrated human PPI network with increased accuracy, confidence and coverage. IHP-PING retrieves human PPI datasets from eight different resources (Refer to the [PDF manual](http://web.cbio.uct.ac.za/ITGOM/post-analysis-tools/ihp-ping-dev/IHP-PING_Manual_2020.pdf) for details).

## Tool administration and usage
textsf{IHP-PING} contains one main module and one main folder containing modules required for generating an integrated human PPI network and formatting results to be written into a file.

#### 1. Requirements
IHP-PING v2.4.1 requires Python &ge; 2.7), requiring the installation of the NCBI BLAST software locally when retrieving interactions predicted from sequence data. It also requires the selenium Python package, as well as chromedriver and chromium-browser, for retrieving the DIP dataset. These needs to be installed prior to the use of IHP-PING. It is worth noting that IHP-PING is an adaptable, portable and expandable Python resource that can run on any computer and any operating system provided that the computer runs Python and satisfies the requirements, but it has been tested only on Linux (ubuntu).

#### 2. Quick start guide
Please refer to the [PDF reference manual](http://web.cbio.uct.ac.za/ITGOM/post-analysis-tools/ihp-ping-dev/IHP-PING_Manual_2020.pdf) on how to download or clone the IHP-PING package. **IHP-PING** is run using the following one line command:

<span class="background-color:#58D3F7;">  python ihppinbuilder.py -r resources -o outputfolder -i outputProtID -f outputfileformat</span>


#### 3. Specific license
These tools are freely downloadable under [GNU General Public License (GPL)](https://www.gnu.org/licenses/gpl-3.0.en.html), precompiled for Linux version and protected by copyright laws, a free software and comes with ABSOLUTELY NO WARRANTY.


### Contact
Please use the [pysml-dev](http://web.cbio.uct.ac.za/ITGOM/post-analysis-tools/mysml-dev/) specific website link to contact the maintainer for related suggestions or for reporting potential bugs and errors or specific concerns related to the post-analysis-tools related codes. 

### Acknowledgements
Any work dependent on open-source software owes debt to those who developed these tools. The authors thank everyone involved with free software, from the core developers to those who contributed to the documentation. Many thanks to the authors of the freely available libraries for making this work possible. This study is supported by the National Institutes of Health (NIH), USA, under Common Fund under H3ABioNet (U24HG006941) and SADaCC (1U01HG007459-01).
