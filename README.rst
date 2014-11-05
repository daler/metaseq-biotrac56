Overview
--------
This repository is used as part of the Bio-Trac 56 RNA-seq course.

Quick view
----------
Don't want to download the data or install anything?  You can view the rendered
notebooks here:

* http://nbviewer.ipython.org/github/daler/metaseq-biotrac56/blob/master/de-example.ipynb
* http://nbviewer.ipython.org/github/daler/metaseq-biotrac56/blob/master/heatmap-example.ipynb

Alternatively, you can download the PDF versions in this repository for offline
viewing.

Files in the repository
-----------------------

- `de-example.ipynb`: an IPython Notebook to demo some differntial expression analysis
- `heatmap-example.ipynb`: an IPython notebook to demo integration with
   ChIP-seq data
- `processor.py`: for post-processing .py files created from .ipynb files
- `prepare-output.sh`: script to prepare files
- `download_data.bash` to download all data used in the demo
- `RNA-seq.R` R script to perform differential expression on the downloaded
   RNA-seq data


Setup
-----

In order to run the IPython notebooks on your own computer, you'll need to:

* install metaseq and its prerequisites
* download the data used in this analysis

The following instructions assume you do not have the scientific Python stack
installed and do not have genomics tools (BEDTools, samtools, tabix, UCSC
utilities) already installed.  If you do, then see the end of this README for
more advanced instructions.

While this installation process takes a few minutes, the nice thing is that
when it's complete you will have a fully-functioning scientific Python
installation as well as some common genomics tools.

Install metaseq and prerequisites
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Instructions are slightly different for Mac and Linux; follow the directions on
https://pythonhosted.org/metaseq/install.html. This will help you run the
`metaseq` installation script, which will:

- download and install genomics tools (BEDTools, samtools, tabix, UCSC
  utilites)
- download and install an isolated Python environment
- download and install prerequisites for metaseq
- download and install metaseq itself

The isolated Python environment will not affect any other Python versions you
have on your computer.

The page at https://pythonhosted.org/metaseq/install.html has more details and
the link to the installation script.


Activate the environment
~~~~~~~~~~~~~~~~~~~~~~~~
Once the installation completes:

1. Open a new terminal
2. Assuming you've accepted the defauls from the installation script, run::

    source activate metaseq-test

3. Whenever you're done using the test environment and want to go back to
   normal, simply close the terminal.


Download the materials and data
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

1. Download the code to run the analysis from
   https://github.com/daler/metaseq-biotrac56/archive/master.zip
2. Unzip the code somewhere convenient on your machine.  Let's say you unzipped
   it to `~/metaseq-biotrac56`.
3. Download the data (755 MB) from
   http://helix.nih.gov/~dalerr/metaseq-biotrac56-data.zip.
4. Extract the `data` folder and place it in the same directory as where you
   unzipped the code.  So if you had unzipped the code to
   `~/metaseq-biotrac56`, you now have a directory called
   `~/metaseq-biotrac56/data` and files like
   `~/metaseq-biotrac56/data/H1-hESC_1.chr11.bam`.


Install demo prerequisites
~~~~~~~~~~~~~~~~~~~~~~~~~~
There are some additional requirements that we use in the demo that need to be
installed.

1. Make sure the `metaseq-test` environment is activated (see above)

2. Go to the directory where you've unzipped the materials

3. Run::

    pip install -r extra-requirements.txt


Run the IPython Notebook
~~~~~~~~~~~~~~~~~~~~~~~~

1. Make sure the `metaseq-test` environment is activated (see above)

2. Go to the directory where you've unzipped the materials

3. Run::

    ipython notebook

4. Your web browser should open showing a list of `.ipynb` files.  Click on one
   to begin.



Advanced
--------


Already have things installed?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If you have a scientific Python installation along with commonly used genomics
tools (BEDTools, samtools, tabix, UCSC utilities) run::

    pip install metaseq
    pip install -r extra-requirements.txt

Alternatively, if you only want to install a subset of these tools, you can run
the `metaseq` installation script with the `-h` option to see available
options.  See https://pythonhosted.org/metaseq/install.html#customizing.


Downloading and processing data from scratch
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
NOTE: this is only needed if you don't download the already-prepared data.
Furthermore, you'll need to have R and DESeq2 installed to perform the
differential expression.

1. Make sure the `metaseq-test` environment is activated (see above)

2. Go to the directory where you've unzipped the materials, and run::

    bash download_data.bash

(this will take a while, something like 20 minutes depending on your
connection)

3. Assuming you have R and DESeq2 installed, run::

    Rscript RNA-seq.R

