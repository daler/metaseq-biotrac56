Overview
--------
This repository is used as part of the Bio-Trac 56 RNA-seq course.

It contains:

    - IPython Notebooks which are used for the demo of metaseq
    - `download_data.bash` to download all data used in the demo
    - `RNA-seq.R` R script to perform differential expression on the downloaded
      RNA-seq data


If you don't care about running the demo yourself, you can simply look at the
rendered notebook at:

    http://nbviewer.ipython.org/github/daler/metaseq-biotrac56/blob/master/de-example.ipynb

Otherwise, follow the instructions below to get set up.


Setup
-----

Note: see the end of this README if you have the scientific Python stack and
genomics tools (BEDTools, samtools, tabix, UCSC utilities) already installed.


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


Activate the environment
~~~~~~~~~~~~~~~~~~~~~~~~
Once the installation completes:

1. Open a new terminal
2. Assuming you've accepted the defauls from the installation script, run::

    source activate metaseq-test

3. When you're done, simply close the terminal.



Download the materials
~~~~~~~~~~~~~~~~~~~~~~

1. Go to https://github.com/daler/metaseq-biotrac56 and click the *Download
   ZIP* button on the lower right-hand side.

2. Unzip this file somewhere convenient on your machine



Downloading the prepared data
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Data will be available at a to-be-determined location (approx 530GB).


Installing demo prerequisites
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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

