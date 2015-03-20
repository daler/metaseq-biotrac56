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

- ``de-example.ipynb``: an IPython Notebook to demo some differential expression analysis
- ``heatmap-example.ipynb``: an IPython notebook to demo integration with
  ChIP-seq data
- ``processor.py``: for post-processing .py files created from .ipynb files
- ``prepare-output.sh``: script to prepare files
- ``download_data.bash`` to download all data used in the demo
- ``RNA-seq.R``:  R script to perform differential expression on the downloaded
  RNA-seq data


Setup
-----

In order to run the IPython notebooks on your own computer, you will need to
install the software and download the data.

The following instructions assume you do not have the scientific Python stack
installed and do not have genomics tools (BEDTools, samtools, tabix, UCSC
utilities) already installed.  If you do, then see the end of this README for
more advanced instructions.

While this installation process takes a few minutes, the nice thing is that
when it's complete you will have a fully-functioning scientific Python
installation as well as some common genomics tools.


Quick version if using Linux
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
If you're using Linux, the following commands will download and install metaseq
and dependencies, download and unpack the example files and example data, and
start the IPython Notebook server::

    # Install metaseq and dependencies
    wget --no-check-certificate https://raw.githubusercontent.com/daler/metaseq/master/create-metaseq-test-environment.sh
    bash create-metaseq-test-environment.sh -v

    # Download and unpack code and data for this example
    wget https://github.com/daler/metaseq-biotrac56/archive/master.zip
    unzip master.zip
    cd metaseq-biotrac56-master
    wget http://helix.nih.gov/~dalerr/metaseq-biotrac56-data.zip
    unzip metaseq-biotrac56-data.zip

    # Activate metaseq environment
    source activate metaseq-test

    # Install additional prerequisites
    pip install -r extra-requirements.txt

    # Start the notebook
    ipython notebook


A web browser should open, and you should see the IPython notebook starting
page.


Quick version if using Mac
~~~~~~~~~~~~~~~~~~~~~~~~~~
::

    # Install metaseq and dependencies
    curl -O https://raw.githubusercontent.com/daler/metaseq/master/create-metaseq-test-environment.sh
    bash create-metaseq-test-environment.sh -v

    # Download and unpack code and data for this example
    curl -O -L https://github.com/daler/metaseq-biotrac56/archive/master.zip
    unzip master.zip
    cd metaseq-biotrac56-master
    curl -O -L http://helix.nih.gov/~dalerr/metaseq-biotrac56-data.zip
    unzip metaseq-biotrac56-data.zip

    # Activate metaseq environment
    source activate metaseq-test

    # Install additional prerequisites
    pip install -r extra-requirements.txt

    # Start the notebook
    ipython notebook

A web browser should open, and you should see the IPython notebook starting
page.

Usage
~~~~~
Any time you want to use the new environment, open a terminal and type::

    source activate metaseq-test

This lets you use any of the installed Python packages and lets you install any
other packages without needing admin rights.  It is completely separate from
any other Python you might have installed on your machine. This command adjusts
your `$PATH` variable so that the first place it looks for programs is the
`~/miniconda/bin` directory.

Note that this environment can be used for *any* scientific Python use, not
just for `metaseq`.

Deactivate the environment by closing the terminal, or::

    source deactivate



Detailed installation instructions
----------------------------------

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

2. Assuming you've accepted the defaults from the installation script, run::

    source activate metaseq-test

3. Whenever you're done using the test environment and want to go back to
   normal, simply close the terminal.

4. When you want to use the test environment again, you need to run::

    source activate metaseq-test



Download the materials and data
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

1. Download the code to run the analysis from
   https://github.com/daler/metaseq-biotrac56/archive/master.zip
2. Unzip the code somewhere convenient on your machine.  Let's say you unzipped
   it to ``~/metaseq-biotrac56``.
3. Download the data (755 MB) from
   http://helix.nih.gov/~dalerr/metaseq-biotrac56-data.zip.
4. Extract the ``data`` folder and place it in the same directory as where you
   unzipped the code.  So if you had unzipped the code to
   ``~/metaseq-biotrac56``, you now have a directory called
   ``~/metaseq-biotrac56/data`` and lots of data files; for example one of them
   should be ``~/metaseq-biotrac56/data/H1-hESC_1.chr11.bam``.


Install demo prerequisites
~~~~~~~~~~~~~~~~~~~~~~~~~~
There are some additional requirements that we use in the demo that need to be
installed (``mygene`` and ``fisher``, for example).

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


Updates
-------
Learning Python
~~~~~~~~~~~~~~~
After the talk, some people asked about learning Python in general. Here are
some links to get you started:

- Think Like a Computer Scientist: http://interactivepython.org/courselib/static/thinkcspy/toc.html
- Codecademy: http://www.codecademy.com/learn
- Software Carpentry:  http://software-carpentry.org/index.html
- Python for Biologists: http://pythonforbiologists.com/index.php/introduction-to-python-for-biologists/

Virtual Machine
~~~~~~~~~~~~~~~
Some people expressed interest in a virtual machine that they could use to try
the demos.

This virtual machine contains Ubuntu 14.04.1 with everything you need to run
the demos.  You can use it on Mac or Windows by first installing the
`VirtualBox <https://www.virtualbox.org/wiki/Downloads>`_ application and the
VirtualBox Extension Pack for your platform.  This will let you run the image
as a guest operating system on your machine.

Then download the following virtual machine image (4.7 GB):
http://helix.nih.gov/~dalerr/metaseq-vm.ova

When it is done downloading, double-click its icon to import.  Or, if
VirtualBox is already open, choose File -> Import Appliance. **NOTE: you only
have to do this once.**

You have the opportunity to make some tweaks, like how much RAM you'd like the
VM to have.  You can always change this later.

Now, any time you want to start the VM:

1. Open the VirtualBox program
2. Select the metaseq VM from the list (it's probably the only one on the list)
3. Press the "Start" button.

Ubuntu will now start up in a separate window.  See the `VirtualBox manual
<http://www.virtualbox.org/manual/>`_ for more details.

username: ``ubuntu``
password: ``ubuntu``

The "installation-details.txt" file on the desktop shows what was installed.

To run the demos, open a terminal, and use the commands::

    source activate metaseq-test

    cd metaseq-biotrac56-master

    ipython notebook

A Firefox window will pop up, running the demo.

When you're done, either go to the gear icon in the upper right and choose
"Shut Down", or simply close the window to shut down the VM. The current state
of the VM will be saved for the next time you start it up, so any changes you
make will be saved on your computer.
