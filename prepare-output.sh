#!/bin/bash

# Prepares IPython notebooks for presentation.

# pip install runipy
#
# this runs the notebook from start to finish to ensure everything works
runipy -o de-example.ipynb

# convert it to .pdf and .py versions
ipython nbconvert de-example.ipynb --to latex --post PDF
ipython nbconvert de-example --to python
