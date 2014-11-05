#!/bin/bash

# Prepares IPython notebooks for presentation.
#
# Cleans cached analysis data, runs each notebook from scratch, then converts
# it to PDF and .py.  Then runs the processor script to:
#
#   - clean up the comments (line wrap them)
#   - get rid of IPython "In[x}:" prompts
#   - comment out the %matplotlib equivalents (lines that start with
#   get_ipython)
#   - add a "plt.show()" to the end of the script so that plots are shown.
#
# requires runipy to be installed (pip install runipy)
set -e

clean_run_files () {

    for fn in \
        gene.npz \
        gene.features \
        tss_*.features \
        tss_*.npz \
        TSS.gtf \
        genes.gtf \
        filtered-genes.gtf \
        filelist \
        encode_specific \
        ; do
        if [ -e $fn ]; then
            rm $fn
        fi
    done
}

clean_run_files

for x in de-example heatmap-example; do
    runipy -o ${x}.ipynb

    # convert it to .pdf and .py versions
    ipython nbconvert ${x}.ipynb --to latex --post PDF
    ipython nbconvert ${x}.ipynb --to python

    # post-process the .py scripts
    python processor.py ${x}.py > tmp.py

    # clean up
    mv tmp.py ${x}.py
    rm -r ${x}_files
    rm ${x}.tex
done

clean_run_files
