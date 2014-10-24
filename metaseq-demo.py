import webbrowser
import yaml
import os
import mygene
mg = mygene.MyGeneInfo()


import metaseq
from metaseq.results_table import DESeq2Results
import numpy as np
from matplotlib import pyplot as plt
import gffutils
import pybedtools
from pybedtools import featurefuncs as ff
cache_dir = '/tmp'
dbfn = 'data/gencode.v19.annotation.chr11.db'
cached_db = os.path.join(cache_dir, os.path.basename(dbfn))
if not os.path.exists(cached_db):
    os.system('cp -v %s %s' % (dbfn, cached_db))
db = gffutils.FeatureDB(cached_db)

d = DESeq2Results('data/DESeq-results.txt')
d.attach_db(db)

gata1 = pybedtools.BedTool('data/K562_GATA1.narrowPeak')
tal1 = pybedtools.BedTool('data/K562_TAL1.narrowPeak')

with_gata = d.genes_with_peak(
    gata1,
    id_attribute='gene_id',
    transform_func=ff.TSS,
    upstream=1500,
    downstream=0,
)

with_tal = d.genes_with_peak(
    tal1,
    id_attribute='gene_id',
    transform_func=ff.TSS,
    upstream=1500,
    downstream=0,
)

with_both = d.genes_with_peak(
    gata1.intersect(tal1),
    id_attribute='gene_id',
    transform_func=ff.TSS,
    upstream=1500,
    downstream=0,
)
def go_callback(x):
    ms = mg.getgene(x.split('.')[0], fields=['go', 'name', 'symbol', 'ensembl.gene', 'genomic_pos_hg19'])
    print yaml.safe_dump(ms, default_flow_style=False)

def ucsc(x):
    g = db[x]
    url = 'http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&position={0.chrom}%3A{0.start}-{0.stop}'.format(g)
    print url
    webbrowser.open(url)


def do_both(x):
    go_callback(x)
    ucsc(x)

d.scatter(
    xfunc=np.log,
    yfunc=None,
    genes_to_highlight=[
        (
            d.upregulated(.05),
            dict(color='#a40000', alpha=0.5)
        ),
        (
            d.downregulated(0.05),
            dict(color='#204a87', alpha=0.5)
        ),
        #(
        #    with_gata,
        #    dict(s=150, color='r', facecolors='None', linewidth=1, alpha=1.0),
        #),
        #(
        #    with_tal,
        #    dict(s=150, color='b', facecolors='None', linewidth=1, alpha=1.0),
        #),
        (
            with_both,
            dict(s=150, color='#f57900', facecolors='None', linewidth=1, alpha=1.0),
        ),
        (
            d.index.isin(['ENSG00000168004.5']),
            dict(color='g', s=100, alpha=1),
        ),
    ],
    x='baseMean',
    y='log2FoldChange',
    marginal_histograms=False,
    callback=do_both,
)



plt.show()
