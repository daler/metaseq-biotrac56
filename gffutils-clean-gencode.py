import gffutils
gffutils.create_db(
    "data/gencode.v19.annotation.chr11.gtf",
    "data/gencode.v19.annotation.chr11.db",
    keep_order=False,
    sort_attribute_values=False,
    id_spec={'gene': 'gene_id', 'transcript': 'transcript_id'},
    verbose=True, merge_strategy='merge', infer_gene_extent=False,)
