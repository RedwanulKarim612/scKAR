import gffutils

gencode_v45 = gffutils.create_db(
    './AmexT_v47-AmexG_v6.0-DD.gtf',
    'AmexT_v47_v6DD.db',
    merge_strategy='merge',
    keep_order=False,
    force=True,
    sort_attribute_values=True,
    disable_infer_transcripts=True,
    disable_infer_genes=True
    # id_spec={'gene': 'gene_id', 'transcript': 'transcript_id'}
)
