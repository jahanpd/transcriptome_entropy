#!/usr/bin/env python3
import gffutils
import anndata as ad
import pandas as pd
import os

base = "/mnt/RNASeq"

loc_db = os.path.join(base, "data/celegans.db")
data = ad.read_h5ad(os.path.join(base, "data/scRNAseq/ad_worm_aging.h5ad"))

db = gffutils.FeatureDB(loc_db)

data = ad.read_h5ad(os.path.join(base, "data/scRNAseq/ad_worm_aging.h5ad"))


genes = data.var.gene_ids

data_dict = dict(
    gene_id=[],
    contig=[],
    start=[],
    end=[]
)

for gene in genes:
    g = db[gene]

    data_dict["gene_id"].append(gene)
    data_dict["contig"].append(g.seqid)
    data_dict["start"].append(g.start)
    data_dict["end"].append(g.end)

pd.DataFrame(data_dict).to_csv(
    os.path.join(base, "data/genes.csv"),
    index=False
)
