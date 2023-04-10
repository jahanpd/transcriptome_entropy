import os
import anndata as ad
import gffutils
import pysam
import numpy as np
from database import api
from tqdm import tqdm
from tqdm.contrib.concurrent import process_map

base = "/mnt/RNASeq"

loc_db = os.path.join(base, "data/celegans.db")
data = ad.read_h5ad(os.path.join(base, "data/scRNAseq/ad_worm_aging.h5ad"))
bamLoc = os.path.join(base, "data/scRNAseq/")
# create sqlite database to store data
dbepath = os.path.join(base, "data/entropy.db")
_ = api(dbepath)

FILENAMES = []
for path in os.listdir(bamLoc):
    if os.path.isfile(os.path.join(bamLoc, path)):
        if ".bam.1" in path and "bai" not in path:
            FILENAMES.append(os.path.join(bamLoc, path))
print(data)
print(data.obs)
print(data.var)

# store calculated entropy for each gene and celltype per timepoint
entropy = np.zeros((len(data.obs), len(data.var)), dtype=np.float32)
# store the raw sequence counts for each base pair which is a flatted array of [[g,c,a,t,del,refskip], ...]
distData = np.zeros((len(data.obs), len(data.var)), dtype="object")

db = gffutils.FeatureDB(loc_db)

# ratio = data.var.n_counts / data.var.n_cells
# genes = [x for x, y in zip(data.var.gene_ids, ratio) if y > 1.]

sub = data[data.obs.annotate_name == "41_1:NA pharyngeal epithelium"]
print(sub.obs)


asd
# for each bam file
# get the time and series info
# subset the barcode markers, and indexes in the h5ad that match that time and series info
# then for each gene in the h5ad
# then for each read in gene space check that read barcode is in subset and create entry in hashmap/dict
#
#
print(FILENAMES)


def process_gene(inputs):
    gene, g, bamfile = inputs
    dbe = api(dbepath)
    samfile = pysam.AlignmentFile(filename, "rb")
    # TODO a check to see if already a an entry and skip

    # store keeps track of all the unique barcodes (aka cell types)
    # for each cell the counts for each g c a t del are recorded for each base
    # entropy is calculated from these counts
    store = dict()

    for pileupcolumn in samfile.pileup(g.seqid, g.start, g.end, truncate=True):
        for pileupread in pileupcolumn.pileups:
            tags = pileupread.alignment.get_tags()
            CR = [x for x in tags if x[0] == "CR"][0][1]
            if CR not in store:
                store[CR] = [np.zeros(6, dtype=np.float32)]
            else:
                store[CR].append(np.zeros(6, dtype=np.float32))
            store[CR][-1][0] += pileupcolumn.pos
            if not pileupread.is_del and not pileupread.is_refskip:
                # query position is None if is_del or is_refskip is set.
                b = pileupread.alignment.query_sequence[
                    pileupread.query_position]
                if b == "G":
                    store[CR][-1][1] += 1.
                elif b == "C":
                    store[CR][-1][2] += 1.
                elif b == "A":
                    store[CR][-1][3] += 1.
                elif b == "T":
                    store[CR][-1][4] += 1.
                else:
                    print(b)
                    raise Exception("sequence token no handled")
            else:               #
                if pileupread.is_del:
                    store[CR][-1][5] += 1.

    # iterate over store and calculate entropy
    # save store and entropy to disk

    inserts = []
    for key in store.keys():
        entropy = 0.
        for arr in store[key]:
            total = arr[1:].sum()
            if total == 0:
                # this condition is met when no reads encompass base
                break
            else:
                part = lambda x, t: 0. if x == 0 else -np.log(x/t)*(x/t)
                entropy += sum([part(x, total) for x in arr[1:]])

        # filt = [x for x in store[key] if x[1:].sum() > 0]
        inserts.append(
            (
                ":".join([key, gene]),
                entropy,
                np.stack(store[key])
            )
        )
    dbe.bulk_insert(
        inserts
    )
    samfile.close()


print("STARTING")
# fpbar = tqdm(FILENAMES, desc="BAM files", leave=True, position=0)
for filename in FILENAMES:
    print(filename)
    time = "_".join(filename[:-5].split("_")[3:6]).replace(".", "")

    gs = [db[gene] for gene in genes]
    inputs = list(zip(genes, gs, [filename]*len(genes)))

    process_map(process_gene, inputs, max_workers=13)
