import os
import requests
from tqdm import tqdm

# a script to download the data needed for the analysis
# cwd = os.getcwd()
cwd = "/mnt/RNASeq"
print(cwd)

# check folders exist and create
if not os.path.exists(os.path.join(cwd, "data/scRNAseq")):
    os.makedirs(os.path.join(cwd, "data/scRNAseq"))

if not os.path.exists(os.path.join(cwd, "data/ecvRNAseq")):
    os.makedirs(os.path.join(cwd, "data/ecvRNAseq"))


# helper function to download files
def download(url: str, fname: str):
    resp = requests.get(url, stream=True)
    total = int(resp.headers.get('content-length', 0))
    # Can also replace 'file' with a io.BytesIO object
    with open(fname, 'wb') as file, tqdm(
        desc=fname,
        total=total,
        unit='iB',
        unit_scale=True,
        unit_divisor=1024,
    ) as bar:
        for data in resp.iter_content(chunk_size=1024):
            size = file.write(data)
            bar.update(size)


# get scRNAseq datasets
# https://www.biorxiv.org/content/10.1101/2022.06.15.496201v1.full
# https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA858484&o=acc_s%3Aa
bamFiles = [
    "https://sra-pub-src-2.s3.amazonaws.com/SRR20172077/possorted_genome_bam_TC2_d1_1.bam.1",
    "https://sra-pub-src-2.s3.amazonaws.com/SRR20172076/possorted_genome_bam_TC2_d3_1.bam.1",
    "https://sra-pub-src-1.s3.amazonaws.com/SRR20172075/possorted_genome_bam_TC2_d5_1.bam.1",
    "https://sra-pub-src-2.s3.amazonaws.com/SRR20172074/possorted_genome_bam_TC2_d8_1.bam.1",
    "https://sra-pub-src-1.s3.amazonaws.com/SRR20172073/possorted_genome_bam_TC2_d11_1.bam.1",
    "https://sra-pub-src-1.s3.amazonaws.com/SRR20172072/possorted_genome_bam_TC2_d1_2.bam.1",
    "https://sra-pub-src-2.s3.amazonaws.com/SRR20172071/possorted_genome_bam_TC2_d3_2.bam.1",
    "https://sra-pub-src-2.s3.amazonaws.com/SRR20172070/possorted_genome_bam_TC2_d5_2.bam.1",
    "https://sra-pub-src-2.s3.amazonaws.com/SRR20172069/possorted_genome_bam_TC2_d8_2.bam.1",
    "https://sra-pub-src-1.s3.amazonaws.com/SRR20172068/possorted_genome_bam_TC2_d11_2.bam.1",
    "https://sra-pub-src-2.s3.amazonaws.com/SRR20172067/possorted_genome_bam_TC2_d15_1.bam.1",
    "https://storage.googleapis.com/worm_public/ad_worm_aging.h5ad",
]

for bf in bamFiles:
    bfn = os.path.join(
        os.path.join(cwd, "data/scRNAseq"),
        bf.split("/")[-1]
    )
    print(bfn)
    if not os.path.isfile(bfn):
        download(bf, bfn)
    if not os.path.isfile(bfn + ".bai") and "bam" in bfn:
        # routine to make index for bams
        # samtools index possorted_genome_bam_TC2_d1_2.bam.1
        os.system("samtools index {} {}".format(
            bfn,
            bfn + ".bai"
        ))

# get annotations for WS268
if not os.path.isfile(os.path.join(cwd, "data/c_elegans.PRJNA13758.WS268.annotations.gff3.gz")):
    annot = "https://downloads.wormbase.org/species/c_elegans/PRJNA13758/gff/c_elegans.PRJNA13758.WS268.annotations.gff3.gz"
    download(
        annot,
        os.path.join(cwd, "data/c_elegans.PRJNA13758.WS268.annotations.gff3.gz")
    )

# gunzip file
if not os.path.isfile(os.path.join(cwd, "data/c_elegans.PRJNA13758.WS268.annotations.gff3")):
    os.system("gunzip {}".format(
        os.path.join(cwd, "data/c_elegans.PRJNA13758.WS268.annotations.gff3.gz"),
    ))
