import pysam
import numpy as np
import pandas as pd

FILENAME = "./data/possorted_genome_bam_TC2_d1_1.bam.1"

samfile = pysam.AlignmentFile(FILENAME, "rb")

# for read in samfile.head(10):
#     print(read)
#     print(read.reference_id)

data_dict = dict(
    chromosome=[],
    position=[],
    coverage=[],
    G=[],
    C=[],
    A=[],
    T=[],
    entropy=[],
    isDel=[],
    isRefSkip=[]
)


def entCalc(arr):
    total = sum(arr) + 1e-8
    ent = [-np.log(x/total + 1e-8)*(x/total) for x in arr]
    return sum(ent)


for pileupcolumn in samfile.pileup("X"):
    print("\ncoverage at base %s = %s" % (
        pileupcolumn.pos, pileupcolumn.n,
    ))
    g = 0.
    c = 0.
    a = 0.
    t = 0.
    nDel = 0.
    nRefSkip = 0.

    for pileupread in pileupcolumn.pileups:
        if not pileupread.is_del and not pileupread.is_refskip:
            # query position is None if is_del or is_refskip is set.
            # tags = pileupread.alignment.get_tags()
            # cr = [x for x in tags if x[0] == "CR"][0][1]
            # print(cr)
            base = pileupread.alignment.query_sequence[
                pileupread.query_position]
            if base == "G":
                g += 1
            elif base == "C":
                c += 1
            elif base == "A":
                a += 1
            elif base == "T":
                t += 1
            else:
                print(base)
                raise Exception("sequence token no handled")

        else:
            if pileupread.is_del:
                nDel += 1
            if pileupread.is_refskip:
                nRefSkip += 1

    data_dict["chromosome"].append("X")
    data_dict["position"].append(pileupcolumn.pos)
    data_dict["coverage"].append(pileupcolumn.n)
    data_dict["G"].append(g)
    data_dict["C"].append(c)
    data_dict["A"].append(a)
    data_dict["T"].append(t)
    data_dict["entropy"].append(entCalc([g, c, a, t]))
    data_dict["isDel"].append(nDel)
    data_dict["isRefSkip"].append(nRefSkip)

samfile.close()

df = pd.DataFrame(data_dict)
print(df)

df.to_parquet("ChrIWholeOrganismDay1.parquet")
