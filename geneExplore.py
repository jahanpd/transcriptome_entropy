import pprint
from BCBio.GFF import GFFExaminer
from BCBio import GFF
import gffutils

# in_file = "/home/jahan/RNASeq/data/ncbi_dataset/data/GCF_000002985.6/genomic.gff"
# in_file = "/home/jahan/RNASeq/data/ncbi_dataset/data/GCA_000002985.3/genomic.gff"
# in_file = "/home/jahan/RNASeq/data/c_elegans.PRJNA13758.current.annotations.gff2"
in_file = "/home/jahan/RNASeq/data/c_elegans.PRJNA13758.WS268.annotations.gff3"
loc_db = "/home/jahan/RNASeq/data/celegans.db"

# examiner = GFFExaminer()
# in_handle = open(in_file)
# pprint.pprint(examiner.parent_child_map(in_handle))
# pprint.pprint(examiner.available_limits(in_handle))
# in_handle.close()


def id_spec(feat):
    attr = feat.attributes
    if feat.featuretype == "gene" and "Name" in attr:
        print(attr)
        return attr["Name"][0]
    elif "ID" in attr:
        return attr["ID"][0]
    else:
        return None


# db = gffutils.create_db(in_file, loc_db, force=True, id_spec=id_spec, merge_strategy="warning")
db = gffutils.FeatureDB(loc_db)
g = db['y74c9a.9']
print(g.start)
print(g.attributes)
print(g.attributes['gene'])
