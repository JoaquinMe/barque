import gzip
import sys

import pandas as pd


# Defining classes
class Fasta(object):
    """Fasta object with name and sequence"""

    def __init__(self, name, sequence):
        self.name = name
        self.sequence = sequence

    def write_to_file(self, handle):
        handle.write(">" + self.name + "\n")
        handle.write(self.sequence + "\n")


# Defining functions
def myopen(infile, mode="rt"):
    if infile.endswith(".gz"):
        return gzip.open(infile, mode=mode)
    else:
        return open(infile, mode=mode)


def fasta_iterator(input_file):
    """Takes a fasta file input_file and returns a fasta iterator"""
    with myopen(input_file) as f:
        sequence = ""
        name = ""
        begun = False
        for line in f:
            line = line.strip()
            line = str(line)
            if line.startswith(">"):
                if begun:
                    yield Fasta(name, sequence)
                name = line.replace(">", "")
                sequence = ""
                begun = True
            else:
                sequence += line

        if name != "":
            yield Fasta(name, sequence)


# Parsing user input
try:
    centroid_file = sys.argv[1]  # Input fasta file
    clusters_file = sys.argv[2]  # clusters tsv
    result_file = sys.argv[3]  # result_097.csv
except:
    print(__doc__)
    sys.exit(0)
colnames = [
    "record_type",
    "cluster_number",
    "sequence_length",
    "pid",
    "strand",
    "legacy",
    "legacy2",
    "alignment",
    "qlabel",
    "tlabel",
]
centroids_iterator = fasta_iterator(centroid_file)
clusters_table = pd.read_csv(clusters_file, sep="\t", names=colnames)
result_table = pd.read_csv(result_file, sep=",")
# armé centroids
# grep '^C' test/clusters.uc para ver como están los clusters
# grep '^S' test/clusters.uc para ver los centroides
# grep '^H' test/clusters.uc para ver los hits
centroid_dict = {}
for seq in centroids_iterator:
    centroid_dict[seq.name] = {}

for k in centroid_dict.keys():
    lista_matches = list(
        clusters_table.loc[
            (clusters_table["record_type"] == "H") & (clusters_table["tlabel"] == k)
        ]["qlabel"]
    )
    centroid_dict[k] = {ele: "" for ele in lista_matches}


# centroid_dict[centroid][child]="asignación taxonomica"
# asigno taxonomia al diccionario
for centroid in centroid_dict.keys():
    for child in centroid_dict[centroid].keys():
        # encontrar en la tabla la asignacion taxonómica del child
        tax = result_table.loc[result_table["name"] == child]["tax"].iloc[0]
        centroid_dict[centroid][child] = tax


# chequeo como da la asignación
tax_check = {k: set() for k in centroid_dict.keys()}
for centroid in centroid_dict.keys():
    for child in centroid_dict[centroid].keys():
        tax_check[centroid].add(centroid_dict[centroid][child])

for k in tax_check:
    print(tax_check[k])
    print(len(tax_check[k]))
    print()
