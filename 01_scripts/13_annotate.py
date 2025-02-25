import gzip
import sys

import pandas as pd
from Bio import Entrez
from Bio.Blast import NCBIWWW, NCBIXML

Entrez.email = "joaquin.messano@gmail.com"


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


def calculate_percentage_identity(hsps):
    """
    Calculate percentage identity between the query and subject sequences based on HSPs.
    """
    total_identities = 0
    total_aligned_length = 0

    for hsp in hsps:
        total_identities += hsp.identities  # Number of identical residues
        total_aligned_length += hsp.align_length  # Total length of the alignment

    # Calculate percentage identity
    if total_aligned_length > 0:
        percentage_identity = (total_identities / total_aligned_length) * 100
    else:
        percentage_identity = 0

    return percentage_identity


def entrez_tax_assign(accession):
    handle = Entrez.efetch(db="nucleotide", id=accession, rettype="gb", retmode="xml")
    record = Entrez.read(handle)
    taxonomy = record[0]["GBSeq_taxonomy"]
    organism = record[0]["GBSeq_organism"]
    return taxonomy, organism


def blast_search(sequence, pid_cutoff):
    cols = ["accession", "pid", "tax", "org"]
    df = pd.DataFrame(columns=cols)
    # Perform a BLAST search using the NCBIWWW.qblast function
    result_handle = NCBIWWW.qblast("blastn", "nt", sequence)

    # Parse the BLAST result
    blast_records = NCBIXML.parse(result_handle)

    # Iterate through the BLAST records
    for blast_record in blast_records:
        for i, alignment in enumerate(blast_record.alignments):
            accession = alignment.accession
            percentage_identity = calculate_percentage_identity(alignment.hsps)
            tax, org = entrez_tax_assign(accession)
            newrow = {
                "accession": accession,
                "pid": percentage_identity,
                "tax": tax,
                "org": org,
            }
            df.loc[i] = newrow

    df = df.loc[df["pid"] >= pid_cutoff]
    taxtable = df["tax"].str.split("; ", expand=True)
    keep = True
    taxlist = []
    i = 0
    while keep and taxtable.shape[1] > i:
        col = taxtable.iloc[:, i]
        if col.unique().size == 1:
            keep = True
            taxlist.append(col.unique()[0])
            i += 1
        else:
            keep = False
    return taxlist


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
    fasta_file = sys.argv[1]  # Input fasta file
    perc_id = float(sys.argv[2]) * 100  # porcentaje de similaridad para cortar
except:
    print(__doc__)
    sys.exit(0)

fasta_sequences = fasta_iterator(fasta_file)

count = 0
result = []
# TODO:
# clusterizar antes de mandar a blast
# agarrar los centroides y tirar esos a blast
# anotar esos en la bd
for seq in fasta_sequences:
    print(f"Sequence {count}")
    # Devuelve una lista de taxones(["Eukaryota","Metazoa",...,"Homo"])
    data = blast_search(seq.sequence, perc_id)
    tax_joined = "_".join(data)
    result.append([seq.name, seq.sequence, tax_joined])
    count += 1
df = pd.DataFrame(result, columns=["name", "sequence", "tax"])
df.to_csv("result_097.csv", index=False)
