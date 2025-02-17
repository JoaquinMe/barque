import gzip
import sys


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
    fasta_file = sys.argv[1]  # Input fasta file
except:
    print(__doc__)
    sys.exit(0)

fasta_sequences = fasta_iterator(fasta_file)

for seq in fasta_sequences:
    print(seq.name)
