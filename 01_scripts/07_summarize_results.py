#!/usr/bin/env python3
"""Summarize vsearch results

Usage:
    <program> input_folder output_folder primer_file min_length min_coverage min_coverage_experiment

Where:
    input_folder is '09_vsearch'
    output_folder is '12_results'
    primer_file is '02_info/primers.csv'
    min_length is the minimum length of the hits to keep (typically >= 100)
    min_coverage is the minimun number of hits a taxon (species, genus or phylum)
        must have in *at least* one sample in order for the taxon to be kept
    min_coverage_experiment is the minimun number of hits a taxon (species, genus or phylum)
        must have in *all the combined samples* in order for the taxon to be kept
"""

# Modules

# Silencing distutils deprecation warning
import warnings

warnings.filterwarnings("ignore", category=DeprecationWarning)

import gzip
import os
import sys
from collections import defaultdict
from distutils.version import LooseVersion


# Defining functions
def myopen(infile, mode="rt"):
    if infile.endswith(".gz"):
        return gzip.open(infile, mode=mode)
    else:
        return open(infile, mode=mode)


# Parsing user input
try:
    input_folder = sys.argv[1]
    output_folder = sys.argv[2]
    primer_file = sys.argv[3]
    min_length = int(sys.argv[4])
    min_coverage = int(sys.argv[5])
    min_experiment_coverage = int(sys.argv[6])
except:
    print(__doc__)
    sys.exit(1)

# Read primer_file
primers = dict()
with open(primer_file) as pfile:
    for line in pfile:
        if line.startswith("#"):
            continue

        l = line.strip().split(",")
        primers[l[0]] = l[5]

# Read vsearch results from input_folder
try:
    result_files = sorted(os.listdir(input_folder), key=LooseVersion)
except TypeError:
    result_files = sorted(os.listdir(input_folder))

species_dictionary = {}
genus_dictionary = {}
phylum_dictionary = {}
similarity_dict = defaultdict(int)

# Iterate through primers, gather taxon counts
for primer in primers:
    multiple_hits_species = defaultdict(int)
    multiple_hits_genus = defaultdict(int)
    multiple_hits_global_infos = defaultdict(list)

    # Get minimum similarity for the primer
    primer_info = [
        x.strip().split(",")
        for x in open(primer_file).readlines()
        if x.startswith(primer + ",")
    ][0]
    min_species_similarity = 100 * float(primer_info[6])
    min_genus_similarity = 100 * float(primer_info[7])
    min_phylum_similarity = 100 * float(primer_info[8])

    # Create summary dictionary
    species_dictionary[primer] = {}
    genus_dictionary[primer] = {}
    phylum_dictionary[primer] = {}

    # List result files for primer
    primer_results = [
        filename for filename in result_files if "_" + primer + "_" in filename
    ]
    primer_results = [
        filename
        for filename in primer_results
        if not (
            filename.endswith("_matched.fasta.gz")
            or filename.endswith("_matched.fasta")
        )
    ]

    # Iterate through result files for primer
    for result_file in primer_results:
        sample = result_file.split("_")[0]
        sys.stdout.write("  " + sample)
        sys.stdout.flush()
        species_dictionary[primer][sample] = defaultdict(int)
        genus_dictionary[primer][sample] = defaultdict(int)
        phylum_dictionary[primer][sample] = defaultdict(int)

        # Get infos from result file
        seen = set()
        sequence_dict = defaultdict(list)
        sequence_list = []

        # Read all hits for each sequence into a dictionary
        with myopen(os.path.join(input_folder, result_file)) as rfile:
            for line in rfile:
                sequence_name = line.split()[0]
                sequence_dict[sequence_name].append(line.strip().split()[1:])
                if not sequence_name in seen:
                    seen.add(sequence_name)
                    sequence_list.append(sequence_name)

        # Treat each sequence
        for seq in sequence_list:

            count = int(seq.split("_")[3])
            best_score = max([float(x[1]) for x in sequence_dict[seq]])

            # Find best species, genus or phylum
            best_species = set(
                [
                    "_".join(x[0].split("_")[:3])
                    for x in sequence_dict[seq]
                    if (
                        float(x[1]) == best_score
                        and float(x[1]) >= min_species_similarity
                        and int(x[2]) >= min_length
                    )
                ]
            )

            best_genus = set(
                [
                    "_".join(x[0].split("_")[:2])
                    for x in sequence_dict[seq]
                    if (
                        float(x[1]) == best_score
                        and float(x[1]) >= min_genus_similarity
                        and int(x[2]) >= min_length
                    )
                ]
            )

            best_phylum = set(
                [
                    x[0].split("_")[0]
                    for x in sequence_dict[seq]
                    if (
                        float(x[1]) == best_score
                        and float(x[1]) >= min_phylum_similarity
                        and int(x[2]) >= min_length
                    )
                ]
            )

            # Species level identification
            if len(best_species) == 1:
                species = list(best_species)[0]
                species_dictionary[primer][sample][species] += count

                genus = list(best_genus)[0]
                genus_dictionary[primer][sample][genus] += count

                phylum = list(best_phylum)[0]
                phylum_dictionary[primer][sample][phylum] += count

                # Collect similarity per species and write to file at the end
                similarity_dict[
                    (result_file.split("_")[0], primer, species, str(best_score))
                ] += count

            # Summaryze multiple hits
            elif len(best_species) > 1:

                species = "zMultiple_Hits_" + " : ".join(
                    sorted(list(best_species))
                ).replace("_", "^")
                species_dictionary[primer][sample][species] += count
                multiple_hits_species[";".join(sorted(list(best_species)))] += count
                # similarity_dict[(result_file.split("_")[0], primer, species.replace(" ", "").replace("^", "_"), str(best_score))] += count

                # Gather multiple hit infos
                multiple_hits_global_infos[";".join(list(best_species))].append(
                    (sample, seq)
                )

                # Genus level identification
                if len(best_genus) == 1:
                    genus = list(best_genus)[0]
                    genus_dictionary[primer][sample][genus] += count

                    phylum = list(best_phylum)[0]
                    phylum_dictionary[primer][sample][phylum] += count

                if len(best_genus) > 1:
                    multiple_hits_genus[":".join(sorted(list(best_genus)))] += count

                # Phylum level identification
                else:
                    phylum = list(best_phylum)[0]
                    phylum_dictionary[primer][sample][phylum] += count

    # Write new line after all the samples of a primer pair are finished
    print("")

    # Write multiple hit summary for species
    lines_species = []
    for group in sorted(multiple_hits_species.keys()):
        lines_species.append((multiple_hits_species[group], group))

    lines_species = sorted(lines_species, reverse=True)
    lines_species = [str(x[0]) + "," + x[1] for x in lines_species]

    # Write multiple hit summary for genus
    lines_genus = []
    for group in sorted(multiple_hits_genus.keys()):
        lines_genus.append((multiple_hits_genus[group], group))

    lines_genus = sorted(lines_genus, reverse=True)
    lines_genus = [str(x[0]) + "," + x[1] for x in lines_genus]

# Get represented taxons
species_found = {}
for primer in species_dictionary:
    species_found[primer] = set()

    for sample in species_dictionary[primer]:
        for species in species_dictionary[primer][sample]:
            count = species_dictionary[primer][sample][species]
            if count > 0:
                # Species
                species_found[primer].add(species)

    species_found[primer] = sorted(list(species_found[primer]))

genus_found = {}
for primer in genus_dictionary:
    genus_found[primer] = set()

    for sample in genus_dictionary[primer]:
        for genus in genus_dictionary[primer][sample]:
            count = genus_dictionary[primer][sample][genus]
            if count > 0:

                # Genus
                genus_found[primer].add(genus)

    genus_found[primer] = sorted(list(genus_found[primer]))

phylum_found = {}
for primer in phylum_dictionary:
    phylum_found[primer] = set()

    for sample in phylum_dictionary[primer]:
        for phylum in phylum_dictionary[primer][sample]:
            count = phylum_dictionary[primer][sample][phylum]
            if count > 0:

                # Phylum
                phylum_found[primer].add(phylum)

    phylum_found[primer] = sorted(list(phylum_found[primer]))

# Summarize results
for primer in sorted(species_dictionary):

    # Create header line
    species_table = [["Group,Genus,Species,TaxonName,Total"]]
    genus_table = [["Group,Genus,Total"]]
    phylum_table = [["Group,Total"]]

    for sample in sorted(species_dictionary[primer]):
        species_table[0].append(sample)
        genus_table[0].append(sample)
        phylum_table[0].append(sample)

    # Add rows to table
    # Species
    for species in species_found[primer]:
        if not "Multiple_Hits" in species:
            if species.endswith("_"):
                print(species)
                print(",".join(species.split("_")).replace("^", "_").split(","))
                print()
        species_table.append(",".join(species.split("_")).replace("^", "_").split(","))
        count_sum = 0
        new_line = []
        for sample in sorted(species_dictionary[primer]):
            count = species_dictionary[primer][sample][species]
            count_sum += count
            new_line.append(str(count))
        species_table[-1].append(str(count_sum))
        species_table[-1] += new_line

    # Genus
    for genus in genus_found[primer]:
        genus_table.append([",".join(genus.split("_"))])
        count_sum = 0
        new_line = []
        for sample in sorted(genus_dictionary[primer]):
            count = genus_dictionary[primer][sample][genus]
            count_sum += count
            new_line.append(str(count))
        genus_table[-1].append(str(count_sum))
        genus_table[-1] += new_line

    # Phylum
    for phylum in phylum_found[primer]:
        phylum_table.append([",".join(phylum.split("_"))])
        count_sum = 0
        new_line = []
        for sample in sorted(phylum_dictionary[primer]):
            count = phylum_dictionary[primer][sample][phylum]
            count_sum += count
            new_line.append(str(count))
        phylum_table[-1].append(str(count_sum))
        phylum_table[-1] += new_line

    # Print results to file
    # Species
    with open(
        os.path.join(output_folder, primer + "_species_table.csv"), "wt"
    ) as outfile:
        header = [line for line in species_table if line[0].startswith("Group")][0]
        outfile.write(",".join(header) + "\n")
        species_table = [
            line for line in species_table if not line[0].startswith("Group")
        ]
        counts = []

        for line in sorted(
            sorted(species_table, key=lambda x: int(x[3]), reverse=True),
            key=lambda x: x[0],
        ):

            # Add full taxon name
            if line[0] == "zMultiple":
                line[3:3] = ["MultipleHits"]
            else:
                line[3:3] = ["_".join(line[:3])]

            counts.append([int(x) for x in line[4:]])

            prepared_line = ",".join(line) + "\n"

            if max([int(x) for x in line[5:]]) >= min_coverage:
                if int(line[4]) >= min_experiment_coverage:
                    outfile.write(prepared_line)

        counts_by_sample = [sum(x) for x in zip(*counts)]
        prepared_line = (
            "Total,reads,by,sample,"
            + ",".join([str(x) for x in counts_by_sample])
            + "\n"
        )
        outfile.write(prepared_line)

    # Genus
    with open(
        os.path.join(output_folder, primer + "_genus_table.csv"), "wt"
    ) as outfile:
        for line in genus_table:
            prepared_line = ",".join(line) + "\n"

            if prepared_line.startswith("Group"):
                outfile.write(prepared_line)

            elif max([int(x) for x in line[2:]]) >= min_coverage:
                outfile.write(prepared_line)

    # Phylum
    with open(
        os.path.join(output_folder, primer + "_phylum_table.csv"), "wt"
    ) as outfile:
        for line in phylum_table:
            prepared_line = ",".join(line) + "\n"

            if prepared_line.startswith("Group"):
                outfile.write(prepared_line)

            elif max([int(x) for x in line[1:]]) >= min_coverage:
                outfile.write(prepared_line)

    # Export multiple hits infos
    with open(
        os.path.join(output_folder, primer + "_multiple_hit_infos.csv"), "wt"
    ) as outfile:
        for multiple_hit in multiple_hits_global_infos:
            for sequence in multiple_hits_global_infos[multiple_hit]:
                sample, seq = sequence
                outfile.write(",".join([multiple_hit, sample, seq]) + "\n")

# Output similarity values per species and site
with open(
    os.path.join(output_folder, "similarity_by_species_and_site.tsv"), "wt"
) as outfile:
    outfile.write("Sample\tPrimer\tSpecies\tSimilarity\tNumSequences\n")

    for s in sorted(similarity_dict):
        outfile.write("\t".join(list(s) + [str(similarity_dict[s])]) + "\n")
