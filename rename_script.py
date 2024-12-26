import os
import sys

cwd = os.getcwd()
targetfolder = sys.argv[1]
project_name = sys.argv[2]

folder = cwd + "/" + targetfolder

files = os.listdir(folder)
files = [file for file in files if file.endswith(".fastq.gz")]

for file in files:
    sampleID = file.split("_")[0]
    R = (file.split(".")[0]).split("_")[1]
    os.rename(
        folder + "/" + file, f"{folder}/{sampleID}_{project_name}_R{R}_001.fastq.gz"
    )
