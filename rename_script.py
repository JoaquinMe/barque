import os
import sys

cwd = os.getcwd()
targetfolder = sys.argv[1]

folder = cwd + "/" + targetfolder

project_name = "stream"

files = os.listdir(folder)
for file in files:
    sampleID = file.split("_")[0]
    R = (file.split(".")[0]).split("_")[1]
    os.rename(
        folder + "/" + file, f"{folder}/{sampleID}_{project_name}_R{R}_001.fastq.gz"
    )
