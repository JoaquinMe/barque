import os
import re
import subprocess

import pandas as pd

index = pd.read_csv("14_tests/tests.csv")

for i, row in index.iterrows():
    bicho = row["bicho"]
    marker = row["marker"]

    result_rm_config = subprocess.run(
        ["rm", "02_info/barque_config.sh", "02_info/primers.csv"]
    )

    result_cp_config= subprocess.run(["cp",f"14_tests/{marker}/{marker}_barque_config.sh","02_info/barque_config.sh"])
    result_cp_primers=subprocess.run(["cp",f"14_tests/{marker}/{marker}_primers.csv","02_info/primers.csv"])

    # descargar todos los accessions
    result_descarga = subprocess.run(
        [
            "./descargar_multiple.sh",
            f"14_tests/{marker}/{marker}_accessions.txt",
            "04_data"
        ],
        env=os.environ,
    )
    # rename
    result_rename = subprocess.run(
        ["python", "rename_script.py", "04_data", bicho], env=os.environ
    )
    # arrancar barque
    result_barque = subprocess.run(
        ["./barque", f"02_info/barque_config.sh"], env=os.environ
    )

    # Cleanup
    print("Cleanup")
    result_zip = subprocess.run(
        [
            "tar",
            "-zcvf",
            f"14_tests/test_results/results_{marker}.tar.gz",
            "12_results/"
        ],env=os.environ
    )
    result_logs = subprocess.run(
        ["tar", "-zcvf", f"14_tests/test_results/logs_{marker}.tar.gz", "99_logfiles/"],env=os.environ
    )
    result_rm = subprocess.run(
        [
            "rm",
            "-rf",
            "04_data/*",
            "05_trimmed/*",
            "06_merged/*",
            "07_split_amplicons/*",
            "08_chimeras/*",
            "09_vsearch/*",
            "10_read_dropout/*",
            "11_non_annotated/*",
            "12_results/*",
            "13_otu_database/*",
            "99_logfiles/*",
        ],env=os.environ
    )

    result_hash_results = subprocess.run(
        ["shasum", f"14_tests/test_results/results_{marker}.tar.gz"],
        capture_output=True,env=os.environ
    )
    checksum_results = re.match(r"^(\w*) ", result_hash_results.stdout.decode("utf-8"))
    result_hash_logs = subprocess.run(
        ["shasum", f"14_tests/test_results/logs_{marker}.tar.gz"], capture_output=True,env=os.environ
    )
    checksum_logs = re.match(r"^(\w*) ", result_hash_results.stdout.decode("utf-8"))

    index.loc[i, "checksum_results"] = checksum_results
    index.loc[i, "checksum_logs"] = checksum_logs

index.to_csv("14_tests/test_results.csv", index=False, header=True)
