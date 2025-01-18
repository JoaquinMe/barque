import pandas as pd
import os
import subprocess

index = pd.read_csv("14_tests/tests.csv")

for i,row in index.iterrows():
    bicho=row["bicho"]
    marker=row["marker"]
    #descargar todos los accessions
    result_descarga= subprocess.run(["./descargar_multiple.sh",f"14_tests/{marker}/{marker}_accessions.txt","04_data"], env=os.environ)
    #rename
    result_rename= subprocess.run(["python","rename_script.py","04_data",bicho],env= os.environ)
    #arrancar barque
    result_barque= subprocess.run(["./barque",f"14_tests/{marker}/{marker}_barque_config.sh"])
    
