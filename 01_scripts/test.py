import gzip
import os
import sys

import pandas as pd

pid_cutoff = float(sys.argv[1]) * 100  # porcentaje de similaridad para cortar
df = pd.read_csv("01_scripts/test.csv")
df = df.loc[df["pid"] >= pid_cutoff]
taxtable = df["tax"].str.split("; ", expand=True)
taxtable_const = taxtable.copy()
keep = True
taxlist = []
i = 0
print(taxtable.head())
while keep and taxtable.shape[1]>i:
    col = taxtable.iloc[:, i]
    if col.unique().size == 1:
        keep = True
        taxlist.append(col.unique()[0])
        i+=1
    else:
        keep = False

print(taxlist)
