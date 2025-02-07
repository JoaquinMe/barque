https://github.com/JoaquinMe/barque

Barque is a fast eDNA metabarcoding analysis pipeline that first denoises and then annotates ASVs or OTUs, using high-quality barcoding databases.
Barque can produce denoised OTUs and annotate them using a custom database.
These annotated OTUs can then be used as a database themselves to find read counts per OTU per sample,
effectively annotating the reads with the OTUs that were previously found. In this process, some of the OTUs are annotated to the species level,
some to the genus or higher levels.

# Instalación

```bash
git clone https://github.com/JoaquinMe/barque
```

# Requerimientos

**Barque solo funciona en sistemas GNU-Linux o MacOS**

Usando conda(https://docs.conda.io/projects/conda/en/stable/user-guide/install/index.html)

```
cd barque/
conda env create -f environment.yaml

```

# Archivos de configuración

## 02_info/primers.csv

Este es un archivo tipo csv que contiene la información de los primers que se van a usar para el procedimiento.
Cada fila es un primer diferente. Las filas comienzan con un "#" si no van a ser usadas para el análisis.

En 14_tests hay modelos de estos archivos para cada primer.

### Ejemplo

| PrimerName            | ForwardSeq                 | ReverseSeq                  | MinAmpliconSize | MaxAmpliconSize | DatabaseName | SimilSpecies | SimilGenus | SimilPhylum |
| --------------------- | -------------------------- | --------------------------- | --------------- | --------------- | ------------ | ------------ | ---------- | ----------- |
| 12s200pb_MiFishU      | GTCGGTAAAACTCGTGCCAGC      | CATAGTGGGGTATCTAATCCCAGTTTG | 150             | 250             | 12S          | 0.98         | 0.9        | 0.85        |
| #COICOIintF_jgHCO2198 | GGWACWGGWTGAACWGTWTAYCCYCC | TAIACYTCIGGRTGICCRAARAAYCA  | 300             | 325             | bold         | 0.97         | 0.9        | 0.85        |

- PrimerName: El nombre del primer, los resultados van a tener este nombre como prefijo (un archivo fasta por nombre de primer)
- ForwardSeq: La secuencia del primer F
- ReverseSeq: La secuencia del primer R
- Min/Max AmpliconSize: El tamaño mínimo y maximo del amplicón. Es recomendable ajustar este parámetro si hay caída de lecturas en el paso de extracción de amplicones
- DatabaseName: El nombre de la base de datos como está nombrada en 03_databases. Si en este campo pongo "16S", el programa va a buscar un archivo en 03_databases/ llamado 16S.fasta.gz
- SimilSpecies,SimilGenus,SimilPhylum: Estos son los porcentajes de identidad que el programa usa para asignación de taxonomía.

## 02_info/barque_config.sh

- NCPUS: número de CPUs para usar. Muchos de los procesos de barque se paralelizan.
- PRIMER_FILE: Archivo con la información de los primers. Recomendable no cambiar esta opción.
- SKIP_DATA_PREP: 1 para saltear preparación de datos. 0 para correr todo el pipeline.
- CROP_LENGTH: cortar reads a esta longitud despues de filtrar.
- MIN_OVERLAP:
- MAX_OVERLAP:
- MAX_PRIMER_DIFF:
- SKIP_CHIMERA_DETECTION:
- MAX_ACCEPTS:
- MAX_REJECTS:
- QUERY_COV:
- MIN_HIT_LENGTH:
- MIN_HITS_SAMPLE:
- MIN_HITS_EXPERIMENT:
- NUM_NON_ANNOTATED_SEQ:
- MIN_DEPTH_MULTI:
- SKIP_OTUS:
- MIN_SIZE_FOR_OTU:
