# Barque

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

```bash
cd barque/
conda env create --name barque --file environment.yaml
conda activate barque
```

# Correr Tests

```bash
git clone https://github.com/JoaquinMe/barque_tests
cd  barque_tests/
conda activate barque
python run_tests.py
```

# Correr Barque

- Copiar las muestras demultiplexadas a 04_data.
- Se necesitan un par de archivos por muestra. (F y R)
- Las secuencias contenidas en estos archivos deben contener los primers usados en el PCR.
- Es posible que se tenga que demultiplexar antes de correr Barque.

_El nombre de los archivos debe cumplir con este formato:_

```
MuestraID_*_R1_001.fastq.gz
MuestraID_*_R2_001.fastq.gz
```

- Cambiar los archivos de configuración a preferencia.
- Ejecutar barque:

```bash
./barque 02_info/barque_config.sh
```

## Nota:

- Cada nombre de muestra, o MuestraID no tiene que contener ningún guión bajo (\_).
- El asterisco (\*) puede ser cualquier texto alfanumérico

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

### Parámetros globales

- NCPUS: número de CPUs para usar. Muchos de los procesos de barque se paralelizan.
- PRIMER_FILE: Archivo con la información de los primers. Recomendable no cambiar esta opción.

### Saltear preparación de datos

- SKIP_DATA_PREP: 1 para saltear preparación de datos. 0 para correr todo el pipeline.

### Filtrado con Trimmomatic

- CROP_LENGTH: cortar reads a esta longitud despues de filtrar.

### Merge reads con flash

- MIN_OVERLAP: número mínimo de nucleótidos que se pueden superponer al mergear (FLASH)
- MAX_OVERLAP: número máximo de nucleótidos que se pueden superponer al mergear (FLASH)

### Extracción de amplicones

- MAX_PRIMER_DIFF: número máximo de diferencias entre el primer y la secuencia

### Quimeras

- SKIP_CHIMERA_DETECTION: 1 para saltear detección de quimeras. 0 para detectar quimeras.

### vsearch

- MAX_ACCEPTS: controla cuantos matches son aceptados antes de parar la búsqueda.
- MAX_REJECTS: controla cuantos no-matches son permitidos antes de parar la búsqueda.
- QUERY_COV: determina cuanta cobertura del query sobre el match para considerarlo accept o reject.
- MIN_HIT_LENGTH: limita la cantidad de matches reportados por query

### Filtros

- MIN_HITS_SAMPLE: es el número mínimo de coincidencias que un taxón debe tener en al menos una muestra para que el taxón se conserve
- MIN_HITS_EXPERIMENT: es el número mínimo de coincidencias que un taxón debe tener en todas las muestras combinadas para que el taxón se conserve.

### Non-annotated

- NUM_NON_ANNOTATED_SEQ: cantidad de reads no anotados mas frecuentas que se deben guardar.

### Multiple_hits

- MIN_DEPTH_MULTI: profundidad mínima con la cual reportar reads únicos por sample

### OTUs

- SKIP_OTUS: 1 para saltearse OTUs, 0 para usarla
- MIN_SIZE_FOR_OTU: tamaño mínimo del OTU para definirlo como tal
