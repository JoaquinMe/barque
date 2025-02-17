# Trimmomatic

Hace un filtrado de calidad de los reads

```bash
java -XX:ParallelGCThreads=1 -cp "$TRIMMOMATIC_JAR" org.usadellab.trimmomatic.TrimmomaticPE \
-phred33 \
"$DATA_FOLDER"/"$BASE"R1_001.fastq.gz \
"$DATA_FOLDER"/"$BASE"R2_001.fastq.gz \
"$TRIMMED_FOLDER"/"$BASE"R1_001.fastq.gz \
"$TRIMMED_FOLDER"/"$BASE"R1_001.single.fastq.gz \
"$TRIMMED_FOLDER"/"$BASE"R2_001.fastq.gz \
"$TRIMMED_FOLDER"/"$BASE"R2_001.single.fastq.gz \
LEADING:20 \ # Cortar bases al principio del read si están abajo de 20 de calidad
TRAILING:20 \ # Cortar bases al final del read si están abajo de 20 de calidad
SLIDINGWINDOW:20:20 \ # Hacer un análisis de ventana deslizante de 20 de largo, cortando el read si el promedio de calidad de la ventana baja de 20
MINLEN:"$MIN_HIT_LENGTH" \ # Remover reads que tengan una longitud menor que #MIN_HIT_LENGTH
CROP:"$CROP_LENGTH" \ # Cortar reads a esta longitud, sacando bases del final
-trimlog "$TRIMMED_FOLDER"/"${BASE%_}".log # hacer un archivo log
```

# Merge reads (FLASH)

"Mergea" reads, alineando el read F y el R. Formando una lectura mas larga.
Al configurar MIN_OVERLAP y MAX_OVERLAP, tener en cuenta el tamaño del amplicón y el de las lecturas.

```
F:            ----->
R:                <-----
merged:       ----------
```

```bash
flash -t 1 -z -m "$MIN_OVERLAP" -M "$MAX_OVERLAP" \
    {}R1_001.fastq.gz {}R2_001.fastq.gz \
    --to-stdout \> "$MERGED_FOLDER"/{/}merged.fastq.gz
#-t 1 -> 1 thread
#-m $MIN_OVERLAP overlap mínimo para mergear
#-M overlap máximo para mergear
```

# Split Amplicons

Usando el archivo primers.csv, extraer amplicones de los reads mergeados con posibilidad de utilizar primers degenerados.

- MAX_PRIMER_DIFF: cantidad de mismatches para permitir entre el primer y la secuencia

# Chimeras

a
