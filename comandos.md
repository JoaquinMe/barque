# Trimmomatic

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

SLIDINGWINDOW:20:20 \
MINLEN:"$MIN_HIT_LENGTH" \
CROP:"$CROP_LENGTH" \
-trimlog "$TRIMMED_FOLDER"/"${BASE%_}".log
```
