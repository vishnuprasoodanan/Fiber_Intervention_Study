#Analysis using 16S Greengenes database provided by Kraken2
for R1 in *_R1.fastq;
do
    # We cut the sample name on the first _ and store it in $sample 
    sample=$(basename $R1 | cut -f1 -d.)

    # We infer the second pair replacing _R1 with _R2
     R2=${R1/_R1/_R2}
    
     echo "running trimmomatic on $sample"
     echo "Reads: $R1 / $R2"
     kraken2 --threads 16 --db ../KRAKEN2_DATABASES/16S_Greengenes_k2db --report "${sample}_kraken.rpt" --paired "$R1" "$R2" > "${sample}_kraken.tsv"
done

#Analysis using 16S SILVA database provided by Kraken2
for R1 in *_R1.fastq;
do
    # We cut the sample name on the first _ and store it in $sample 
    sample=$(basename $R1 | cut -f1 -d.)

    # We infer the second pair replacing _R1 with _R2
     R2=${R1/_R1/_R2}
    
     echo "running trimmomatic on $sample"
     echo "Reads: $R1 / $R2"
     kraken2 --threads 16 --db ../KRAKEN2_DATABASES/16S_SILVA138_k2db --report "${sample}_kraken.rpt" --paired "$R1" "$R2" > "${sample}_kraken.tsv"
done

