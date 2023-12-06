for file in *_L001_R1_001.fastq
do
  stub=${file%_L001_R1_001.fastq}
  echo $stub
  file2=${stub}_L001_R2_001.fastq

  fastq_quality_filter -q 25 -p 70 -i $file -o ${stub}_QC_R1.fastq
  fastq_quality_filter -q 25 -p 70 -i $file2 -o ${stub}_QC_R2.fastq
done
