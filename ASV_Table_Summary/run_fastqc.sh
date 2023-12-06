for file in *_QC_R1.fastq
do
  stub=${file%_QC_R1.fastq}
  echo $stub
  file2=${stub}_QC_R2.fastq

  fastqc $file -o .
  fastqc $file2 -o .
done
