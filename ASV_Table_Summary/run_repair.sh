for file in *_QC_R1.fastq
do
  stub=${file%_QC_R1.fastq}
  echo $stub
  file2=${stub}_QC_R2.fastq

  repair.sh in1=$file in2=$file2 out1=OUT/${stub}_QC_R1.fastq  out2=OUT/${stub}_QC_R2.fastq outs=OUT/${stub}_single.fastq
done
