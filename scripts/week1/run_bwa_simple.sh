#!/bin/bash
module load BWA/0.7.17-GCCcore-11.2.0
input_fastq=${1}

# strip path and suffix
base_file_name="${input_fastq%.fastq}"
base_file_name=${base_file_name##*/}

echo "running $input_fastq"

ref_fasta_local="/shared/biodata/reference/iGenomes/Homo_sapiens/UCSC/hg19/Sequence/BWAIndex/genome.fa"

bwa mem \
      -p -v 3 -M \
      -R "@RG\t test" \
      "${ref_fasta_local}" "${input_fastq}" > \
      "${base_file_name}.sam"

module purge
