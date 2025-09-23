#/bin/bash
input_fastq=${1}
base_file_name="${input_fastq\.fastq}"
sample_name="SM:${base_file_name}"
read_group_id="ID:${base_file_name}"
platform_info="PL:Illumina"
ref_fasta_local="/shared/biodata/reference/iGenomes/Homo_sapiens/UCSC/hg19/Sequence/BWAIndex/genome.fa"

bwa mem \
      -p -v 3 -M \
      -R '@RG\t~{read_group_id}\t~{sample_name}\t~{platform_info}' \
      ~{ref_fasta_local} ~{input_fastq} > \
      ~{base_file_name}.sam 