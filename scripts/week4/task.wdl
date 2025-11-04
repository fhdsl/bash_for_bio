version 1.0
  input {
    File ref_fasta
    File ref_bwt
    File ref_pac
    File ref_amb
    File ref_ann
    File ref_sa
    Array[File] fastq_r1
  }

  scatter (fq in fastq_r1){
    call bwa_align {
        input:
            reference_bwt = ref_bwt,
            reference_ann = ref_ann,
            reference_pac = ref_pac,
            reference_amb = ref_amb,
            reference_sa = ref_sa,
            reference_fasta = ref_fasta,
            reads = fq
        }

  }

output {
    Array[File] sam_file = bwa_align.bwa_output
  }


task bwa_align {
  meta {
    description: "Aligns paired-end reads to a reference using BWA-MEM"
    outputs: {
      bwa_output: "Output SAM file containing aligned reads"
    }
  }

  parameter_meta {
    reference_fasta: "Reference genome FASTA file",
    reads: "FASTQ of forward (R1) reads"
  }

  input {
    File reference_fasta
    File reference_bwt
    File reference_ann
    File reference_pac
    File reference_amb
    File reference_sa
    File reads
  }

  # Get name of reference fasta file within 'bwa_index' folder
  String ref_name = reference_fasta
  String out_name = basename(reads)
  command <<<
    # BWA alignment using 8 threads
    bwa mem -t 8 "~{ref_name}" "~{reads}" > "~{out_name}.sam"
  >>>
  output {
    File bwa_output = "~{out_name}.sam"
  }

  runtime {
    docker: "getwilds/bwa:0.7.17"
    cpu: 8
    memory: "16 GB"
  }
}