version 1.0
workflow bwa_samtools_workflow {
  meta {
    description: "Workflow for aligning, sorting, and index paired-end reads"
    outputs: {
      sorted_bam: "Final sorted BAM alignment file",
      sorted_bai: "BAM index file for the sorted BAM"
    }
  }
  parameter_meta {
    ref_fasta: "Reference genome FASTA file (used to create the BWA index)"
    fastq_r1: "FASTQ file for forward (R1) reads"
    name: "Sample name for output files"
  }
  input {
    File ref_fasta
    File ref_bwt
    File ref_pac
    File ref_amb
    File ref_ann
    File ref_sa
    File fastq_r1
    String name
  }
  call bwa_align {
    input:
      reference_bwt = ref_bwt,
      reference_ann = ref_ann,
      reference_pac = ref_pac,
      reference_amb = ref_amb,
      reference_sa = ref_sa,
      reference_fasta = ref_fasta,
      reads = fastq_r1,
      name = name
  }
  call samtools_process {
    input:
      sam_file = bwa_align.bwa_output,
      name = name
  }
  output {
    File sorted_bam = samtools_process.sorted_bam
    File sorted_bai = samtools_process.sorted_bai
  }
}
task bwa_align {
  meta {
    description: "Aligns paired-end reads to a reference using BWA-MEM"
    outputs: {
      bwa_output: "Output SAM file containing aligned reads"
    }
  }
  parameter_meta {
    reference_fasta: "Reference genome FASTA file"
    reads: "FASTQ of forward (R1) reads"
    name: "Sample name for output SAM file"
  }
  input {
    File reference_fasta
    File reference_bwt
    File reference_ann
    File reference_pac
    File reference_amb
    File reference_sa
    File reads
    String name
  }
  # Get name of reference fasta file within 'bwa_index' folder
  String ref_name = reference_fasta
  command <<<
    # BWA alignment using 8 threads
    bwa mem -t 8 "~{ref_name}" "~{reads}" > "~{name}.sam"
  >>>
  output {
    File bwa_output = "~{name}.sam"
  }
  runtime {
    docker: "getwilds/bwa:0.7.17"
    cpu: 8
    memory: "16 GB"
  }
}
task samtools_process {
  meta {
    description: "Sorts SAM alignment file by coordinate and creates an index"
    outputs: {
      sorted_bam: "Coordinate-sorted BAM file",
      sorted_bai: "BAM index file (.bai) for the sorted BAM"
    }
  }
  parameter_meta {
    sam_file: "SAM alignment file to be sorted and indexed"
    name: "Sample name for output files"
  }
  input {
    File sam_file
    String name
  }
  command <<<
    # Sort SAM file and convert to BAM format using 8 threads
    samtools sort -@ 8 -o "~{name}_sorted.bam" "~{sam_file}"
    # Create index for sorted BAM file
    samtools index "~{name}_sorted.bam"
  >>>
  output {
    File sorted_bam = "~{name}_sorted.bam"
    File sorted_bai = "~{name}_sorted.bam.bai"
  }
  runtime {
    docker: "getwilds/samtools:1.19"
    cpu: 8
    memory: "16 GB"
  }
}
