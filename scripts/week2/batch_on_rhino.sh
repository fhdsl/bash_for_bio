#!/bin/bash
for file in *.fastq
do
  wc $file
done
