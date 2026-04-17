#!/bin/bash
for file in ./data/*.fastq
do
  wc $file
done