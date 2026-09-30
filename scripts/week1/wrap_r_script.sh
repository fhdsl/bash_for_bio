#!/bin/bash
module load fhR
Rscript process_data.R ${1}
module purge
