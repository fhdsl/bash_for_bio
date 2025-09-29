#!/bin/bash
module load fhPython
python3 process_file.py ${1}
module purge