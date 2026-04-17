library(tidyverse)
library(janitor)

args <- commandArgs()

# Use arg$CSVFILE in read.csv
csv_file <- read.csv(file=args$input_file)

# Do some work with csv_file
csv_filtered <- csv_file |>  summary()

# Write output
write.csv(csv_filtered, file = paste0(args$input_file, "_summary.csv"))