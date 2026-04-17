library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)

# Use the first argument to read the file
csv_file <- read.csv(file=args[1])

# Do some work with csv_file
csv_filtered <- csv_file |>  count(tumor_stage)

# Write output using the first argument
write.csv(csv_filtered, file = paste0(args[1], "_summary.csv"))