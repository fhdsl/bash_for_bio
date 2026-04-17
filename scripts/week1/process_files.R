library(tidyverse)
library(janitor)

args <- commandArgs()

# Use arg$CSVFILE in read.csv
csv_file <- read.csv(file=args$input_file)

# Do some work with csv_file
csv_filtered <- csv_file |> janitor::clean_names()

# Write output
write.csv(csv_filtered, file = paste0(args$CSVFILE, "_cleaned.csv"))