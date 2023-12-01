library(tidyverse)

file <- "~/Documents/SeaFlow/curation/SeaFlow_curation_parameters.tsv"
ctable <- read.table(file, sep = '\t', header = TRUE)

cruise_list <- unique(ctable$cruise)

# Export individual files

for (cruise in cruise_list){
  print(cruise)
  dir.create(paste0("~/Documents/Seaflow/curation/cruises/", cruise, "/"))
  ind <- which(ctable$cruise == cruise)
  this_line <- ctable[ind, c(1, 5:29)]
  file_name <- paste0("~/Documents/Seaflow/curation/cruises/", cruise, "/", cruise, '.curation_params.tsv')
  write.table(this_line, file_name, sep = "\t", row.names = FALSE)
}