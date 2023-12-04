# seaflow-curation

Contains the function and parameters to flag outliers in SeaFlow data.

### Run 'curateSF()

INPUT:
- `db`: Full path to SeaFlow SQL .db file, typically *.vct.db in the shared Google Drive snakemake directory
- `save_path`: Path to save *.outlier.tsv files 
- `show_plots`: Logical to show plots or not.  Default = FALSE. You will be prompted to press **ENTER** after each plot if `show_plots` = TRUE.

USAGE:
```
for (cruise in cruise_list){
  print(cruise)
  db <- paste0(path, cruise, '/', cruise, ".vct.db")
  curateSF(db, paste(local_path, "cruises/", cruise), show_plots = TRUE)}
```

Annette Hynes, Baker Van Buren, Chris Berthiaume, and François Ribalet contributed to this project.
