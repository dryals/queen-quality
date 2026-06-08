# queen-quality
GS and GWAS on queen size and fertility traits

Dylan Ryals, Bradley Metz, Luiz Brito, David Tarpy, Brock Harpur

Sample and sequence data are currently available on NCBI under BioProject accession PRJNA1469566

- The script `1_data_prep.sh` is designed for SLURM protocol on an HPC. It filters raw data and performs initial analyses:
    - The `pheno` directory contains raw phenotypes for queens.
    - The `scripts` directory contains scripts for cleaning, filtering, and adjusting data.
    - Processed data are output to the `data` directory.
- The R markdown file `2_analysis.Rmd` runs analyses and generates figures, with `2_analysis.html` displaying all outputs.
    - Publication figures are output to the `fig` directory.
- Cleaned phenotypic data are temporarily available as `cleaned_pheno.csv` until permanant archiving upon publication.


