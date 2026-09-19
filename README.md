# CKM Stage 2 phenotypes and mortality (NHANES 2007–2018)

Code for: *Phenotypic heterogeneity and mortality in cardiovascular-kidney-metabolic syndrome Stage 2: a survey-weighted cluster analysis of NHANES 2007–2018.*

## Files
- `CKM_Stage2_FINAL_pipeline.R` — downloads NHANES and NCHS linked mortality data, assigns CKM stages, derives phenotypes, and runs all analyses.
- `Make_Supplement_and_Package.R` — builds the figures, Table 1, and the Supplementary Material.

## Requirements
R ≥ 4.6 and an internet connection on first run. Required packages install automatically.

## Run
```r
source("CKM_Stage2_FINAL_pipeline.R")
source("Make_Supplement_and_Package.R")
```
The random seed is fixed, so results are reproducible. No data are stored in this repository; all data are public and are downloaded from NCHS by the script.

## License
MIT
