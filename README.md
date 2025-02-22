# :brain: ASiDentify

This repository contains the data and scripts required to run ASiDentify.

## Instructions
1. Download this repository to your local machine.
2. Some files are too large to store in this repository. Download the following files to the /data directory:
* Download Supplementary Data 3 from Gandal et al. 2022 [https://doi.org/10.1038/s41586-022-05377-7](https://doi.org/10.1038/s41586-022-05377-7). Save to the /data directory, named 'Gandal_2022.xlsx'.
* Download 'Gene Expression by Cluster, trimmed means' from [https://portal.brain-map.org/atlases-and-data/rnaseq/human-multiple-cortical-areas-smart-seq](https://portal.brain-map.org/atlases-and-data/rnaseq/human-multiple-cortical-areas-smart-seq). Save to the /data directory, named 'multiple_cortical_trimmed_means.csv'.

## Prepare input features and gene labels. 
3. Run `getRPKM_Burke.R` to generate 'Burke_metadata.csv' and 'Burke_timecourse.csv'.
4. Run `make_ASiD_dataset_PC.R` to generate the input data for ASiDentify. 
5. Run `generating_RBP_lists.R` to generate the 4 input dataframes for the ASiD model of RBPs (comprehensive, high-confidence, canonical and non-canonical). 
6. Run `make_ASiD_dataset_RBP.R` to generate the input data for ASiDentify - modified for RBPs. 
7. Run `make_ASiD_dataset_CR.R` to generate the input data for ASiDentify - modified for chromatin regulators. 

## Run ASiD.
8. Run `analyzing_PC_output.R` to analyze the output of ASiD. This script converts beta coefficients to odds ratios, calculates confidence intervals and AUROC and AUPRC values. The output data, including gene predictions, and lambda values used in the model are available in 'PC_genes_output.xlsx'. This script also produces Figure 3.
* Prerequisite: This script requires running section 1 from `nested_cv.R`.

9. Run `analyzing_RBP_output.R` to analyze the output of the modified ASiD model for RBPs. This script converts beta coefficients to odds ratios, calculates confidence intervals and AUROC and AUPRC values. The output data, including gene predictions, and lambda values used in the model are available in 'RBP_output.xlsx'. This script also produces Figures S3 and S4.
* Prerequisite: This script requires running section 2 from `nested_cv.R`.

10. Run `analyzing_CR_output.R` to analyze the output of the modified ASiD model for chromatin regulators. This script converts beta coefficients to odds ratios, calculates confidence intervals and AUROC and AUPRC values. The output data, including gene predictions, and lambda values used in the model are available in 'CR_output.xlsx'. This script also produces Figure S5.
* Prerequisite: This script requires running section 3 from `nested_cv.R`.

11. Run `analyzing_random_output.R` to analyze the ten random sets of 801 genes. This script converts beta coefficients to odds ratios, and calculates confidence intervals. The output data are available in 'Random_output.xlsx'. 
* Prerequisite: This script requires running `nested_cv_random.R`.

## Further analysis figures, and tables.
12. Run `making_figs.R` to produce Figures 2,4-8, Figures S1-2 and Tables S3-6.



### Notes
* Running `nested_cv.R` uses `nested_cv_funcs.R`, which contains the script to run the ASiDentify model. 
* All .xlsx output files are also saved as .csv files in the /tables directory. 
* All code was run on R version 4.3.3.
* [![DOI](https://zenodo.org/badge/864926068.svg)](https://doi.org/10.5281/zenodo.14901532)

