# PRS-CS
`PRS-CS` is Python based command line tool that infers posterior SNP effect sizes under continuous shrinkage (CS) priors
using GWAS summary statistics and an external LD reference panel. Details of the method are described in the bioRxiv preprint:

T Ge, CY Chen, Y Ni, YCA Feng, JW Smoller. Polygenic Prediction via Bayesian Regression and Continuous Shrinkage Priors. bioRxiv preprint, doi: https://doi.org/10.1101/416859, 2018.
 

## Getting Started


### Requirements


### Installing PRS-CS


## Using PRS-CS

`
python PRScs.py --ref_dir=PATH_TO_REFERENCE --bim_prefix=VALIDATION_BIM_PREFIX --sst_file=SUM_STATS_FILE --n_gwas=GWAS_SAMPLE_SIZE --out_dir=OUTPUT_DIR [--a=PARAM_A --b=PARAM_B --phi=PARAM_PHI --n_iter=MCMC_ITERATIONS --n_burnin=MCMC_BURNIN --thin=MCMC_THINNING_FACTOR --CHROM=CHROM]
`
 - PATH_TO_REFERENCE: Full path to the folder ldblk_1kg that contains information on the LD reference panel (snpinfo_1kg_hm3 and ldblk_1kg_chr*.hdf5).

 - VALIDATION_BIM_PREFIX: Full path and the prefix of the bim file for the validation set. 

 - SUM_STATS_FILE: Full path and the file name of the GWAS summary statistics.
                   Summary statistics file must have the following format:



## Support
Please report any problems or questions to Tian Ge (tge1@mgh.harvard.edu).
