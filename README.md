# PRS-CS

**PRS-CS** is a Python based command line tool that infers posterior SNP effect sizes under continuous shrinkage (CS) priors
using GWAS summary statistics and an external LD reference panel. Details of the method are described in the bioRxiv preprint:

T Ge, CY Chen, Y Ni, YCA Feng, JW Smoller. Polygenic Prediction via Bayesian Regression and Continuous Shrinkage Priors. bioRxiv preprint, doi: https://doi.org/10.1101/416859, 2018.
 

## Getting Started

- Clone this repository using the following git command:
   
    `git clone https://github.com/getian107/PRScs.git`

    Alternatively, download the source files from the github website (`https://github.com/getian107/PRScs`)

- Download the [LD reference](https://www.dropbox.com/s/v2n1jjwnihk3nwb/ldblk_1kg.tar.gz?dl=0 "LD reference") (~4.5G) and extract the files using e.g. `tar -zxvf ldblk_1kg.tar.gz`

- PRScs requires Python packages **scipy** (https://www.scipy.org/) and **h5py** (https://www.h5py.org/) installed.
 
- Once Python and its dependencies have been installed, you can run

    `./PRScs.py --help` or `./PRScs.py -h`

    to print a list of command-line options.


## Using PRS-CS

`
python PRScs.py --ref_dir=PATH_TO_REFERENCE --bim_prefix=VALIDATION_BIM_PREFIX --sst_file=SUM_STATS_FILE --n_gwas=GWAS_SAMPLE_SIZE --out_dir=OUTPUT_DIR [--a=PARAM_A --b=PARAM_B --phi=PARAM_PHI --n_iter=MCMC_ITERATIONS --n_burnin=MCMC_BURNIN --thin=MCMC_THINNING_FACTOR --CHROM=CHROM]
`
 - PATH_TO_REFERENCE (required): Full path to the folder `ldblk_1kg` that contains information on the LD reference panel (`snpinfo_1kg_hm3` and `ldblk_1kg_chr*.hdf5`).

 - VALIDATION_BIM_PREFIX (required): Full path and the prefix of the bim file for the validation set. 

 - SUM_STATS_FILE (required): Full path and the file name of the GWAS summary statistics. The summary statistics file must have the following format (including the header line):

```
    SNP          A1   A2   BETA      P
    rs4970383    C    A    -0.0064   4.7780e-01
    rs4475691    C    T    -0.0145   1.2450e-01
    rs13302982   A    G    -0.0232   2.4290e-01
    ...
```
Or:
```
    SNP          A1   A2   OR        P
    rs4970383    A    C    0.9825    0.5737                 
    rs4475691    T    C    0.9436    0.0691
    rs13302982   A    G    1.1337    0.0209
    ...
```
where SNP is the rs ID, A1 is the reference/effect allele, A2 is the alternative allele, BETA/OR is the effect/odds ratio of the reference allele, P is the p-value of the effect.

 - GWAS_SAMPLE_SIZE (required): Sample size of the GWAS.

 - OUTPUT_DIR (required): Output directory of the posterior effect size estimates.

 - PARAM_A (optional): Parameter a in the gamma-gamma prior. Default is 1.

 - PARAM_B (optional): Parameter b in the gamma-gamma prior. Default is 0.5.

 - PARAM_PHI (optional): Global shrinkage parameter phi. If phi is not specified, it will be learnt from the data using a full Bayesian approach. This usually works well for polygenic traits in real applications. For ultra-sparse genetic architectures, fixing phi to a small number may have improved prediction accuracy. Alternatively, a small-scale grid search (e.g., phi=1e-6, 1e-4, 1e-2, 1) can be used to find the phi value that produces the most accurate prediction.

 - MCMC_ITERATIONS (optional): Total number of MCMC iterations. Default is 1,000.

 - MCMC_BURNIN (optional): Number of burnin iterations. Default is 500.

 - MCMC_THINNING_FACTOR (optional): Thinning of the Markov chain. Default is 5.

 - CHROM (optional): The chromosome on which the model is fitted. Parallel computation for the 22 autosomes is recommended. Default is iterating through 22 autosomes (can be time-consuming).


## Output

PRS-CS writes posterior SNP effect size estimates for each chromosome to the user-specified directory. The output file contains chromosome, rs ID, base pair, A1, A2 and posterior effect size estimate for each SNP. Polygenic prediction can be produced by using `PLINK`'s `--score` command.


## Test Data

The test data contains GWAS summary statistics and a bim file for 1,000 SNPs on chromosome 22.
An example to use the test data:

`
python PRScs.py --ref_dir=path_to_ref --bim_prefix=path_to_bim/test --sst_file=path_to_sumstats/sumstats.txt --n_gwas=200000 --chrom=22 --out_dir=path_to_output
`


## Support

Please direct any problems or questions to Tian Ge (tge1@mgh.harvard.edu).


