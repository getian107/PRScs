#!/usr/bin/env python

"""
PRS-CS: a polygenic prediction method that infers posterior SNP effect sizes under continuous shrinkage (CS) priors
using GWAS summary statistics and an external LD reference panel.

Reference: T Ge, CY Chen, Y Ni, YCA Feng, JW Smoller. Polygenic Prediction via Bayesian Regression and Continuous Shrinkage Priors.
           Nature Communications, 10:1776, 2019.


Usage:
python PRScs.py --ref_dir=PATH_TO_REFERENCE --bim_prefix=VALIDATION_BIM_PREFIX --sst_file=SUM_STATS_FILE --n_gwas=GWAS_SAMPLE_SIZE --out_dir=OUTPUT_DIR
                [--a=PARAM_A --b=PARAM_B --phi=PARAM_PHI --n_iter=MCMC_ITERATIONS --n_burnin=MCMC_BURNIN --thin=MCMC_THINNING_FACTOR
                 --chrom=CHROM --write_psi=WRITE_PSI --write_pst=WRITE_POSTERIOR_SAMPLES --seed=SEED]

"""


import os
import sys
import getopt

import parse_genet
import mcmc_gtb
import gigrnd


def parse_param():
    long_opts_list = ['ref_dir=', 'bim_prefix=', 'sst_file=', 'a=', 'b=', 'phi=', 'n_gwas=',
                      'n_iter=', 'n_burnin=', 'thin=', 'out_dir=', 'chrom=', 'beta_std=', 'write_psi=', 'write_pst=', 'seed=', 'help']

    param_dict = {'ref_dir': None, 'bim_prefix': None, 'sst_file': None, 'a': 1, 'b': 0.5, 'phi': None, 'n_gwas': None,
                  'n_iter': 1000, 'n_burnin': 500, 'thin': 5, 'out_dir': None, 'chrom': range(1,23),
                  'beta_std': 'FALSE', 'write_psi': 'FALSE', 'write_pst': 'FALSE', 'seed': None}

    print('\n')

    if len(sys.argv) > 1:
        try:
            opts, args = getopt.getopt(sys.argv[1:], "h", long_opts_list)          
        except:
            print('Option not recognized.')
            print('Use --help for usage information.\n')
            sys.exit(2)

        for opt, arg in opts:
            if opt == "-h" or opt == "--help":
                print(__doc__)
                sys.exit(0)
            elif opt == "--ref_dir": param_dict['ref_dir'] = arg
            elif opt == "--bim_prefix": param_dict['bim_prefix'] = arg
            elif opt == "--sst_file": param_dict['sst_file'] = arg
            elif opt == "--a": param_dict['a'] = float(arg)
            elif opt == "--b": param_dict['b'] = float(arg)
            elif opt == "--phi": param_dict['phi'] = float(arg)
            elif opt == "--n_gwas": param_dict['n_gwas'] = int(arg)
            elif opt == "--n_iter": param_dict['n_iter'] = int(arg)
            elif opt == "--n_burnin": param_dict['n_burnin'] = int(arg)
            elif opt == "--thin": param_dict['thin'] = int(arg)
            elif opt == "--out_dir": param_dict['out_dir'] = arg
            elif opt == "--chrom": param_dict['chrom'] = arg.split(',')
            elif opt == "--beta_std": param_dict['beta_std'] = arg.upper()
            elif opt == "--write_psi": param_dict['write_psi'] = arg.upper()
            elif opt == "--write_pst": param_dict['write_pst'] = arg.upper()
            elif opt == "--seed": param_dict['seed'] = int(arg)
    else:
        print(__doc__)
        sys.exit(0)

    if param_dict['ref_dir'] == None:
        print('* Please specify the directory to the reference panel using --ref_dir\n')
        sys.exit(2)
    elif param_dict['bim_prefix'] == None:
        print('* Please specify the directory and prefix of the bim file for the target dataset using --bim_prefix\n')
        sys.exit(2)
    elif param_dict['sst_file'] == None:
        print('* Please specify the summary statistics file using --sst_file\n')
        sys.exit(2)
    elif param_dict['n_gwas'] == None:
        print('* Please specify the sample size of the GWAS using --n_gwas\n')
        sys.exit(2)
    elif param_dict['out_dir'] == None:
        print('* Please specify the output directory using --out_dir\n')
        sys.exit(2)

    for key in param_dict:
        print('--%s=%s' % (key, param_dict[key]))

    print('\n')
    return param_dict


def main():
    param_dict = parse_param()

    for chrom in param_dict['chrom']:
        print('##### process chromosome %d #####' % int(chrom))

        if '1kg' in os.path.basename(param_dict['ref_dir']):
            ref_dict = parse_genet.parse_ref(param_dict['ref_dir'] + '/snpinfo_1kg_hm3', int(chrom))
        elif 'ukbb' in os.path.basename(param_dict['ref_dir']):
            ref_dict = parse_genet.parse_ref(param_dict['ref_dir'] + '/snpinfo_ukbb_hm3', int(chrom))

        vld_dict = parse_genet.parse_bim(param_dict['bim_prefix'], int(chrom))

        sst_dict = parse_genet.parse_sumstats(ref_dict, vld_dict, param_dict['sst_file'], param_dict['n_gwas'])

        ld_blk, blk_size = parse_genet.parse_ldblk(param_dict['ref_dir'], sst_dict, int(chrom))

        mcmc_gtb.mcmc(param_dict['a'], param_dict['b'], param_dict['phi'], sst_dict, param_dict['n_gwas'], ld_blk, blk_size,
            param_dict['n_iter'], param_dict['n_burnin'], param_dict['thin'], int(chrom), param_dict['out_dir'], param_dict['beta_std'],
	    param_dict['write_psi'], param_dict['write_pst'], param_dict['seed'])

        print('\n')


if __name__ == '__main__':
    main()


