#!/usr/bin/env python

"""
Markov Chain Monte Carlo (MCMC) sampler for polygenic prediction with continuous shrinkage (CS) priors.

"""


import scipy as sp
from scipy import linalg 
from scipy import random
import gigrnd


def mcmc(a, b, phi, sst_dict, n, ld_blk, blk_size, n_iter, n_burnin, thin, chrom, out_dir, beta_std, seed):
    print('... MCMC ...')

    # seed
    if seed != None:
        random.seed(seed)

    # derived stats
    beta_mrg = sp.array(sst_dict['BETA'], ndmin=2).T
    maf = sp.array(sst_dict['MAF'], ndmin=2).T
    n_pst = (n_iter-n_burnin)/thin
    p = len(sst_dict['SNP'])
    n_blk = len(ld_blk)

    # initialization
    beta = sp.zeros((p,1))
    psi = sp.ones((p,1))
    sigma = 1.0
    if phi == None:
        phi = 1.0; phi_updt = True
    else:
        phi_updt = False

    beta_est = sp.zeros((p,1))
    psi_est = sp.zeros((p,1))
    sigma_est = 0.0
    phi_est = 0.0

    # MCMC
    for itr in range(1,n_iter+1):
        if itr % 100 == 0:
            print('--- iter-' + str(itr) + ' ---')

        mm = 0; quad = 0.0
        for kk in range(n_blk):
            if blk_size[kk] == 0:
                continue
            else:
                idx_blk = range(mm,mm+blk_size[kk])
                dinvt = ld_blk[kk]+sp.diag(1.0/psi[idx_blk].T[0])
                dinvt_chol = linalg.cholesky(dinvt)
                beta_tmp = linalg.solve_triangular(dinvt_chol, beta_mrg[idx_blk], trans='T') + sp.sqrt(sigma/n)*random.randn(len(idx_blk),1)
                beta[idx_blk] = linalg.solve_triangular(dinvt_chol, beta_tmp, trans='N')
                quad += sp.dot(sp.dot(beta[idx_blk].T, dinvt), beta[idx_blk])
                mm += blk_size[kk]

        err = max(n/2.0*(1.0-2.0*sum(beta*beta_mrg)+quad), n/2.0*sum(beta**2/psi))
        sigma = 1.0/random.gamma((n+p)/2.0, 1.0/err)

        delta = random.gamma(a+b, 1.0/(psi+phi))

        for jj in range(p):
            psi[jj] = gigrnd.gigrnd(a-0.5, 2.0*delta[jj], n*beta[jj]**2/sigma)
        psi[psi>1] = 1.0

        if phi_updt == True:
            w = random.gamma(1.0, 1.0/(phi+1.0))
            phi = random.gamma(p*b+0.5, 1.0/(sum(delta)+w))

        # posterior
        if (itr>n_burnin) and (itr % thin == 0):
            beta_est = beta_est + beta/n_pst
            psi_est = psi_est + psi/n_pst
            sigma_est = sigma_est + sigma/n_pst
            phi_est = phi_est + phi/n_pst

    # convert standardized beta to per-allele beta
    if beta_std == 'False':
        beta_est /= sp.sqrt(2.0*maf*(1.0-maf))

    # write posterior effect sizes
    if phi_updt == True:
        eff_file = out_dir + '_pst_eff_a%d_b%.1f_phiauto_chr%d.txt' % (a, b, chrom)
    else:
        eff_file = out_dir + '_pst_eff_a%d_b%.1f_phi%1.0e_chr%d.txt' % (a, b, phi, chrom)

    with open(eff_file, 'w') as ff:
        for snp, bp, a1, a2, beta in zip(sst_dict['SNP'], sst_dict['BP'], sst_dict['A1'], sst_dict['A2'], beta_est):
            ff.write('%d\t%s\t%d\t%s\t%s\t%.6e\n' % (chrom, snp, bp, a1, a2, beta))

    # print estimated phi
    if phi_updt == True:
        print('... Estimated global shrinkage parameter: %1.2e ...' % phi_est )

    print('... Done ...')


