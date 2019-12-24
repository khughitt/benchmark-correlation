#!/bin/env python
"""
Efficient correlation implementation in python

Source: https://github.com/ikizhvatov/efficient-columnwise-correlation/blob/master/columnwise_corrcoef_perf.py
"""
import numpy as np
import timeit

def corrcoeff_einsum_optimized(O, P):
    (n, t) = O.shape      # n traces of t samples
    (n_bis, m) = P.shape  # n predictions for each of m candidates

    DO = O - (np.einsum("nt->t", O, optimize='optimal') / np.double(n)) # compute O - mean(O)
    DP = P - (np.einsum("nm->m", P, optimize='optimal') / np.double(n)) # compute P - mean(P)

    cov = np.einsum("nm,nt->mt", DP, DO, optimize='optimal')

    varP = np.einsum("nm,nm->m", DP, DP, optimize='optimal')
    varO = np.einsum("nt,nt->t", DO, DO, optimize='optimal')
    tmp = np.einsum("m,t->mt", varP, varO, optimize='optimal')

    return cov / np.sqrt(tmp)

# main
mat = np.load(snakemake.input[0])

# numpy
num_times = int(snakemake.config['benchmark']['num_times'])

# benchmark
times = timeit.repeat("corrcoeff_einsum_optimized(mat.T, mat.T)", 
                      globals=globals(), number=1, repeat=num_times)

# save timings
with open(snakemake.output['timings'], 'w') as fp:
    entry = ['Pearson', 'Python', 'numpy.einsum', snakemake.wildcards['nrows'], str(min(times))]
    fp.write(", ".join(entry) + "\n")

# store correlation matrix result for comparison
cor_mat = corrcoeff_einsum_optimized(mat.T, mat.T)
cor_mat[np.tril_indices_from(cor_mat)] = np.nan

np.save(snakemake.output['cor_mat'], cor_mat)
