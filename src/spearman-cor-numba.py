"""
spearman correlation using numba
https://stackoverflow.com/questions/52371329/fast-spearman-correlation-between-two-pandas-dataframes
"""
from numba import njit
import pandas as pd
import numpy as np
import timeit

@njit
def mean1(a):
  n = len(a)
  b = np.empty(n)
  for i in range(n):
    b[i] = a[i].mean()
  return b

@njit
def std1(a):
  n = len(a)
  b = np.empty(n)
  for i in range(n):
    b[i] = a[i].std()
  return b

@njit
def spearman_cor(a, b):
    ''' Correlation '''
    n, k = a.shape
    m, k = b.shape

    mu_a = mean1(a)
    mu_b = mean1(b)
    sig_a = std1(a)
    sig_b = std1(b)

    out = np.empty((n, m))

    for i in range(n):
        for j in range(m):
            out[i, j] = (a[i] - mu_a[i]) @ (b[j] - mu_b[j]) / k / sig_a[i] / sig_b[j]

    return out

# load test data
mat = np.load(snakemake.input[0])

# benchmark method
num_times = int(snakemake.config['benchmark']['num_times'])
times = timeit.repeat("spearman_cor(mat, mat)", globals=globals(), number=1, repeat=num_times)

# save timings
with open(snakemake.output['timings'], 'w') as fp:
    entry = ['Spearman', 'Python', 'numba', snakemake.wildcards['nrows'], str(min(times))]
    fp.write(", ".join(entry) + "\n")

# store correlation matrix result for comparison
cor_mat = spearman_cor(mat, mat)
cor_mat[np.tril_indices_from(cor_mat)] = np.nan

np.save(snakemake.output['cor_mat'], cor_mat)
