"""
Pearson correlation cupy benchmark
KH Dec 16, 2019
"""
import numpy as np
import cupy as cp
import timeit

# load test data
mat = np.load(snakemake.input[0])

# benchmark method
num_times = int(snakemake.config['benchmark']['num_times'])

try:
    times = timeit.repeat("cp.asnumpy(cp.corrcoef(cp.array(mat)))", 
                        globals=globals(), number=1, repeat=num_times)
except:
    # if graphics card runs out of memory, return empty result
    times = [np.nan]

# save timings
with open(snakemake.output['timings'], 'w') as fp:
    entry = ['Pearson', 'Python', 'cupy.corrcoef', snakemake.wildcards['num_rows'], str(min(times))]
    fp.write(", ".join(entry) + "\n")

# store correlation matrix result for comparison
try:
    cor_mat = cp.asnumpy(cp.corrcoef(cp.array(mat)))
    cor_mat[np.tril_indices_from(cor_mat)] = np.nan
except:
    num_rows = int(snakemake.wildcards['num_rows'])
    cor_mat = np.repeat(np.nan, num_rows**2).reshape((num_rows, num_rows))

np.save(snakemake.output['cor_mat'], cor_mat)
