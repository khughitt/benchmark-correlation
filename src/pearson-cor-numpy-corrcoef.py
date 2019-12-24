"""
Pearson correlation numpy.corrcoef benchmark
KH Dec 17, 2019
"""
import numpy as np
import timeit

# load test data
mat = np.load(snakemake.input[0])

# benchmark method
num_times = int(snakemake.config['benchmark']['num_times'])
times = timeit.repeat("np.corrcoef(mat)", globals=globals(), number=1, repeat=num_times)

# save timings
with open(snakemake.output['timings'], 'w') as fp:
    entry = ['Pearson', 'Python', 'numpy.corrcoef', snakemake.wildcards['num_rows'], str(min(times))]
    fp.write(", ".join(entry) + "\n")

# store correlation matrix result for comparison
cor_mat = np.corrcoef(mat)
cor_mat[np.tril_indices_from(cor_mat)] = np.nan

np.save(snakemake.output['cor_mat'], cor_mat)
