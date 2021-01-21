"""
Correlation speed benchmark
KH Dec 16, 2019
"""
import os
import numpy as np
import pandas as pd

# number of columns to use for test data (static)
num_cols = int(config['test_data']['num_cols'])

# determine expected output files; some methods cannot be run on all dataset
# sizes due to time/memory constraints
num_rows = {}

cor_methods = ['pearson_cor_numpy_cor_coef', 'pearson_cor_numpy_einsum',
               'pearson_cor_base_r', 'pearson_cor_julia', 'spearman_cor_numba']
               #'pearson_cor_cupy', 'pearson_cor_coop_pcor'] 
               

for method in cor_methods:
    num_rows[method] = config['test_data']['num_rows']

    if method in config['method_row_limits']:
        row_limit = config['method_row_limits'][method]
        num_rows[method] = [x for x in num_rows[method] if x <= row_limit]

out_dir = config['output_dir']

rule combine_timings:
    input:
        expand(os.path.join(out_dir, 'timings', 'pearson_cor_numpy_cor_coef', 'timings_{nrows}.csv'), nrows=num_rows['pearson_cor_numpy_cor_coef']),
        expand(os.path.join(out_dir, 'timings', 'pearson_cor_numpy_einsum', 'timings_{nrows}.csv'), nrows=num_rows['pearson_cor_numpy_einsum']),
        expand(os.path.join(out_dir, 'timings', 'pearson_cor_base_r', 'timings_{nrows}.csv'), nrows=num_rows['pearson_cor_base_r']),
        expand(os.path.join(out_dir, 'timings', 'pearson_cor_julia', 'timings_{nrows}.csv'), nrows=num_rows['pearson_cor_julia']),
        expand(os.path.join(out_dir, 'timings', 'spearman_cor_numba', 'timings_{nrows}.csv'), nrows=num_rows['spearman_cor_numba'])
        # expand(os.path.join(out_dir, 'timings', 'pearson_cor_coop_pcor', 'timings_{nrows}.csv'), nrows=num_rows['pearson_cor_coop_pcor']),
        # expand(os.path.join(out_dir, 'timings', 'pearson_cor_cupy', 'timings_{nrows}.csv'), nrows=num_rows['pearson_cor_cupy']),
    output: 
        os.path.join(out_dir, 'timings', 'all_timings.csv')
    run:
        dfs = []

        for infile in input:
            dfs.append(pd.read_csv(infile, header=None))

        dat = pd.concat(dfs)

        dat.columns = ['Method', 'Language', 'Implementation', 'Num Rows', 'Time (Secs)']
        dat.to_csv(output[0], index=False)

rule measure_cor_mat_differences:
    input:
        os.path.join(out_dir, 'cor_mats', 'pearson_cor_numpy_cor_coef', 'cor_mat_{nrows}.npy'),
        os.path.join(out_dir, 'cor_mats', 'pearson_cor_numpy_einsum', 'cor_mat_{nrows}.npy'),
        os.path.join(out_dir, 'cor_mats', 'pearson_cor_base_r', 'cor_mat_{nrows}.npy'),
        os.path.join(out_dir, 'cor_mats', 'spearman_cor_numba', 'cor_mat_{nrows}.npy')
        # os.path.join(out_dir, 'cor_mats', 'pearson_cor_cupy', 'cor_mat_{nrows}.npy'),
        # os.path.join(out_dir, 'cor_mats', 'pearson_cor_coop_pcor', 'cor_mat_{nrows}.npy'),
    output:
        os.path.join(out_dir, 'cor_mats', 'cor_mat_{nrows}_differences.tsv')
    run:
        # create an empty matrix to store differences in
        res = np.zeros((len(input), len(input)))

        # iterate over correlation matrices and compute diffs
        for ind1 in range(len(input)):
            m1 = np.load(input[ind1])
            m1[np.isnan(m1)] = 0

            for ind2 in range(len(input)):
                if (ind1 == ind2) or (res[ind1, ind2] != 0):
                    continue

                m2 = np.load(input[ind2])
                m2[np.isnan(m2)] = 0
                diffs = abs(m1 - m2)
                res[ind1, ind2] = diffs.sum()

        # convert to a dataframe and add column and rownames
        ind_names = [os.path.basename(os.path.dirname(x)) for x in input]
        df = pd.DataFrame(res, columns=ind_names, index=ind_names)

        # save result
        df.to_csv(output[0], sep='\t')

rule benchmark_spearman_cor_numpy_corrcoef:
    input:
        os.path.join(out_dir, 'data', 'mat_{nrows}.npy')
    output:
        cor_mat=os.path.join(out_dir, 'cor_mats', 'spearman_cor_numba', 'cor_mat_{nrows}.npy'),
        timings=os.path.join(out_dir, 'timings', 'spearman_cor_numba', 'timings_{nrows}.csv')
    threads: config['benchmark']['num_threads']
    script: 'src/spearman-cor-numba.py'

rule benchmark_pearson_cor_numpy_corrcoef:
    input:
        os.path.join(out_dir, 'data', 'mat_{nrows}.npy')
    output:
        cor_mat=os.path.join(out_dir, 'cor_mats', 'pearson_cor_numpy_cor_coef', 'cor_mat_{nrows}.npy'),
        timings=os.path.join(out_dir, 'timings', 'pearson_cor_numpy_cor_coef', 'timings_{nrows}.csv')
    threads: config['benchmark']['num_threads']
    script: 'src/pearson-cor-numpy-corrcoef.py'

rule benchmark_pearson_cor_numpy_einsum:
    input:
        os.path.join(out_dir, 'data', 'mat_{nrows}.npy')
    output:
        cor_mat=os.path.join(out_dir, 'cor_mats', 'pearson_cor_numpy_einsum', 'cor_mat_{nrows}.npy'),
        timings=os.path.join(out_dir, 'timings', 'pearson_cor_numpy_einsum', 'timings_{nrows}.csv')
    threads: config['benchmark']['num_threads']
    script: 'src/pearson-cor-numpy-einsum.py'

rule benchmark_pearson_cor_coop_pcor:
    input:
        os.path.join(out_dir, 'data', 'mat_{nrows}.npy')
    output:
        cor_mat=os.path.join(out_dir, 'cor_mats', 'pearson_cor_coop_pcor', 'cor_mat_{nrows}.npy'),
        timings=os.path.join(out_dir, 'timings', 'pearson_cor_coop_pcor', 'timings_{nrows}.csv')
    threads: config['benchmark']['num_threads']
    script: 'src/pearson-cor-coop-pcor.R'

rule benchmark_pearson_base_r:
    input:
        os.path.join(out_dir, 'data', 'mat_{nrows}.npy')
    output:
        cor_mat=os.path.join(out_dir, 'cor_mats', 'pearson_cor_base_r', 'cor_mat_{nrows}.npy'),
        timings=os.path.join(out_dir, 'timings', 'pearson_cor_base_r', 'timings_{nrows}.csv')
    threads: config['benchmark']['num_threads']
    script: 'src/pearson-cor-base-r.R'

rule benchmark_pearson_cor_julia:
    input:
        os.path.join(out_dir, 'data', 'mat_{nrows}.npy')
    output:
        cor_mat=os.path.join(out_dir, 'cor_mats', 'pearson_cor_julia', 'cor_mat_{nrows}.npy'),
        timings=os.path.join(out_dir, 'timings', 'pearson_cor_julia', 'timings_{nrows}.csv')
    threads: config['benchmark']['num_threads']
    script: 'src/pearson-cor.jl'

rule benchmark_pearson_cor_cupy:
    input:
        os.path.join(out_dir, 'data', 'mat_{nrows}.npy')
    output:
        cor_mat=os.path.join(out_dir, 'cor_mats', 'pearson_cor_cupy', 'cor_mat_{nrows}.npy'),
        timings=os.path.join(out_dir, 'timings', 'pearson_cor_cupy', 'timings_{nrows}.csv')
    threads: config['benchmark']['num_threads']
    script: 'src/pearson-cor-cupy.py'

rule build_test_data:
    output:
        os.path.join(out_dir, 'data', 'mat_{nrows}.npy')
    run:
        np.random.seed(0)
        np.save(output[0], np.random.randn(int(wildcards['nrows']), num_cols))
        
