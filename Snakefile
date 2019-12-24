"""
Correlation speed benchmark
KH Dec 16, 2019
"""
import os
import numpy as np
import pandas as pd

# number of columns to use for test data (static)
num_cols = int(config['test_data']['num_cols'])

out_dir = config['output_dir']

rule all:
    input:
        expand(os.path.join(out_dir, 'data', 'mat_{num_rows}.npy'), num_rows=config['test_data']['num_rows']),
        expand(os.path.join(out_dir, 'cor_mats', 'cor_mat_{num_rows}_differences.tsv'), num_rows=config['test_data']['num_rows']),
        os.path.join(out_dir, 'timings', 'all_timings.csv')

rule combine_timings:
    input:
        expand(os.path.join(out_dir, 'timings', 'pearson_cor_numpy_cor_coef', 'timings_{num_rows}.csv'), num_rows=config['test_data']['num_rows']),
        expand(os.path.join(out_dir, 'timings', 'pearson_cor_numpy_einsum', 'timings_{num_rows}.csv'), num_rows=config['test_data']['num_rows']),
        expand(os.path.join(out_dir, 'timings', 'pearson_cor_base_r', 'timings_{num_rows}.csv'), num_rows=config['test_data']['num_rows']),
        expand(os.path.join(out_dir, 'timings', 'pearson_cor_coop_pcor', 'timings_{num_rows}.csv'), num_rows=config['test_data']['num_rows']),
        expand(os.path.join(out_dir, 'timings', 'pearson_cor_cupy', 'timings_{num_rows}.csv'), num_rows=config['test_data']['num_rows'])
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
        os.path.join(out_dir, 'cor_mats', 'pearson_cor_numpy_cor_coef', 'cor_mat_{num_rows}.npy'),
        os.path.join(out_dir, 'cor_mats', 'pearson_cor_numpy_einsum', 'cor_mat_{num_rows}.npy'),
        os.path.join(out_dir, 'cor_mats', 'pearson_cor_coop_pcor', 'cor_mat_{num_rows}.npy'),
        os.path.join(out_dir, 'cor_mats', 'pearson_cor_base_r', 'cor_mat_{num_rows}.npy'),
        os.path.join(out_dir, 'cor_mats', 'pearson_cor_cupy', 'cor_mat_{num_rows}.npy')
    output:
        os.path.join(out_dir, 'cor_mats', 'cor_mat_{num_rows}_differences.tsv')
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

rule benchmark_pearson_cor_numpy_corrcoef:
    input:
        os.path.join(out_dir, 'data', 'mat_{num_rows}.npy')
    output:
        cor_mat=os.path.join(out_dir, 'cor_mats', 'pearson_cor_numpy_cor_coef', 'cor_mat_{num_rows}.npy'),
        timings=os.path.join(out_dir, 'timings', 'pearson_cor_numpy_cor_coef', 'timings_{num_rows}.csv')
    threads: config['benchmark']['num_threads']
    script: 'src/pearson-cor-numpy-corrcoef.py'

rule benchmark_pearson_cor_numpy_einsum:
    input:
        os.path.join(out_dir, 'data', 'mat_{num_rows}.npy')
    output:
        cor_mat=os.path.join(out_dir, 'cor_mats', 'pearson_cor_numpy_einsum', 'cor_mat_{num_rows}.npy'),
        timings=os.path.join(out_dir, 'timings', 'pearson_cor_numpy_einsum', 'timings_{num_rows}.csv')
    threads: config['benchmark']['num_threads']
    script: 'src/pearson-cor-numpy-einsum.py'

rule benchmark_pearson_cor_coop_pcor:
    input:
        os.path.join(out_dir, 'data', 'mat_{num_rows}.npy')
    output:
        cor_mat=os.path.join(out_dir, 'cor_mats', 'pearson_cor_coop_pcor', 'cor_mat_{num_rows}.npy'),
        timings=os.path.join(out_dir, 'timings', 'pearson_cor_coop_pcor', 'timings_{num_rows}.csv')
    threads: config['benchmark']['num_threads']
    script: 'src/pearson-cor-coop-pcor.R'

rule benchmark_pearson_base_r:
    input:
        os.path.join(out_dir, 'data', 'mat_{num_rows}.npy')
    output:
        cor_mat=os.path.join(out_dir, 'cor_mats', 'pearson_cor_base_r', 'cor_mat_{num_rows}.npy'),
        timings=os.path.join(out_dir, 'timings', 'pearson_cor_base_r', 'timings_{num_rows}.csv')
    threads: config['benchmark']['num_threads']
    script: 'src/pearson-cor-base-r.R'

rule benchmark_pearson_cor_cupy:
    input:
        os.path.join(out_dir, 'data', 'mat_{num_rows}.npy')
    output:
        cor_mat=os.path.join(out_dir, 'cor_mats', 'pearson_cor_cupy', 'cor_mat_{num_rows}.npy'),
        timings=os.path.join(out_dir, 'timings', 'pearson_cor_cupy', 'timings_{num_rows}.csv')
    threads: config['benchmark']['num_threads']
    script: 'src/pearson-cor-cupy.py'

rule build_test_data:
    output:
        os.path.join(out_dir, 'data', 'mat_{num_rows}.npy')
    run:
        np.save(output[0], np.random.randn(int(wildcards['num_rows']), num_cols))
        
