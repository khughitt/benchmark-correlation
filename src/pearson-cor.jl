#
# Pearson correlation using Statistics.cor()
#
using Statistics
using BenchmarkTools
using LinearAlgebra
using Printf
using NPZ

mat = npzread(snakemake.input[1])

# fix type for config param
nrows = parse(Int64, snakemake.wildcards["nrows"])

# load input data
res = @benchmark Statistics.cor(mat') samples=snakemake.config["benchmark"]["num_times"] evals=1

# get minimum time in seconds
min_time = minimum(res.times) / 1e9

# save timings
open(snakemake.output["timings"], "w") do io
    entry = @sprintf "Pearson, Julia, Statistics.cor, %d, %f" nrows min_time;
    write(io, entry)
end;

# store result for comparison
cor_mat = Statistics.cor(mat')

# set lower triangular matrix to NaN
for k in 1:size(cor_mat, 1) cor_mat[diagind(cor_mat, k)] .= NaN end
cor_mat = cor_mat'

# set diagonal to NaN 
cor_mat[diagind(cor_mat)] .= NaN

# save correlation matrix
npzwrite(snakemake.output["cor_mat"], cor_mat)
