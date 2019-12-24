#
# Pearson correlation using coop::pcor
#
library(coop)
library(RcppCNPy)
library(microbenchmark)

# load input data
mat <- npyLoad(snakemake@input[[1]])

# benchmark performance
times <- microbenchmark(
  coop::pcor(t(mat)),
  times = snakemake@config$benchmark$num_times
)

# get minimum time in seconds
min_time <- min(times$time) / 1e9

# store result for comparison
cor_mat <- coop::pcor(t(mat))

cor_mat[lower.tri(cor_mat)] <- NA
diag(cor_mat) <- NA

# save timings
entry <- c('Pearson', 'R', 'coop::pcor', snakemake@wildcards$num_rows, min_time)
entry <- paste0(entry, collapse = ', ')

fp <- file(snakemake@output[['timings']], 'w')
write(entry, file = fp)
close(fp)

# save correlation matrix
npySave(snakemake@output[['cor_mat']], cor_mat)
