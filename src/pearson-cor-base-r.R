#
# Pearson correlation using stats::cor
#
library(RcppCNPy)
library(microbenchmark)

# load input data
mat <- npyLoad(snakemake@input[[1]])

# benchmark performance
times <- microbenchmark(
  cor(t(mat)),
  times = 5 #snakemake@config$benchmark$num_times
)

# get minimum time in seconds
min_time <- min(times$time) / 1e9

# store result for comparison
cor_mat <- cor(t(mat))

cor_mat[lower.tri(cor_mat)] <- NA
diag(cor_mat) <- NA

# save timings
entry <- c('Pearson', 'R', 'stats::cor', snakemake@wildcards$num_rows, min_time)
entry <- paste0(entry, collapse = ', ')

fp <- file(snakemake@output[['timings']], 'w')
write(entry, file = fp)
close(fp)

# save correlation matrix
npySave(snakemake@output[['cor_mat']], cor_mat)
