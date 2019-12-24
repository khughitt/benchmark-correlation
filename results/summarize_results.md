R/Python Correlation Performance Benchmark Results Summary
================
2019-12-24

``` r
library(tidyverse)
```

# Overview

``` r
res <- read_csv('/data/benchmark/correlation/timings/all_timings.csv')
```

    ## Parsed with column specification:
    ## cols(
    ##   Method = col_character(),
    ##   Language = col_character(),
    ##   Implementation = col_character(),
    ##   `Num Rows` = col_double(),
    ##   `Time (Secs)` = col_double()
    ## )

``` r
knitr::kable(res, digits = 2)
```

| Method  | Language | Implementation | Num Rows | Time (Secs) |
| :------ | :------- | :------------- | -------: | ----------: |
| Pearson | Python   | numpy.corrcoef |     1000 |        0.05 |
| Pearson | Python   | numpy.corrcoef |     5000 |        0.29 |
| Pearson | Python   | numpy.corrcoef |    10000 |        1.02 |
| Pearson | Python   | numpy.corrcoef |    15000 |        2.52 |
| Pearson | Python   | numpy.corrcoef |    20000 |        4.42 |
| Pearson | Python   | numpy.corrcoef |    25000 |        6.49 |
| Pearson | Python   | numpy.corrcoef |    37500 |       14.04 |
| Pearson | Python   | numpy.corrcoef |    50000 |       25.36 |
| Pearson | Python   | numpy.einsum   |     1000 |        0.03 |
| Pearson | Python   | numpy.einsum   |     5000 |        0.46 |
| Pearson | Python   | numpy.einsum   |    10000 |        1.63 |
| Pearson | Python   | numpy.einsum   |    15000 |        3.65 |
| Pearson | Python   | numpy.einsum   |    20000 |        6.67 |
| Pearson | Python   | numpy.einsum   |    25000 |        9.14 |
| Pearson | Python   | numpy.einsum   |    37500 |       21.58 |
| Pearson | Python   | numpy.einsum   |    50000 |       42.10 |
| Pearson | R        | stats::cor     |     1000 |        0.66 |
| Pearson | R        | stats::cor     |     5000 |       17.16 |
| Pearson | R        | stats::cor     |    10000 |       71.99 |
| Pearson | R        | stats::cor     |    15000 |      167.05 |
| Pearson | R        | coop::pcor     |     1000 |        0.02 |
| Pearson | R        | coop::pcor     |     5000 |        0.28 |
| Pearson | R        | coop::pcor     |    10000 |        0.87 |
| Pearson | R        | coop::pcor     |    15000 |        1.76 |
| Pearson | R        | coop::pcor     |    20000 |        3.11 |
| Pearson | R        | coop::pcor     |    25000 |        4.41 |
| Pearson | R        | coop::pcor     |    37500 |       11.21 |
| Pearson | Python   | cupy.corrcoef  |     1000 |        0.02 |
| Pearson | Python   | cupy.corrcoef  |     5000 |        0.30 |
| Pearson | Python   | cupy.corrcoef  |    10000 |        1.17 |
| Pearson | Python   | cupy.corrcoef  |    15000 |        2.61 |

# Time (seconds)

``` r
ggplot(res, aes(x = Implementation, y = `Time (Secs)`, fill = Implementation)) +
  geom_bar(stat = 'identity') + 
  facet_wrap(~`Num Rows`, scales = 'free_y') +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90)) +
  ggtitle("Correlation performance comparison") +
  ylab("Time (seconds)")
```

![](summarize_results_files/figure-gfm/correlation_benchmark_barplot-1.png)<!-- -->

# Time (seconds, excluding stats::cor)

``` r
res_subset <- res %>% 
  filter(Implementation != 'stats::cor')

ggplot(res_subset, aes(x = Implementation, y = `Time (Secs)`, fill = Implementation)) +
  geom_bar(stat = 'identity') + 
  facet_wrap(~`Num Rows`, scales = 'free_y') +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90)) +
  ggtitle("Correlation performance comparison (Excluding stats::cor)") +
  ylab("Time (seconds)")
```

![](summarize_results_files/figure-gfm/correlation_benchmark_barplot-3-1.png)<!-- -->

# Time (log-seconds)

``` r
ggplot(res, aes(x = Implementation, y = `Time (Secs)`, fill = Implementation)) +
  geom_bar(stat = 'identity') + 
  scale_y_continuous(trans='log1p') +
  facet_wrap(~`Num Rows`) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90)) +
  ggtitle("Correlation performance comparison (Log-scaled)") +
  ylab("Time (log-seconds)")
```

![](summarize_results_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->
