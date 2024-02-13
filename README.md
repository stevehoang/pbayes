# pbayes: Bayesian Probability Mapping for P-Values

pbayes is an R package designed to map p-values to posterior probabilities based on a uniform-beta mixture model. This approach allows users to interpret p-values from hypothesis tests within a Bayesian framework, offering a nuanced understanding of statistical results.

## Installation

You can install pbayes directly from GitHub using the `remotes` package:

```r
# Install devtools if you haven't already
if (!requireNamespace("remotes", quietly = TRUE))
  install.packages("remotes")

# Install pbayes from GitHub
devtools::install_github("stevehoang/pbayes")
```

## Usage

This packages exposes several functions, but most users will only need to interface with `pbayes`. This function takes as input many independent p-values (at least hundreds), and outputs 1) uniform-beta mixture distribution fitted to the data, and 2) a posterior probability for each p-value that it corresponds to a non-uniform component of the mixture. The latter value is interpreted as the posterior probability of the alternative hypothesis being true.

