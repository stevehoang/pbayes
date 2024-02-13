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

To demonstrate usage, we will simulate a distribution p-values and then use `pbayes` to fit a uniform-beta distribution to the data:

```r
# Simulate p-values sampled from a mixture distribution -------------------

# set a random seed
set.seed(42)

# number of p-values to simulate
sample_size <- 10000

# set parameters 
l0 <- 0.8 # mixing fraction of the uniform component of the mixture
l1 <- 0.2 # mixing fraction of the non-uniform beta component
r1 <- 0.5 # shape parameter 1 of the non-uniform beta component
s1 <- 3 # shape parameter 2 of the non-uniform beta component

# simulate p-values based on the parameters above
n_uniform <- sample_size * l0
n_nonuniform <- sample_size * l1
comp_uniform <- rbeta(n_uniform, 1, 1) # sample from a uniform distribtion
comp_nonuniform <- rbeta(n_nonuniform, r1, s1) # sample from a non-uniform beta
p_sim <- c(comp_uniform, comp_nonuniform) # simulated p-values
hist(p_sim)
```
![pvalue_hist](https://github.com/stevehoang/pbayes/assets/3991279/6aeefa18-b2a4-4dae-b2a5-b46ad5627972)

This type of p-value distribution is typical in biological "omic" experiments where hypothesis tests are performed across thousands of entities (e.g., genes) simultaneously.

Now we'll use `pbayes` to fit the mixure model and calculate posterior probabilities.
