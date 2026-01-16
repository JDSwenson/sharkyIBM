
<!-- README.md is generated from README.Rmd. Please edit that file -->

# sharkyIBM

<!-- badges: start -->

<!-- badges: end -->

sharkyIBM includes code to simulate and sample a population of sharks,
with individual “simsharks” growing, breeding and dying based on age-
and/or sex-specific probability distributions governing growth,
movement, fecundity, and mortality. The simulations track first and
second-order relatives, and allow combinations of gene flow and sampling
design that resemble a real population of interest e.g., via simulating
multiple genetically distinct subpopulations that are sampled during an
aggregation.

The primary purposes of this package are to 1) assist with close-kin
mark-recapture (CKMR) sample design, and 2) develop/hone CKMR models.

## Installation

You can install the development version of sharkyIBM from
[GitHub](https://github.com/) with:

``` r
# install.packages("pdevtools")
devtools::install_github("JDSwenson/sharkyIBM")
```

## Usage

Explain a potential use case …

## Example

This is a basic example showing how to prepare input data for two
populations with different life histories and biological attributes,
simulate those populations forward in time, and then draw 100 samples
from each population over the final 5 years of the simulation.

``` r
library(sharkyIBM)
## basic example code
```

You’ll need to render `README.Rmd` regularly, to keep `README.md`
up-to-date. `devtools::build_readme()` is handy for this.

You can also embed plots, for example:

<img src="man/figures/README-pressure-1.png" width="100%" />
