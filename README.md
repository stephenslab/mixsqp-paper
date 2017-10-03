# mixopt

Algorithms for solving the "mixture distribution" optimization
problem. See [here](https://stephenslab.github.io/mixopt) for an
overview of the problem and results of our experiments.

## License

All source code and software in this repository are made available
under the terms of the
[MIT license](https://opensource.org/licenses/mit-license.html).

## Quick Start

1. Install [Julia](http://julialang.org).

2. Install the LowRankApprox package in Julia:

   ```julia
   Pkg.add("LowRankApprox")
   ```

3. in R, install and load the mixopt package from this repository:

   ```R
   library(devtools)
   install_github("stephenslab/mixopt")
   library(mixopt)
   ```

4. Add step here.

## How to build the webpages

*Add instructions here.*

## How to update the package documentation

Within a local copy of this repository, run the following R commands:

```R
library(devtools)
document()
```

## How to check the package for errors and warnings

```bash
R CMD build mixopt --resave-data
R CMD check --as-cran mixopt_0.1-0.tar.gz
```

## Credits

This project was developed by
[Youngseok Kim](https://github.com/youngseok-kim)
at the [University of Chicago](https://www.uchicago.edu),
with contributions from
[Peter Carbonetto](https://pcarbo.github.io),
[Matthew Stephens](http://stephenslab.uchicago.edu) and
[Mihai Anitescu](http://www.mcs.anl.gov/~anitescu).

