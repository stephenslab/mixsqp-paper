# mixopt

Algorithms for solving the "mixture distribution" optimization
problem. See [here](https://stephenslab.github.io/mixopt) for an
overview of the problem and results of our experiments.

## License

All source code and software in this repository are made available
under the terms of the
[MIT license](https://opensource.org/licenses/mit-license.html).

## Quick Start

1. Install [Julia](http://julialang.org) version 0.6. 

2. Install the Plots and LowRankApprox packages in Julia:

   ```julia
   Pkg.add("Plots")
   Pkg.add("DataFrames")
   Pkg.add("RCall")
   Pkg.add("LowRankApprox")
   ```

   The Plots package is not needed to run the core mixopt functions,
   but is useful for running the examples.

3. In R, install and load the mixopt package from a local copy of the
   repository downloaded (or "cloned") onto your computer.

   ```R
   library(devtools)
   install_local("mixopt")
   library(mixopt)
   ```

   This command should automatically retrieve and install the
   appropriate version of the `rjulia` package from Github. (Note that
   cannot use the `install_github` command here while the repository
   is private.)

4.

   ```
   library(rjulia)
   julia_init()
   ```

4. Run this small example demonstrating the SQP algorithm for fitting
   a mixture model to a data set with 5,000 samples.

   ```R
   example("mixsqp")
   ```

5. Check out the package documentation for more details on the methods.

   ```R
   help(package = mixopt)
   ```

## How to build the webpages

*Add instructions here.*

## How to update the package documentation

Within a local copy of this repository, run the following R commands:

```R
library(devtools)
document()
```

## How to test the package build

```bash
R CMD build mixopt --resave-data
R CMD check --as-cran mixopt_0.1-5.tar.gz
```

## Credits

This project was developed by
[Youngseok Kim](https://github.com/youngseok-kim)
at the [University of Chicago](https://www.uchicago.edu),
with contributions from
[Peter Carbonetto](https://pcarbo.github.io),
[Matthew Stephens](http://stephenslab.uchicago.edu) and
[Mihai Anitescu](http://www.mcs.anl.gov/~anitescu).

