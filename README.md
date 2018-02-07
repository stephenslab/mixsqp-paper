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

2. Install the Plots and LowRankApprox packages in Julia, then
   precompile them:

   ```julia
   Pkg.add("Plots")
   Pkg.add("DataFrames")
   Pkg.add("RCall")
   Pkg.add("LowRankApprox")
   using Plots
   using DataFrames
   using RCall
   using LowRankApprox
   ```

   The Plots package is not needed to run the core mixopt functions,
   but is useful for running the examples.

## Credits

This project was developed by
[Youngseok Kim](https://github.com/youngseok-kim)
at the [University of Chicago](https://www.uchicago.edu),
with contributions from
[Peter Carbonetto](https://pcarbo.github.io),
[Matthew Stephens](http://stephenslab.uchicago.edu) and
[Mihai Anitescu](http://www.mcs.anl.gov/~anitescu).

