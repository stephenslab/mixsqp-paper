# mixopt

Algorithms for solving the "mixture distribution" optimization
problem. See [here](https://stephenslab.github.io/mixopt) for an
overview of the problem and results of our experiments.

## License

Copyright (c) 2017, Youngseok Kim.

All source code and software in this repository are made available
under the terms of the
[MIT license](https://opensource.org/licenses/mit-license.html).

## Quick Start

1. Install [Julia](http://julialang.org).

2. Install Python 3.x and Jupyter bundled together with
   [Anaconda](https://www.anaconda.com). Alternatively, install
   Python 3.x and Jupyter separately; see
   [here](https://jupyter.org/install) and [here](http://python.org)
   for other options to install Jupyter and Python 3.x.

3. Install the [Julia kernel for
   Jupyter](https://github.com/JuliaLang/IJulia.jl).

4. Download and install [MOSEK](https://www.mosek.com).

5. Install the DataFrames, PyCall and PyPlot packages in Julia:

   ```julia
   Pkg.add("DataFrames")
   Pkg.add("CSV")
   Pkg.add("PyCall")
   Pkg.build("PyCall")
   Pkg.add("PyPlot")
   Pkg.add("Mosek")
   ```

   The `Pkg.build` step is important for installing the `PyPlot`
   package.

6. Restart Julia, then precompile the packages and load them to make
   sure that they work:

   ```julia
   using DataFrames
   using CSV
   using PyCall
   using PyPlot
   using Mosek
   ```

   Precompiling the PyPlot package may take some time since it may
   download and install additional Python packages.

7. *Add next step here.*

## What's included

*Give overview of important files here.*

## Credits

This project was developed by
[Youngseok Kim](https://github.com/youngseok-kim)
at the [University of Chicago](https://www.uchicago.edu),
with contributions from
[Peter Carbonetto](https://pcarbo.github.io),
[Matthew Stephens](http://stephenslab.uchicago.edu) and
[Mihai Anitescu](http://www.mcs.anl.gov/~anitescu).

