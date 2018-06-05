# mix-SQP: sequentional quadratic programming for maximum-likelihood estimation of mixture proportions

This repository contains code resources to accompany our research
paper,

> Youngseok Kim, Peter Carbonetto, Matthew Stephens and Mihai Anitescu
> (2018). *A fast algorithm for maximum likelihood estimation of
> mixture proportions using sequential quadratic programming.*

We provide code implementating optimization algorithms for
maximum-likelihood estimation of mixture proportions, including a fast
approach based on sequential quadratic programming, which we call
**mix-SQP**.

All methods and examples are implemented in Julia, aside from the
interior point solver based on [MOSEK][mosek], which is implemented in
R. For an R package implementing some of these methods, see the
[mixsqp repository][mixsqp]


and data to reproduce some of the
empirical results presented in
[Xing & Stephens (2016)][smashr-preprint]. For additional source code,
see the [dscr-smash][dscr-smash]
repository. The new methods can be found in the
[smashr package][smashr].

If you find a bug, please post an [issue][issues].


Algorithms for solving the "mixture distribution" optimization
problem. See [here](https://stephenslab.github.io/mixopt) for an
overview of the problem and results of our experiments.

## License

Copyright (c) 2017-2018, Youngseok Kim, Peter Carbonetto, Matthew
Stephens and Mihai Anitescu.

All source code and software in this repository are made available
under the terms of the
[MIT license](https://opensource.org/licenses/mit-license.html). See
the [LICENSE](LICENSE) file for the full text of the license.

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

5. Install several Julia packages, including the DataFrames, PyCall
   and PyPlot packages:

   ```julia
   Pkg.add("Distributions")
   Pkg.add("PyCall")
   Pkg.build("PyCall")
   Pkg.add("PyPlot")
   Pkg.add("JuMP")
   Pkg.add("Mosek")
   ```

   The `Pkg.build` step is important for installing the `PyPlot`
   package.

6. Restart Julia, then precompile the packages and load them to make
   sure that they work:

   ```julia
   using Distributions
   using PyCall
   using PyPlot
   using JuMP
   using Mosek
   ```

   Precompiling the PyPlot package may take some time since it may
   download and install additional Python packages.

7. *Add next step here.*

## What's included

*Give overview of important files here.*

## Developer notes

Run these commands to re-build the website from the Jupyter
notebooks.

```bash
cd analysis
jupyter nbconvert --template=../docs/nbconvert.tpl \
  index.ipynb mixSQP_example.ipynb low_rank_approximation.ipynb \
  --output-dir=../docs
```

## Credits

This project was developed by [Youngseok Kim][youngseok],
[Peter Carbonetto][peter], [Matthew Stephens][matthew] and
[Mihai Anitescu][mihai] at the
[University of Chicago](https://www.uchicago.edu).

[youngseok]: https://github.com/youngseok-kim
[peter]: https://pcarbo.github.io
[matthew]: http://stephenslab.uchicago.edu
[mihai]: http://www.mcs.anl.gov/~anitescu
[mosek]: http://mosek.com
