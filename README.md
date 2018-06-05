# mix-SQP: Sequentional Quadratic Programming for fast maximum-likelihood estimation of mixture proportions

This repository contains code resources to accompany our research
paper,

> Youngseok Kim, Peter Carbonetto, Matthew Stephens and Mihai Anitescu
> (2018). *A fast algorithm for maximum likelihood estimation of
> mixture proportions using sequential quadratic programming.*
> (Submitted for review.)

We provide code implementing optimization method for
maximum-likelihood estimation of mixture proportions, including a fast
algorithm based on sequential quadratic programming, which we call
"mix-SQP".

All methods and demonsstrations are implemented in [Julia][julia],
aside from the interior point solver based on [MOSEK][mosek], which is
implemented in R. For an R package implementing some of these methods,
please see the [mixsqp repository][mixsqp]

For demos of the mix-SQP algorithm and comparisons with other
optimization approaches, see [here][github-site].

If you find a bug, or you have a question or feedback on our work,
please post an [issue][issues].

## License

Copyright (c) 2017-2018, Youngseok Kim, Peter Carbonetto, Matthew
Stephens and Mihai Anitescu.

All source code and software in this repository are made available
under the terms of the [MIT license][mit-license]. See the
[LICENSE](LICENSE) file for the full text of the license.

## Citing this work

If you find any of the source code in this repository useful for your
work, please cite our manuscript, Kim *et al* (2018). The full
citation is given above.

## Setup instructions

To work through the [demos][github-site], follow these steps.

**Note:** If you encounter a problem running the demos, or installing
any of the programs, you want to try replicating our computing setup
as closely as possible. At the bottom of each Jupyter notebook and
webpage included [here][github-site], we have recorded information
about the setup used to run the Julia code.

1. Install [Julia](http://julialang.org).

2. Install the [Distributions][distributions-julia] and
   [LowRankApprox][lowrankapprox-julia] Julia packages:

   ```
   Pkg.add("Distributions")
   Pkg.add("LowRankApprox")
   ```

3. Restart Julia, then precompile the packages and load them to make
   sure that they work:

   ```julia
   using Distributions
   using LowRankApprox
   ```

4. Clone or download this git repository on your computer.

5. Work through the [first demo][mixsqp-example]. To do this, start up
   Julia, making sure that the working directory is set to the
   `analysis` directory inside the git repository. To verify, run
   `pwd()` in Julia; it should return a directory with
   `mixsqp-paper/analysis` at the end.

   Alternatively, you may open the
   [mixsqp_example notebook](analysis/mixSQP_example.ipynb) in
   [Jupyter][jupyter] (e.g., Jupyter notebook), and run the code
   inside Jupyter. This requires that you have already installed the
   [Jupyter kernel for Julia][ijulia].

6. *Add additional steps here.*

Download and install [MOSEK](https://www.mosek.com).

Install several Julia packages, including the DataFrames, PyCall
and PyPlot packages:

   ```julia
   Pkg.add("PyCall")
   Pkg.build("PyCall")
   Pkg.add("PyPlot")
   Pkg.add("JuMP")
   Pkg.add("Mosek")
   ```

The `Pkg.build` step is important for installing the `PyPlot`
package.

Precompiling the PyPlot package may take some time since it may
download and install additional Python packages.

## What's included

```
├── analysis
├── code
├── data
├── docs
├── misc
└── output
```

+ The `analysis` directory contains Jupyter notebooks illustrating
  the optimization methods in simulated data sets, and R code to
  generate the plots for the manuscript.

+ The `code` directory contains Julia code implementing the
  optimization algorithms, as well as other code used in the numerical
  experiments.

+ The `data` directory is where all raw and processed data are stored.

+ The `docs` directory contains webpages generated from the Jupyter
  notebooks.

+ Results of the demos and numerical experiments are stored in the
  `output` directory.

## Developer notes

Run these commands from the `analysis` directory to rebuild the
website from the Jupyter notebooks.

```bash
jupyter nbconvert --template=../docs/nbconvert.tpl index.ipynb \
  mixSQP_example.ipynb low_rank_approximation.ipynb --output-dir=../docs
```

## Credits

This project was developed by [Youngseok Kim][youngseok],
[Peter Carbonetto][peter], [Matthew Stephens][matthew] and
[Mihai Anitescu][mihai] at the
[University of Chicago](https://www.uchicago.edu).

[github-site]: https://stephenslab.github.io/mixsqp-paper
[mixsqp]: https://github.com/youngseok-kim/mixsqp
[issues]: https://github.com/stephenslab/mixsqp-paper/issues
[mit-license]: https://opensource.org/licenses/mit-license.html
[youngseok]: https://github.com/youngseok-kim
[peter]: https://pcarbo.github.io
[matthew]: http://stephenslab.uchicago.edu
[mihai]: http://www.mcs.anl.gov/~anitescu
[julia]: http://julialang.org
[mosek]: http://mosek.com
[jupyter]: http://jupyter.org
[ijulia]: https://github.com/JuliaLang/IJulia.jl
[distributions-julia]: https://github.com/JuliaStats/Distributions.jl
[lowrankapprox-julia]: https://github.com/klho/LowRankApprox.jl
[mixsqp-example]: https://stephenslab.github.io/mixsqp-paper/mixSQP_example.html
