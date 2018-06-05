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

## Setup instructions (simple)

To try the first and second demos [here][github-site], follow these
steps.

**Note:** If you encounter a problem running the demos, or installing
any of the programs, you want to try replicating our computing setup
as closely as possible. At the bottom of each Jupyter notebook and
webpage included [here][github-site], we have recorded information
about the setup used to run the Julia code.

1. Install [Julia](http://julialang.org).

2. Install the [Distributions][distributions-julia] and
   [LowRankApprox][lowrankapprox-julia] Julia packages:

   ```julia
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

5. Work through the first and second demos [here][github-site]. To do
   this, start up Julia, making sure that the working directory is set
   to the `analysis` directory inside the git repository. To verify,
   run `pwd()` in Julia; it should return a directory with
   `mixsqp-paper/analysis` at the end.

   Alternatively, you may open the
   [mixsqp_example notebook](analysis/mixSQP_example.ipynb) in
   [Jupyter][jupyter] (e.g., Jupyter notebook), and run the code
   inside Jupyter. This requires that you have already installed the
   [Jupyter kernel for Julia][ijulia].

## Setup instructions (more complicated)

The simple setup instructions above will allow you to run the first
two examples. If you have successfully run the first two examples, and
you would work through more examples, follow these additional setup
instructions. The other examples compare performance of different
methods, and require additional programs and Julia packages to
run.

1. Install [R][R]. Note that version 3.4.0 or later is required.

2. Install [MOSEK][mosek].

3. Install the [Rmosek][rmosek] package in R (see [here][mosek-docs] for
   installation instructions).

4. Run [this example](code/test.rmosek.R) in R to verify that MOSEK
   and Rmosek are correctly installed.

5. Install the R package [REBayes][rebayes] from CRAN:

   ```R
   install.packages("REBayes")
   ```

2. Install the [Mosek][mosek-julia], [JuMP][jump-julia] and
   [RCall][rcall-julia] Julia packages:

   ```julia
   Pkg.add("Mosek")
   Pkg.add("JuMP")
   Pkg.add("RCall")
   ```

(*Old instructions:* The `Pkg.build` step is important for installing
the `PyPlot` package. Precompiling the PyPlot package may take some
time since it may download and install additional Python packages.)

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
  mixSQP_example.ipynb low_rank_approximation.ipynb \
  mixSQP_REBayes_EM.ipynb --output-dir=../docs
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
[R]: https://www.r-project.org
[julia]: http://julialang.org
[mosek]: http://mosek.com
[mosek-docs]: https://www.mosek.com/documentation
[jupyter]: http://jupyter.org
[ijulia]: https://github.com/JuliaLang/IJulia.jl
[rmosek]: https://CRAN.R-project.org/package=Rmosek
[rebayes]: https://cran.r-project.org/package=REBayes
[distributions-julia]: https://github.com/JuliaStats/Distributions.jl
[lowrankapprox-julia]: https://github.com/klho/LowRankApprox.jl
[jump-julia]: https://github.com/JuliaOpt/JuMP.jl
[mosek-julia]: https://github.com/JuliaOpt/Mosek.jl
[rcall-julia]: https://github.com/JuliaInterop/RCall.jl
