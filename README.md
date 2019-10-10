# mix-SQP: sequential quadratic programming for fast maximum-likelihood estimation of mixture proportions

This repository contains code resources to accompany our research
paper,

> Youngseok Kim, Peter Carbonetto, Matthew Stephens and Mihai Anitescu
> (2019). *A fast algorithm for maximum likelihood estimation of
> mixture proportions using sequential quadratic programming.*
> [arXiv:1806.01412][arxiv].

We provide code implementing optimization methods for
maximum-likelihood estimation of mixture proportions, in particular a
fast algorithm based on sequential quadratic programming, which we
call "mix-SQP". This method is implemented [Julia][julia]. For
demonstrations of the mix-SQP algorithm and comparisons with other
optimization approaches, see [here][github-site].

For an R implementation of the mix-SQP algorithm, see the
[mixsqp R package][mixsqp].

If you find a bug, or you have a question or feedback on our work,
please post an [issue][issues].

## License

Copyright (c) 2017-2019, Youngseok Kim, Peter Carbonetto, Matthew
Stephens and Mihai Anitescu.

All source code and software in this repository are made available
under the terms of the [MIT license][mit-license]. See the
[LICENSE](LICENSE) file for the full text of the license.

## Citing this work

If you find any of the source code in this repository useful for your
work, please cite our paper. The full citation is given above.

## What's included

```
├── data
├── analysis
├── scripts
├── code
├── output
└── docs
```

+ The `data` directory is where raw and simulated data are stored.

+ The `analysis` directory contains Jupyter notebooks illustrating the
  optimization methods on simulated data sets. The Julia code in these
  notebooks has been tested in Julia 1.1.1.

+ The `scripts` directory contains various Julia and MATLAB scripts
  implementing some of the numerical experiments described in the
  paper. It also includes the R scripts to generate plots for the
  paper. The Julia and MATLAB scripts were tested in Julia 1.1.1 and
  MATLAB 9.5.0 (2018b).

+ The `code` directory contains Julia and MATLAB code implementing the
  optimization algorithms, as well as other code used in the numerical
  experiments.

+ Results of the demos and numerical experiments are stored in the
  `output` directory.

+ The `docs` directory contains webpages generated from the Jupyter
  notebooks using the [nbconvert][nbconvert] program.

## Setup instructions

To run the first two Julia demos [here][github-site], follow these
steps.  *Note:* If you encounter a problem running the demos, or
installing any of the programs, we have recorded information at the
bottom of each Jupyter notebook [here][github-site] about the exact
setup used to run the Julia code.

1. Install [Julia][julia]. *Note that all the code has been tested in
   Julia 1.1.1.*

2. Install the Distributions, LowRankApprox packages in Julia:

   ```julia
   Pkg.add("Distributions")
   Pkg.add("LowRankApprox")
   ```

3. Precompile the packages and load them to make sure that they work:

   ```julia
   using Distributions
   using LowRankApprox
   ```

4. Clone or download this git repository on your computer.

5. Work through the first and second demos [here][github-site]. To do
   this, start up Julia, making sure that the working directory is set
   to the `analysis` directory inside the git repository. To verify,
   run `readdir()` in Julia; it should return a list of the Jupyter
   notebooks such as `mixSQP_example.ipynb`. If it does not, then you
   are in the wrong directory.

   Alternatively, you may open the
   [mixsqp_example notebook](analysis/mixSQP_example.ipynb) in
   [Jupyter][jupyter] (e.g., Jupyter notebook), and run the code
   inside Jupyter. This requires Jupyter and the
   [Julia kernel for Jupyter][ijulia].

## More setup instructions

The simple setup instructions above will allow you to run the first
two examples. If you have successfully run the first two examples, and
you would like to try the other demos, follow these additional setup
instructions. These additional examples compare performance of
different methods, and require additional programs and Julia packages
to run.

1. Install [R][R].

2. Install [MOSEK][mosek].

3. Install the Rmosek package in R (see [here][mosek-docs] for
   installation instructions).

4. Run [this example](scripts/test_rmosek.R) in R to verify that MOSEK
   and Rmosek are correctly installed.
   
5. Install the REBayes package in R:

   ```R
   install.packages("REBayes")
   ```

2. Install the Mosek, JuMP and RCall packages in Julia, and precompile
   the packages to make sure they work:

   ```julia
   Pkg.add("Mosek")
   Pkg.add("JuMP")
   Pkg.add("RCall")
   using Mosek
   using JuMP
   using RCall
   ```

Now try running the other demos listed [here][github-site].

## Developer notes

Run these commands from the `analysis` directory to rebuild the
website from the Jupyter notebooks.

```bash
jupyter nbconvert --template=../docs/nbconvert.tpl index.ipynb \
  mixSQP_example.ipynb low_rank_approximation.ipynb \
  mixSQP_REBayes_EM.ipynb ash_computation_time.ipynb \
  QPsubproblem.ipynb --output-dir=../docs
```

## Credits

This project was developed by [Youngseok Kim][youngseok],
[Peter Carbonetto][peter], [Matthew Stephens][matthew] and
[Mihai Anitescu][mihai] at the [University of Chicago][uchicago].

[github-site]: https://stephenslab.github.io/mixsqp-paper
[mixsqp]: https://github.com/stephenslab/mixsqp
[issues]: https://github.com/stephenslab/mixsqp-paper/issues
[arxiv]: https://arxiv.org/abs/1806.01412
[mit-license]: https://opensource.org/licenses/mit-license.html
[youngseok]: https://github.com/youngseok-kim
[peter]: https://pcarbo.github.io
[matthew]: http://stephenslab.uchicago.edu
[mihai]: http://www.mcs.anl.gov/~anitescu
[uchicago]: https://www.uchicago.edu
[R]: https://www.r-project.org
[julia]: http://julialang.org
[mosek]: http://mosek.com
[mosek-docs]: https://www.mosek.com/documentation
[jupyter]: http://jupyter.org
[nbconvert]: https://nbconvert.readthedocs.io
[ijulia]: https://github.com/JuliaLang/IJulia.jl

