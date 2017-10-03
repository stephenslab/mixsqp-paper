# mixopt

Algorithms for solving the "mixture distribution" optimization
problem. See [here](https://stephenslab.github.io/mixopt) for an
overview of the problem and results of our experiments.

## License

All source code and software in this repository are made available
under the terms of the
[MIT license](https://opensource.org/licenses/mit-license.html).

## Setup

+ Install [Julia](http://julialang.org).

+ Install the LowRankApprox package in Julia:

   ```julia
   Pkg.add("LowRankApprox")
   ```

+ in R, install the mixopt package from this repository:

  ```R
  library(devtools)
  install_github("stephenslab/mixopt")
  ```

## How to build the webpages

*Add instructions here.*

## How to update the package documentation.

Within a local copy of this repository, run the following R commands:

```R
library(devtools)
document()
```

## Credits

This project was developed by
[Youngseok Kim](https://github.com/youngseok-kim)
at the [University of Chicago](https://www.uchicago.edu),
with contributions from
[Peter Carbonetto](https://pcarbo.github.io),
[Matthew Stephens](http://stephenslab.uchicago.edu) and
[Mihai Anitescu](http://www.mcs.anl.gov/~anitescu).

