Linear Regresion with PDE regularization for data distributed on 2D Manifolds
================================================================================

The output is an extension of the `R` Package [fdaPDE](https://cran.r-project.org/web/packages/fdaPDE/index.html)

The source code is written in `C++` and linked to `R` throught the API `RcppEigen` and `.Call`

Subfolder structure
--------------------------
- `src` contains all `C++` code and a special file named `Makevars` necessary to build and install the `R` package
- `R` contains the `R` functions that wrap the `C++` calls
- `data` contains all .rda and .RData files useful for testings
- `tests` contains basic `R` script to run tests

To install the package, please make sure that you have the package [devtools](https://cran.r-project.org/web/packages/devtools/index.html) already intalled. If using a Linux machine, it is also advisable to install [rgl](https://cran.r-project.org/web/packages/rgl/index.html) before `fdaPDE`

From the root folder then type

    R -e "library(devtools); install()" --silent
    R -e "library(devtools); document()" --silent
