# Contributing

We welcome contributions to immunedeconv!

## Filing an issue

Bug reports and feature requests are indispensible for improving immunedeconv. To make them as useful as possible: 

* Search the repository to see if someone has already reported the same issue.
  This allows contributors to spend less time responding to issues, and more time adding new features!
* Please provide a minimal complete verifiable example for any bug.
  If you're not sure what this means, check out
  [this blog post](http://matthewrocklin.com/blog/work/2018/02/28/minimal-bug-reports)
  by Matthew Rocklin or [this definition](https://stackoverflow.com/help/mcve) from StackOverflow.
* Let us know about your environment. Environment information is available via `sessionInfo()`.

## Contributing code

We are absolutely enthusiastic about code contributions! If you prepare a PR, we'd like you to follow 
the following guidelines: 

### Tests
Please write tests! We use [`testthat`](https://testthat.r-lib.org/) to ensure the package works correctly. 
You can refer to the [existing test suite](https://github.com/icbi-lab/immunedeconv/tree/master/tests)
and the [Testing chapter from the R packages book](http://r-pkgs.had.co.nz/tests.html) when adding new tests. 

### Documentation
When adding a new function, make sure you write appropriate documentation. We make use 
of [roxygen2](https://cran.r-project.org/web/packages/roxygen2/vignettes/roxygen2.html) to 
document functions and [pkgdown](https://pkgdown.r-lib.org/) to render the documentation website. 

 * The roxygen2 docstring should describe what the function is doing (try to think about it from 
   the prespective of the user).
 * The docstring should describe all parameters and return values.
 * Consider adding an example to one of the [vignettes](https://github.com/icbi-lab/immunedeconv/tree/master/vignettes). 
