# Contributing to immunedeconv

We welcome contributions to immunedeconv! For contributions of any kind, please adhere to our
[Code of Conduct](https://github.com/omnideconv/.github/blob/main/CODE_OF_CONDUCT.md). 

## Answering questions

TODO

The easiest way to contribute to immunedeconv is by answering questions in our
community forum. 

## Filing an issue

Bug reports and feature requests are indispensible for improving immunedeconv. To make them as useful as possible:

- Search the repository to see if someone has already reported the same issue.
  This allows contributors to spend less time responding to issues, and more time adding new features!
- Please provide a minimal complete verifiable example for any bug.
  If you're not sure what this means, check out
  [this blog post](http://matthewrocklin.com/blog/work/2018/02/28/minimal-bug-reports)
  by Matthew Rocklin or [this definition](https://stackoverflow.com/help/mcve) from StackOverflow.
- Let us know about your environment. Environment information is available via `sessionInfo()`.

## Contributing code

We are absolutely enthusiastic about code contributions! If you prepare a pull-request (PR), the 
following guidelines help you to get started quickly and to meet our coding
standards. 

### Getting set-up

We assume, you are familar with the basic fork-and-pull-request workflow of
GitHub. This workflow is not specific to immunedeconv and there are plento of
tutorials available online, e.g. TODO. 

### Tests

Please write tests! We use [`testthat`](https://testthat.r-lib.org/) to ensure the package works correctly.
You can refer to the [existing test suite](https://github.com/icbi-lab/immunedeconv/tree/master/tests)
and the [Testing chapter from the R packages book](http://r-pkgs.had.co.nz/tests.html) when adding new tests.

TODO continuous integration
TODO test locally
TODO integration with rstudio

### Conda

TODO adding a dependency
TODO testing to build the recipe locally
TODO testing to build the recipe on the CI
The conda recipe in [bioconda-recipes](TODO) needs to be updated when a new
release of immunedeconv is created. For more details, see the
[release](#making-a-release) section. 

### Documentation

 * The roxygen2 docstring should describe what the function is doing (try to think about it from 
   the prespective of the user).
 * The docstring should describe all parameters and return values.
 * Consider adding an example to one of the [vignettes](https://github.com/icbi-lab/immunedeconv/tree/master/vignettes). 

#### Previewing documentation

TODO build locally 
TODO preview on CI

#### Publishing documentation

The documentation is automatically built by the continuous integration, and
published on the omnideconv website every time a pull-request is merged into the
main branch. 

## Making a release

The package maintainers are in charge of creating releases whenever bugfixes or
features were merged into the main branch. 

This project adheres to [semantic versioning](https://semver.org) for choosing a
version number. 

To make a release: 

1. Create a release in the GitHub interface. Both the tag and the release title
   should have the format `vX.X.X`. 
2. Provide a summary of the changes since the last release. Each change should
   link to the corresponding issue or pull-request. 
3. Make sure to give appropriate credits to contributors
4. Once the release is created, update the recipe on Bioconda. The Bioconda bot
   should create a pull-request automatically within a day. Update the recipe if
   necessary and approve and merge the pull-request. `
