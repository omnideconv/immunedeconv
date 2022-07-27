<!-- TODO: rewrite for generic omnideconv package and put on the website; link it from `.github/CONTRIBUTING.md` -->

# Contributing to the omnideconv communitiy

We welcome contributions to omnideconv! For contributions of any kind, please adhere to our
[Code of Conduct](https://github.com/omnideconv/.github/blob/main/CODE_OF_CONDUCT.md).

<!-- this is for the case we have a global community forum at some point -->
<!-- ## Answering questions

The easiest way to contribute to omnideconv is by answering questions in our community forum.

-->

## Filing an issue

Bug reports and feature requests are indispensible for improving the quality of omnideconv tools.
To make them as useful as possible:

- Search the repository to see if someone has already reported the same issue.
  This allows contributors to spend less time responding to issues, and more time adding new features!
- Please provide a minimal complete verifiable example for any bug.
  If you're not sure what this means, check out
  [this blog post](http://matthewrocklin.com/blog/work/2018/02/28/minimal-bug-reports)
  by Matthew Rocklin or [this definition](https://stackoverflow.com/help/mcve) from StackOverflow.
- Let us know about your environment. Environment information is available via `sessionInfo()`.

## Contributing code

If you want to propose code changes to omnideconv repositories (e.g. to add a feature or fix a bug),
you'll need to prepare a pull-request (PR) that will be reviewed by core-team members. We are very excited
about code contributions, however, before adding a feature, consider opening an issue for discussion.

This section gives an overview of our coding standards and provides you all necessary information
to get started quickly!

### Getting set-up

We assume that you are already familiar with gith and with making pull requests on GitHub. If not, there are plenty
of tutorials available online, e.g. the ["first contributions" tutorial][first-contributions]

[first-contributions]: (https://github.com/firstcontributions/first-contributions/blob/main/README.md).

### Installing dev dependencies

In addition to the package you are developing and all its runtime dependencies, you may need additional packages
to run tests and build vignettes. These dependencies are usually included in the `Suggests:` section of the
`DESCRIPTION` file.

We recommend installing all dependencies through the [r-lib/pak][pak] package manager, a faster and more reliable
alternative to `install_github()`:

```r
# install the `r-lib/pak` package manager, if you haven't already
install.packages("pak")

# install package including all dev dependencies, e.g. immunedeconv
pak::pkg_install("omnideconv/immunedeconv", dependencies = TRUE)
```

### Code-style

We use [pre-commit][] to enforce consistent code-styles across projects. On every commit, pre-commit checks will either
automatically fix issues with the code, or raise an error message. See [pre-commit checks](#pre-commit-checks) for
a full list of checks included in omnideconv repositories per default. Individual project maintainers may decide
to add or remove checks, though.

To enable pre-commit locally, install the `pre-commit` Python tool:

```bash
pip install pre-commit
```

Then, simply run

```bash
pre-commit install
```

to activate the checks.

Alternatively, you can rely on the [pre-commit.ci][] service enabled on GitHub. If you didn't run `pre-commit`
before pushing changes to GitHub, it will automatically commit fixes to your pull request, or show an error message.

If pre-commit.ci added a commit on a branch you still have been working on locally, simply use

```bash
git pull --rebase
```

to integrate the changes into yours.

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

- The roxygen2 docstring should describe what the function is doing (try to think about it from
  the prespective of the user).
- The docstring should describe all parameters and return values.
- Consider adding an example to one of the [vignettes](https://github.com/icbi-lab/immunedeconv/tree/master/vignettes).

TODO run roxygen locally
TODO run roxygen on CI

#### Previewing documentation

TODO build locally
TODO preview on CI

#### Publishing documentation

The documentation is automatically built by the continuous integration, and
published on the omnideconv website every time a pull-request is merged into the
main branch.

### Checks

bioconductor vs CRAN vs freestyle

### Making a release

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

### Repository set-up

This section gives an overview which services integrate with omnideconv repositories and how they were set-up.
It is meant primarily as a reference for (future) mainainers of omnideconv repositories.

#### Documentation using `pkgdown` hosted on GitHub Pages

#### Coverage tests with CodeCov

#### Pre-commit checks

<!-- links -->

[pak]: https://github.com/r-lib/pak
[pre-commit]: https://pre-commit.com
[pre-commit.ci]: https://pre-commit.ci
