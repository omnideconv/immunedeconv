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

### Continuous integration

Continuous integration are automated services that run on every pull request such as

- code consistency checks
- build documentation
- deploy website.

We have continuous integration set-up for all repositories. What the different services
are doing and how you can benefit from them is detailed in the following sections.

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

### Checks

Before a pull request can get merged, the code must pass `R CMD check`.

You can test your packages passes the checks locally by running

```r
devtools::check()
```

in an R console.

Alternatively, continuous integration will run the checks on every pull request.

### Tests

We use [`testthat`](https://testthat.r-lib.org/) for automated testing. Automated tests are small pieces of code that
run the package with small example data to ensure that everything works as expected. If you add new functionality
to one of the packages, please add the corresponding tests.

To learn more about tests, please refer to the [testing chapter from the R packages book](http://r-pkgs.had.co.nz/tests.html),
or take a look at one of our existing test suites, e.g. [immunedeconv](https://github.com/omnideconv/immunedeconv/tree/master/tests).

If you are working with RStudio, you can simply press _Cmd/Ctrl + Shift + T_ to execute tests. It is also possible
to execute tests from the R console:

```r
devtools::test()
```

Finally, our continuous integration will automatically run the tests on all pull requests.

### Conda

For some of our packages (omnideconv, immunedeconv) we cannot distribute the package on the usual repositories
(CRAN or Bioconductor), since several of our dependencies are not available from there. To provide an alternative
way of installation in addition to `install_github`, we decided to provide [bioconda packages][bioconda] for those
repositories.

Conda is a platform-independent package manager that automatically resolves version conflicts and provides
pre-compiled binaries of the packages. Therefore, conda tends to be faster and more reliable than a direct installation
from GitHub.

Each conda package has a YAML "recipe", which is hosted on [bioconda-recipes][] (e.g. [immunedeconv][immunedeconv-bioconda]).
All dependencies must be available either on bioconda or [conda-forge][]. If a dependency is not available from
there, you'll need to add it. To this end, please follow the instructions from [conda-forge][conda-forge-add-package]
for general-purpose packages and [bioconda][bioconda-add-package] for packages specific to the biological sciences.

The conda recipe in [bioconda-recipes][] needs to be updated when a new
release of an omnideconv package is created. For more details, see the [release](#making-a-release) section.

For omnideconv packages that provide a bioconda version, there should be a copy of the `meta.yml` from bioconda-recipes
in the `.conda` directory. This file is used by the continuous integration to check that the conda recipe
can be built and used at any time. This file should be kept in sync with the file from bioconda-recipes (i.e.
when you have to add a dependency to the local `meta.yml` to make the tests pass, you will also need to
add that depependency to the version from bioconda-recipes on the next release).

It is also possible to test building the package locally. This might be useful for debugging:

```bash
cd .conda  # or whatever directory your `meta.yml` is in
conda build . --no-anaconda-upload
```

### Documentation

Documentation is essential for the users to correctly use our package, so please take a lot of care when writing
documentation.

We use [roxygen2][] to build documentation from code comments. It automatically generates manual files (`.Rd`) and the
`NAMESPACE` file that declares which external functions are used by the package. To learn more about roxygen2, please
refer to their [getting started vignette][roxygen-get-started].

To build the documentation, run

```r
devtools::document()
```

in an R console. If you are using Rstudio, you can simply press _Crtl/Cmd + Shift + D_.

Alternatively, you can rely on our continuous integration to build the documentation of every push to a pull request.
If the documentation was out of date, it will automatically add a commit with the updated documentation. If the CI
added a commit to a branch you still have been working on locally, simply use

```bash
git pull --rebase
```

to integrate the changes into yours.

#### Writing documentation

To write good docstrings, it's good to keep the following points in mind

- The roxygen2 docstring should describe what the function is doing (try to think about it from
  the prespective of the user).
- The docstring should describe all parameters and return values.
- Consider adding an example to one of the vignettes

To learn more about how to write good documentation for R packages, we recommend the "Documentation" chapters of the
[R packages book][r-pkgs-doc].

### pkgdown website

We use [pkgdown][] to automatically create a documentation website for each package. The websites are served
by GitHub pages from the `gh-pages` branch of each repository.

The website is built automatically by continuous integration of every pull request. If a pull request gets merged
into the main branch, the website gets updated automatically. For each pull request a "preview link" is automatically
generated and added as a comment. From this link, you can inspect how the website will look like given the changes
from your pull request.

Alternatively, you can build and preview the website locally by running in an R console

```r
pkgdown::build_site()
```

### Making a release

The package maintainers are in charge of creating releases whenever bugfixes or
features were merged into the main branch.

Please adhere to [semantic versioning](https://semver.org) for choosing a
version number, in brief:

> Given a version number MAJOR.MINOR.PATCH, increment the:
>
> - MAJOR version when you make incompatible API changes,
> - MINOR version when you add functionality in a backwards compatible manner, and
> - PATCH version when you make backwards compatible bug fixes.
>
> Additional labels for pre-release and build metadata are available as extensions to the MAJOR.MINOR.PATCH format.

#### Release checklist

- [ ] Update the version number in `DESCRIPTION`
- [ ] Create a release in the GitHub interface. Both the tag and the release title should have the format `vX.X.X`.
- [ ] Write [release notes](#writing-release-notes) as described below.
- [ ] Update the bioconda recipe, if applicable. The BiocondaBot should create a pull-request automatically within a day and
      send a notification. Update the recipe if necessary, then approve and merge the pull-request.
- [ ] Consider advertising the release, e.g. on twitter

#### Writing release notes

The release notes or "changelog" should provide an overview of the changes introduced since the last release.
It is good practice to separate it in sections that correspond to the categories of changes used for semantic versioning, i.e.

- backwards incompatible changes,
- new features, and
- bug fixes.

Each change should link to the corresponding issue or pull-request. Please make sure to give appropriate credits
to contributors.

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
[bioconda]: https://bioconda.github.io/
[bioconda-recipes]: https://github.com/bioconda/bioconda-recipes/
[immunedeconv-bioconda]: https://github.com/bioconda/bioconda-recipes/blob/master/recipes/r-immunedeconv/meta.yaml
[conda-forge]: https://conda-forge.org/
[conda-forge-add-package]: https://conda-forge.org/docs/maintainer/adding_pkgs.html
[bioconda-add-package]: https://bioconda.github.io/contributor/index.html
[roxygen2]: https://roxygen2.r-lib.org/
[roxygen-get-started]: https://roxygen2.r-lib.org/articles/roxygen2.html
[r-pkgs-doc]: https://r-pkgs.org/man.html
[pkgdown]: https://pkgdown.r-lib.org/
