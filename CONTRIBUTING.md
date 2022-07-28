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

[first-contributions]: https://github.com/firstcontributions/first-contributions/blob/main/README.md

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
- automated tests
- building documentation

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

You can test if your packages passes the checks locally by running

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
for general-purpose packages, and [bioconda][bioconda-add-package] for packages specific to the biological sciences.

The conda recipe in [bioconda-recipes][] needs to be updated when a new
release of an omnideconv package is created. For more details, see the [making a release](#making-a-release) section.

For omnideconv packages that provide a bioconda package, there should be a copy of the `meta.yml` from bioconda-recipes
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

Documentation is essential for the users to correctly use our packages, so please take a lot of care writing it.

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

#### omnideconv bot

We have an `omnideconv-bot` github account that is used for single-sign on (SSO) for external services and
to generate github tokens for automated commits. Unlike the `GITHUB_TOKEN` used by github actions per default,
access tokens tied to a github account can be used to trigger github actions, which is required for the e.g. the
`roxygenize` action that makes commits to a pull request. Omnideconv repos have access to the bot token
as `secrets.BOT_GH_TOKEN`. The password of the omnideconv-bot account is guarded by the steering committee.

#### Documentation using `pkgdown` hosted on GitHub Pages

The documentation website is built on every pull request by the `pkgdown.yml` github action.
If the check passes, it deploys the website to a temporary website hosted by [netlify][] that allows to inspect
how the website will look like after merging the PR.

If the PR gets merged into the main branch, the website will be updated (i.e. the changes are pushed to the `gh-pages`)
branch of the repository.

Setting up netlify for a new repository:

- login to netlify using "Sign in with GitHub" when logged in with the omnideconv-bot account.
- Create a new site
- in "Site Settings" copy the "Site ID" token and store it as a repository secret named `NETLIFY_SITE_ID`.

An access token to authorize against netlify is available to all repositories as `NETLIFY_AUTH_TOKEN`.

The action uses [actions-netlify][] to deploy to the temporary site, see e.g. [immunedeconv][immunedeconv-netlify-deploy]

#### Coverage tests with CodeCov

Coverage tells what fraction of the code is “covered” by unit tests, thereby encouraging contributors to [write tests](#tests).
To enable coverage checks, head over to [codecov][] and sign in with your GitHub account. You’ll find more information
in “getting started” section of the [codecov docs][].

In brief, you need to

1. Generate a Codecov Token by clicking _setup repo_ in the codecov dashboard.
2. Go to the _Settings_ of your repository on GitHub.
3. Go to _Security > Secrets > Actions_
4. Create new repository secret with name `CODECOV_TOKEN` and paste the token generated by codecov

#### Pre-commit checks

[Pre-commit][] checks are fast programs that check code for errors, inconsistencies and code styles, before the code is committed.
The `.pre-commit-config.yaml` at the root of each repository is the standardized place to describe which
checks to run, and how.

Here's an overview of pre-commit checks included in omnideconv repositories by default. Individual maintainers
may decide to add or remove checks as they prefer.

- [R-specific checks][]

  - **style files**: Autoformat R files with [styler][]
  - **readme-rmd-rendered**: Make sure README.Rmd hasn't been edited more recently than README.md
  - **parseable-r**: Make sure there are not syntax errors in R files
  - **no-browser-statement**: Don't allow to commit code with a `browser()` statement in it
  - **no-debug-statement**: Don't allow to commit code with a `debug()` statement in it.
  - **deps-in-desc**: Ensure that all packages used with `pkgname::fun()` are declared in the `DESCRIPTION` file.

- [pre-commit hooks][]: Genric pre-commit hooks
  - end-of-file-fixer: Makes sure files end in a newline and only a newline.
  - forbid-to-commit: Forbid to commit `.Rhistory`, `.RData`, `.Rds`, or `.rds` files, except for the `data` directory.
- [prettier][]: General purpose autoformatter for many file types (not R)

[r-specific checks]: https://lorenzwalthert.github.io/precommit/articles/available-hooks.html
[pre-commit hooks]: https://github.com/pre-commit/pre-commit-hooks
[prettier]: https://prettier.io/

#### Pre-commit CI

We recommend setting up [pre-commit.ci][] to enforce consistency checks on every commit and pull-request.

Pre-commit.ci was installed into the omnideconv organization. To activate it for a certain repository, head
to the organization settings on GitHub. Choose _Integrations > GitHubApps > pre-commit ci > Configure_.
There, you can choose in the section _Repository access_ for which repos to activate the service.

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
[netlify]: https://www.netlify.com/
[actions-netlify]: https://github.com/nwtgck/actions-netlify
[immunedeconv-netlify-deploy]: https://github.com/omnideconv/immunedeconv/blob/e5fdf46e5c01a8ac7db011051c58ca9d1d219db8/.github/workflows/pkgdown.yml#L51-L61
[codecov]: https://app.codecov.io/login/gh
[codecov docs]: https://docs.codecov.com/docs
