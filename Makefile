check:
	Rscript -e 'install.packages(c("devtools", "roxygen2", "covr"), repos="https://cran.rstudio.com/")'
	Rscript -e 'devtools::install_github("hadley/devtools")'
	Rscript -e 'devtools::install_dev_deps()'
	Rscript -e 'devtools::check()'

.PHONY: docs
docs:
	Rscript -e 'pkgdown::build_site()'

.PHONY: deploy_docs
deploy_docs: docs
	git add docs && git commit -m "update docs" && git push


