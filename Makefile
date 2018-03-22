help:
	@echo "The following commands are available:"
	@echo "    roxygenize        Build roxygen2 documentation"
	@echo "    check             install dependencies and check if package can be loaded"
	@echo "    install           install package and dependencies"
	@echo "    docs              build docs with pkgdown"
	@echo "    deploy_docs       build docs and commit them to the git repository (-> github pages)"
	@echo "    clean             clean repository (docs, man and all dev files)"

.PHONY: roxygenize
roxygenize:
	Rscript -e "library(methods); library(devtools); document()"

.PHONY: check
check:
	Rscript -e 'install.packages(c("devtools", "roxygen2", "covr"), repos="https://cran.rstudio.com/")'
	Rscript -e 'devtools::install_github("hadley/devtools")'
	Rscript -e 'devtools::install_dev_deps()'
	Rscript -e 'devtools::check()'

.PHONY: install
install: | check roxygenize
	R CMD INSTALL .

.PHONY: docs
docs: roxygenize
	Rscript -e 'pkgdown::build_site()'

.PHONY: deploy_docs
deploy_docs: docs
	git add docs && git commit -m "update docs" && git push

.PHONY: clean
clean:
	rm -rfv docs/*
	rm -rfv builds
	rm -rfv man
	rm -rfv vignettes/*.html
	find . -type f -name "*~" -exec rm '{}' \;
	find . -type f -name ".Rhistory" -exec rm '{}' \;
