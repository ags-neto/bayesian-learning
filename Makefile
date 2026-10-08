.PHONY: all test install-deps report clean

R ?= Rscript

# Debian packages used by the analysis and by the renderer. They are system
# packages, not project-local ones; there is no project library to restore.
install-deps:
	apt-get install -y r-base-core jags r-cran-knitr r-cran-rmarkdown r-cran-ggplot2 r-cran-dplyr pandoc
	Rscript -e 'install.packages("runjags", repos="https://cloud.r-project.org")'

all: test

test:
	@$(R) -e 'for (p in c("runjags","dplyr")) if (!requireNamespace(p, quietly=TRUE)) { cat("missing R package:", p, "\n"); quit(status=1) }'
	@$(R) tests/run_tests.R .

report:
	@$(R) -e 'rmarkdown::render("Bayesian Learning and Monte Carlo Simulation - Final Project.Rmd", output_format="pdf_document")'

clean:
	@rm -f tests/run_tests.log
