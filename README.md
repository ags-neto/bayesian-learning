# bayesian-learning

> Bayesian analysis of the Forest Fires dataset (Montesinho, Portugal) fitted in JAGS from R Markdown, and the report it produces.

## What it is

A university final project in Bayesian learning and Monte Carlo simulation. A single R
Markdown document, `Bayesian Learning and Monte Carlo Simulation - Final Project.Rmd`,
loads `forestfire.txt` and fits three Bayesian normal linear models of the burnt forest
area with JAGS through the `runjags` package, two chains of 5 000 samples after 2 000
burn-in iterations each:

- model 1 (all variables): area regressed on FFMC, DMC, DC, ISI, temperature, relative
  humidity, wind, rain and month, with weakly informative normal priors on the
  coefficients and a Gamma(0.001, 0.001) prior on the precision. The prior standard
  deviations are the variances of the regressors.
- model 2 (weather only): area regressed on temperature, relative humidity, wind and rain.
- model 3: written as `y[i] ~ dnorm(beta0 + beta5, invsigma2)` — that is, a model with two
  intercepts and no temperature term; the `temp` variable passed in the data is not used
  (JAGS reports `Unused variable "temp" in data`).

The document then compares the three models by DIC and closes with the posterior plots.
The finished report is committed as
`Bayesian Learning and Monte Carlo Simulation - Final Project.pdf` (18 pages). Rendering
the `.Rmd` as it stands today on R 4.5.0 and TeX Live 2024 produces 6 pages without the
appendix figures, so the committed PDF does not correspond to the current text: its
wording ("modelled the burned area") differs from the `.Rmd` ("modeled the burnt area"),
and its figure pages come from the posterior diagnostic plots that the last chunk of the
`.Rmd` no longer draws. The data tables and the DIC values (5768 / 5763 / 5765) are
identical in both. Which of the two is the earlier one is not settled here.

`forestfire.txt` is the UCI *Forest Fires* dataset: 517 observations of 13 variables
recorded in Montesinho Natural Park, Portugal, published by P. Cortez and A. Morais
(2007), "A data mining approach to predict forest fires using meteorological data", and
distributed by the UCI Machine Learning Repository under CC BY 4.0
(<https://doi.org/10.24432/C5D88D>). The file in this repository is a space-separated copy
of `forestfires.csv` with an added row-number column; the attribute names and ranges match
the published ones exactly.

The response is hard to model: 247 of the 517 observations have `area` equal to zero, the
maximum is 1090.84 ha and the mean is 12.85 ha (sd 63.66). Fitting the models again
reproduces the published tables but also shows how little the regressors explain: sigma
settles at about 63.5 in all three models, against a raw sd of the response of 63.66, and
the intervals are wide. The tests below reproduce the statistics and do not change them;
they record, without correcting, two things
about how the numbers are read: the report calls temperature the only significant predictor
and cites the interval (0.199, 0.916) for it, while the table printed above that text has
an interval that contains zero, and the DIC spread it uses to pick a favourite model is a
few units wide. The numbers themselves are untouched.

## Requirements

- R 4.5 (tested on Debian 13, arm64) with `runjags` and `dplyr`.
- JAGS 4.3.2.
- `pandoc` 3.1 and a LaTeX installation (`texlive-latex-base`, `texlive-latex-recommended`,
  `texlive-latex-extra`, `texlive-fonts-recommended`, `lmodern`) to render the report.
- GNU make, for the `test` target.

## Install / Build

```sh
sudo apt-get install -y r-base-core jags r-cran-knitr r-cran-rmarkdown r-cran-ggplot2 \
    r-cran-dplyr pandoc texlive-latex-base texlive-latex-recommended texlive-latex-extra \
    lmodern
sudo Rscript -e 'install.packages("runjags", repos="https://cloud.r-project.org")'
```

`make install-deps` runs the same commands.

## Usage

Run the analysis headless, or render the report to PDF:

```sh
Rscript tests/run_tests.R .          # data checks, three models, DIC, negative check
make test                            # the same, with a package pre-check

Rscript -e 'rmarkdown::render("Bayesian Learning and Monte Carlo Simulation - Final Project.Rmd", output_format = "pdf_document")'
```

To replay one model directly, the posterior means of model 2 (about 18 s on a Raspberry Pi
5) are:

```sh
Rscript -e 'library(runjags); d <- read.table("forestfire.txt");
  m <- "model { for (i in 1:N){ y[i] ~ dnorm(beta0 + beta5*temp[i], invsigma2) }
       beta0 ~ dnorm(0, 0.01); beta5 ~ dnorm(0, 1/var(d$temp)); invsigma2 ~ dgamma(0.001, 0.001) }"
  print(run.jags(m, n.chains = 2, data = list(y = d$area, temp = d$temp, N = nrow(d)),
       monitor = c("beta0", "beta5"), burnin = 5000, sample = 5000), digits = 3)'
```

## Tests

`make test` runs `tests/run_tests.R`, which loads the data, fits the three models and
checks 20 facts against values known before the run. It exits non-zero when any check
fails. Two checks are negative: with the `temp` column removed the analysis must abort
instead of returning numbers, and the suite itself must exit non-zero when it is pointed
at that corrupted copy. Real output on a Raspberry Pi 5 (R 4.5.0, JAGS 4.3.2):

```
== bayesian-learning | run_tests.R ==
R R version 4.5.0 (2025-04-11) | 2026-10-08 15:31:45 | root=/home/dsh/bl-work

-- artefacts --
Rmd lines: 346 | PDF bytes: 570625

-- data (forestfire.txt) --
dim: 517 x 13 | columns: X,Y,month,day,FFMC,DMC,DC,ISI,temp,RH,wind,rain,area

-- models (JAGS) -- burnin=5000 sample=5000 chains=2
Warning message:
In rjags::jags.model(model, data = dataenv, inits = inits, n.chains = length(runjags.object$end.state),  :
  Unused variable "temp" in data
elapsed: 18.7 s
posterior means  model1(all):  beta0=-0.000 beta1=-0.044 beta2=0.112 beta3=-0.019 beta4=-0.506 beta5=0.899 beta6=-0.173 beta7=0.847 beta8=-0.004 beta9=0.956 sigma=63.568
posterior means  model2(wx):   beta0=-0.613 beta5=0.895 beta6=-0.138 beta7=0.667 beta8=-0.009 sigma=63.481
posterior means  model3(temp): beta0=9.213 beta5=2.964 sigma=63.763
DIC: model1_all=5759.9 model2_weather=5758.9 model3_temp=5762.9
DIC ordering obtained: model2_weather < model1_all < model3_temp

-- negative check: corrupted data --
Warning message:
Unknown or uninitialised column: `temp`. 
aborted as expected: 'x' is NULL
Warning message:
In file(file, "rt") :
  cannot open file '/tmp/RtmpbzJzse/does-not-exist.txt': No such file or directory

== verification table ==
| exists: forestfire.txt                                                    | file                  | present                                                         | yes      |
| exists: Bayesian Learning and Monte Carlo Simulation - Final Project.Rmd  | file                  | present                                                         | yes      |
| exists: Bayesian Learning and Monte Carlo Simulation - Final Project.pdf  | file                  | present                                                         | yes      |
| Rmd line count                                                            | > 300                 | 346                                                             | yes      |
| data rows                                                                 | 517                   | 517                                                             | yes      |
| data columns                                                              | 13                    | 13                                                              | yes      |
| required columns present                                                  | 13 named columns      | X,Y,month,day,FFMC,DMC,DC,ISI,temp,RH,wind,rain,area            | yes      |
| no missing values                                                         | 0 NA                  | 0                                                               | yes      |
| area is right-skewed and zero-inflated                                    | max>1000 & >200 zeros | max=1090.84 zeros=247                                           | yes      |
| posterior sigma (model 3) within 5 of the raw sd of area (63.66)          | |mean-63.66|<=5       | 63.76                                                           | yes      |
| posterior interval for beta5 (model 2) covers the report's mean (0.918)   | 0.918 inside interval | 0.196 .. 1.629                                                  | yes      |
| posterior mean of beta5 (model 3) is near the report's 3.2                | 3.2 +/- 1.0           | 2.964                                                           | yes      |
| model 1 interval for beta1 (FFMC) covers zero, as the report states       | lower<0<upper         | -0.518 .. 0.499                                                 | yes      |
| all three models converge (psrf <= 1.05)                                  | <=1.05                | 1.0193                                                          | yes      |
| model 1 exposes all 11 monitored parameters                               | 11                    | 11                                                              | yes      |
| DIC of the three models is finite and near the report's 5763-5768         | within 5700-5800      | 5759.9/5758.9/5762.9                                            | yes      |
| DIC finite and model 2 (weather) has the lowest DIC, as the report argues | min at model2_weather | model2_weather 5758.9                                           | yes      |
| corrupted data aborts the analysis                                        | error raised          | error raised                                                    | yes      |
| missing data file aborts the load                                         | error raised          | error raised                                                    | yes      |
| the test suite fails on a corrupted repository (exit != 0)                | non-zero exit         | suite exited 1 and aborted on the corrupted column, as expected | yes      |
checks: 20 | failed: 0
RESULT: PASS
log written to /home/dsh/bl-work/tests/run_tests.log
```

This is the whole output except the JAGS progress lines; `tests/run_tests.log` holds the
same transcript and is rewritten on every run. The three models as published are what the
run fits — the tests rehearse their results, they do not reshape them.

## Structure

```
Bayesian Learning and Monte Carlo Simulation - Final Project.Rmd   346 lines; the report and the analysis
Bayesian Learning and Monte Carlo Simulation - Final Project.R     178 lines; Windows draft of the same analysis
Bayesian Learning and Monte Carlo Simulation - Final Project.pdf   18 pages; the report as submitted
forestfire.txt                                                     518 lines; 517 observations x 13 variables
tests/run_tests.R                                                  20 checks, three models, two negative checks
Makefile                                                           test, install-deps, report, clean
```

The `.R` file is an earlier draft of the same three models with CRLF line endings and a
hard-coded path (`D:/Ann Wanjiru/Bayesian Analysis/forestfire.txt`) to a machine that is
not this repository's. Its authorship is not stated anywhere in it or in the git history —
the two 2022 commits push both the `.R` and the `.Rmd` in the same upload — so it is left
unattributed here.

## License

Por definir — co-autoria (Inês Jorge da Silva e Ferreira).

The report's YAML header and the PDF metadata name two authors: Inês Jorge da Silva e
Ferreira and André Guilherme dos Santos Neto. No `LICENSE` file exists, so no licence is
granted until both authors agree.

The dataset is not covered by that decision: `forestfire.txt` is the UCI *Forest Fires*
dataset by Cortez and Morais (2007), CC BY 4.0, and redistribution requires the
attribution above.

Os dados não são deste trabalho: ver `DATA-SOURCES.md` para a origem, a autoria e a licença do `forestfire.txt`.
