# A Gamma-Frailty Laguerre–Weibull Mortality Model with Actuarial Applications

Code and manuscript for the paper by Fredy Castellares (Universidade Federal de São Paulo) and Silvio C. Patricio (University of Southern Denmark).

The paper derives the Laguerre–Weibull (LW) distribution from the generating function of the Laguerre polynomials and adds gamma frailty to obtain the gamma–Laguerre–Weibull (GLW) mortality model, which contains the gamma-Gompertz model as the case b = 1. It gives the statistical and actuarial properties of the model, estimation by Poisson maximum likelihood, a likelihood-ratio test of b = 1 with a bootstrap calibration, a Monte Carlo simulation study, and an application to female mortality in six countries from the Human Mortality Database, with a comparison against the Gompertz, Kannisto and gamma-Gompertz laws.

## Preprint and archived version

- Preprint: arXiv:XXXX.XXXXX (https://arxiv.org/abs/XXXX.XXXXX)
- Archived code and results: Zenodo, https://doi.org/10.5281/zenodo.XXXXXXX

## Repository structure

```
manuscript.pdf                 compiled manuscript
LaTeX/
  main.tex                     manuscript
  references.bib               bibliography
  figures/                     figures read by the manuscript
R/
  0_glw_functions.R            GLW hazard and its derivatives, Poisson log-likelihood,
                               constrained maximum likelihood fit, Fisher information
  0_simulation_summary.R       summary measures and Table 1 from the simulation results
  1_laguerre_weibull.R         Figure 1
  2_gamma_laguerre_weibull.R   Figure 2
  3_application_estimation.R   fits of the GLW and gamma-Gompertz models to the HMD data
  4_application_results.R      parameter figures, life expectancy, annuities and assurances
  5_application_benchmarks.R   comparison with the Gompertz, Kannisto and gamma-Gompertz laws
  6_application_bootstrap.R    bootstrap calibration of the test of b = 1
  7_simulation_study.R         Monte Carlo simulation study
  figures/                     copies of the figures produced by the scripts
results/                       CSV files and LaTeX tables produced by the scripts
data/                          HMD input files (see Data)
```

## Running the code

The scripts need R with the packages `ggplot2`, `patchwork` and `scales`, and the base package `parallel`. Set the variable `path` at the top of each script to the location of the repository.

Run the numbered scripts in order. The two scripts numbered 0 are not run directly: `0_glw_functions.R` is loaded by the others, and `0_simulation_summary.R` is called at the end of `7_simulation_study.R`.

| Script | Output | Use in the manuscript |
|---|---|---|
| `1_laguerre_weibull.R` | `laguerre_weibull.pdf` | Figure 1 |
| `2_gamma_laguerre_weibull.R` | `gamma_laguerre_weibull.pdf` | Figure 2 |
| `3_application_estimation.R` | `application_estimates.csv` | input to scripts 4 to 6 |
| `4_application_results.R` | `estimates_*.pdf`, `life_expectancy_difference.pdf`, `application_actuarial.csv`, `application_annuity_assurance.csv`, `application_actuarial_table.tex` | parameter and life expectancy figures (Section 8), appendix table |
| `5_application_benchmarks.R` | `benchmark_*.pdf`, `application_benchmarks.csv`, `application_benchmarks_summary.csv`, `application_benchmarks_table.tex` | comparison table (Section 8), appendix figures |
| `6_application_bootstrap.R` | `application_bootstrap.csv` | bootstrap results (Section 8) |
| `7_simulation_study.R` | `simulation_raw.csv`, `simulation_summary.csv`, `simulation_table.tex` | Section 7 and Table 1 |

Figures are written to `R/figures/` and `LaTeX/figures/`, and all other output to `results/`. Scripts 3 to 6 need the HMD files; scripts 4 to 6 also need the output of script 3. Running `0_simulation_summary.R` on its own rebuilds Table 1 from the stored simulation results without rerunning the simulation.

The three LaTeX tables in `results/` are copied into `main.tex`, so after rerunning the scripts, check that the tables in the manuscript still match them.

Scripts 3, 5, 6 and 7 run in parallel on all but one core through `mclapply`, which uses a single core on Windows. The simulation study (1,000 replications for each of 32 settings) takes considerably longer than the other scripts.

## Data

The application uses death counts (`Deaths_1x1.txt`), exposures (`Exposures_1x1.txt`) and female life tables (`fltper_1x1.txt`) from the Human Mortality Database (https://www.mortality.org) for Denmark (`DNK`), France (`FRATNP`), England and Wales (`GBRTENW`), Italy (`ITA`), Sweden (`SWE`) and the United States (`USA`). The scripts read them from `data/<country code>/STATS/`. The data are free after registration, but the HMD user agreement does not allow their redistribution, so they should be downloaded from the HMD rather than shared with this repository.

## Compiling the manuscript

From the `LaTeX/` folder, run `latexmk -pdf main.tex`. This writes `LaTeX/main.pdf`; the compiled version in the repository is `manuscript.pdf`.

## License

See `LICENSE`.
