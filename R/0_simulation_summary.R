# Monte Carlo simulation study: summary measures and Table 1 of the manuscript
# Reads results/simulation_raw.csv (written by R/7_simulation_study.R) and writes
#   results/simulation_summary.csv  absolute bias, RMSE, empirical SD, average asymptotic
#                                   standard error (ASE), coverage of 95% Wald intervals and
#                                   rejection rate of the likelihood-ratio test of b = 1
#   results/simulation_table.tex    Table 1 of the manuscript, built from simulation_summary.csv
# Set path below to the repository folder

if (!exists("path")) {
  path = "~/Desktop/GitHub/laguerre_weibull_distribution"
}

raw = read.csv(file.path(path, "results", "simulation_raw.csv"))

# Scenarios: same 2^4 factorial design as in R/7_simulation_study.R
scenarios = expand.grid(
  theta  = c(0.02, 0.04),
  a      = c(0.075, 0.10),
  b      = c(1, 1.1),
  sigma2 = c(0.05, 0.1)
)
scenarios = cbind(scenario = paste0("S", seq_len(nrow(scenarios))), scenarios)
populations = c("small", "large")
cohort_labels = c(small = "20{,}000", large = "50{,}000")
parameters = c("theta", "a", "b", "sigma2")

# ------------------------------------------------------------------
# Summary measures
# ------------------------------------------------------------------
z = qnorm(0.975)
summary_rows = list()
for (k in seq_len(nrow(scenarios))) {
  for (pop in populations) {
    r = raw[raw$scenario == scenarios$scenario[k] & raw$population == pop, ]
    for (p in parameters) {
      true = scenarios[k, p]
      est = r[[p]]
      se = r[[paste0("se_", p)]]
      summary_rows[[length(summary_rows) + 1]] = data.frame(
        scenario = scenarios$scenario[k], population = pop, parameter = p, true = true,
        mean = mean(est), abs_bias = abs(mean(est) - true), rmse = sqrt(mean((est - true)^2)),
        sd = sd(est), ase = mean(se, na.rm = TRUE),
        coverage = mean(abs(est - true) <= z * se, na.rm = TRUE),
        reject_b1 = mean(r$lr_b1 > qchisq(0.95, df = 1)),
        non_converged = sum(r$convergence != 0)
      )
    }
  }
}
summary_table = do.call(rbind, summary_rows)
write.csv(summary_table, file.path(path, "results", "simulation_summary.csv"), row.names = FALSE)

# ------------------------------------------------------------------
# Table 1 (LaTeX), one row per scenario and cohort size
# ------------------------------------------------------------------
idx = expand.grid(i = 1:2, j = 1:2, k = 1:2, l = 1:2)
labels = sprintf("$(\\theta_%d,a_%d,b_%d,\\sigma^2_%d)$", idx$i, idx$j, idx$k, idx$l)
f4 = function(v) sprintf("%.4f", v)
f3 = function(v) sprintf("%.3f", v)
tex = c("\\begin{sidewaystable}", "\\centering", "\\scriptsize", "\\setlength{\\tabcolsep}{2pt}",
  "\\caption{Simulation results: absolute bias, RMSE, empirical standard deviation (SD), average asymptotic standard error (ASE) and coverage probability (CP) of the 95\\% Wald intervals for each parameter, and rejection rate (RR) of the likelihood-ratio test of $H_0\\colon b=1$ at the 5\\% level. 1{,}000 replications.}",
  "\\label{tab:sim}",
  "\\begin{tabular*}{\\textheight}{@{\\extracolsep{\\fill}}l*{20}{r}r}", "\\toprule",
  "& \\multicolumn{5}{c}{$\\theta$} & \\multicolumn{5}{c}{$a$} & \\multicolumn{5}{c}{$b$} & \\multicolumn{5}{c}{$\\sigma^2$} & \\\\",
  "\\cmidrule(lr){2-6}\\cmidrule(lr){7-11}\\cmidrule(lr){12-16}\\cmidrule(lr){17-21}",
  paste0("Scenario", strrep(" & Bias & RMSE & SD & ASE & CP", 4), " & RR\\\\"), "\\midrule")
for (pop in populations) {
  tex = c(tex, paste0("\\multicolumn{22}{l}{\\emph{$n=", cohort_labels[[pop]], "$}}\\\\[1pt]"))
  for (k in seq_len(nrow(scenarios))) {
    cells = c()
    for (p in parameters) {
      s = summary_table[summary_table$scenario == scenarios$scenario[k] &
                          summary_table$population == pop & summary_table$parameter == p, ]
      cells = c(cells, f4(s$abs_bias), f4(s$rmse), f4(s$sd), f4(s$ase), f3(s$coverage))
    }
    cells = c(cells, f3(s$reject_b1))
    tex = c(tex, paste(labels[k], "&", paste(cells, collapse = " & "), "\\\\"))
  }
  if (pop == populations[1]) {
    tex = c(tex, "\\midrule")
  }
}
tex = c(tex, "\\bottomrule", "\\end{tabular*}", "\\end{sidewaystable}")
writeLines(tex, file.path(path, "results", "simulation_table.tex"))
