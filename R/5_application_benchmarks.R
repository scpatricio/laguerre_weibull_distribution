# Application: comparison of the GLW model with standard mortality laws
# Models fitted by Poisson maximum likelihood to the same data as R/3_application_estimation.R
# (females, six countries, every year 1950-2020, ages 50-109, x = age - 50 + 0.5):
#   Gompertz           mu(x) = alpha exp(beta x)
#   Kannisto           mu(x) = alpha exp(beta x) / {1 + alpha exp(beta x)}
#   gamma-Gompertz     the GLW model with b = 1 (Beard's law), taken from application_estimates.csv
#   GLW                taken from application_estimates.csv
# Outputs: results/application_benchmarks.csv (all population-years)
#          results/application_benchmarks_summary.csv and application_benchmarks_table.tex
#            (table in Section 8 of the manuscript)
#          figures/application/benchmark_ex50.pdf, benchmark_ex65.pdf, benchmark_ex80.pdf
#            (appendix figures)
#          figures/application/benchmark_daic.pdf and benchmark_hazards.pdf (appendix figures)
# Run R/3_application_estimation.R first. Set path below to the repository folder

library(ggplot2)
library(patchwork)
library(parallel)

path = "~/Desktop/GitHub/laguerre_weibull_distribution"
source(file.path(path, "R", "0_glw_functions.R"))

countries = c("DNK", "FRATNP", "GBRTENW", "ITA", "SWE", "USA")
country_names = c(DNK = "Denmark", FRATNP = "France", GBRTENW = "England and Wales",
                  ITA = "Italy", SWE = "Sweden", USA = "United States")
sex = "Female"
age_min = 50
age_max = 109
ages_report = c(50, 65, 80)
figure_year = 2020
n_cores = max(1, detectCores() - 1)
out_dirs = file.path(path, c("R/figures/application", "LaTeX/figures/application"))

# ------------------------------------------------------------------
# Data (same selection as R/3_application_estimation.R)
# ------------------------------------------------------------------
read_hmd = function(file) {
  out = read.table(file, skip = 2, header = TRUE, na.strings = ".", stringsAsFactors = FALSE)
  out$Age = as.numeric(sub("\\+", "", out$Age))
  out
}

hmd = list()
for (country in countries) {
  hmd[[country]] = list(
    deaths = read_hmd(file.path(path, "data", country, "STATS", "Deaths_1x1.txt")),
    exposures = read_hmd(file.path(path, "data", country, "STATS", "Exposures_1x1.txt")),
    lt = read_hmd(file.path(path, "data", country, "STATS", "fltper_1x1.txt"))
  )
}

get_data = function(country, year) {
  deaths = hmd[[country]]$deaths
  exposures = hmd[[country]]$exposures
  keep_d = deaths$Year == year & deaths$Age >= age_min & deaths$Age <= age_max
  keep_e = exposures$Year == year & exposures$Age >= age_min & exposures$Age <= age_max
  stopifnot(identical(deaths$Age[keep_d], exposures$Age[keep_e]))
  out = data.frame(age = deaths$Age[keep_d], D = round(deaths[[sex]][keep_d]),
                   E = exposures[[sex]][keep_e])
  out = out[out$E > 0 & !is.na(out$D) & !is.na(out$E), ]
  out$x = out$age - age_min + 0.5
  out
}

# ------------------------------------------------------------------
# Hazards and survival functions from exact age 50 (x = 0)
# ------------------------------------------------------------------
gompertz_hazard = function(x, par) par[1] * exp(par[2] * x)
gompertz_survival = function(x, par) exp(-par[1] / par[2] * (exp(par[2] * x) - 1))
kannisto_hazard = function(x, par) par[1] * exp(par[2] * x) / (1 + par[1] * exp(par[2] * x))
kannisto_survival = function(x, par) ((1 + par[1]) / (1 + par[1] * exp(par[2] * x)))^(1 / par[2])
glw_surv = function(x, par) {
  if (par[4] < boundary_tol) {
    return(exp(-par[1] * (exp(par[2] * x^par[3]) - 1)))
  }
  (1 + par[4] * par[1] * (exp(par[2] * x^par[3]) - 1))^(-1 / par[4])
}

# Remaining life expectancy at exact age y
life_expectancy = function(y, survival, par) {
  x = y - age_min
  s_x = survival(x, par)
  integrate(function(t) survival(x + t, par) / s_x, 0, 120, rel.tol = 1e-10)$value
}

# Poisson maximum likelihood for a hazard with parameters bounded below
fit_law = function(dat, hazard, starts, lower) {
  nll = function(par) {
    value = -sum(dpois(dat$D, lambda = dat$E * hazard(dat$x, par), log = TRUE))
    if (!is.finite(value)) {
      value = 1e100
    }
    value
  }
  best = NULL
  for (st in starts) {
    fit = tryCatch(nlminb(pmax(st, lower), nll, lower = lower,
                          control = list(eval.max = 5000, iter.max = 5000, rel.tol = 1e-12)),
                   error = function(e) NULL)
    if (!is.null(fit) && (is.null(best) || fit$objective < best$objective)) {
      best = fit
    }
  }
  list(par = unname(best$par), loglik = -best$objective,
       convergence = as.integer(grepl("limit", best$message)))
}

# ------------------------------------------------------------------
# Fits
# ------------------------------------------------------------------
est = read.csv(file.path(path, "results", "application_estimates.csv"))
est = est[est$sex == sex, ]

fit_one = function(i) {
  e = est[i, ]
  dat = get_data(e$country, e$year)
  n = nrow(dat)
  gomp = coef(lm(log(D / E) ~ x, data = dat[dat$D > 0, ]))
  alpha = exp(gomp[1])
  beta = gomp[2]
  gom = fit_law(dat, gompertz_hazard,
                starts = list(c(alpha, beta), c(alpha, 1.2 * beta), c(alpha / 2, beta)),
                lower = c(1e-10, 1e-6))
  kan = fit_law(dat, kannisto_hazard,
                starts = list(c(alpha, beta), c(alpha, 1.2 * beta), c(alpha / 2, beta)),
                lower = c(1e-10, 1e-6))
  models = list(
    "Gompertz" = list(par = gom$par, loglik = gom$loglik, k = 2, survival = gompertz_survival,
                      hazard = gompertz_hazard, convergence = gom$convergence),
    "Kannisto" = list(par = kan$par, loglik = kan$loglik, k = 2, survival = kannisto_survival,
                      hazard = kannisto_hazard, convergence = kan$convergence),
    "Gamma-Gompertz" = list(par = c(e$theta_b1, e$a_b1, 1, e$sigma2_b1), loglik = e$loglik_b1, k = 3,
                            survival = glw_surv, hazard = glw_hazard, convergence = e$convergence_b1),
    "GLW" = list(par = c(e$theta, e$a, e$b, e$sigma2), loglik = e$loglik, k = 4,
                 survival = glw_surv, hazard = glw_hazard, convergence = e$convergence)
  )
  lt = hmd[[e$country]]$lt
  out = list()
  curves = list()
  for (m in names(models)) {
    mod = models[[m]]
    ex_hat = sapply(ages_report, life_expectancy, survival = mod$survival, par = mod$par)
    ex = sapply(ages_report, function(y) lt$ex[lt$Year == e$year & lt$Age == y])
    out[[m]] = data.frame(country = e$country, year = e$year, model = m, n_ages = n, k = mod$k,
                          loglik = mod$loglik, aic = 2 * mod$k - 2 * mod$loglik,
                          bic = log(n) * mod$k - 2 * mod$loglik, convergence = mod$convergence,
                          ex_50 = ex[1], ex_65 = ex[2], ex_80 = ex[3],
                          ex_err_50 = ex_hat[1] - ex[1], ex_err_65 = ex_hat[2] - ex[2],
                          ex_err_80 = ex_hat[3] - ex[3])
    if (e$year == figure_year) {
      curves[[m]] = data.frame(country = e$country, age = dat$age, model = m,
                               mu = mod$hazard(dat$x, mod$par), observed = dat$D / dat$E)
    }
  }
  list(fit = do.call(rbind, out), curves = do.call(rbind, curves))
}

fits = mclapply(seq_len(nrow(est)), fit_one, mc.cores = n_cores)
bench = do.call(rbind, lapply(fits, function(f) f$fit))
rownames(bench) = NULL
model_levels = c("Gompertz", "Kannisto", "Gamma-Gompertz", "GLW")
bench$model = factor(bench$model, levels = model_levels)

# ------------------------------------------------------------------
# Summary over all population-years: Delta AIC_m = AIC_m - min_m AIC_m within each
# population-year, share of population-years in which each model has the lowest AIC,
# and mean absolute error of the remaining life expectancy at ages 50, 65 and 80
# ------------------------------------------------------------------
key = paste(bench$country, bench$year)
bench$daic = bench$aic - tapply(bench$aic, key, min)[key]
bench$best_aic = bench$daic == 0
write.csv(bench, file.path(path, "results", "application_benchmarks.csv"), row.names = FALSE)
summary_table = do.call(rbind, lapply(rev(model_levels), function(m) {
  d = bench[bench$model == m, ]
  data.frame(model = m, k = d$k[1], median_daic = median(d$daic), share_best_aic = mean(d$best_aic),
             mae_ex_50 = mean(abs(d$ex_err_50)), mae_ex_65 = mean(abs(d$ex_err_65)),
             mae_ex_80 = mean(abs(d$ex_err_80)), non_converged = sum(d$convergence != 0))
}))
write.csv(summary_table, file.path(path, "results", "application_benchmarks_summary.csv"), row.names = FALSE)
print(summary_table)

# ------------------------------------------------------------------
# Supplementary figure: Delta AIC by country and year
# ------------------------------------------------------------------
bench$country_name = factor(country_names[bench$country], levels = country_names)
fig = ggplot(bench, aes(x = year, y = daic, colour = model)) +
  geom_line(linewidth = 0.6) +
  scale_y_continuous(trans = "log1p", breaks = c(0, 10, 100, 1000, 10000)) +
  facet_wrap(~country_name, ncol = 3) +
  guides(colour = guide_legend(nrow = 1)) +
  labs(x = "Year", y = expression(Delta * "AIC (log scale)"), colour = NULL,
       title = expression(Delta * "AIC of four mortality laws, females")) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom")
for (out_dir in out_dirs) {
  ggsave(file.path(out_dir, "benchmark_daic.pdf"), fig, width = 11, height = 5)
}

# ------------------------------------------------------------------
# Figure: observed and fitted hazards in the last year, by country
# ------------------------------------------------------------------
curves = do.call(rbind, lapply(fits, function(f) f$curves))
curves$country_name = factor(country_names[curves$country], levels = country_names)
curves$model = factor(curves$model, levels = model_levels)
obs = unique(curves[curves$observed > 0, c("country_name", "age", "observed")])
fig = ggplot(curves, aes(x = age)) +
  geom_point(data = obs, aes(y = observed), colour = "grey60", size = 0.8) +
  geom_line(aes(y = mu, colour = model), linewidth = 0.6) +
  scale_y_log10(labels = scales::label_log()) +
  facet_wrap(~country_name, ncol = 3) +
  guides(colour = guide_legend(nrow = 1)) +
  labs(x = "Age", y = "Hazard (log scale)", colour = NULL,
       title = paste0("Observed and fitted hazards, females, ", figure_year)) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom")
for (out_dir in out_dirs) {
  ggsave(file.path(out_dir, "benchmark_hazards.pdf"), fig, width = 11, height = 5)
}

# ------------------------------------------------------------------
# Appendix figures: estimated minus actual remaining life expectancy at ages 50, 65 and 80,
# one figure per age; top two rows: difference in years, bottom two rows: relative
# difference in percent of the actual value; one panel per country, one line per model
# ------------------------------------------------------------------
measures = c("Difference (years)", "Relative difference (%)")
for (y in ages_report) {
  error = bench[[paste0("ex_err_", y)]]
  d = rbind(
    data.frame(country_name = bench$country_name, year = bench$year, model = bench$model,
               measure = measures[1], value = error),
    data.frame(country_name = bench$country_name, year = bench$year, model = bench$model,
               measure = measures[2], value = 100 * error / bench[[paste0("ex_", y)]])
  )
  d$measure = factor(d$measure, levels = measures)
  panels = lapply(measures, function(m) {
    ggplot(d[d$measure == m, ], aes(x = year, y = value, colour = model)) +
      geom_hline(yintercept = 0, linewidth = 0.4) +
      geom_line(linewidth = 0.6) +
      facet_wrap(~country_name, ncol = 3) +
      labs(x = "Year", y = m, colour = NULL) +
      theme_minimal(base_size = 11)
  })
  fig = panels[[1]] / panels[[2]] +
    plot_layout(guides = "collect") +
    plot_annotation(title = paste0("Estimated minus actual remaining life expectancy at age ", y, ", females")) &
    guides(colour = guide_legend(nrow = 1)) &
    theme(legend.position = "bottom")
  for (out_dir in out_dirs) {
    ggsave(file.path(out_dir, paste0("benchmark_ex", y, ".pdf")), fig, width = 11, height = 10)
  }
}

# ------------------------------------------------------------------
# Table of the manuscript (LaTeX), built from the summary above
# ------------------------------------------------------------------
tab = summary_table
tex = c("\\begin{table}[htb!]", "\\centering", "\\small",
  "\\caption{Comparison of the GLW model with three mortality laws over the 426 population-years (females, ages 50--109, 1950--2020). $\\Delta\\mathrm{AIC}_m=\\mathrm{AIC}_m-\\min_m\\mathrm{AIC}_m$ within each population-year; Best AIC is the percentage of population-years in which the model has the lowest AIC; MAE is the mean absolute error, in years, of the estimated remaining life expectancy.}",
  "\\label{tab:benchmark}",
  "\\begin{tabular*}{\\textwidth}{@{\\extracolsep{\\fill}}lrrrrrr}", "\\toprule",
  "Model & Parameters & Median $\\Delta\\mathrm{AIC}$ & Best AIC (\\%) & MAE $\\mathring{e}_{50}$ & MAE $\\mathring{e}_{65}$ & MAE $\\mathring{e}_{80}$\\\\",
  "\\midrule")
for (i in seq_len(nrow(tab))) {
  tex = c(tex, paste0(tab$model[i], " & ", tab$k[i], " & ", sprintf("%.1f", tab$median_daic[i]), " & ",
                      sprintf("%.1f", 100 * tab$share_best_aic[i]), " & ", sprintf("%.3f", tab$mae_ex_50[i]), " & ",
                      sprintf("%.3f", tab$mae_ex_65[i]), " & ", sprintf("%.3f", tab$mae_ex_80[i]), "\\\\"))
}
tex = c(tex, "\\bottomrule", "\\end{tabular*}", "\\end{table}")
writeLines(tex, file.path(path, "results", "application_benchmarks_table.tex"))
