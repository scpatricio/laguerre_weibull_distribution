# Application: parametric bootstrap of the likelihood-ratio test of b = 1
# For the fits from 1990 onwards in which the estimate of sigma2 is on the boundary
# (sigma2 = 0) under the GLW model or under the gamma-Gompertz model, the chi-squared(1)
# reference distribution may not hold. For each such population-year we simulate death
# counts D* ~ Poisson(E * mu0(x)) from the fitted gamma-Gompertz model (the null model),
# refit both models, and compute the bootstrap p-value (1 + #{LR* >= LR}) / (B + 1).
# Run R/3_application_estimation.R first. Set path below to the repository folder

library(parallel)

path = "~/Desktop/GitHub/laguerre_weibull_distribution"
source(file.path(path, "R", "0_glw_functions.R"))

n_boot = 200
first_year = 1990
age_min = 50
age_max = 109
sex = "Female"
n_cores = max(1, detectCores() - 1)

read_hmd = function(file) {
  out = read.table(file, skip = 2, header = TRUE, na.strings = ".", stringsAsFactors = FALSE)
  out$Age = as.numeric(sub("\\+", "", out$Age))
  out
}

get_data = function(country, year) {
  deaths = read_hmd(file.path(path, "data", country, "STATS", "Deaths_1x1.txt"))
  exposures = read_hmd(file.path(path, "data", country, "STATS", "Exposures_1x1.txt"))
  keep_d = deaths$Year == year & deaths$Age >= age_min & deaths$Age <= age_max
  keep_e = exposures$Year == year & exposures$Age >= age_min & exposures$Age <= age_max
  out = data.frame(age = deaths$Age[keep_d], D = round(deaths[[sex]][keep_d]),
                   E = exposures[[sex]][keep_e])
  out = out[out$E > 0 & !is.na(out$D) & !is.na(out$E), ]
  out$x = out$age - age_min + 0.5
  out
}

est = read.csv(file.path(path, "results", "application_estimates.csv"))
affected = est[est$sex == sex & est$year >= first_year & (est$boundary | est$boundary_b1), ]
cat(nrow(affected), "population-years to bootstrap\n")

datasets = lapply(seq_len(nrow(affected)), function(i) get_data(affected$country[i], affected$year[i]))

tasks = expand.grid(r = seq_len(n_boot), i = seq_len(nrow(affected)))
one_boot = function(t) {
  i = tasks$i[t]
  set.seed(2026 + t)
  dat = datasets[[i]]
  par0 = c(affected$theta_b1[i], affected$a_b1[i], 1, affected$sigma2_b1[i])
  dat$D = rpois(nrow(dat), dat$E * glw_hazard(dat$x, par0))
  full = fit_glw(dat, s2_grid = c(0.02, 0.1), b_grid = c(0.9, 1, 1.2))
  null = fit_glw(dat, fixed_b = 1, s2_grid = c(0.02, 0.1), b_grid = 1)
  2 * (full$loglik - null$loglik)
}
start_time = Sys.time()
lr_boot = unlist(mclapply(seq_len(nrow(tasks)), one_boot, mc.cores = n_cores))
cat("Bootstrap done in", round(as.numeric(difftime(Sys.time(), start_time, units = "mins")), 1), "min\n")

out = affected[, c("country", "year", "sigma2", "sigma2_b1", "boundary", "boundary_b1", "lr_b1", "p_b1")]
out$lr_boot_95 = tapply(lr_boot, tasks$i, quantile, probs = 0.95)
out$lr_boot_max = tapply(lr_boot, tasks$i, max)
out$p_boot = sapply(seq_len(nrow(affected)), function(i) {
  (1 + sum(lr_boot[tasks$i == i] >= affected$lr_b1[i])) / (n_boot + 1)
})
out$n_boot = n_boot
write.csv(out, file.path(path, "results", "application_bootstrap.csv"), row.names = FALSE)
print(out)
