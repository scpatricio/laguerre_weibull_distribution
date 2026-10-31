# Monte Carlo simulation study for the GLW model fitted by Poisson maximum likelihood
# Individual lifetimes are simulated from the GLW model (gamma frailty, then LW lifetime
# given the frailty) for a cohort followed from x = 0 until the last death. Deaths and
# exposures are then aggregated in one-year intervals (midpoints x = 0.5, 1.5, ...) and
# the model is fitted with the Poisson likelihood
# Set path below to the repository folder

library(parallel)

path = "~/Desktop/GitHub/laguerre_weibull_distribution"
source(file.path(path, "R", "0_glw_functions.R"))

n_rep = 1000
n_cores = max(1, detectCores() - 1)

# ------------------------------------------------------------------
# Scenarios: full factorial design, two values per parameter (2^4 = 16)
# b = 1 (gamma-Gompertz) in half of the scenarios, so the size of the
# likelihood-ratio test is assessed in 8 scenarios and its power in the other 8
# ------------------------------------------------------------------
scenarios = expand.grid(
  theta  = c(0.02, 0.04),
  a      = c(0.075, 0.10),
  b      = c(1, 1.1),
  sigma2 = c(0.05, 0.1)
)
scenarios = cbind(scenario = paste0("S", seq_len(nrow(scenarios))), scenarios)

# ------------------------------------------------------------------
# Cohort sizes (number of simulated lives)
# ------------------------------------------------------------------
cohort_sizes = c(small = 20000, large = 50000)

# ------------------------------------------------------------------
# Simulate one cohort followed from x = 0 until the last death
# 1. Frailty Z ~ gamma with mean 1 and variance sigma2
# 2. Lifetime given Z by inversion of the LW survival function
#    S(x | z) = exp[-z theta {exp(a x^b) - 1}], that is,
#    X = [log{1 - log(U) / (z theta)} / a]^(1/b), U ~ uniform(0, 1)
# 3. Deaths D and exposures E (person-years lived) in each interval [k, k + 1),
#    k = 0, 1, ..., up to the interval of the last death, with x = k + 0.5
# ------------------------------------------------------------------
simulate_cohort = function(n, par) {
  theta = par[1]
  a = par[2]
  b = par[3]
  s2 = par[4]
  z = rgamma(n, shape = 1 / s2, scale = s2)
  u = runif(n)
  lifetime = (log(1 - log(u) / (z * theta)) / a)^(1 / b)
  interval = floor(lifetime) + 1
  n_int = max(interval)
  D = tabulate(interval, nbins = n_int)
  alive_start = n - c(0, cumsum(D))[1:n_int]
  alive_end = alive_start - D
  partial = as.numeric(tapply(lifetime - (interval - 1), factor(interval, levels = 1:n_int), sum))
  partial[is.na(partial)] = 0
  E = alive_end + partial
  data.frame(x = 1:n_int - 0.5, D = D, E = E)
}

# ------------------------------------------------------------------
# One replication: simulate, fit GLW and gamma-Gompertz (b = 1), store results
# Each task gets its own seed, so results do not depend on the number of cores
# ------------------------------------------------------------------
one_rep = function(task) {
  set.seed(2026 + task$id)
  par = unlist(scenarios[task$k, c("theta", "a", "b", "sigma2")])
  dat = simulate_cohort(cohort_sizes[[task$pop]], par)
  full = fit_glw(dat, s2_grid = c(0.05, 0.2), b_grid = c(0.8, 1, 1.2))
  null = fit_glw(dat, fixed_b = 1, s2_grid = c(0.05, 0.2), b_grid = 1)
  info = expected_information(dat, full$par)
  se = tryCatch(sqrt(diag(solve(info))), error = function(e) rep(NA, 4))
  c(full$par, se, full$convergence, 2 * (full$loglik - null$loglik))
}

# ------------------------------------------------------------------
# Run: all scenario x population x replication tasks in one parallel pool
# ------------------------------------------------------------------
grid = expand.grid(rep = seq_len(n_rep), pop = names(cohort_sizes), k = seq_len(nrow(scenarios)),
                   stringsAsFactors = FALSE)
grid$id = seq_len(nrow(grid))
tasks = split(grid, grid$id)

chunks = split(tasks, ceiling(seq_along(tasks) / n_rep))
out = list()
start_time = Sys.time()
for (j in seq_along(chunks)) {
  out = c(out, mclapply(chunks[[j]], one_rep, mc.cores = n_cores, mc.preschedule = TRUE))
  elapsed = round(as.numeric(difftime(Sys.time(), start_time, units = "secs")), 1)
  cat(length(out), "of", length(tasks), "replications done,", elapsed, "secs\n")
}
out = do.call(rbind, out)
colnames(out) = c("theta", "a", "b", "sigma2", "se_theta", "se_a", "se_b", "se_sigma2",
                  "convergence", "lr_b1")
raw = data.frame(scenario = scenarios$scenario[grid$k], population = grid$pop, rep = grid$rep, out)
write.csv(raw, file.path(path, "results", "simulation_raw.csv"), row.names = FALSE)

# ------------------------------------------------------------------
# Summary measures and Table 1 of the manuscript
# ------------------------------------------------------------------
source(file.path(path, "R", "0_simulation_summary.R"))
