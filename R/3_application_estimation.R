# Application: GLW model fitted to HMD death counts and exposures
# Six countries, females, every year from 1950 to 2020, ages 50-109
# Age is truncated at 50 and shifted to mid-interval: x = age - 50 + 0.5
# The open interval 110+ is excluded
# Death counts are rounded to integers (HMD counts can be fractional)
# Set path below to the repository folder

library(parallel)

path = "~/Desktop/GitHub/laguerre_weibull_distribution"
source(file.path(path, "R", "0_glw_functions.R"))

countries = c("DNK", "FRATNP", "GBRTENW", "ITA", "SWE", "USA")
years = 1950:2020
n_cores = max(1, detectCores() - 1)
sexes = c("Female")
age_min = 50
age_max = 109

# ------------------------------------------------------------------
# Data
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
    exposures = read_hmd(file.path(path, "data", country, "STATS", "Exposures_1x1.txt"))
  )
}

get_data = function(country, year, sex) {
  deaths = hmd[[country]]$deaths
  exposures = hmd[[country]]$exposures
  keep_d = deaths$Year == year & deaths$Age >= age_min & deaths$Age <= age_max
  keep_e = exposures$Year == year & exposures$Age >= age_min & exposures$Age <= age_max
  stopifnot(identical(deaths$Age[keep_d], exposures$Age[keep_e]))
  out = data.frame(
    age = deaths$Age[keep_d],
    D = round(deaths[[sex]][keep_d]),
    E = exposures[[sex]][keep_e]
  )
  out = out[out$E > 0 & !is.na(out$D) & !is.na(out$E), ]
  out$x = out$age - age_min + 0.5
  out
}

# ------------------------------------------------------------------
# Estimation: one fit of the GLW model and one of the gamma-Gompertz model (b = 1)
# for each country, sex and year, run in parallel
# ------------------------------------------------------------------
fit_one = function(task) {
  dat = get_data(task$country, task$year, task$sex)
  full = fit_glw(dat)
  gg = fit_glw(dat, fixed_b = 1)
  se = sqrt(diag(solve(expected_information(dat, full$par))))
  lr = 2 * (full$loglik - gg$loglik)
  est = data.frame(
    country = task$country, sex = task$sex, year = task$year,
    theta = full$par[1], a = full$par[2], b = full$par[3], sigma2 = full$par[4],
    se_theta = se[1], se_a = se[2], se_b = se[3], se_sigma2 = se[4],
    loglik = full$loglik, convergence = full$convergence, boundary = full$boundary,
    theta_b1 = gg$par[1], a_b1 = gg$par[2], sigma2_b1 = gg$par[4],
    loglik_b1 = gg$loglik, convergence_b1 = gg$convergence, boundary_b1 = gg$boundary,
    lr_b1 = lr, p_b1 = pchisq(lr, df = 1, lower.tail = FALSE), n_ages = nrow(dat)
  )
  est
}

tasks = expand.grid(year = years, sex = sexes, country = countries, stringsAsFactors = FALSE)
tasks = split(tasks, seq_len(nrow(tasks)))
start_time = Sys.time()
fits = mclapply(tasks, fit_one, mc.cores = n_cores)
cat("Estimation done in", round(as.numeric(difftime(Sys.time(), start_time, units = "mins")), 1), "min\n")
results = do.call(rbind, fits)
rownames(results) = NULL
write.csv(results, file.path(path, "results", "application_estimates.csv"), row.names = FALSE)
