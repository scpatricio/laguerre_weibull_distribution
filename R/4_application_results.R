# Application: figures and actuarial values from the GLW fits
# Run R/3_application_estimation.R first; this script reads results/application_estimates.csv
# 1. Parameter estimates (black line) with 95% Wald confidence intervals (grey band),
#    one figure per parameter, one panel per country
# 2. Estimated minus actual (HMD) remaining life expectancy at ages 50, 65 and 80:
#    difference in years (top row) and relative difference in % (bottom row),
#    columns are ages, colours are countries
# 3. Single life annuity and single assurance at ages 50, 65 and 80, force of interest 0.05
# Model time x is measured from exact age 50, so exact age y corresponds to x = y - 50
# Set path below to the repository folder

library(ggplot2)

path = "~/Desktop/GitHub/laguerre_weibull_distribution"
source(file.path(path, "R", "0_glw_functions.R"))

country_names = c(DNK = "Denmark", FRATNP = "France", GBRTENW = "England and Wales",
                  ITA = "Italy", SWE = "Sweden", USA = "United States")
age_min = 50
ages_report = c(50, 65, 80)
delta = 0.05
table_years = seq(1950, 2020, by = 10)
out_dirs = file.path(path, c("R/figures/application", "LaTeX/figures/application"))

est = read.csv(file.path(path, "results", "application_estimates.csv"))
est$country_name = factor(country_names[est$country], levels = country_names)

# ------------------------------------------------------------------
# 1. Parameter estimates with 95% confidence intervals
# ------------------------------------------------------------------
parameters = c(theta = "theta", a = "a", b = "b", sigma2 = "sigma^2")
z = qnorm(0.975)
for (p in names(parameters)) {
  d = data.frame(country_name = est$country_name, sex = est$sex, year = est$year,
                 value = est[[p]], se = est[[paste0("se_", p)]])
  d$lower = pmax(d$value - z * d$se, 0)
  d$upper = d$value + z * d$se
  fig = ggplot(d, aes(x = year, y = value)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), fill = "grey75") +
    geom_line(colour = "black", linewidth = 0.6) +
    facet_wrap(~country_name, ncol = 3) +
    # scale_y_log10() +
    labs(x = "Year", y = parse(text = parameters[[p]]),
         title = parse(text = paste0("'Estimates of '*", parameters[[p]], "*' and 95% confidence intervals, females'"))) +
    theme_minimal(base_size = 11)
  for (out_dir in out_dirs) {
    ggsave(file.path(out_dir, paste0("estimates_", p, ".pdf")), fig, width = 11, height = 5)
  }
}

# ------------------------------------------------------------------
# GLW survival from exact age 50, and actuarial functions at exact age y
# ------------------------------------------------------------------
glw_survival = function(x, par) {
  theta = par[1]
  a = par[2]
  b = par[3]
  s2 = par[4]
  if (s2 < 1e-10) {
    return(exp(-theta * (exp(a * x^b) - 1)))
  }
  (1 + s2 * theta * (exp(a * x^b) - 1))^(-1 / s2)
}

# Integral of exp(-force * t) * S(x + t) / S(x) over t > 0
discounted_life = function(y, par, force) {
  x = y - age_min
  s_x = glw_survival(x, par)
  integrand = function(t) exp(-force * t) * glw_survival(x + t, par) / s_x
  integrate(integrand, 0, 120, rel.tol = 1e-10)$value
}

life_expectancy = function(y, par) discounted_life(y, par, 0)
annuity = function(y, par) discounted_life(y, par, delta)
assurance = function(y, par) 1 - delta * annuity(y, par)

actuarial = list()
for (i in seq_len(nrow(est))) {
  par = unlist(est[i, c("theta", "a", "b", "sigma2")])
  for (y in ages_report) {
    actuarial[[length(actuarial) + 1]] = data.frame(
      country = est$country[i], sex = est$sex[i], year = est$year[i], age = y,
      ex_hat = life_expectancy(y, par), annuity = annuity(y, par), assurance = assurance(y, par)
    )
  }
}
actuarial = do.call(rbind, actuarial)

# ------------------------------------------------------------------
# 2. Estimated minus actual (HMD life tables) remaining life expectancy
# ------------------------------------------------------------------
read_lt = function(country, sex) {
  file = paste0(ifelse(sex == "Female", "f", "m"), "ltper_1x1.txt")
  lt = read.table(file.path(path, "data", country, "STATS", file), skip = 2, header = TRUE,
                  na.strings = ".", stringsAsFactors = FALSE)
  lt$Age = as.numeric(sub("\\+", "", lt$Age))
  data.frame(country = country, sex = sex, year = lt$Year, age = lt$Age, ex = lt$ex)
}
lt = do.call(rbind, lapply(names(country_names), function(cn) {
  do.call(rbind, lapply(unique(est$sex), function(s) read_lt(cn, s)))
}))
actuarial = merge(actuarial, lt, by = c("country", "sex", "year", "age"), all.x = TRUE)
actuarial$ex_diff = actuarial$ex_hat - actuarial$ex
actuarial$ex_reldiff = 100 * (actuarial$ex_hat - actuarial$ex) / actuarial$ex
actuarial = actuarial[order(actuarial$country, actuarial$sex, actuarial$year, actuarial$age), ]
write.csv(actuarial, file.path(path, "results", "application_actuarial.csv"), row.names = FALSE)

measures = c("Difference (years)", "Relative difference (%)")
d = rbind(
  data.frame(actuarial[, c("country", "year", "age")], measure = measures[1], value = actuarial$ex_diff),
  data.frame(actuarial[, c("country", "year", "age")], measure = measures[2], value = actuarial$ex_reldiff)
)
d$measure = factor(d$measure, levels = measures)
d$country_name = factor(country_names[d$country], levels = country_names)
d$age_label = factor(paste("Age", d$age), levels = paste("Age", ages_report))
fig = ggplot(d, aes(x = year, y = value, colour = country_name)) +
  geom_hline(yintercept = 0, linewidth = 0.4) +
  geom_line(linewidth = 0.6) +
  facet_grid(measure ~ age_label, scales = "free_y", switch = "y") +
  labs(x = "Year", y = NULL, colour = NULL,
       title = "Estimated minus actual remaining life expectancy, GLW model, females") +
  guides(colour = guide_legend(nrow = 1)) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom", strip.placement = "outside")
for (out_dir in out_dirs) {
  ggsave(file.path(out_dir, "life_expectancy_difference.pdf"), fig, width = 11, height = 5)
}

# ------------------------------------------------------------------
# 3. Single life annuity and single assurance for the appendix tables
# ------------------------------------------------------------------
tab = actuarial[actuarial$year %in% table_years, c("country", "sex", "year", "age", "annuity", "assurance")]
tab = reshape(tab, idvar = c("country", "sex", "year"), timevar = "age", direction = "wide")
tab$country = country_names[tab$country]
tab = tab[order(tab$sex, tab$country, tab$year), ]
write.csv(tab, file.path(path, "results", "application_annuity_assurance.csv"), row.names = FALSE)

# ------------------------------------------------------------------
# Appendix table (LaTeX): annuities and assurances every 10 years
# ------------------------------------------------------------------
tex = c("\\begin{table}[h]", "\\centering", "\\footnotesize", "\\renewcommand{\\arraystretch}{0.92}",
  "\\caption{Estimated single life annuity $\\bar{a}_x$ and single assurance $\\bar{A}_x$ at ages $x=50$, 65 and 80 for females, with force of interest $\\delta=0.05$.}",
  "\\label{tab:actuarial}",
  "\\begin{tabular*}{\\textwidth}{@{\\extracolsep{\\fill}}llrrrrrr}", "\\toprule",
  "& & \\multicolumn{3}{c}{Single life annuity} & \\multicolumn{3}{c}{Single assurance}\\\\",
  "\\cmidrule(lr){3-5}\\cmidrule(lr){6-8}",
  "Country & Year & $\\bar{a}_{50}$ & $\\bar{a}_{65}$ & $\\bar{a}_{80}$ & $\\bar{A}_{50}$ & $\\bar{A}_{65}$ & $\\bar{A}_{80}$\\\\",
  "\\midrule")
tab = tab[order(tab$country, tab$year), ]
previous = ""
for (i in seq_len(nrow(tab))) {
  label = ifelse(tab$country[i] == previous, "", tab$country[i])
  if (label != "" && previous != "") {
    tex = c(tex, "\\addlinespace[2pt]")
  }
  previous = tab$country[i]
  values = unlist(tab[i, c(paste0("annuity.", ages_report), paste0("assurance.", ages_report))])
  tex = c(tex, paste0(label, " & ", tab$year[i], " & ", paste(sprintf("%.4f", values), collapse = " & "), "\\\\"))
}
tex = c(tex, "\\bottomrule", "\\end{tabular*}", "\\end{table}")
writeLines(tex, file.path(path, "results", "application_actuarial_table.tex"))
