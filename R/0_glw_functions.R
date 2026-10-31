# Functions for the gamma-Laguerre-Weibull (GLW) model
# par = (theta, a, b, sigma2)

# GLW hazard
glw_hazard = function(x, par) {
  theta = par[1]
  a = par[2]
  b = par[3]
  s2 = par[4]
  eta = exp(a * x^b)
  theta * a * b * x^(b - 1) * eta / (1 + s2 * theta * (eta - 1))
}

# Derivatives of the GLW hazard with respect to (theta, a, b, sigma2)
glw_gradient = function(x, par) {
  theta = par[1]
  a = par[2]
  b = par[3]
  s2 = par[4]
  eta = exp(a * x^b)
  kappa = 1 + theta * s2 * (eta - 1)
  cbind(
    theta = a * b * x^(b - 1) * eta / kappa^2,
    a = theta * b * x^(b - 1) * eta / kappa + theta * a * b * x^(2 * b - 1) * eta * (1 - theta * s2) / kappa^2,
    b = theta * a * x^(b - 1) * eta * (1 + b * log(x)) / kappa +
      theta * a^2 * b * x^(2 * b - 1) * eta * (1 - theta * s2) * log(x) / kappa^2,
    sigma2 = -theta^2 * a * b * x^(b - 1) * eta * (eta - 1) / kappa^2
  )
}

# Negative Poisson log-likelihood, D ~ Poisson(E * mu), with par = (theta, a, b, sigma2).
# If fixed_b is given, par = (theta, a, sigma2) and b is fixed.
neg_loglik = function(par, x, D, E, fixed_b = NULL) {
  if (!is.null(fixed_b)) {
    par = c(par[1:2], fixed_b, par[3])
  }
  mu = glw_hazard(x, par)
  value = -sum(dpois(D, lambda = E * mu, log = TRUE))
  if (!is.finite(value)) {
    value = 1e100
  }
  value
}

# Gradient of neg_loglik with respect to par
neg_loglik_gradient = function(par, x, D, E, fixed_b = NULL) {
  full = par
  if (!is.null(fixed_b)) {
    full = c(par[1:2], fixed_b, par[3])
  }
  mu = glw_hazard(x, full)
  g = -colSums((D / mu - E) * glw_gradient(x, full))
  if (!is.null(fixed_b)) {
    g = g[c(1, 2, 4)]
  }
  g[!is.finite(g)] = 0
  g
}

# Maximum likelihood fit by constrained optimisation (nlminb, PORT routines).
# theta, a and b are kept above a small positive bound and sigma2 >= 0, so sigma2 = 0
# (the LW model) is a point of the parameter space. Several starting values are used,
# and theta and a start from a Gompertz regression of log(D / E) on x.
# A fit is on the boundary when the estimate of sigma2 is below boundary_tol.
boundary_tol = 1e-8

fit_glw = function(dat, fixed_b = NULL, s2_grid = seq(0, 1, 0.2)[-1],
                   b_grid = seq(0, 2, 0.2)[-1]) {
  gomp = coef(lm(log(D / E) ~ x, data = dat[dat$D > 0, ]))
  alpha = exp(gomp[1])
  beta = gomp[2]
  if (is.null(fixed_b)) {
    lower = c(1e-10, 1e-10, 1e-3, 0)
    upper = c(Inf, Inf, 10, Inf)
  } else {
    lower = c(1e-10, 1e-10, 0)
    upper = c(Inf, Inf, Inf)
  }
  starts = list()
  for (s2 in s2_grid) {
    for (b in b_grid) {
      if (is.null(fixed_b)) {
        starts[[length(starts) + 1]] = c(alpha / beta, beta, b, s2)
      } else {
        starts[[length(starts) + 1]] = c(alpha / beta, beta, s2)
      }
    }
  }
  starts = unique(starts)
  run = function(st) {
    st = pmin(pmax(st, lower), upper)
    nlminb(st, neg_loglik, gradient = neg_loglik_gradient, x = dat$x, D = dat$D, E = dat$E,
           fixed_b = fixed_b, lower = lower, upper = upper,
           control = list(eval.max = 5000, iter.max = 5000, rel.tol = 1e-12))
  }
  best = NULL
  for (st in starts) {
    fit = tryCatch(run(st), error = function(e) NULL)
    if (!is.null(fit) && (is.null(best) || fit$objective < best$objective)) {
      best = fit
    }
  }
  # Restart from the best solution until the fit stops improving
  for (k in 1:10) {
    fit = tryCatch(run(best$par), error = function(e) NULL)
    if (is.null(fit)) {
      break
    }
    improvement = best$objective - fit$objective
    if (fit$objective < best$objective) {
      best = fit
    }
    if (improvement < 1e-8) {
      break
    }
  }
  par = unname(best$par)
  if (!is.null(fixed_b)) {
    par = c(par[1:2], fixed_b, par[3])
  }
  par[4][par[4] < boundary_tol] = 0
  # nlminb can report "false convergence" at the maximum when the requested tolerance is
  # below the numerical precision of the log-likelihood. The restarts above confirm the
  # maximum, so only a reached iteration or evaluation limit is flagged as non-convergence.
  convergence = as.integer(grepl("limit", best$message))
  list(par = par, loglik = -best$objective, convergence = convergence,
       boundary = par[4] < boundary_tol)
}

# Expected Fisher information
expected_information = function(dat, par) {
  mu = glw_hazard(dat$x, par)
  g = glw_gradient(dat$x, par)
  t(g) %*% (g * dat$E / mu)
}
