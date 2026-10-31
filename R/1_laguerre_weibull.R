# Laguerre-Weibull (LW) distribution: density, survival and hazard functions
# F(x) = 1 - exp{-theta (exp(a x^b) - 1)},  mu(x) = theta a b x^(b-1) exp(a x^b)
# Set path below to the repository folder

library(ggplot2)
library(patchwork)

path = "~/Desktop/GitHub/laguerre_weibull_distribution"

# Colour scheme for the LW distribution, used throughout the manuscript
# Parameter sets (theta, a, b)
pars = data.frame(
  theta = c(0.5, 0.3, 0.2, 0.05),
  a     = c(1.0, 1.0, 0.6, 1.0),
  b     = c(0.5, 1.0, 1.5, 2.0)
)
pars$label = sprintf("theta == %s*','~a == %s*','~b == %s", pars$theta, pars$a, pars$b)

lw_hazard = function(x, theta, a, b) {
  theta * a * b * x^(b - 1) * exp(a * x^b)
}

lw_survival = function(x, theta, a, b) {
  exp(-theta * (exp(a * x^b) - 1))
}

lw_density = function(x, theta, a, b) {
  lw_hazard(x, theta, a, b) * lw_survival(x, theta, a, b)
}

x = seq(0.005, 3.5, length.out = 1000)

curves = do.call(rbind, lapply(seq_len(nrow(pars)), function(i) {
  p = pars[i, ]
  data.frame(
    x = x,
    density = lw_density(x, p$theta, p$a, p$b),
    survival = lw_survival(x, p$theta, p$a, p$b),
    hazard = lw_hazard(x, p$theta, p$a, p$b),
    set = p$label
  )
}))
curves$set = factor(curves$set, levels = pars$label)

base_theme = theme_classic(base_size = 15) +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 13),
    plot.tag = element_text(face = "bold")
  )

panel = function(variable, ylab, ylim = NULL, log_y = FALSE) {
  ggplot(curves, aes(x = x, y = .data[[variable]], colour = set)) +
    geom_line(linewidth = 1) +
    scale_colour_discrete(labels = scales::label_parse()) +
    coord_cartesian(ylim = ylim) +
    guides(colour = guide_legend(nrow = 1, byrow = TRUE)) +
    labs(x = "x", y = ylab) +
    (if (log_y) scale_y_log10(labels = scales::label_log()) else NULL) +
    base_theme
}

p_density = panel("density", "Density", ylim = c(0, 2))
p_survival = panel("survival", "Survival function")
p_hazard = panel("hazard", "Hazard function (log scale)", log_y = TRUE)

fig = (p_density | p_survival | p_hazard) +
  plot_layout(guides = "collect") +
  plot_annotation(title = "Laguerre-Weibull distribution", tag_levels = "A") &
  theme(legend.position = "bottom")

# Save to the R figures folder and to the LaTeX figures folder used by the manuscript
for (out_dir in file.path(path, c("R/figures", "LaTeX/figures"))) {
  ggsave(file.path(out_dir, "laguerre_weibull.pdf"), fig, width = 11, height = 5)
}
