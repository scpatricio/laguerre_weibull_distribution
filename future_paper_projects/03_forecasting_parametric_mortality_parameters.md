# Forecasting the parameters of parametric mortality laws: interpretable mortality projections and their accuracy for annuity pricing

Suggested repository name: `forecasting-parametric-mortality-laws`

## Research question

Can forecasts of the parameters of a parametric mortality law (for example gamma-Gompertz–Makeham or GLW) give annuity prices as accurate as standard stochastic mortality models and machine-learning forecasts, while keeping parameters that actuaries can interpret?

## Pipeline

1. **Data.** HMD period death counts and exposures, ages 50 and above, several countries, both sexes, 1950 to the latest year.
   Why: long series are needed for both fitting and backtesting.

2. **Parameter series.** Fit the parametric law to every calendar year (Poisson likelihood, or Bayesian in Stan) to obtain time series of parameters, for example the level, the rate of ageing, the frailty variance and the shape parameter b.
   Why: the parameters summarise the mortality schedule in a few interpretable numbers.

3. **Parameter forecasting.** Model the parameter series jointly, for example with a multivariate random walk with drift, a vector autoregression, or a Bayesian dynamic model in Stan, with constraints that keep the parameters valid.
   Why: joint models keep the correlation between parameters, which matters for the forecast survival curves.

4. **Benchmarks.** Lee–Carter, Cairns–Blake–Dowd (CBD) and one machine-learning forecast (for example a neural network or LocalGLMnet-type model).
   Why: an editor will expect a comparison with the standard models and with the machine-learning methods now in the journals.

5. **Backtesting.** Rolling-origin out-of-sample forecasts (for example fit to 1950–1999, forecast 2000–2019) of death rates, life expectancy at 65 and cohort annuity values.
   Why: the comparison must be on prediction, and annuity values are the quantity that matters for actuaries.

6. **Prediction intervals.** Compare the coverage of the prediction intervals of all methods.
   Why: uncertainty is as important as the point forecast for reserving and capital.

7. **Interpretation.** Show how the forecast parameters (for example a rising b or a falling frailty variance) explain the projected changes in old-age mortality.
   Why: interpretability is the advantage over black-box methods, so it must be shown explicitly.

8. **Write-up and code.**

## Why an editor would take it

- Forecasting is the dominant mortality topic in actuarial journals.
- Interpretable parametric forecasts compared with black-box methods is a recognisable and debated angle.
- The evaluation is in annuity values, not only in death rates.

## Target journals (in order)

1. Scandinavian Actuarial Journal
2. Annals of Actuarial Science
3. North American Actuarial Journal
4. Insurance: Mathematics and Economics
5. European Actuarial Journal
