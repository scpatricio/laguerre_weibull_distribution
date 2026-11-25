# Model risk from the choice of mortality law in annuity and assurance valuation: Makeham term, frailty and the old-age plateau

Suggested repository name: `mortality-law-model-risk-annuities`

## Research question

How much do life annuity, deferred annuity and assurance values change with the components of the mortality law (Makeham term, gamma frailty, plateau or no plateau at old ages), and how should an actuary account for this model risk alongside parameter uncertainty?

## Pipeline

1. **Data.** HMD death counts and exposures (single ages and years), from age 30 and from age 50, several countries, both sexes, selected years (for example every 10 years from 1950).
   Why: starting at 30 identifies the Makeham term; starting at 50 matches pension and annuity ages. Comparing the two shows how the starting age changes the conclusions.

2. **Candidate models.** Gompertz, Gompertz–Makeham, gamma-Gompertz (Beard), gamma-Gompertz–Makeham, GLW, and optionally Kannisto/Perks.
   Why: the set is nested by components, so each difference in annuity values can be traced to one component (Makeham term, frailty, shape parameter b).

3. **Bayesian estimation in Stan** with a Poisson likelihood for each model and population-year; check convergence and posterior predictive fit.
   Why: the posterior carries parameter uncertainty, including the trade-off between the Makeham term c and the frailty variance, directly into the actuarial values.

4. **Actuarial values from the posterior.** For each posterior draw, compute the remaining life expectancy, the whole life annuity, deferred annuities (for example a 20-year deferred annuity at 65, payments from 85), the whole life assurance and a term assurance at working ages (for example 40), at several forces of interest (for example 0%, 2%, 4%).
   Why: deferred annuities depend most on old-age survival, term assurances on working-age mortality, so the two product types expose different components.

5. **Decomposition of uncertainty.** Split the total variance of each actuarial value into parameter uncertainty (within model) and model uncertainty (between models). Report it by issue age, deferral period, product and interest rate.
   Why: this is the central result: where the choice of law matters more than the estimation error, and where it does not.

6. **Model averaging.** Combine the models by Bayesian stacking or model averaging and compare the averaged actuarial values with those of the single best model.
   Why: gives actuaries a practical way to price under model risk instead of picking one law.

7. **Interpretation of where each component matters.**
   - The Makeham term acts as an increase in the force of interest in annuity values (ā_x at δ equals the Gompertz annuity at δ + c), so it matters little for annuities and more for assurances at working ages.
   - Frailty matters most for deferred annuities and old issue ages.
   - Use the mixture representation (Patricio and Missov 2024, Demographic Research) to show the share of premature deaths conditional on survival to the issue age.
   Why: turns the numbers into a rule an actuary can remember.

8. **Old-age plateau check.** Refit on ages 50–90 and 50–109 and compare the shape parameter b and the deferred-annuity values.
   Why: shows whether the plateau question is decided by the data at the oldest ages or extrapolated from middle ages.

9. **Analytical support.** Closed forms for the gamma-Gompertz(–Makeham) annuity (Castellares et al. 2024) and elasticities of ā_x with respect to c and σ² as functions of age and δ.
   Why: formal results that explain the simulation and numerical findings.

10. **Write-up and code.** Reproducible R and Stan code, figures of uncertainty decomposition by age and product.

## Why an editor would take it

- It answers a pricing and reserving question: how much the choice of mortality law matters for liabilities, and how to allow for it.
- It sits in the current model-risk and mis-estimation literature (Richards, British Actuarial Journal) and the current interest in deferred and longevity annuities.
- It combines the classical closed-form results with Bayesian model averaging, which is a modern and accepted approach to model selection in actuarial journals.
- It gives practical guidance on which component matters for which product, age and interest rate.

## Target journals (in order)

1. Scandinavian Actuarial Journal
2. Annals of Actuarial Science
3. European Actuarial Journal
4. North American Actuarial Journal
5. Insurance: Mathematics and Economics (only if formal results and portfolio or capital consequences are added)
