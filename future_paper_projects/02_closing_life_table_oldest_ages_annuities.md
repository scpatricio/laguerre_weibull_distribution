# Closing the life table at the oldest ages: out-of-sample validation of old-age mortality models and the cost for deferred annuity pricing

Suggested repository name: `closing-life-table-oldest-ages`

## Research question

Which method of extending the life table above the ages with reliable data predicts old-age mortality best out of sample, and how much does the choice change the values of annuities that depend on survival at the oldest ages?

## Pipeline

1. **Data.** HMD cohort data for extinct (or almost extinct) cohorts, several countries; optionally the International Database on Longevity (IDL) for ages above 105.
   Why: extinct cohorts give the true survival at the oldest ages, so predictions can be checked against what actually happened.

2. **Validation design.** For each cohort, fit the models on ages up to a cut-off (for example 90 or 95) and predict mortality and survival above it.
   Why: an out-of-sample test is what this literature lacks; in-sample AIC comparisons say little about extrapolation.

3. **Candidate methods.**
   - Kannisto, gamma-Gompertz, GLW, Coale–Kisker and the HMD method.
   - An extreme-value approach, a generalized Pareto distribution above a threshold (Gbari et al. 2017).
   Why: covers the plateau, no-plateau and finite-endpoint views, and the methods used in practice.

4. **Prediction accuracy.** Compare predicted and observed death rates, survival probabilities and remaining life expectancy above the cut-off, with proper scoring rules where possible.
   Why: gives a ranking based on prediction, not fit.

5. **Actuarial cost.** For each method, compute deferred annuity values (for example payments starting at 85 or 90), remaining life expectancy and whole life annuities at 80 and 90, at several interest rates. Compare with the values implied by the observed cohort survival.
   Why: translates prediction errors into money, which is what the actuary needs.

6. **Prudence.** Check which methods are conservative (overstate survival) and by how much.
   Why: links to the regulatory and prudential view that a constant ultimate force of mortality is conservative for annuities.

7. **Write-up and code.** Reproducible pipeline from HMD cohort files to tables of annuity errors by method, country and cut-off.

## Why an editor would take it

- Closing the table is a recurring practical problem in annuity and pension work; the CMI has working papers on high-age mortality.
- The out-of-sample design with extinct cohorts gives a clean answer instead of another in-sample comparison.
- The cost is expressed in annuity values, which makes the result directly usable.

## Target journals (in order)

1. North American Actuarial Journal
2. Annals of Actuarial Science
3. Scandinavian Actuarial Journal
4. European Actuarial Journal
5. British Actuarial Journal
