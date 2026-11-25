# A note on the Makeham term and gamma frailty in annuity values: interest-rate equivalence and closed-form elasticities under the gamma-Gompertz–Makeham law

Suggested repository name: `makeham-frailty-annuity-elasticities`

## Research question

How do the Makeham term and the frailty variance each change the values of life annuities and assurances under the gamma-Gompertz–Makeham law, as functions of issue age and force of interest?

## Pipeline

1. **Interest-rate equivalence.** Show that, under Gompertz–Makeham and gamma-Gompertz–Makeham, the Makeham factor exp(−ct) combines with the discount factor, so the annuity at force of interest δ equals the annuity without the Makeham term at δ + c, and derive the corresponding results for assurances.
   Why: gives a simple and exact way to read the effect of the Makeham term as an interest-rate shift.

2. **Closed forms.** Write the annuity and assurance values in closed form (hypergeometric functions, Castellares et al. 2024; Makeham–Beard, Bowie 2021).
   Why: the elasticities below need explicit expressions.

3. **Elasticities.** Derive the elasticities of ā_x and Ā_x with respect to c, σ², a and b, and study their sign and monotonicity in issue age and δ (for example whether the frailty effect falls as δ rises and grows with issue age).
   Why: formal statements about where each component matters are what makes a note publishable.

4. **Numerical illustration.** Evaluate the elasticities with parameter values estimated from HMD data (for example France, Italy, Japan and Sweden, as in Castellares et al. 2024) for annuities at 50, 65 and 80 and term assurances at working ages.
   Why: shows that the effects are of practical size.

5. **Connection to the mixture representation.** Use the share of premature deaths conditional on survival to the issue age (Patricio and Missov 2024) to explain why the Makeham term matters for assurances at working ages and frailty for annuities at old ages.
   Why: a short, intuitive explanation of the formal results.

6. **Write-up.** A short note (around 12 to 15 pages).

## Why an editor would take it

- Short, exact and immediately usable, in the style of Bowie (2021).
- Gives practitioners a rule of thumb: the Makeham term is an interest-rate shift, frailty is an old-age effect, with closed-form sizes for both.
- Extends the life-expectancy elasticities of Missov and Lenart (2013) to discounted actuarial values.

## Target journals (in order)

1. Annals of Actuarial Science
2. European Actuarial Journal
3. Scandinavian Actuarial Journal
4. North American Actuarial Journal
