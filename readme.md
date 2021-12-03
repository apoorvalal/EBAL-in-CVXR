# `CVXRebal`: Entropy balancing in CVXR

estimation functions in `R/CVXRebal.R`. The `CVXR:solve` call
currently uses the commercial MOSEK solver, is free for academics for
[personal use](https://www.mosek.com/license/request/). Otherwise,
[other
solvers](https://cvxr.rbind.io/cvxr_examples/cvxr_using-other-solvers/)
can be used by passing it as the `solver` argument in the
`ebalRegAdjust` function call.

```{r}
rm(list = ls())
library(data.table)
library(causalsens)

source("R/CVXRebal.R")

# lalonde PSID sample
data(lalonde.psid)
y = 're78'; w =  'treat'
x = setdiff(colnames(lalonde.psid), c(y, w))
ebalRegAdjust(lalonde.psid, y, w, x) |> summary()
```

```
Call:
   felm(formula = y ~ treatment, weights = ω_all)

Weighted Residuals:
   Min     1Q Median     3Q    Max
 -1104      3     16     49   4116

Coefficients:
            Estimate Std. Error t value     Pr(>|t|)
(Intercept)     3924        660    5.95 0.0000000031 ***
treatment       2425        876    2.77       0.0057 **
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 144 on 2673 degrees of freedom
Multiple R-squared(full model): 0.00682   Adjusted R-squared: 0.00645
Multiple R-squared(proj model): 0.00682   Adjusted R-squared: 0.00645
F-statistic(full model):18.4 on 1 and 2673 DF, p-value: 0.0000189
F-statistic(proj model): 18.4 on 1 and 2673 DF, p-value: 0.0000189
```
