# %%
rm(list = ls())
library(LalRUtils)
libreq(CVXR, tidyverse, data.table, ebal, tictoc, knitr)

# %% # create a toy dataset (50 controls, 30 treated)
d <- c(rep(0,50),rep(1,30))
# 3 covariates - normal, treated have higher means
X           <- rbind(replicate(3,rnorm(50, 0)),
                     replicate(3,rnorm(30, .5)))
colnames(X) <- paste("x",1:3,sep="")
df = data.table(d, X)

# %% # entropy balancing - original
eb.out <- ebalance(Treatment=df$d, X=df[, 2:4])
# %% # means in reweighted control group data match treated
rbind(
  # means in T and C groups
  df[, lapply(.SD, mean), by = d],
  # means in reweighted control group data
  data.table(d = -1, t(apply(df[d==0, 2:4], 2, weighted.mean, w=eb.out$w)))
)

# %% manual balance with treatment means
# target
b = df[d == 1, lapply(.SD, mean)][1, 2:4] |> as.numeric()
# control observations
X = df[d == 0, 2:4] |> as.matrix()
# %% manual entropy balancing
ω = Variable(df[d == 0, .N])                  # weights are of length n_0
objective <- Maximize(sum(entr(ω)))           # entropy objective fn
constraints <- list(
    ω >= 0, sum(ω) == 1,  # proper weights
    t(X) %*% ω == b       # balance
  )
prob <- Problem(objective, constraints)
result <- solve(prob, solver = "MOSEK")
ω_hat <- result$getValue(ω)
# %%
rbind(
  # means in T and C groups
  df[, lapply(.SD, mean), by = d],
  # means in reweighted control group data
  data.table(d = -1, t(apply(df[d==0, 2:4], 2, weighted.mean, w=ω_hat)))
)

# %% compare solutions - identical
cbind(eb.out$w / sum(eb.out$w), ω_hat)  |> head()

# %%
libreq(causalsens, glue, lfe)
data(lalonde.exp); data(lalonde.psid);
dtlalonde = data.table(lalonde.exp)
dtpsid = data.table(lalonde.psid)
y = 're78'
w =  'treat'
x = setdiff(colnames(dtpsid), c(y, w))

# %% ebal reg
eb_reg = function(df, y, w, ctrls){
  setorderv(df, w)
  X = df[, ..ctrls]; w = df[[w]]; y = df[[y]]
  out_eb <- ebalance(Treatment= w, X= X, print.level = -1)
  wt1 = c(out_eb$w, w[w==1])
  eff = felm(y ~ w, weights = wt1) %>% robustify()
  summary(eff)$coefficients[2,1:2]
}

# %%
eb_reg2 = function(df, y, w, x){
  setorderv(df, w)
  # treated and control obs
  ctrl  = df[eval(parse(text = glue("{w} == 0")))]
  treat = df[eval(parse(text = glue("{w} == 1")))]
  # balance constraints - means of all Xs in treated group
  b = treat[, lapply(.SD, mean), .SDcols = x]  |> as.numeric()
  X = ctrl[, ..x] |> as.matrix() # matrix of control obs
  # solve optimisation problem in CVXR
  ω = Variable(nrow(ctrl))                      # weights are of length n_0
  objective <- Maximize(sum(entr(ω)))           # entropy objective fn
  constraints <- list(
      ω >= 0, sum(ω) == 1,  # proper weights
      t(X) %*% ω     == b       # balance
  )
  prob <- Problem(objective, constraints)
  result <- solve(prob, solver = "MOSEK")
  ω_hat <- result$getValue(ω)
  # line up weights (solved weights for control units and 1/n1 for treat units)
  ω_all = c(ω_hat, rep(1/nrow(df), nrow(treat) ) )
  treatment = df[[w]]; y = df[[y]]
  eff = felm(y ~ treatment, weights = ω_all)  |> robustify()
  summary(eff)$coefficients[2,1:2]
}


# %%
rbind(
  eb_reg(dtlalonde, y, w, x),
  eb_reg2(dtlalonde, y, w, x)
)


# %%
rbind(
  eb_reg(dtpsid, y, w, x),
  eb_reg2(dtpsid, y, w, x)
)
