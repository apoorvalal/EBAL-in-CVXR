#' function to solve for entropy balancing weights
#' @param b p-vector of treatment means
#' @param X n_0 X p matrix of control observations
#' @keywords entropy balancing primal
#' @export
ebalSolve = function(b, X, solver = "MOSEK"){
  require(CVXR)
  #####################################################################
  # initialise problem weights are of length n_0
  ω = Variable(nrow(X))
  # entropy objective fn - uses entr() atom
  objective = Maximize(sum(entr(ω)))
  constraints = list(
      t(X) %*% ω     == b,                   # balance
      ω >= 0, sum(ω) == 1                    # proper weights
  )
  # solve using commercial solver MOSEK by default - use another solver if unavailable
  result = solve(Problem(objective, constraints), solver = solver)
  ω_hat = result$getValue(ω)
  ω_hat
}

# %%
#' Regression adjustment using entropy balancing weights
#' @param df data.table with treatment, outcome, and covariates
#' @param y string with outcome name
#' @param w string with treatment name
#' @param x character vector with control names
#' @export
ebalRegAdjust = function(df, y, w, x){
  require(data.table); require(glue); require(lfe)
  # move control obs to the top
  setDT(df)
  setorderv(df, w)
  # treated and control obs in separate dataframes
  ctrl  = df[eval(parse(text = glue::glue("{w} == 0")))]
  treat = df[eval(parse(text = glue::glue("{w} == 1")))]
  # vector of treated means for balance constraints
  b = treat[, lapply(.SD, mean), .SDcols = x]  |> as.numeric()
  # matrix of control obs
  X = ctrl[, ..x] |> as.matrix()
  # solve optimisation problem in CVXR
  ω_hat = ebalSolve(b, X)
  # line up weights (solved weights for control units and 1/n for treat units)
  treatment = df[[w]]; y = df[[y]]
  ω_all = c(ω_hat, rep(1/nrow(df), sum(treatment) ) )
  # fit linear model and return LFE object
  eff = felm(y ~ treatment, weights = ω_all) |> robustify()
}

# %%
