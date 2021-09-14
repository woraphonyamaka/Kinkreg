#### Positive and Negative funcitons

pos.part <- function(x) x*(x>0)  # Regime 1
neg.part <- function(x) x*(x<0)  # Regime 2


##################################################
Kinkreg1 = function(y,x,Call="optim"){
  LSkink <- function(par, y, x, Call) {
    K1 = (ncol(x) * 2) + 1
    k = length(par)
    n = length(y)
    r = par[c((K1 + 1):k)]
    T = length(r)
    X = matrix(0, n, T * 2)
    for (j in 1:1) {
      X[, c(j:(2 * j))] = cbind(neg.part(x[, j] - r[j]),
                                pos.part(x[, j] - r[j]))
    }
    xmat = cbind(X)
    KK1 = ncol(xmat) + 1
    alpha = par[1:KK1]
    e1 = y - (cbind(1, xmat)) %*% (alpha)
    if (Call == "optim")
      return(mean(e1^2))
    if (Call == "residuals")
      return(e1)
  }
  Xmat = cbind(x)
  r = colMeans(Xmat)
  K1 = ncol(Xmat) + ncol(Xmat)
  par = c(b0 = 0, x = rep(0.1, K1), kink = r)
  KK = length(par)
  model <- optim(par, LSkink, y = y, x = Xmat, Call = "optim",
                 control = list(maxit = 1e+05, fnscale = 1), method = "BFGS",
                 hessian = TRUE)
  coef <- model$par
  k = length(coef)
  model$se <- sqrt(-diag(solve(model$hessian)))
  for (i in 1:k) {
    if (is.nan(model$se[i]))
      model$se[i] <- sqrt(-diag(solve(-model$hessian)))[i]
  }
  S.E. = model$se
  (paramsWithTs = cbind(model$par, coef/S.E.))
  stat = coef/S.E.
  pvalue <- 2 * (1 - pnorm(abs(stat)))
  result <- cbind(coef, S.E., stat, pvalue)
  RSS = model$value
  output = list(result = result, RSS = model$value)
  output

}



#### Example Simulation data
#set.seed(111)  # Set seed for reproducibility
#k = 1                         #number of dependent variable
#n = 500                       #number of observation
#r1=1.5                        # Kink parameter
#x = matrix(rnorm(n*k,r1,2),ncol=k)
#e = rnorm(n,0,sd=1)
#x1 = cbind(neg.part(x-r1),pos.part(x-r1))

#y=0.5+(0.5*x1[,1])-(1*x1[,2])+e
#plot(x,y, col="red", lwd=2)

#Kinkreg1(y,x,Call="optim")




