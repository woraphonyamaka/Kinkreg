# This is an  function named 'Copula based Stochastic frontier'
# Written by Dr.Worphon Yamaka,
# Center of excellence in Econometrics, Faculty of Economics,
# Chiang Mai University
#### Positive and Negative funcitons

pos.part <- function(x) x*(x>0)  # Regime 1
neg.part <- function(x) x*(x<0)  # Regime 2


##################################################3
Kinkreg = function(y,x,z,Call="optim"){
LSkink <- function(par,y,x,z,Call){
<<<<<<< HEAD
  k=length(par)
=======
 k=length(par)
>>>>>>> 200382cd9d05316d9fbd1663d78217b517425c1a
  K=ncol(x)
  T=length(r)
  U=K-T
  n=length(y)
  r=par[c((k-T+1):k)]
  exo=z
  n <- length(y)
  X=matrix(0,n,T*2)
  for (j in 1:T){
    X[,c(((2*j)-1):(2*j))] = cbind(neg.part(x[,j]-r[j]),pos.part(x[,j]-r[j]))
  }
<<<<<<< HEAD

  xmat=cbind(X,exo)

=======
xmat=cbind(X,exo)
>>>>>>> 200382cd9d05316d9fbd1663d78217b517425c1a
K1=ncol(xmat)+1
alpha=par[1:K1]
e1 = y - (cbind(1,xmat))%*%(alpha)

 if (Call=="optim")
return(sum(e1^2))

  if (Call=="residuals")
  return(e1)

  }
##=================================
Xmat=cbind(x)
Zmat=cbind(z)
r=colMeans(Xmat)
K1=ncol(Xmat)+ncol(Xmat)
K2=ncol(Zmat)
par=c(b0=1,x=rep(0.1,K1),z=rep(0.1,K2), kink=r)
KK=length(par)
model<- optim(par,LSkink,y=y,x=Xmat,z=Zmat,Call="optim",control = list(maxit=100000,fnscale=1)
          ,method="BFGS",hessian=TRUE)


# table of results
coef<- model$par
k=length(coef)
model$se <- sqrt(-diag(solve(model$hessian)))

for(i in 1:k){
  if (is.nan(model$se[i]))  # control for optimization
    model$se[i] <- sqrt(-diag(solve(-model$hessian)))[i]
}


S.E.= model$se
(paramsWithTs = cbind (model$par , coef/S.E. ) )
stat=coef/S.E.
pvalue <- 2*(1 - pnorm(abs(stat)))
result <- cbind(coef,S.E.,stat,pvalue)
RSS= model$value

output=list(
  result=result,
  RSS=model$value
)
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
#z = cbind(rnorm(n,0,sd=1))
#y=0.5+(0.5*x1[,1])-(1*x1[,2])+0.5*z+e
#plot(x,y, col="red", lwd=2)

#Kinkreg(y,x,z,Call="optim")




