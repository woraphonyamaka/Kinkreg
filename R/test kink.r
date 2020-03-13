#### Kink test
# Ho: Linear
# Ha: Kink

kinktest=function(y,x,level=0.90, boot=100,search=0.01) {
# Controls
Mean=mean(x)
lower=range(x)[1]+abs(Mean)
upper=range(x)[2]-abs(Mean)
gammas = seq(lower,upper,by=search) # Grid on Threshold parameter for estimation specify follow the data x
dx = seq(min(x),max(x),by=0.1)              # Grid on regression function for display
level = level			              # For confidence sets
boot = boot   			              # Number of bootstrap replications
Ceps = c(0.5,1,2,4)		              # For numerical delta method bootstrap
n=length(y)
############################################

# Some useful functions
# expectile regression using least square method
pos.part <- function(x) x*(x>0)  # Regime 1
neg.part <- function(x) x*(x<0)  # Regime 2

m0=lm(y~x)
e0=y-fitted(m0)
sse0 = sum(e0^2)

kink <- function(theta,y,x) {
r1=theta[1]
x11 = cbind(neg.part(x-r1),pos.part(x-r1))

exp1 <- lm(y~x11[,1]+x11[,2])
res=y-fitted(exp1)
e2=sum(res^2)
e2}
par=mean(x)
lower =-Inf
upper =Inf
fit1 <- optim(par,kink ,y=y,x=x,
          control = list(maxit=100000,fnscale=1),method="L-BFGS-B",
           lower =lower,upper =upper, hessian=TRUE )


x11 = cbind(neg.part(x-fit1$par),pos.part(x-fit1$par))
XX=as.matrix(cbind(1,x11))
model<- lm(y~x11[,1]+x11[,2])

thres=cbind(fit1$par,sqrt(diag(solve(fit1$hessian))))


#---------- Threshold Model Grid Search----------------------
  grid = length(gammas)  # specify grid search
  rd = length(dx)
  sse = matrix(0,grid,1) # store error term to find lowest error
  k = 4
  # find threshold value
  for (j in 1:grid) {
  gamj=gammas[j]
  x1 = cbind(neg.part(x-gamj),pos.part(x-gamj))
  exp1 <- lm(y~x1[,1]+x1[,2])
  e1=y-fitted(exp1)
  sse[j] = sum(e1^2)
  }
  gi = which.min(sse)              # find location of sum square error
  gammahat = gammas[gi]            # define of threshold value
  ssemin = sse[gi]                 # minimum sum square of error.

#-----------------------------------------
# setup the kink regression
# Regression Estimate
#----------------------------
  x1 = cbind(neg.part(x-gammahat),pos.part(x-gammahat)) # x
  bt = coef(model)
  et = y - cbind(1,x1)%*% bt
  hg = cbind(neg.part(x-gammahat),pos.part(x-gammahat))
  hg1 = summary(lm(y~hg[,1]+hg[,2]))
  hg2 = hg1[1:3]
  betahat = cbind(bt,gammahat)  # estimated result for Kink regression
  wt = n*(sse0-ssemin)/ssemin
  wg = n*(sse-ssemin)/ssemin

# Plot Grid Search
#plot(gammas,sse,type="l",ylab="Least Squares Criterion",xlab="kink Parameter")



# predict kink regression
G = cbind(1,neg.part(dx-gammahat),pos.part(dx-gammahat))
yf = G%*%bt
#The parameter estimates from this fitted regression kink model
bt # estimated parameter
#yf # expected Yhat




  m21 <- lm(y~x)
  e0 <- y-fitted(m21)
# Bootstrap & Testing
  waldb = matrix(0,grid,boot)
  sseb  = matrix(0,grid,boot)
  betab = array(0,c(grid,k-1,boot))
  u = matrix(rnorm(n*boot),n,boot)
  eb = matrix(e0,n,boot)*u
  yb = matrix(cbind(1,x1)%*%bt,n,boot) + matrix(et,n,boot)*u
  h2= matrix(0,2,boot)
  for (i in 1:boot) {
  h1 = coef(lm(eb[,i]~x))
  h2[,i] = h1[1:2]
   }
  hg2=matrix(0,3,boot)
  bb=matrix(0,3,boot)
  x0=cbind(1,x)
  eb0 = eb - (x0%*%h2)
  bsse0 = colSums(eb0^2)

 # create progress bar
pb <- winProgressBar(title = "Bootrapping %", min = 0,
                     max = grid, width = 300)


for (j in 1:grid) {
    gamj = gammas[j]
    x2 = cbind(neg.part(x-gamj),pos.part(x-gamj))
    x3=cbind(1,x2)
   for (i in 1:boot) {
    hg1 = coef(lm(eb[,i]~x2[,1]+x2[,2]))
    hg2[,i] = hg1[1:3]
    }
    eb0 = eb - (x3%*%hg2)
    bsse = colSums(eb0^2)
    waldb[j,] = n*(bsse0-bsse)/bsse
    for (i in 1:boot) {
    b0 = coef(lm(yb[,i]~x2[,1]+x2[,2]))
    bb[,i] = b0[1:3]
    }
    eb1 = yb - x3%*%bb
    sseb[j,] = colSums(eb1^2)
    betab[j,,] = bb
   Sys.sleep(0.1)
   setWinProgressBar(pb, j, title=paste( round(j/grid*100, 0),
                                        "% done"))
  }
close(pb)
# Multiplier Bootstrap test for Threshold
  wb = apply(waldb,2,max)
  pv = mean(wb > matrix(wt,boot,1))
  crit = quantile(wb,probs=level)
windows()
plot(gammas,sse,type="l",main="Grid Search",ylab="Least Squares Criterion",xlab="kink Parameter")
windows()
plot(x,y, main= "Kink fitted line")
lines(dx,yf, col="red", lwd=3)
output=list(
  Wald.test=crit,pvalue=pv
)
output

}



#### Example Simulation data
#set.seed(111)  # Set seed for reproducibility
#k = 1                         #number of dependent variable
#n = 500                       #number of observation
#r1=1.5                        # Kink parameter
#x = rnorm(n,r1,sd=1)
#e = rnorm(n,0,sd=1)
#x1 = cbind(neg.part(x-r1),pos.part(x-r1))
#y=0.5+(0.5*x1[,1])-(1*x1[,2])+e

# Ho: Linear
# Ha: Kink
#kinktest(y,x,level=0.90, boot=10,search=0.01)
