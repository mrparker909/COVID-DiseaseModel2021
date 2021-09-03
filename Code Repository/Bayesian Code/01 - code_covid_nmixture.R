##########
# Disease Analytics (COVID19) Model: alpha = 1
#
#######

#############################################
# Northern Health Authority
# load data:
dat = read.csv("bccdc_NHA_weekly_non_cumulative.csv")
n = dat$New.Cases
r = dat$New.Recoveries
D = dat$New.Deaths

# fix leading recoveries and deaths to 0, as they belong to previous data
r[1] = 0
D[1] = 0

a<-c(11,15,17,13,15,11,13,11,9,8,6,6,6,6,6,6,6,16,18,24,22,21,25,50,38,39,44,41,22)

sink(file="northdata1.txt")

### The model
nTimes=29
nSites=1

cat("
model {
  #prior 
  lambda ~ dgamma(0.5, 0.05)
  gamma  ~ dunif(0, 30)
  omega  ~ dunif(0, 5)
  p      ~ dunif(0, 1)
  pdeath ~ dunif(0, 1)
  prec   ~ dunif(0, 1-pdeath)

  # model
  N[1] ~ dpois(lambda)
  n[1] ~ dbin(p, N[1])
  S[1] ~ dpois(0)
  G[1] ~ dpois(0)
  R[1] ~ dpois(0)

  for(t in 1:nTimes) {
      D[t+1] ~ dbin(pdeath,N[t])
      R[t+1] ~ dbin(prec/(1-pdeath), N[t]-D[t+1])
      A[t+1] <- N[t]-D[t+1]-R[t+1]
      
      G[t+1] ~ dpois(gamma)        # new importations
      S[t+1] ~ dpois(omega*(N[t])) # domestic spread
     
      N[t+1]<-A[t+1]+G[t+1]+S[t+1]
      
      n[t+1]~ dbin(p, N[t+1]-a[t]+r[t]+D[t])
      r[t+1]~ dbin(prec, a[t])
  }
}
", fill=TRUE)
sink()


library(rjags)


# Bundle data
dat.const <- list(nTimes=nTimes, n=n, D=D, r=r, a=a)

# Initial values
N <- rep(280000,30) # Initial values for total cases N
R <- rep(28000,30)  # Initial values for total recoveries R
G <- rep(28000,30)  # Initial values for importation
S <- N-G            # Initial values for domestic spread
N[-1] <- NA

pdeath <- runif(1, 0, 1)
prec   <- runif(1, 0, 1-pdeath)

# Bundle inital values
init.const <- function() list(
  lambda=runif(1, 0, 30),
  gamma=runif(1, 0, 30),
  N=N,
  omega=runif(1, 0, 5),
  G=G,
  S=S,
  R=R,
  p=runif(1, 0, 1),
  pdeath=pdeath,
  prec=prec)
pars.const <- c("lambda", "gamma", "omega", "p", "pdeath", "prec")

# Compile model
jm.const10 <- jags.model("northdata1.txt", dat.const, init.const,
                         n.chains=1,n.adapt=1200000)

# Posterior samples
update(jm.const10, n.iter=1200000)
ps.const1 <- coda.samples(jm.const10, pars.const, n.iter=1400000,thin=200)

summary(ps.const1)
saveRDS(ps.const1, file = "bayes_results1.RDS")
