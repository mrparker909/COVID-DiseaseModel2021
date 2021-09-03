library(doParallel) # parallelization
library(foreach)    # parallelization
library(pracma)     # hessian() function
source("loglikelihood.R") # the likelihood function

# Setup 20 processor cluster
cluster <- makeCluster(20)
registerDoParallel(cluster)
foreach::getDoParWorkers()

calcAIC <- function(nll, par) {
  return(2*nll+2*length(par))
}

# k = num par, n = sampling occasions
calcAICc <- function(nll, k, n) {
  return(2*(nll+k) + 2*(k^2 + k)/(n-k-1))
}

#############################################
# Northern Health Authority
# load data:
dat = read.csv("bccdc_NHA_weekly_non_cumulative.csv")
nit = matrix(dat$New.Cases, nrow=1)
rit = matrix(dat$New.Recoveries, nrow=1)
Dit = matrix(dat$New.Deaths, nrow=1)

# fix leading recoveries and deaths to 0, as they belong to previous data
rit[,1] = 0
Dit[,1] = 0

# setting up covariates:
#############################################
# 4 phases: omega
#############################################
phase1 = as.numeric(dat$Phase==1) # phase 1
phase2 = as.numeric(dat$Phase==2) # phase 2
phase3 = as.numeric(dat$Phase==3) # phase 3a
phase4 = as.numeric(dat$Phase==4) # phase 3b

phase_cov = list(p2=phase2, p3=phase3, p4=phase4) # note that baseline B0 is phase1


#############################################
# mobility: omega = B0
#############################################
mobility1 = as.numeric(dat$retail_and_recreation_weekly_mean_percent_change_from_baseline)/100
mobility2 = as.numeric(dat$grocery_and_pharmacy_weekly_mean_percent_change_from_baseline)/100
mobility3 = as.numeric(dat$parks_weekly_mean_percent_change_from_baseline)/100
mobility4 = as.numeric(dat$transit_stations_weekly_mean_percent_change_from_baseline)/100
mobility5 = as.numeric(dat$workplaces_weekly_mean_percent_change_from_baseline)/100
mobility6 = as.numeric(dat$residential_weekly_mean_percent_change_from_baseline)/100

mobility_cov = list(mob_retail=mobility1,
                    mob_grocery=mobility2,
                    mob_parks=mobility3,
                    mob_transit=mobility4,
                    mob_work=mobility5,
                    mob_residential=mobility6)

#############################################
# test volume: pdet
#############################################
testVolume = list(testVolume = as.numeric(scale(dat$New.Tests..weekly.)))


#############################################
# Model1: no covariates
#############################################
# par := starting parameter values (lambda, gamma, omega, pdet, pmor, prec)
#par = c(1, -1, -1, 0, -1, -1)
par = c(4.00578029, -19.17485603, 0.03656662, -1.00718013, -6.06948179, -0.07205155)

# K := upper bound on summations
K = matrix(200, nrow=nrow(nit), ncol=ncol(nit))

# check for finite loglik at initial values:
loglik(par = par, nit = nit, K = K, Dit = Dit, rit = rit)
# [1] 324.3285
#############################################


##
# Testing: K=200
##
K = matrix(200, nrow=nrow(nit), ncol=ncol(nit))
tictoc::tic()
op1 = optim(par = par,
            fn  = loglik, method = "BFGS", hessian = F,
            nit = nit,
            Dit = Dit, 
            rit = rit, 
            K   = K, 
            HPDfile = NULL,
            control = list(trace=1, REPORT=1, reltol=1e-10))
tictoc::toc()
# 1174.95 sec elapsed
saveRDS(op1, "results/mod1.RDS")
# > op1
# $par
# [1] 3.96235656 -19.17485601   0.01473496  -0.83720519  -6.03339132  -0.06979559
# 
# $value
# [1] 298.0685
# 
# $counts
# function gradient 
# 36        9 
# 
# $convergence
# [1] 0
# 
# $message
# NULL


calcAIC(op1$value, op1$par)
# > calcAIC(op1$value, op1$par)
# [1] 608.1371

# Hessian:
tictoc::tic()
op1_H = pracma::hessian(f = loglik, 
                        x0 = op1$par,
                        nit = nit,
                        Dit = Dit, 
                        rit = rit, 
                        K   = K)
tictoc::toc()
# 668.36 sec elapsed
saveRDS(op1_H, "results/mod1_hessian.RDS")
# > op1_H
# [,1]          [,2]          [,3]       [,4]          [,5]       [,6]
# [1,] 4.389144e+01  2.861023e-06  8.921143e+00  4.7113543  1.167356e+00   9.507561
# [2,] 2.861023e-06  7.629395e-06 -9.536743e-07  0.0000000 -9.536743e-07   0.000000
# [3,] 8.921143e+00 -9.536743e-07  1.179393e+03 45.3710108 -1.198456e+00  -9.760872
# [4,] 4.711354e+00  0.000000e+00  4.537101e+01 20.4982185 -6.329145e-01  -5.155686
# [5,] 1.167356e+00 -9.536743e-07 -1.198456e+00 -0.6329145  7.841908e+00   2.591887
# [6,] 9.507561e+00  0.000000e+00 -9.760872e+00 -5.1556864  2.591887e+00 571.558308

sd_err1 = sqrt(abs(diag(solve(op1_H))))
# [1] 0.15357335 362.03867983   0.03044454   0.23440207   0.35865004   0.04199539

# exp(op1$par[1])
# 95% CI: exp(op1$par[1]-1.96*sd_err1[1]), exp(op1$par[1]+1.96*sd_err1[1]) 
# lambda: 52.58109 (38.91395, 71.04833)
#
# exp(op1$par[2])
# 95% CI: exp(op1$par[2]-1.96*sd_err1[2]), exp(op1$par[2]+1.96*sd_err1[2])
# gamma: 4.703984e-09 (0, Inf) <- border estimate, CI not useable.
#
# exp(op1$par[3])
# 95% CI: exp(op1$par[3]-1.96*sd_err1[3]), exp(op1$par[3]+1.96*sd_err1[3])
# omega: 1.014844 (0.9560583, 1.077244)
#
# plogis(op1$par[4])
# 95% CI: plogis(op1$par[4]-1.96*sd_err1[4]), plogis(op1$par[4]+1.96*sd_err1[4])
# pdet:  0.3021237 (0.2147322, 0.4066631)
#
# plogis(op1$par[5])
# 95% CI: plogis(op1$par[5]-1.96*sd_err1[5]), plogis(op1$par[5]+1.96*sd_err1[5])
# pmort: 0.002391616 (0.00118557, 0.004818621)
#
# plogis(op1$par[6])
# 95% CI: plogis(op1$par[6]-1.96*sd_err1[6]), plogis(op1$par[6]+1.96*sd_err1[6])
# prec:  0.4825582 (0.4620465, 0.5031288)




#############################################
# Model2: pdet ~ test volume
#############################################

# par := starting parameter values (lambda, gamma, omega, 
#                                   pdet, pdet_testVolume, 
#                                   pmor, prec)
par = c(4.17025131, -19.17485600, 0.01003416,
        -1.38620295, 0.75147471,
        -6.21141387, -0.08374649)
          
# check for finite loglik at initial values:
K = matrix(200, nrow=nrow(nit), ncol=ncol(nit))
loglik(par = par, nit = nit, K = K, Dit = Dit, rit = rit, p_t_c = testVolume)
# [1] 296.5104

##
# Testing: K=200
##
K = matrix(200, nrow=nrow(nit), ncol=ncol(nit))
tictoc::tic()
op2 = optim(par = par,
            fn  = loglik, method = "BFGS", hessian = F,
            nit = nit,
            Dit = Dit, 
            rit = rit, 
            K   = K, 
            p_t_c = testVolume,
            HPDfile = NULL,
            control = list(trace=1, REPORT=1, reltol=1e-10))
tictoc::toc()
# 17412.95 sec elapsed
saveRDS(op2, "results/mod2.RDS")
#
# > op2
# $par
# [1]   3.97792960 -19.17485595  -0.02525427  -0.65758255   0.92695731  -6.04646873  -0.07122534
# 
# $value
# [1] 290.9164
# 
# $counts
# function gradient 
# 45       12 
# 
# $convergence
# [1] 0
# 
# $message
# NULL


calcAIC(op2$value, op2$par)
# > calcAIC(op2$value, op2$par)
# [1] 595.8327

# Hessian: (NOT RUN)
tictoc::tic()
op2_H = pracma::hessian(f = loglik, 
                        x0 = op2$par,
                        nit = nit,
                        Dit = Dit, 
                        rit = rit, 
                        p_t_c = testVolume,
                        K   = K)
tictoc::toc()
#  sec elapsed
# saveRDS(op2_H, "results/mod2_hessian.RDS")
# > op2_H

# sd_err2 = sqrt(abs(diag(solve(op2_H))))





#############################################
# Model3a: omega ~ phases
#############################################

# par := starting parameter values (lambda, gamma, 
#                                   omega_phase1, omega_phase2, omega_phase3, omega_phase4
#                                   pdet, pmor, prec)
par = c(3.9, -19.17485603, 
        -0.06828611, -0.23926542, 0.32845523, -0.04211158,
        -0.83720519, -6.03339132, -0.06979559)

# check for finite loglik at initial values:
K = matrix(200, nrow=nrow(nit), ncol=ncol(nit))
loglik(par = par, nit = nit, K = K, Dit = Dit, rit = rit, o_t_c = phase_cov)
# [1] 277.9212

##
# Testing: K=200
##
K = matrix(200, nrow=nrow(nit), ncol=ncol(nit))
tictoc::tic()
op3 = optim(par = par,
            fn  = loglik, method = "BFGS", hessian = F,
            nit = nit,
            Dit = Dit, 
            rit = rit, 
            K   = K, 
            o_t_c = phase_cov,
            HPDfile = NULL,
            control = list(trace=1, REPORT=1, reltol=1e-10))
tictoc::toc()
# 27818.55 sec elapsed
saveRDS(op3, "results/mod3.RDS")
#
# > op3
# $par
# [1]   4.09559076 -19.17485600  -0.11257949  -0.21033570   0.39102319  -0.03625018  -1.27514369  -6.14672935  -0.08059150
# 
# $value
# [1] 274.1925
# 
# $counts
# function gradient 
# 44       15 
# 
# $convergence
# [1] 0
# 
# $message
# NULL



calcAIC(op3$value, op3$par)
# > calcAIC(op3$value, op3$par)
# [1] 

# Hessian: (NOT RUN)
tictoc::tic()
op3_H = pracma::hessian(f = loglik, 
                        x0 = op3$par,
                        nit = nit,
                        Dit = Dit, 
                        rit = rit, 
                        o_t_c = phase_cov,
                        K   = K)
tictoc::toc()
#  sec elapsed
# > op3_H

# sd_err3 = sqrt(abs(diag(solve(op3_H))))




#############################################
# Model3b: omega ~ mobility
#############################################

# par := starting parameter values (lambda, gamma, 
#                                   omega, omega_mob1, omega_mob2, omega_mob3, omega_mob4, omega_mob5, omega_mob6
#                                   pdet, pmor, prec)
par = c(4.01100794, -19.17485603, 
        0.32681933, 5.33101709, -5.71155838, 0.46610098, 0.44765711, -0.36431596, 5.31133289,
        -0.87870088, -6.07408587, -0.07249772)
      
# check for finite loglik at initial values:
K = matrix(200, nrow=nrow(nit), ncol=ncol(nit))
loglik(par = par, nit = nit, K = K, Dit = Dit, rit = rit, o_t_c = mobility_cov)
# [1] 

##
# Testing: K=200
##
K = matrix(200, nrow=nrow(nit), ncol=ncol(nit))
tictoc::tic()
op4 = optim(par = par,
            fn  = loglik, method = "BFGS", hessian = F,
            nit = nit,
            Dit = Dit, 
            rit = rit, 
            K   = K, 
            o_t_c = mobility_cov,
            HPDfile = NULL,
            control = list(trace=1, REPORT=1, reltol=1e-10))
tictoc::toc()
# 127748.06 sec elapsed
saveRDS(op4, "results/mod4.RDS")
# > op4
# $par
# [1]   4.10657717 -19.17485474   0.68732418   3.95057399  -5.22459321   0.49893297   1.98892378  -0.12195271   5.58128977
# [10]  -1.29974289  -6.15622128  -0.08133562
# 
# $value
# [1] 273.3618
# 
# $counts
# function gradient 
# 95       51 
# 
# $convergence
# [1] 0
# 
# $message
# NULL



calcAIC(op4$value, op4$par)
# > calcAIC(op4$value, op4$par)
# [1] 570.7235

# Hessian: (NOT RUN)
tictoc::tic()
op4_H = pracma::hessian(f = loglik, 
                        x0 = op4$par,
                        nit = nit,
                        Dit = Dit, 
                        rit = rit, 
                        o_t_c = mobility_cov,
                        K   = K)
tictoc::toc()
#  sec elapsed
#> op4_H

# sd_err4 = sqrt(abs(diag(solve(op4_H))))




#############################################
# Model5: pdet ~ testVolume, omega ~ mobility + phases
#############################################

# par := starting parameter values (lambda, gamma, 
#                                   omega, omega_mob1, omega_mob2, omega_mob3, omega_mob4, omega_mob5, omega_mob6,
#                                   omega_phase2, omega_phase3, omega_phase4,
#                                   pdet, pdet_testVolume
#                                   pmor, prec)
par = c(4.37553225, -19.17485598,  
        -0.41203458,  14.42090562,   1.68131988,  -0.65582151,   0.17686817,  -6.79683923,  22.78451110,
        -0.74427820,  -0.50164094,  -0.48658199,  
        -1.50654176, 1.10157056,  
        -6.39465062, -0.09215878)

# check for finite loglik at initial values:
K = matrix(200, nrow=nrow(nit), ncol=ncol(nit))
loglik(par = par, nit = nit, K = K, Dit = Dit, rit = rit, o_t_c = c(mobility_cov, phase_cov), p_t_c=testVolume)
# [1] 274.2649

##
# Testing: K=200
##
K = matrix(200, nrow=nrow(nit), ncol=ncol(nit))
tictoc::tic()
op5 = optim(par = par,
            fn  = loglik, method = "BFGS", hessian = F,
            nit = nit,
            Dit = Dit, 
            rit = rit, 
            K   = K, 
            o_t_c = c(mobility_cov, phase_cov), 
            p_t_c = testVolume,
            HPDfile = NULL,
            control = list(trace=1, REPORT=1, reltol=1e-10))
tictoc::toc()
# 230485.87 sec elapsed (used 20 cores)
saveRDS(op5, "results/mod5.RDS")

# > op5
# $par
# [1]   4.31198857 -19.17485465   1.81566931   6.50186129   2.66400892  -0.74187487  11.47352810  -4.59811073  33.68423379  -0.47287583
# [11]   0.33545348   0.22596618  -1.56863908   0.45872592  -6.33746653  -0.09154571  ##-0.07205155## <- accidental extra input par, remove 1 from params for AIC
# 
# $value
# [1] 268.307
# 
# $counts
# function gradient 
# 127       71 
# 
# $convergence
# [1] 0
# 
# $message
# NULL


calcAIC(op5$value, op5$par[-1])
# > calcAIC(op5$value, op5$par)
# [1] 568.614


# Hessian: (NOT RUN)
tictoc::tic()
op5_H = pracma::hessian(f = loglik, 
                        x0 = op5$par,
                        nit = nit,
                        Dit = Dit, 
                        rit = rit,  
                        o_t_c = c(mobility_cov, phase_cov), 
                        p_t_c = testVolume,
                        K   = K)
tictoc::toc()
# 44089.33 sec elapsed
#> op5_H

# sd_err5 = sqrt(abs(diag(solve(op5_H[1:16,1:16]))))






#############################################
# Model6: pdet ~ testVolume, omega ~ phases
#############################################

# par := starting parameter values (lambda, gamma, 
#                                   omega_phase1, omega_phase2, omega_phase3, omega_phase4
#                                   pdet, pdet_testVolume
#                                   pmor, prec)
par = c(4.15703041, -19.16999998, 
        -0.04281684, -0.29448922, 0.26672537, -0.12424797,
        -1.24660654,  0.61661025,
        -6.33746653, -0.09154571)

# check for finite loglik at initial values:
K = matrix(200, nrow=nrow(nit), ncol=ncol(nit))
loglik(par = par, nit = nit, K = K, Dit = Dit, rit = rit, o_t_c = phase_cov, p_t_c=testVolume)
# [1] 274.3148

##
# Testing: K=200
##
K = matrix(200, nrow=nrow(nit), ncol=ncol(nit))
tictoc::tic()
op6 = optim(par = par,
            fn  = loglik, method = "BFGS", hessian = F,
            nit = nit,
            Dit = Dit, 
            rit = rit, 
            K   = K, 
            o_t_c = phase_cov, 
            p_t_c = testVolume,
            HPDfile = NULL,
            control = list(trace=1, REPORT=1, reltol=1e-10))
tictoc::toc()
# 37042.17 sec elapsed (used 20 cores)
saveRDS(op6, "results/mod6.RDS")

# > op6
# $par
# [1]   4.17232608 -19.16999995  -0.09110563  -0.24562316   0.34021945  -0.13665910  -1.42683222   0.30790731  -6.21347856
# [10]  -0.08535224
# 
# $value
# [1] 270.4497
# 
# $counts
# function gradient 
# 54       16 
# 
# $convergence
# [1] 0
# 
# $message
# NULL


calcAIC(op6$value, op6$par)
# > calcAIC(op6$value, op6$par)
# [1] 560.8995

# Hessian:
tictoc::tic()
op6_H = pracma::hessian(f = loglik, 
                        x0 = op6$par,
                        nit = nit,
                        Dit = Dit, 
                        rit = rit,  
                        o_t_c = phase_cov, 
                        p_t_c = testVolume,
                        K   = K)
tictoc::toc()
saveRDS(op6_H, "results/mod6_Hessian.RDS")
# 20723.75 sec elapsed
#> op6_H
# [,1]          [,2]        [,3]          [,4]          [,5]          [,6]          [,7]          [,8]               [,9]         [,10]
# [1,]  3.388739e+01 -9.536743e-07   24.418430  3.363552e-01  1.570797e-02  6.580353e-05  1.302226e+01 -3.446047e+00  3.449878e+00  1.839298e+01
# [2,] -9.536743e-07  3.814697e-06    0.000000 -1.907349e-06 -1.907349e-06 -1.907349e-06 -1.907349e-06  9.536743e-07  0.000000e+00 -9.536743e-07
# [3,]  2.441843e+01  0.000000e+00 1397.038971  7.934465e+01  6.661450e+02  3.840279e+02  1.972869e+01  1.476714e+02 -2.719311e+00 -1.449407e+01
# [4,]  3.363552e-01 -1.907349e-06   79.344651  4.895654e+01  1.044937e+01  4.458427e-02 -2.469660e+00 -5.807304e-02 -3.746223e-02 -2.002287e-01
# [5,]  1.570797e-02 -1.907349e-06  666.144991  1.044937e+01  6.002946e+02  5.447216e+01 -1.890773e+01  6.261608e+01 -1.749992e-03 -9.346008e-03
# [6,]  6.580353e-05 -1.907349e-06  384.027949  4.458427e-02  5.447216e+01  3.295073e+02  1.730482e+01  9.657480e+01 -5.722046e-06 -3.814697e-05
# [7,]  1.302226e+01 -1.907349e-06   19.728693 -2.469660e+00 -1.890773e+01  1.730482e+01  4.476342e+01  1.219426e+01 -1.450192e+00 -7.738833e+00
# [8,] -3.446047e+00  9.536743e-07  147.671373 -5.807304e-02  6.261608e+01  9.657480e+01  1.219426e+01  9.812044e+01  3.837614e-01  2.047958e+00
# [9,]  3.449878e+00  0.000000e+00   -2.719311 -3.746223e-02 -1.749992e-03 -5.722046e-06 -1.450192e+00  3.837614e-01  7.614578e+00  1.788172e+00
# [10,]  1.839298e+01 -9.536743e-07  -14.494065 -2.002287e-01 -9.346008e-03 -3.814697e-05 -7.738833e+00  2.047958e+00  1.788172e+00  6.246219e+02

# NOTE: gamma standard error removed due to singular matrix, eigenvalue 0
sd_err6 = sqrt(diag(solve(op6_H)))
# [1]   0.19824309 512.00003732   0.06897970   0.17353814   0.08276447   0.09656522   0.17170612   0.12699960   0.37912277
# [10]   0.04063373 

# exp(op6$par[1])
# 95% CI: exp(op6$par[1]-1.96*sd_err6[1]), exp(op6$par[1]+1.96*sd_err6[1])
# lambda: 64.86616 (43.98152, 95.66787)
#
# exp(op6$par[2])
# 95% CI: 
# gamma: 4.726882e-09 (0, Inf)
#
# exp(op6$par[3:6])
# # 95% CI: exp(op6$par[3]-1.96*sqrt(sd_err6[3]^2)), exp(op6$par[3]+1.96*sqrt(sd_err6[3]^2))
#           exp(sum(op6$par[3:4])-1.96*sqrt(sum(sd_err6[3:4]^2))), exp(sum(op6$par[3:4])+1.96*sqrt(sum(sd_err6[3:4]^2)))
#           exp(sum(op6$par[c(3,5)])-1.96*sqrt(sum(sd_err6[c(3,5)]^2))), exp(sum(op6$par[c(3,5)])+1.96*sqrt(sum(sd_err6[c(3,5)]^2)))
#           exp(sum(op6$par[c(3,6)])-1.96*sqrt(sum(sd_err6[c(3,6)]^2))), exp(sum(op6$par[c(3,6)])+1.96*sqrt(sum(sd_err6[c(3,6)]^2)))
#           
#          
# omega_phase: 
#      Phase = 1:   0.9129213 (0.7974742, 1.045081)
#      Phase = 2:   0.7141025 (0.4952221, 1.029725)
#      Phase = 3:   1.282888 (1.03867, 1.584528)
#      Phase = 4:   0.7963116 (0.6310553, 1.004844)
#
# plogis(op6$par[7] + op6$par[8]*testVolume$testVolume)
# # 95% CI: plogis(op6$par[7] + op6$par[8]*testVolume$testVolume-1.96*sqrt(sd_err6[7:8]^2)), 
#           plogis(op6$par[7] + op6$par[8]*testVolume$testVolume+1.96*sqrt(sd_err6[7:8]^2)) 
# pdet_t: 
#      t = 1:  0.1709495  (0.1283690, 0.2240239)
#      t = 2:  0.2174460  (0.1780622, 0.2627559)
#      t = 3:  0.1581698  (0.1183182, 0.2082732)
#      t = 4:  0.1578993  (0.1275428, 0.1938754)
#      t = 5:  0.1665303  (0.1248848, 0.2185945)
#      t = 6:  0.1901051  (0.1546944, 0.2314026)
#      t = 7:  0.1856085  (0.1399933, 0.2419061)
#      t = 8:  0.1677329  (0.1357907, 0.2054029)
#      t = 9:  0.1640066  (0.1228991, 0.2154857)
#      t = 10: 0.1683723  (0.1363283, 0.2061504)
#      t = 11: 0.1680168  (0.1260558, 0.2204228)
#      t = 12: 0.1692279  (0.1370479, 0.2071502)
#      t = 13: 0.1637281  (0.1226802, 0.2151423)
#      t = 14: 0.1610325  (0.1301668, 0.1975549)
#      t = 15: 0.1589155  (0.1189026, 0.2091965)
#      t = 16: 0.1594597  (0.1288491, 0.1957086)
#      t = 17: 0.1681589  (0.1261678, 0.2205975)
#      t = 18: 0.1726124  (0.1398972, 0.2111003)
#      t = 19: 0.1830109  (0.1379260, 0.2387516)
#      t = 20: 0.1668833  (0.1350767, 0.2044094)
#      t = 21: 0.1758292  (0.1322271, 0.2299983)
#      t = 22: 0.1942865  (0.1582492, 0.2362275)
#      t = 23: 0.1917536  (0.1448969, 0.2493444)
#      t = 24: 0.2321421  (0.1907454, 0.2794209)
#      t = 25: 0.2218884  (0.1692096, 0.2853347)
#      t = 26: 0.2410445  (0.1984705, 0.2894529)
#      t = 27: 0.2938805  (0.2291432, 0.3681719)
#      t = 28: 0.3976693  (0.3398181, 0.4585288)
#      t = 29: 0.3320896  (0.2620592, 0.4104259)
#      t = 30: 0.2371614  (0.1950971, 0.2850830)
#
# plogis(op6$par[9])
# 95% CI: plogis(op6$par[9]-1.96*sd_err6[9]), plogis(op6$par[9]+1.96*sd_err6[9])
# pmort: 0.001998259 (0.0009514617, 0.004191911)
#
# plogis(op6$par[10])
# 95% CI: plogis(op6$par[10]-1.96*sd_err6[10]), plogis(op6$par[10]+1.96*sd_err6[10])
# prec:  0.4786749 (0.4588447, 0.4985725)




#############################################
# Model7: pdet ~ testVolume, omega ~ mobility
#############################################

# par := starting parameter values (lambda, gamma, 
#                                   omega, omega_mob1, omega_mob2, omega_mob3, omega_mob4, omega_mob5, omega_mob6,
#                                   pdet, pdet_testVolume,
#                                   pmor, prec)
par = c(4.25396516, -19.17, 
        -2.59983434, 12.03382905, -0.82133631, -0.22355791, -5.29845423, -6.90907070, 7.49651192,
        -1.42124022, 1.05652925,  
        -6.28536982, -0.08789352)

# check for finite loglik at initial values:
K = matrix(200, nrow=nrow(nit), ncol=ncol(nit))
loglik(par = par, nit = nit, K = K, Dit = Dit, rit = rit, o_t_c = mobility_cov, p_t_c=testVolume)
# [1] 

##
# Testing: K=200
##
K = matrix(200, nrow=nrow(nit), ncol=ncol(nit))
tictoc::tic()
op7 = optim(par = par,
            fn  = loglik, method = "BFGS", hessian = F,
            nit = nit,
            Dit = Dit, 
            rit = rit, 
            K   = K, 
            o_t_c = mobility_cov, 
            p_t_c = testVolume,
            HPDfile = NULL,
            control = list(trace=1, REPORT=1, reltol=1e-10))
tictoc::toc()
# 142374.08 sec elapsed (used 20 cores)
saveRDS(op7, "results/mod7.RDS")
# > op7
# $par
# [1]   4.20930804 -19.16999902   0.37672452   6.45479483  -2.04517220  -0.09622287   3.54834014  -4.51801891  10.46568076
# [10]  -1.44277740   0.31872481  -6.24600987  -0.08728158
# 
# $value
# [1] 270.7433
# 
# $counts
# function gradient 
# 105       55 
# 
# $convergence
# [1] 0
# 
# $message
# NULL


calcAIC(op7$value, op7$par)
# > calcAIC(op7$value, op7$par)
# [1] 567.4866

# Hessian:
tictoc::tic()
op7_H = pracma::hessian(f = loglik, 
                        x0 = op7$par,
                        nit = nit,
                        Dit = Dit, 
                        rit = rit,  
                        o_t_c = mobility_cov, 
                        p_t_c = testVolume,
                        K   = K)
tictoc::toc()
# 28370.72 sec elapsed
#> op7_H

sd_err7 = sqrt(abs(diag(solve(op7_H))))
# > sd_err7
# [1]   

# exp(op7$par[1])
# 95% CI: exp(op7$par[1]-1.96*sd_err7[1]), exp(op7$par[1]+1.96*sd_err7[1])
# lambda: 70.38394 (45.57406, 108.7)
#
# exp(op7$par[2])
# 95% CI: exp(op7$par[2]-1.96*sd_err7[2]), exp(op7$par[2]+1.96*sd_err7[2])
# gamma: 4.726883e-09 (0, Inf)
#
# exp(op7$par[3])
# 95% CI: exp(op7$par[3]-1.96*sd_err7[3]), exp(op7$par[3]+1.96*sd_err7[3])
# omega base: 0.07428588 (0.0004220821, 13.07421)
#
# exp(op7$par[3] + op7$par[4]*mobility_cov$mob_retail + op7$par[5]*mobility_cov$mob_grocery + op7$par[6]*mobility_cov$mob_parks + op7$par[7]*mobility_cov$mob_transit + op7$par[8]*mobility_cov$mob_work + op7$par[9]*mobility_cov$mob_residential)
# # 95% CI: exp(op7$par[3] + op7$par[4]*mobility_cov$mob_retail + op7$par[5]*mobility_cov$mob_grocery + op7$par[6]*mobility_cov$mob_parks + op7$par[7]*mobility_cov$mob_transit + op7$par[8]*mobility_cov$mob_work + op7$par[9]*mobility_cov$mob_residential - 1.96*sqrt(sum(sd_err7[3:9]^2))), 
#           exp(op7$par[3] + op7$par[4]*mobility_cov$mob_retail + op7$par[5]*mobility_cov$mob_grocery + op7$par[6]*mobility_cov$mob_parks + op7$par[7]*mobility_cov$mob_transit + op7$par[8]*mobility_cov$mob_work + op7$par[9]*mobility_cov$mob_residential + 1.96*sqrt(sum(sd_err7[3:9]^2)))
#          
# omega_t: 
#      t = 1:   0.4712659 (3.921088e-13, 5.664029e+11)
#      t = 2:   0.8567722 (7.128627e-13, 1.029733e+12)
#      t = 3:   1.0974019 (9.130747e-13, 1.318940e+12)
#      t = 4:   1.5546843 (1.293549e-12, 1.868537e+12)
#      t = 5:   0.8906733 (7.410696e-13, 1.070478e+12)
#      t = 6:   0.7302301 (6.075755e-13, 8.776455e+11)
#      t = 7:   0.7528903 (6.264296e-13, 9.048804e+11)
#      t = 8:   0.6612080 (5.501469e-13, 7.946896e+11)
#      t = 9:   1.0476663 (8.716929e-13, 1.259164e+12)
#      t = 10:  0.8352844 (6.949842e-13, 1.003908e+12)
#      t = 11:  0.7796727 (6.487134e-13, 9.370694e+11)
#      t = 12:  0.8793501 (7.316483e-13, 1.056869e+12)
#      t = 13:  0.9342407 (7.773191e-13, 1.122841e+12)
#      t = 14:  0.9759618 (8.120324e-13, 1.172984e+12)
#      t = 15:  1.8001265 (1.497765e-12, 2.163528e+12)
#      t = 16:  1.4819273 (1.233012e-12, 1.781092e+12)
#      t = 17:  1.2618902 (1.049934e-12, 1.516635e+12)
#      t = 18:  1.1809266 (9.825699e-13, 1.419327e+12)
#      t = 19:  1.1513111 (9.579288e-13, 1.383732e+12)
#      t = 20:  1.8578401 (1.545784e-12, 2.232892e+12)
#      t = 21:  1.1629449 (9.676086e-13, 1.397715e+12)
#      t = 22:  0.9687528 (8.060344e-13, 1.164320e+12)
#      t = 23:  1.0017535 (8.334920e-13, 1.203983e+12)
#      t = 24:  0.9792830 (8.147958e-13, 1.176976e+12)
#      t = 25:  1.5058656 (1.252930e-12, 1.809863e+12)
#      t = 26:  0.8605862 (7.160361e-13, 1.034317e+12)
#      t = 27:  0.8549088 (7.113124e-13, 1.027494e+12)
#      t = 28:  0.8497677 (7.070348e-13, 1.021315e+12)
#      t = 29:  0.6548574 (5.448629e-13, 7.870569e+11)
#      t = 30:  0.9619383 (8.003644e-13, 1.156130e+12)
#
# plogis(op7$par[10] + op7$par[11]*testVolume$testVolume)
# # 95% CI: plogis(op7$par[10] + op7$par[11]*testVolume$testVolume-1.96*sqrt(sd_err7[10:11]^2)), 
#           plogis(op7$par[10] + op7$par[11]*testVolume$testVolume+1.96*sqrt(sd_err7[10:11]^2)) 
#           
# pdet_t: 
#      t = 1:  0.12531010 (0.09059036, 0.1708368)
#      t = 2:  0.28505705 (0.19005317, 0.4038709)
#      t = 3:  0.09430409 (0.06751213, 0.1302433)
#      t = 4:  0.09371005 (0.05736173, 0.1494397)
#      t = 5:  0.11393233 (0.08206925, 0.1560629)
#      t = 6:  0.18265996 (0.11623468, 0.2752239)
#      t = 7:  0.16804245 (0.12314980, 0.2250990)
#      t = 8:  0.11695928 (0.07231258, 0.1837126)
#      t = 9:  0.10774711 (0.07746261, 0.1479724)
#      t = 10: 0.11858987 (0.07337244, 0.1860777)
#      t = 11: 0.11768158 (0.08487040, 0.1609467)
#      t = 12: 0.12079457 (0.07480785, 0.1892676)
#      t = 13: 0.10707839 (0.07696564, 0.1470952)
#      t = 14: 0.10074784 (0.06185602, 0.1599243)
#      t = 15: 0.09595515 (0.06872972, 0.1324316)
#      t = 16: 0.09717210 (0.05956919, 0.1546094)
#      t = 17: 0.11804417 (0.08514165, 0.1614182)
#      t = 18: 0.12977243 (0.08068147, 0.2021645)
#      t = 19: 0.15991625 (0.11688951, 0.2149265)
#      t = 20: 0.11481545 (0.07092138, 0.1805954)
#      t = 21: 0.13868578 (0.10068664, 0.1880274)
#      t = 22: 0.19686388 (0.12606999, 0.2940363)
#      t = 23: 0.18819052 (0.13881340, 0.2500321)
#      t = 24: 0.34749793 (0.23862968, 0.4750447)
#      t = 25: 0.30352168 (0.23255261, 0.3852759)
#      t = 26: 0.38676887 (0.27070158, 0.5173023)
#      t = 27: 0.61461311 (0.52582086, 0.6963793)
#      t = 28: 0.88595211 (0.82052273, 0.9295761)
#      t = 29: 0.74593490 (0.67121334, 0.8085193)
#      t = 30: 0.36954927 (0.25648814, 0.4990001)
#
# plogis(op7$par[12])
# 95% CI: plogis(op7$par[12]-1.96*sd_err7[12]), plogis(op7$par[12]+1.96*sd_err7[12])
# pmort: 0.001859902 (0.0008647486, 0.003995699)
#
# plogis(op7$par[13])
# 95% CI: plogis(op7$par[13]-1.96*sd_err7[13]), plogis(op7$par[13]+1.96*sd_err7[13])
# prec:  0.4780408 (0.45865, 0.4974979)







