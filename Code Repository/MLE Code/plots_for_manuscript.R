library(magrittr)
library(ggplot2)
library(dplyr)
library(viridis)
library(patchwork)
library(wesanderson)

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

###############################################
# Data Figure
ggplot(data=data.frame(nit=nit, rit=rit, Dit=Dit, Week=1:30)) +
  geom_point(aes(x=Week, y=nit), color="dodgerblue") +
  theme_classic() +
  scale_x_continuous(breaks=c(1,seq(0,30,5)[-1])) +
  ylab("Newly Observed\nCase Counts") -> p_data1
ggplot(data=data.frame(nit=nit, rit=rit, Dit=Dit, Week=1:30)) +
  geom_point(aes(x=Week, y=rit), color="green") +
  theme_classic() +
  scale_x_continuous(breaks=c(1,seq(0,30,5)[-1])) +
  ylab("Newly Observed\nRecoveries") -> p_data2
ggplot(data=data.frame(nit=nit, rit=rit, Dit=Dit, Week=1:30)) +
  geom_point(aes(x=Week, y=Dit), color="firebrick") +
  theme_classic() +
  scale_y_continuous(breaks=c(0,1,2)) +
  scale_x_continuous(breaks=c(1,seq(0,30,5)[-1])) +
  ylab("Newly Observed\nDeaths") -> p_data3

data_plot = p_data1/p_data2/p_data3
ggsave(
  "./results/plots/data_plot.png",
  plot = data_plot,
  device = "png",
  path = NULL,
  scale = 1,
  width = 8,
  height = 5,
  units = "in",
  dpi = 300
)
###############################################

# Calculate ait from data:
Calc_at2 <- function(t) {
  if(t==0) { return(0) }
  return(sum(nit[1:t])-sum(rit[1:t])-sum(Dit[1:t]))
}

a = 0
for(t in 1:29) {
  a = c(a, Calc_at2(t))
}

# Covariates:
testVolume = list(testVolume = as.numeric(scale(dat$New.Tests..weekly.)))
testVolume2 = list(testVolume = dat$New.Tests..weekly.)

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


# Fitted model data:
op6 = readRDS(file = "results/mod6.RDS")

op6_H = readRDS(file = "results/mod6_Hessian.RDS")
sd_err6 = sqrt(diag(solve(op6_H)))

# lambda_est
lambda_est = exp(op6$par[1])
# 95% CI: 
lambda_low = exp(op6$par[1]-1.96*sd_err6[1])
lambda_upp = exp(op6$par[1]+1.96*sd_err6[1])
# lambda: 64.86616 (43.98152, 95.66787)

# gamma_est
gamma_est =  exp(op6$par[2])
# 95% CI: 
gamma_low = exp(op6$par[2]-1.96*sd_err6[2])
gamma_upp = exp(op6$par[2]+1.96*sd_err6[2])
# gamma: 4.726883e-09 (0, Inf)

# omega_est
omega_est = exp(c(op6$par[3], op6$par[3] + op6$par[4], op6$par[3] + op6$par[5], op6$par[3] + op6$par[6]))
# 95% CI: 
omega_low = exp(c(op6$par[3] - 1.96*sqrt(sd_err6[3]^2), op6$par[3] + op6$par[4] - 1.96*sqrt(sum(sd_err6[c(3,4)]^2)), op6$par[3] + op6$par[5] - 1.96*sqrt(sum(sd_err6[c(3,5)]^2)), op6$par[3] + op6$par[6] - 1.96*sqrt(sum(sd_err6[c(3,6)]^2))))
omega_upp = exp(c(op6$par[3] + 1.96*sqrt(sd_err6[3]^2), op6$par[3] + op6$par[4] + 1.96*sqrt(sum(sd_err6[c(3,4)]^2)), op6$par[3] + op6$par[5] + 1.96*sqrt(sum(sd_err6[c(3,5)]^2)), op6$par[3] + op6$par[6] + 1.96*sqrt(sum(sd_err6[c(3,6)]^2))))

# pdet_est
pdet_est = plogis(op6$par[7] + op6$par[8]*testVolume$testVolume)
# 95% CI: 
pdet_low = plogis(op6$par[7] + op6$par[8]*testVolume$testVolume-1.96*sqrt(sum(sd_err6[7:8]^2)))
pdet_upp = plogis(op6$par[7] + op6$par[8]*testVolume$testVolume+1.96*sqrt(sum(sd_err6[6:8]^2)))

# pmort_est      
pmort_est = plogis(op6$par[9])
# 95% CI: 
pmort_low = plogis(op6$par[9]-1.96*sd_err6[9])
pmort_upp = plogis(op6$par[9]+1.96*sd_err6[9])
# pmort: 0.001998259 (0.0009514617, 0.004191911)

# prec_est
prec_est = plogis(op6$par[10])
# 95% CI: 
prec_low = plogis(op6$par[10]-1.96*sd_err6[10])
prec_upp = plogis(op6$par[10]+1.96*sd_err6[10])
# prec:  0.4786749 (0.4588447, 0.4985725)


# PLOTTING
df = data.frame(p = pdet_est,
                nit = nit,
                TestVolume = testVolume2$testVolume,
                Lower = pdet_low, 
                Upper = pdet_upp,
                Week = 1:30)

ggplot(data=df) + 
  geom_ribbon(aes(ymin=Lower, ymax=Upper, x=Week), alpha = 0.2) +
  geom_errorbar(aes(ymin=Lower, ymax=Upper, x=Week), width=0.6) +
  geom_point(aes(x=Week, y=p)) +
  theme_classic() +
  scale_y_continuous(breaks=seq(-0.1,0.95,0.15)) + 
  scale_x_continuous(breaks=c(1,seq(0,30,5)[-1])) +
  ylab("Estimated Probability\nof Detection") -> p1

ggplot(data=df) +
  geom_line(aes(x=Week, y=TestVolume)) +
  geom_point(aes(x=Week, y=TestVolume)) +
  theme_classic() +
  scale_y_continuous(breaks=seq(0,3000,500)) + 
  scale_x_continuous(breaks=c(1,seq(0,30,5)[-1])) +
  ylab("Weekly Tests\nAdministered") -> p2

pdet_plot = p1/p2
ggsave(
  "./results/plots/pdet_estimates.png",
  plot = pdet_plot,
  device = "png",
  path = NULL,
  scale = 1,
  width = 8,
  height = 4,
  units = "in",
  dpi = 300
)


###########################################
# omega ~ phases, pdet ~ testVolume
###########################################
omega_est2 = rep(omega_est, times=c(8,5,12,5))
omega_upp2 = rep(omega_upp, times=c(8,5,12,5))
omega_low2 = rep(omega_low, times=c(8,5,12,5))

df4 = data.frame(omega = omega_est2,
                 nit = nit,
                 TestVolume = testVolume2$testVolume,
                 Lower = omega_low2, 
                 Upper = omega_upp2,
                 Week = 1:30,
                 Phase = ordered(rep(c("Phase 1", "Phase 2", "Phase 3a", "Phase3b"), times=c(8,5,12,5)), 
                                 levels=c("Phase 1", "Phase 2", "Phase 3a", "Phase3b")))

ggplot(data=df4) +
  geom_ribbon(aes(ymin=0.45, ymax=1.65, x=Week, fill=Phase), alpha = 0.5, show.legend = F) +
  geom_ribbon(aes(ymin=Lower, ymax=Upper, x=Week), alpha = 0.2) +
  geom_errorbar(aes(ymin=Lower, ymax=Upper, x=Week), width=0.6) +
  geom_line(aes(x=Week, y=omega)) +
  geom_point(aes(x=Week, y=omega)) +
  annotate("text", x = c(4,11,19,28), y = 1.7, label = unique(df4$Phase)) +
  theme_classic() +
  scale_x_continuous(breaks=c(1,seq(0,30,5)[-1])) +
  scale_fill_manual(values=wes_palette(n=4, name="Darjeeling1", type = "continuous")) +
  ylab("Estimated Domestic\nSpread") -> p5

ggsave(
  "./results/plots/omega_estimates2.png",
  plot = p5,
  device = "png",
  path = NULL,
  scale = 1,
  width = 8,
  height = 3,
  units = "in",
  dpi = 300
)

N_est = (nit)/pdet_est
N_upp = (nit)/pdet_low
N_low = (nit)/pdet_upp
df5 = data.frame(N = c(N_est),
                 n = c(nit),
                 Upper = c(N_upp),
                 Lower = c(N_low),
                 Week = 1:30)

ggplot(data=df5) +
  geom_ribbon(aes(ymin=Lower, ymax=Upper, x=Week), alpha = 0.2) +
  geom_errorbar(aes(ymin=Lower, ymax=Upper, x=Week), width=0.6) +
  geom_line(aes(x=Week, y=N), color="dodgerblue") +
  geom_point(aes(x=Week, y=N), color="blue") +
  geom_line(aes(x=Week, y=n), color="red") +
  geom_point(aes(x=Week, y=n), color="firebrick") +
  theme_classic() +
  scale_x_continuous(breaks=c(1,seq(0,30,5)[-1])) +
  ggtitle(label = "", subtitle="Blue (top line): estimated active cases\nRed (bottom line): newly observed active cases") + 
  ylab("Estimated Active\nCases") -> p6

ggsave(
  "./results/plots/N_estimates.png",
  plot = p6,
  device = "png",
  path = NULL,
  scale = 1,
  width = 8,
  height = 3,
  units = "in",
  dpi = 300
)
  
  