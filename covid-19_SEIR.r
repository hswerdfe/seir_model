rm(list = ls())

#title: "Deterministic compartmental model for 2019-nCov 2.0 "
#author: "Antoinette Ludwig, Erin Rees'

#date: "March, 16th, 2020"

#################################################################################################

### Resources:

#http://sherrytowers.com/2015/10/01/difference-between-markov-chain-monte-carlo-stochastic-differential-equations-and-agent-based-models/
#http://sherrytowers.com/2016/02/06/stochastic-compartmental-modelling-with-stochastic-differential-equations-2/
#https://cran.r-project.org/web/packages/MultiBD/vignettes/SIR-MCMC.pdf
#https://cran.r-project.org/web/packages/GillespieSSA/index.html
#https://math.unm.edu/~sulsky/mathcamp/SimpleStochModel.pdf


# Gillispie method for random noise added to demographics
# Includes stochasticity for the timing of the next event, and which event is next
# Continuous time, discrete population changes

# But instead, this approach adds stochasticity via the mechanistic parameters (e.g. contact rate)

#################################################################################################

# install.packages('deSolve')
# install.packages('ggplot2')
# install.packages('tidyr')
# install.packages('dplyr')
# install.packages('tidyverse')
# install.packages('matrixStats')
#install.packages('xlsx')
#install.packages('triangle')
### Libinstall.packages()

library(deSolve)
library(ggplot2)
library(tidyr)
library(dplyr)
library(tidyverse)
library(matrixStats)
library(readxl)
library(janitor)

#################################################################################################
### User input for the initial conditions
time_stuff <- read_excel("~/Projects/covid_19/model/input_sheet.xlsx", 
                          sheet = "time")

time_param_names <- 
  time_stuff %>% 
  select(-tmin, -tmax) %>% 
  colnames()

time_stuff <- time_stuff %>% 
  mutate(t_n = tmax - tmin) %>% 
  mutate(., isim = 1:nrow(.))
  
all_param <- time_stuff %>% pivot_longer(cols = time_param_names)
# 
# with(data = init, S)
# init[names(init) == "S"]
# 
# init[S]
# eval(parse(text="N/5"))

#params[names(params) == "cq"]

input_stuff <- read_excel("~/Projects/covid_19/model/input_sheet.xlsx", 
                        sheet = "input")


init <- input_stuff$VALUE
names(init) <- input_stuff$NAME
nSim <- nrow(time_stuff)
#################################################################################################


#########################################333333
#
# SEIR Model
#
seir <- function(time, state, parms) {
  
  with(as.list(c(state, parms)), {
    
    #dS <- -beta*c*(1-lambda)*tau*S*(I_a + I_aqn + I_sm + I_ss + I_smisn + I_ssisn) -beta*cr*(1-lambda)*(1-tau)*S*(I_ar+I_smr+I_ssr+I_smrisn+I_ssrisn)-beta*cq*lambda*S*I_aq
    dS <- -beta*((c*(1-lambda)*tau)+(cq*lambda)+(cr*(1-lambda)*(1-tau)))*S*(I_a + I_aqn + I_sm + I_ss + I_smisn + I_ssisn+I_ar+I_smr+I_ssr+I_smrisn+I_ssrisn+phi*I_aq)
    dL <- (1-lambda)*beta*c*tau*S*(I_a + I_aqn + I_sm + I_ss + I_smisn + I_ssisn+I_ar+I_smr+I_ssr+I_smrisn+I_ssrisn+phi*I_aq) - sigma*L
    dL_q <- lambda*beta*cq*S*(I_a + I_aqn + I_sm + I_ss + I_smisn + I_ssisn+I_ar+I_smr+I_ssr+I_smrisn+I_ssrisn+phi*I_aq) - sigma*L_q
    dL_r <-  beta*cr*(1-lambda)*(1-tau)*S*(I_a + I_aqn + I_sm + I_ss + I_smisn + I_ssisn+I_ar+I_smr+I_ssr+I_smrisn+I_ssrisn+phi*I_aq)-sigma*L_r
    dI_a <- sigma*L - I_a*delta*epsilon - I_a*(1-delta)*upsilon
    dI_aq <- sigma*rho*L_q - I_aq*delta*epsilon - I_aq*(1-delta)*upsilon
    dI_ar <- sigma*L_r - I_ar*delta*epsilon - I_ar*(1-delta)*upsilon
    dI_aqn <- sigma*(1-rho)*L_q - I_aqn*delta*epsilon - I_aqn*(1-delta)*upsilon
    dI_sm <- I_a*delta*epsilon*alpha - kappa*I_sm
    dI_ss <- I_a*delta*epsilon*(1-alpha) - kappa*I_ss
    dI_smr <- I_ar*delta*epsilon*alpha - kappa*I_smr
    dI_ssr <- I_ar*delta*epsilon*(1-alpha) - kappa*I_ssr
    dI_smis <- kappa*feim*I_sm + kappa*feimr*I_smr + delta*alpha*epsilonq*feimq*I_aq - num*I_smis
    dI_smisn <- kappa*(1-feim)*I_sm - num*I_smisn
    dI_ssis <- kappa*feisi*(I_ss+I_ssr) - I_ssis*((1-mu)*nus + mu*nud)
    dI_ssisn <- kappa*(1-feisi-feish)*(I_ss+I_ssr) - I_ssisn*((1-mu)*nus + mu*nud)
    dI_ssh <- kappa*feish*(I_ss+I_ssr) + delta*(1-alpha)*epsilonq*I_aq - I_ssh*((1-mu)*nus + mu*nud)
    dI_smrisn <- kappa*(1-feimr)*I_smr - num*I_smrisn
    dI_ssrisn <- kappa*(1-feisi-feish)*(I_ssr) - I_ssrisn*((1-mu)*nus + mu*nud)
    dI_smqisn <- I_aq*delta*alpha*epsilonq*(1-feimq) - num*I_smqisn
    
    dR <- (I_a +I_aq+I_aqn+I_ar)*(1-delta)*upsilon + num*(I_smis+I_smisn+I_smqisn+I_smrisn) + (I_ssis + I_ssisn + I_ssh+I_ssrisn)*(1-mu)*nus
    dD <- mu*nud*(I_ssis + I_ssisn + I_ssh + I_ssrisn)
    
    return(list(c(dS, dL, dL_q, dL_r, dI_a, dI_aq, dI_ar, dI_aqn, dI_sm, dI_ss,dI_smr, dI_ssr, dI_smis, dI_smisn, dI_ssis, dI_ssisn, dI_ssh, dI_smrisn, dI_ssrisn, dI_smqisn, dR, dD)))
  })
}# 
# ### Set time frame in days



#################################################################################################
# Simulation
# Create list to store model output
listOut <- NULL
listOut <- list()
i=1
for(i in seq(1, nSim, 1)){

  ### Set time frame in days
  t <- seq(0, time_stuff$t_n[i], by = 1)
  
  curr_param <- all_param %>% filter(isim == i)
  
  params <- curr_param$value
  names(params) <- curr_param$name
  
  
  ##########################################
  #
  #
  # This is an example of how to dynamically set params based on conditions for SEIR state
  # It is not elegant and it will override the input file but it should work
  #
  print(paste0("i = ", i))
  print(paste0("beta First = ", params[names(params) == "beta"]))
  
  # *********************************************
  #    ***** EXAMPLE DYNAMIC VARIABLES ****
  params[names(params) == "beta"]  <- 0.05 / init[names(init) == "S"] 
  
  
  print(paste0("beta after = ", params[names(params) == "beta"]))
  
  
  #####################################
  #all(params == old_params)
  #all(names(params) == names(old_params))
  
  #all(init == old_init)
  #all(names(init) == names(old_init))
  #all(t == old_t)
  ### Run the model
  out <- ode(y = init, 
             times = t, 
             func = seir, 
             parms = params
             )
  out <- out %>% 
    as_tibble()
  
  out$time <- seq(time_stuff$tmin[i], time_stuff$tmax[i], 1)
  
  # Add model output to list
  listOut[[i]] <- out
  
  out_for_init <- 
    out %>% 
    slice(nrow(out)) %>% 
    pivot_longer(-time) 
  
  
    
  init <- out_for_init$value
  names(init) <- out_for_init$name
}

big_out <- bind_rows(listOut) %>% distinct()

big_out$L_tot <- big_out %>% select_at(vars(starts_with("L"))) %>% rowSums()
big_out$I_tot <- big_out %>% select_at(vars(starts_with("I"))) %>% rowSums()


big_out %>% mutate(I_old = lag(I_tot))

#############################################
# MAKE PRETTY PICTURES



# Time lines graph

df_p <- 
  big_out %>% 
  select(time, S, L_tot, I_tot, R, D) %>% 
  rename( susceptible := S, Latent := L_tot, Infected = I_tot, Recovered = R, Dead = D) %>% 
  pivot_longer(-time) 

df_p_l <-df_p %>% group_by(name) %>% summarise(max_value =max(value), 
                                               min_value =min(value),
                                               )  

df_p_l$time <- 750
df_p_l <- df_p_l %>% mutate(value = (max_value - min_value) + min_value)
df_p_l <- df_p_l %>% filter(name == "susceptible")
  #big_out %>% slice(which.max(I_tot)) 

df_p %>% 
  ggplot(aes(x = time, color = name, y = value)) +
  geom_point() + 
#  geom_label(data = df_p_l, mapping = aes(label = as.character(max_value))) +
  facet_grid(rows = vars(name), scales = "free_y") +
  labs(title = "Covid-19 - SEIR Model Results", x = "Time [days]", y = "Counts" )




################################ HAVE NOT TESTED PAST HERE ###############################33


timeout
out$time
out[[, "time"]]
plot(out$S)
#################################################################################################

# Plot the results


## Plot


plotdata <- data.frame(x=t, y=out)

matplot(x = plotdata[,c("y.time")], y = plotdata[,c("y.S","y.L","y.L_q","y.I_a","y.I_aq","y.I_smis","y.I_ssh","y.R",'y.D')], type = "l",
        xlab = "Time", ylab = "individuals", main = "SEIR Model 2.0",
        lwd = 1, lty = 1, bty = "l", col = 2:11)

legend(250, 1000000, c("Susceptible", "Latent","Latent_q",  "Pre-sympt", 'Pre-sympto in quarantine','Symptomatic mild isolated','Symptomatic severe hospitalized','Recovered','Dead'), pch = 1, col = 2:11, bty = "n")



###################################################################################
##Calculate ouputs#######################
#############################################



  ### Number of infected over time
  ################################

plotdata$Infect <- plotdata$y.I_a+plotdata$y.I_aq+plotdata$y.I_ar + plotdata$y.I_aqn+plotdata$y.I_sm+plotdata$y.I_ss+plotdata$y.I_smis+plotdata$y.I_smr + plotdata$y.I_ssr + plotdata$y.I_smisn+plotdata$y.I_ssis+plotdata$y.I_ssisn+plotdata$y.I_ssh+ plotdata$y.I_smqisn+plotdata$y.I_smrisn+plotdata$y.I_ssrisn

matplot(x = plotdata[,c("y.time")], y = plotdata[,c('Infect')], type = "l",
        xlab = "Time", ylab = "individuals", main = "Infected",
        lwd = 2, lty = 1, bty = "l", col = 2)

#legend(250, 460000, c("Infected"), pch = 1, col = 2, bty = "n")

  # Max number of cases and date of peak
maxInfect<-max(plotdata$Infect)
maxInfect
datemaxInfect<-plotdata$y.time[plotdata$Infect==max(plotdata$Infect)]
datemaxInfect



  ### Number of new cases each day, over time (incidence)
  #######################################################

n<-nrow(plotdata)

#Calculate the incidence



plotdata$IncI <-NA
plotdata$IncI[1]<-plotdata$Infect[1]
for(i in 2:n)
  plotdata$IncI[i] <- ((plotdata$y.L[i-1]+plotdata$y.L_q[i-1]+plotdata$y.L_r[i-1]) * rand_sigma ) 


#Calculate the incidence max and the related date
maxInc<-max(plotdata$IncI)
maxInc
datemaxInc<-plotdata$y.time[plotdata$IncI==max(plotdata$IncI)]
datemaxInc

# Total number of cases
sumInfect<-sum(plotdata$IncI)
sumInfect

# Outbreak duration (when no new cases)

plotdata$cumI<-NA
plotdata$cumI[1]<-plotdata$IncI[1]
for(i in 2:n)
  plotdata$cumI[i] <- (plotdata$IncI[i] +plotdata$cumI[i-1]) 

matplot(x = plotdata[,c("y.time")], y = plotdata[,c('cumI')], type = "l",
        xlab = "Time", ylab = "individuals", main = "Cumulative incidence",
        lwd = 2, lty = 1, bty = "l", col = 6)

datemaxCumI<-plotdata$y.time[plotdata$cumI==max(plotdata$cumI)]
min(datemaxCumI)

# Attack rate
AR <- sumInfect*100/N
AR



### Hospitalized

plotdata$Hosp <- plotdata$y.I_ssh
matplot(x = plotdata[,c("y.time")], y = plotdata[,c('Hosp')], type = "l",
       xlab = "Time", ylab = "individuals", main = "Hospitalized",
       lwd = 2, lty = 1, bty = "l", col = 3)

#legend(250, 1000000, c("Hospitalized"), pch = 1, col = 3, bty = "n")


maxHosp<-max(plotdata$Hosp)
maxHosp
datemaxHosp<-plotdata$y.time[plotdata$Hosp==max(plotdata$Hosp)]
datemaxHosp

#### Quarantined

plotdata$Quarant <- plotdata$y.I_aq

matplot(x = plotdata[,c("y.time")], y = plotdata[,c('Quarant')], type = "l",
        xlab = "Time", ylab = "individuals", main = "Quarantined",
        lwd = 2, lty = 1, bty = "l", col = 4)

#legend(250, 1000000, c("Hospitalized"), pch = 1, col = 3, bty = "n")


maxQuarant<-max(plotdata$Quarant)
maxQuarant
datemaxQuarant<-plotdata$y.time[plotdata$Quarant==max(plotdata$Quarant)]
datemaxQuarant

 
##### Isolated
plotdata$Isolat <- plotdata$y.I_smis + plotdata$y.I_ssis

matplot(x = plotdata[,c("y.time")], y = plotdata[,c('Isolat')], type = "l",
        xlab = "Time", ylab = "individuals", main = "Isolated",
        lwd = 2, lty = 1, bty = "l", col = 5)


maxIsolat<-max(plotdata$Isolat)
maxIsolat
datemaxIsolat<-plotdata$y.time[plotdata$Isolat==max(plotdata$Isolat)]
datemaxIsolat
















##################################################################################################

#### Calculate R0 or Rt Reproductive number

# dI <- sigma*rho*L - (delta_I*(1-omega) +  gamma_I*omega)*I + mu*epsilon*A
# dA <- sigma*(1 - rho)*L - (1-epsilon)*gamma_A*A - mu*epsilon*A
###################################################################################################
epsilon  <-params['epsilon']
omega  <-params['omega']
plotdata <- data.frame(x=t, y=out)
n<-nrow(plotdata)

plotdata$Rt <-NA
plotdata$Rt[1]<-plotdata$y.I[1]
for(i in 2:n)
  plotdata$Rt[i] <- (plotdata$y.L[i-1]*rand_sigma*rand_rho -(rand_delta_I*(1-omega)+rand_gamma_I*omega)*plotdata$y.I[i-1]+plotdata$y.L[i-1]*rand_sigma*(1-rand_rho)-rand_gamma_A*(1-epsilon)*plotdata$y.A[i-1])

matplot(x = plotdata[,c("y.time")], y = plotdata[,c("Rt")], type = "l",
        xlab = "Time", ylab = "Rt", main = "Rt",
        lwd = 1, lty = 1, bty = "l", col = 2:11)

###################################################################
### for modelling report March 3rd 2020############################
###################################################################

# number of cases

  #Calculate the max number of cases

maxIinf <-max(plotdata$y.I)
maxIinf
  #Plot the number of cases

matplot(x = plotdata[,c("y.time")], y = plotdata[,c("y.I")], type = "l",
        xlab = "Time", ylab = "individuals", main = "Number of cases",
        lwd = 2, lty = 1, bty = "l", col = 2)


#Incidence
  
#Calculate the incidence
epsilon  <-params['epsilon']
plotdata <- data.frame(x=t, y=out)
n<-nrow(plotdata)

plotdata$IncI <-NA
plotdata$IncI[1]<-plotdata$y.I[1]
for(i in 2:n)
  plotdata$IncI[i] <- (plotdata$y.L[i-1]*rand_sigma*rand_rho + rand_mu*epsilon*plotdata$y.A[i-1]) 

  #Calculate the incidence max and the related date
maxInc<-max(plotdata$IncI)
maxInc
datemaxInc<-plotdata$y.time[plotdata$IncI==max(plotdata$IncI)]
datemaxInc



  #Plot the incidence

matplot(x = plotdata[,c("y.time")], y = plotdata[,c("IncI")], type = "l",
        xlab = "Time", ylab = "individuals", main = "Incidence",
        lwd = 2, lty = 1, bty = "l", col = 3)


#Cumulative incidence

plotdata$cumI<-NA
plotdata$cumI[1]<-plotdata$y.I[1]
for(i in 2:n)
  plotdata$cumI[i] <- (plotdata$y.L[i-1]*rand_sigma*rand_rho + rand_mu*epsilon*plotdata$y.A[i-1]+plotdata$cumI[i-1]) 

  #Calculate the max cumulative incidence and its related date (epidemic duration)
maxCumI<-max(plotdata$cumI)
maxCumI
datemaxCumI<-plotdata$y.time[plotdata$cumI==max(plotdata$cumI)]
min(datemaxCumI)

  #Plot the cumulative incidence

matplot(x = plotdata[,c("y.time")], y = plotdata[,c("cumI")], type = "l",
        xlab = "Time", ylab = "individuals", main = "Cumulative incidence",
        lwd = 2, lty = 1, bty = "l", col = 4)












###################################################################
# Validate the results with Wuhan data
###################################################################


#Comparison with First stages of the epidemic (Li, 2020)

setwd('D:/PED/Coronavirus/SEIR model/validation')
d<-read.table("WuhanForValidation.txt",header = T,fill=TRUE)


plotdata <- data.frame(x=t, y=out)


#New I cases at each time step
#dI <- sigma*rho*L - (delta_I*(1-omega) +  gamma_I*omega)*I + mu*epsilon*A
epsilon  <-params[3]

n<-nrow(plotdata)

plotdata <- data.frame(x=t, y=out)
plotdata$IncI <-NA
plotdata$IncI[1]<-plotdata$y.I[1]
for(i in 2:n)
  plotdata$IncI[i] <- (plotdata$y.L[i-1]*rand_sigma*rand_rho + rand_mu*epsilon*plotdata$y.A[i-1]) 


plotdata <-data.frame(plotdata$y.time,plotdata$IncI)
names(plotdata)[names(plotdata) == "plotdata.y.time"] <- "time"
names(plotdata)[names(plotdata) == "plotdata.IncI"] <- "IncI"


d2 <- data.frame(d)
plotdata2 <-merge (plotdata,d2, by='time', all = TRUE)

matplot(x = plotdata2[,c("time")], y = plotdata2[,c("IncI","casesL")], type = "l",
        xlab = "Time", ylab = "individuals", main = "Model comparison",
        lwd = 1, lty = 1, bty = "l", col = 2:11)

legend(250, 400000, c("Incidence rate of Inf","cases obs"), pch = 1, col = 2:11, bty = "n")

# zoon=m on the first 100 days

plotdata2zoom <-subset(plotdata2,plotdata2$time<45)

matplot(x = plotdata2zoom[,c("time")], y = plotdata2zoom[,c("IncI","casesL")], type = "l",
        xlab = "Time", ylab = "individuals", main = "Model comparison",
        lwd = 1, lty = 1, bty = "l", col = 2:11)












#plot for modeling report Feb27th 2020

datareport12 <- data.frame(x=t, y=out)#change the name of the dataset to reflect the value of the c parameter in the name
datareport12 <-data.frame(datareport12$y.time,datareport12$y.I,datareport12$y.A)
n<-nrow(datareport12)
datareport12$TotInf<-NA
for(i in 1:n)
  datareport12$TotInf[i] <- (datareport12$datareport12.y.I[i] + datareport12$datareport12.y.A[i]) 
names(datareport12)[names(datareport12) == "TotInf"] <- "TotInf12"
names(datareport12)[names(datareport12) == "datareport12.y.I"] <- "Inf12"
names(datareport12)[names(datareport12) == "datareport12.y.A"] <- "Asympt12"


datareport10 <- data.frame(x=t, y=out) 
n<-nrow(datareport10)
datareport10$TotInf<-NA
for(i in 1:n)
  datareport10$TotInf[i] <- (datareport10$y.I[i] + datareport10$y.A[i]) 


datareport8 <- data.frame(x=t, y=out) 
n<-nrow(datareport8)
datareport8$TotInf<-NA
for(i in 1:n)
  datareport8$TotInf[i] <- (datareport8$y.I[i] + datareport8$y.A[i]) 


datareport6 <- data.frame(x=t, y=out) 
n<-nrow(datareport6)
datareport6$TotInf<-NA
for(i in 1:n)
  datareport6$TotInf[i] <- (datareport6$y.I[i] + datareport6$y.A[i]) 

datareporttot <- merge(datareport10, datareport12, by='y.time') 
datareporttot <- merge(datareport8, datareporttot, by='y.time') 
datareporttot <- merge(datareport6, datareporttot, by='y.time') 

head(datareporttot)


matplot(x = datareporttot[,c("y.time")], y = datareporttot[,c("TotInf","TotInf.x","TotInf.y")], type = "l",
        xlab = "Time", ylab = "total number of infected", main = "Number of infected (asymptomatic and symptomatic)",
        lwd = 1, lty = 1, bty = "l", col = 2:11)

legend(450, 4000, c("contact rate=8","contact rate=10", "contact rate = 12"), pch = 1, col = 2:11, bty = "n")


######################### Output statistics

#Maximum number of cases
maxI<-max(plotdata$y.I)
maxI
datemaxI<-plotdata$y.time[plotdata$y.I==max(plotdata$y.I)]
datemaxI

maxA<-max(plotdata$y.A)
maxA
datemaxA<-plotdata$y.time[plotdata$y.A==max(plotdata$y.A)]
datemaxA

maxH<-max(plotdata$y.H)
maxH
datemaxH<-plotdata$y.time[plotdata$y.H==max(plotdata$y.H)]
datemaxH

maxR<-max(plotdata$y.R)
maxR

maxD<-max(plotdata$y.D)
maxD


# Cumulative number of infected non-infectious (A)
  # dA <- sigma*(1 - rho)*L - (1-epsilon)*gamma_A*A - mu*epsilon*A

n<-nrow(plotdata)

plotdata$cumA<-NA
plotdata$cumA[1]<-plotdata$y.A[1]
for(i in 2:n)
  plotdata$cumA[i] <- (plotdata$y.L[i-1]*rand_sigma*(1-rand_rho)) 
finalcumA <-sum(plotdata$cumA) #certainement inutil de faire une somme final: le dernier élémentd eligne devrait déjà donner le chiffre recherché


#Cumulative number of infected infectious (I)
  #dI <- sigma*rho*L - (delta_I*(1-omega) +  gamma_I*omega)*I + mu*epsilon*A
epsilon  <-params[3]
plotdata$cumI<-NA
plotdata$cumI[1]<-plotdata$y.I[1]
for(i in 2:n)
  plotdata$cumI[i] <- (plotdata$y.L[i-1]*rand_sigma*rand_rho + rand_mu*epsilon*plotdata$y.A[i-1]) 
finalcumI <-sum(plotdata$cumI)

#Cumulative number of Latent quarantined L_q)
  # dL_q <-lambda*beta*c*S*(I+A) - L_q*sigma
lambda  <-params[4]
plotdata$cumL_q<-NA
plotdata$cumL_q[1]<-plotdata$y.L_q[1]
for(i in 2:n)
  plotdata$cumL_q[i] <- (lambda*rand_beta*rand_c*plotdata$y.S[i-1]*(plotdata$y.A[i-1]+plotdata$y.I[i-1]))
finalcumL_q <-sum(plotdata$cumL_q)

#Cumulative number of recovered (R)

#Cumulative number of Hospitalized-home (H-H)

#Cumulative number of dead (D)

# Epidemic duration (Between the firts case (I or A) and the last case lreleased by the H-H compartiment

# Date of the maximum number of cases (I)

#Maximum number of cases (I)


# R(t) or R(0)




######################################################################################################

# validation

#calculation of the number of hopitalised that are dying at each time step and then we compare
#to the final number of dead peopple. The difference should be 0.

n<-nrow(plotdata)

plotdata$test<-NA
for(i in 1:n)
  plotdata$test[i] <- (plotdata$y.H[i]*0.046)

finalstat <-sum(plotdata$test)
plotdata$y.D[365]-finalstat





#finalcumI <-sum(plotdata$cumI)

d<-read.table("WuhanFirstPart.txt",header = T)


plotdata <- data.frame(x=t, y=out)

plotdata <-data.frame(plotdata$y.time,plotdata$y.I,plotdata$y.A)
names(plotdata)[names(plotdata) == "plotdata.y.time"] <- "time"
names(plotdata)[names(plotdata) == "plotdata.y.I"] <- "Inf"
names(plotdata)[names(plotdata) == "plotdata.y.A"] <- "Asym"

d2 <- data.frame(d)
plotdata2 <-merge (plotdata,d2, by='time', all = TRUE)

matplot(x = plotdata2[,c("time")], y = plotdata2[,c("Inf","Asym","casesL")], type = "l",
        xlab = "Time", ylab = "individuals", main = "Model comparison",
        lwd = 1, lty = 1, bty = "l", col = 2:11)

legend(250, 400000, c("Inf sim", "Asym sim","cases obs"), pch = 1, col = 2:11, bty = "n")

