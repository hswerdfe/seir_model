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

df_p <- 
  big_out %>% 
  select(time, S, L_tot, I_tot, R, D) %>% 
  rename( susceptible := S, Latent := L_tot, Infected = I_tot, Recovered = R, Dead = D) %>% 
  pivot_longer(-time) 

#df_p_L %>% 


df_p_l <-df_p %>% group_by(name) %>% summarise(max_value =max(value), 
                                               min_value =min(value),
                                               )  

df_p_l$time <- 750
df_p_l$value <- (max_value - min_value) + min_value
df_p_l <- df_p_l %>% filter(name = "susceptible")
  #big_out %>% slice(which.max(I_tot)) 

df_p %>% 
  ggplot(aes(x = time, color = name, y = value)) +
  geom_point() + 
#  geom_label(data = df_p_l, mapping = aes(label = as.character(max_value))) +
  facet_grid(rows = vars(name), scales = "free_y") +
  labs(title = "Covid-19 - SEIR Model Results", x = "Time [days]", y = "Counts" )
