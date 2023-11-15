

### Mathematical modeling to support the fight against ###
### antibiotic-resistant bacteria in hospital settings in Ghana ###





################################################################################
#################### Single resistance model (XSR1 model) ######################
################################################################################

###   Creating the Model

#Run the deSolve package----
library(deSolve)

#ODE of the model will be solved in 11 steps----
#step1-create the model----
single.model <- function(t,x,params){
  
  
  #step2- List the state variables 
  
  X <- x[1]
  S <- x[2] 
  R1 <- x[3]
  
  
  N <- X + S + R1 
  
  
  #step3- list of parameters----
  
  
  beta <- params["beta"]
  gamma <- params["gamma"]
  sigma <- params["sigma"]
  mu <- params["mu"]
  tau1 <- params["tau1"]
  m <- params["m"]
  m1 <- params["m1"]
  c1 <- params["c1"]
  
  
  #step4- code the model equations----
  dX <-  mu*(1-m-m1-X) -beta*S*X + (tau1+gamma)*S -beta*(1-c1)*R1*X +  gamma*R1 
  dS <-  mu*(m-S) + beta*S*X - (tau1+gamma)*S - sigma*beta*(1-c1)*R1*S
  dR1 <- mu*(m1-R1)+ sigma*beta*(1-c1)*R1*S + beta*(1-c1)*R1*X - gamma*R1  
  
  #step5- combine results into a single vector----
  dN <- c(dX, dS, dR1)
  
  #step6- return as a list----
  list(dN)}

#step7- assign values to the parameters----
parameters <- c(beta= 0.0007,# Estimated value
                
                gamma=  1/7,# gamma value from https://doi.org/10.1371/journal.pcbi.1004340
                
                mu=   0.0001,# Estimated value  / #0.003,# mu value from https://doi.org/10.1371/journal.pcbi.1005745
                
                m=   0.07,# m value from https://doi.org/10.1371/journal.pcbi.1005745
                
                m1=   0.001,# m1 value from https://doi.org/10.1371/journal.pcbi.1005745
                
                sigma= 0.25, # sigma value from https://doi.org/10.1080/17513758.2010.488300 AND https://doi.org/10.1073/pnas.0402298101
                
                tau1= seq(0.1, 1), # tau1 value from https://doi.org/10.1371/journal.pcbi.1005745 
                
                c1= 0.05) # c1 value from https://doi.org/10.1080/17513758.2010.488300 AND https://doi.org/10.1073/pnas.97.4.1938

#step8- state the times---- 
times <- seq(from=0, to=150, by=1)

#step9- state the initial conditions----
xstart<- c(X= 720, S= 79, R1= 1)

#step10- compute a model trajectory----
library(tidyverse)

ode(
  func=single.model,
  y=xstart,
  times=times,
  parms=parameters
)%>%
  as.data.frame() -> dat


#step11- plot the results----
plot_dat <- dat %>%
  gather(Compartments, value,-time)%>%
  ggplot(aes(x=time,y=value,color=Compartments, group=Compartments))+
  geom_line(lwd=1.5) +
  # theme_classic()+#
  labs(x='Time(Days)', y='Number of individuals') + theme_bw()


#show normal plot (with ggplot2)----
plot_dat

#Step12- Interactive plot with plotly----
library(plotly)

#Show plot with plotly----
ggplotly(plot_dat)






################################################################################
#########    Model calibration and testing parameter combinations   ############
################################################################################

###
get_R1 <- function(times=times,odefunction=single.model,xstart=xstart,parameters){
  out <- ode(  
    func=odefunction,
    y=xstart,
    times=times,
    parms=parameters)%>%
    as.data.frame()
  R1 <- max(out$R1)
  return(R1)}

###
parms_multi  <- tidyr::expand_grid (beta= 0.0007,
                                    gamma=1/7,
                                    mu= seq(0.0001, 0.1001,0.05),# ok ini used in dynamic model
                                    m=  seq(0.01,0.7,0.05),# ok see report
                                    m1= seq(0.001,0.1,0.02),# A reprendre uniquement 0.001 sur plot
                                    #BUT currently using final value of 0.01 (700 combis) 
                                    #TRY 0.1 later and see (7000 combis)
                                    sigma= 0.25,
                                    tau1= seq(0.1, 1, 0.2),#ok all values from lit
                                    c1= 0.05) 

###
# Result max R1 for a single run
b <- get_R1(times = times,   odefunction = single.model,   xstart = xstart, parameters = parms_multi)


# Result max R1 for multiple combinations
p1 <- parms_multi
p1$R1 <- -999

for (i in 1:nrow(p1)){
  print(i)
  p1$R1[i] <- get_R1(times = times,odefunction = single.model ,xstart = xstart, parameters = p1[i,1:9]) }


#Save results as "p1f.csv"
#Save P1 as a n Excel dataset

write.csv(p1f, file= p1f.csv)
write.table(p1f, file="p1f.csv", sep=",")



#Import P1.csv dataset

library(readxl)
p1f <- read_excel("~/Gnimatin_FLI_15 August/0_Many parametersa/SIngle with dynamics/p1f.xlsx")
View(p1f)




#Plots
library(ggplot2)


#Type 1
ggplot(p1f, aes(mu, m, fill= R1)) + geom_tile()

#Type 2
ggplot(p1f, aes(mu, m, fill= R1)) +  geom_tile()+facet_grid(m1~.)

#Type 3
ggplot(p1f, aes(mu, m, fill= R1)) + geom_tile()+facet_grid(.~tau1)

#Type 4
plot4<-ggplot(p1f, aes(mu, m, fill= R1)) + geom_tile()+facet_grid(m1~tau1)


#
plot4

#Step12- Interactive plot with plotly----
library(plotly)

#Show plot with plotly----
ggplotly(plot4)




################################################################################
#########################   Reproductive numbers   #############################
################################################################################

#parameters

beta <- 0.0007

tau1 <- 0.1

mu <- 0.0001

gamma <- 1/7

c1 <- 0.05 


### Real-time effective reproductive number RR1
################################################################################


RR1 <- (beta* (1 - c1) * dat$X)/ (mu+gamma)


### Plot RR1
# Libraries
library(ggplot2)


# create data
RR1 <- RR1
Time <- times
data <- data.frame(x=Time, y=RR1)

# Plot
plot_RR1 <-ggplot(data, aes(x=Time, y=RR1)) + geom_line() + theme_bw()



# Interactive plot with plotly----
library(plotly)

#Show plot with plotly----
ggplotly(plot_RR1)




### Real-time effective reproductive number RS
###############################################################################

RS <- (beta * dat$X)/ (tau1 + mu + gamma) 

### Plot RS
# Libraries
library(ggplot2)

# create data
RS <- RS
Time <- times
data <- data.frame(x=Time, y=RS)

# Plot
plot_RS <-ggplot(data, aes(x=Time, y=RS)) + geom_line() + theme_bw()

# Interactive plot with plotly----
library(plotly)

#Show plot with plotly----
ggplotly(plot_RS)





## Plot  Reproductive numbers together
################################################################################
databasicrepro <- data.frame(timeaxis=Time, RoS= RS,RoR1= RR1)

#show dataframe
head(databasicrepro)

#Plot
library(ggplot2)
library(tidyverse)

plot_databasicrepro <-databasicrepro %>%  gather(Legend,value,-timeaxis)%>% ggplot(aes(x=timeaxis,y=value,color=Legend))+geom_line(lwd=1.5) +
  labs(x='Time(Days)', y='Reproductive Numbers') + theme_bw()

# Interactive plot with plotly----
library(plotly)

#Show plot with plotly----
ggplotly(plot_databasicrepro)





################################################################################
#########################   Varying parameter beta   ###########################
################################################################################

#Run the deSolve package----
library(deSolve)


#ODE of the model will be solved in 11 steps----
#step1-create the model----
single.modelbeta <- function(t,x,params){
  #step2- List the state Legends 
  X <- x[1]
  S <- x[2] 
  R1 <- x[3]
  N <- X + S + R1 
  #step3- list of parameters----
  beta <- params["beta"]
  gamma <- params["gamma"]
  sigma <- params["sigma"]
  mu <- params["mu"]
  tau1 <- params["tau1"]
  m <- params["m"]
  m1 <- params["m1"]
  c1 <- params["c1"]
  #step4- code the model equations----
  dX <-  mu*(1-m-m1-X) -beta*S*X + (tau1+gamma)*S -beta*(1-c1)*R1*X +  gamma*R1 
  dS <-  mu*(m-S) + beta*S*X - (tau1+gamma)*S - sigma*beta*(1-c1)*R1*S
  dR1 <- mu*(m1-R1)+ sigma*beta*(1-c1)*R1*S + beta*(1-c1)*R1*X - gamma*R1  
  #step5- combine results into a single vector----
  dN <- c(dX, dS, dR1)
  #step6- return as a list----
  list(dN)}
#step7- assign values to the parameters---
parametersb <- c(beta= 0.0003,# Estimated value
                 
                 gamma=  1/7,# gamma value from https://doi.org/10.1371/journal.pcbi.1004340
                 
                 mu=   0.0001,# Estimated value  / #0.003,# mu value from https://doi.org/10.1371/journal.pcbi.1005745
                 
                 m=   0.07,# m value from https://doi.org/10.1371/journal.pcbi.1005745
                 
                 m1=   0.001,# m1 value from https://doi.org/10.1371/journal.pcbi.1005745
                 
                 sigma= 0.25, # sigma value from https://doi.org/10.1080/17513758.2010.488300 AND https://doi.org/10.1073/pnas.0402298101
                 
                 tau1= seq(0.1, 1), # tau1 value from https://doi.org/10.1371/journal.pcbi.1005745 
                 
                 c1= 0.05) # 

###
parametersc <- c(beta= 0.0005,# Estimated value
                 
                 gamma=  1/7,# gamma value from https://doi.org/10.1371/journal.pcbi.1004340
                 
                 mu=   0.0001,# Estimated value  / #0.003,# mu value from https://doi.org/10.1371/journal.pcbi.1005745
                 
                 m=   0.07,# m value from https://doi.org/10.1371/journal.pcbi.1005745
                 
                 m1=   0.001,# m1 value from https://doi.org/10.1371/journal.pcbi.1005745
                 
                 sigma= 0.25, # sigma value from https://doi.org/10.1080/17513758.2010.488300 AND https://doi.org/10.1073/pnas.0402298101
                 
                 tau1= seq(0.1, 1), # tau1 value from https://doi.org/10.1371/journal.pcbi.1005745 
                 
                 c1= 0.05) # 

###
parametersd <- c(beta= 0.0007,# Estimated value
                 
                 gamma=  1/7,# gamma value from https://doi.org/10.1371/journal.pcbi.1004340
                 
                 mu=   0.0001,# Estimated value  / #0.003,# mu value from https://doi.org/10.1371/journal.pcbi.1005745
                 
                 m=   0.07,# m value from https://doi.org/10.1371/journal.pcbi.1005745
                 
                 m1=   0.001,# m1 value from https://doi.org/10.1371/journal.pcbi.1005745
                 
                 sigma= 0.25, # sigma value from https://doi.org/10.1080/17513758.2010.488300 AND https://doi.org/10.1073/pnas.0402298101
                 
                 tau1= seq(0.1, 1), # tau1 value from https://doi.org/10.1371/journal.pcbi.1005745 
                 
                 c1= 0.05) # 


###
parameterse <- c(beta= 0.0009,# Estimated value
                 
                 gamma=  1/7,# gamma value from https://doi.org/10.1371/journal.pcbi.1004340
                 
                 mu=   0.0001,# Estimated value  / #0.003,# mu value from https://doi.org/10.1371/journal.pcbi.1005745
                 
                 m=   0.07,# m value from https://doi.org/10.1371/journal.pcbi.1005745
                 
                 m1=   0.001,# m1 value from https://doi.org/10.1371/journal.pcbi.1005745
                 
                 sigma= 0.25, # sigma value from https://doi.org/10.1080/17513758.2010.488300 AND https://doi.org/10.1073/pnas.0402298101
                 
                 tau1= seq(0.1, 1), # tau1 value from https://doi.org/10.1371/journal.pcbi.1005745 
                 
                 c1= 0.05) # 

###
parametersf <- c(beta= 0.0011,# Estimated value
                 
                 gamma=  1/7,# gamma value from https://doi.org/10.1371/journal.pcbi.1004340
                 
                 mu=   0.0001,# Estimated value  / #0.003,# mu value from https://doi.org/10.1371/journal.pcbi.1005745
                 
                 m=   0.07,# m value from https://doi.org/10.1371/journal.pcbi.1005745
                 
                 m1=   0.001,# m1 value from https://doi.org/10.1371/journal.pcbi.1005745
                 
                 sigma= 0.25, # sigma value from https://doi.org/10.1080/17513758.2010.488300 AND https://doi.org/10.1073/pnas.0402298101
                 
                 tau1= seq(0.1, 1), # tau1 value from https://doi.org/10.1371/journal.pcbi.1005745 
                 
                 c1= 0.05) # 


#step8- state the times---- 
times <- seq(from=0, to=150, by=1)###
#step9- state the initial conditions----
xstart<- c(X= 720, S= 79, R1= 1)


###Plot b

#step10- compute a model trajectory----
library(tidyverse)
ode(
  func=single.modelbeta,
  y=xstart,
  times=times,
  parms=parametersb
)%>%
  as.data.frame() -> datb
#step11- plot the results----
plot_datb <- datb %>%
  gather(Legend, value,-time)%>%
  ggplot(aes(x=time,y=value,color=Legend, group=Legend))+
  geom_line(lwd=1.5) +
  # theme_classic()+#
  labs(x='Time(Days)', y='Number of individuals')
#show normal plot (with ggplot2)----
plot_datb
#Step12- Interactive plot with plotly----
library(plotly)
#Show plot with plotly----
betaplotb <-ggplotly(plot_datb)

###Plot c

#step10- compute a model trajectory----
library(tidyverse)
ode(
  func=single.modelbeta,
  y=xstart,
  times=times,
  parms=parametersc
)%>%
  as.data.frame() -> datc
#step11- plot the results----
plot_datc <- datc %>%
  gather(Legend, value,-time)%>%
  ggplot(aes(x=time,y=value,color=Legend, group=Legend))+
  geom_line(lwd=1.5) +
  # theme_classic()+#
  labs(x='Time(Days)', y='Number of individuals')
#show normal plot (with ggplot2)----
plot_datc
#Step12- Interactive plot with plotly----
library(plotly)
#Show plot with plotly----
betaplotc <-ggplotly(plot_datc)


### Plot d   

#step10- compute a model trajectory----
library(tidyverse)
ode(
  func=single.modelbeta,
  y=xstart,
  times=times,
  parms=parametersd
)%>%
  as.data.frame() -> datd
#step11- plot the results----
plot_datd <- datd %>%
  gather(Legend, value,-time)%>%
  ggplot(aes(x=time,y=value,color=Legend, group=Legend))+
  geom_line(lwd=1.5) +
  # theme_classic()+#
  labs(x='Time(Days)', y='Number of individuals')
#show normal plot (with ggplot2)----
plot_datd
#Step12- Interactive plot with plotly----
library(plotly)
#Show plot with plotly----
betaplotd <-ggplotly(plot_datd)


###   Plot e  

#step10- compute a model trajectory----
library(tidyverse)
ode(
  func=single.modelbeta,
  y=xstart,
  times=times,
  parms=parameterse
)%>%
  as.data.frame() -> date
#step11- plot the results----
plot_date <- date %>%
  gather(Legend, value,-time)%>%
  ggplot(aes(x=time,y=value,color=Legend, group=Legend))+
  geom_line(lwd=1.5) +
  # theme_classic()+#
  labs(x='Time(Days)', y='Number of individuals')
#show normal plot (with ggplot2)----
plot_date
#Step12- Interactive plot with plotly----
library(plotly)
#Show plot with plotly----
betaplote <-ggplotly(plot_date)



### Plot f

#step10- compute a model trajectory----
library(tidyverse)
ode(
  func=single.modelbeta,
  y=xstart,
  times=times,
  parms=parametersf
)%>%
  as.data.frame() -> datf
#step11- plot the results----
plot_datf <- datf %>%
  gather(Legend, value,-time)%>%
  ggplot(aes(x=time,y=value,color=Legend, group=Legend))+
  geom_line(lwd=1.5) +
  # theme_classic()+#
  labs(x='Time(Days)', y='Number of individuals')
#show normal plot (with ggplot2)----
plot_datf
#Step12- Interactive plot with plotly----
library(plotly)
#Show plot with plotly----
betaplotf <- ggplotly(plot_datf)


### Arrange Combined plots 
################################################################################
library(ggpubr)
ggarrange(plot_datb,plot_datc,plot_datd,plot_date,plot_datf, ncol = 3, nrow = 2)







################################################################################
#########################   Varying parameter mu   #############################
################################################################################

#Run the deSolve package----
library(deSolve)


#ODE of the model will be solved in 11 steps----
#step1-create the model----
single.modelmu <- function(t,x,params){
  #step2- List the state Legends 
  X <- x[1]
  S <- x[2] 
  R1 <- x[3]
  N <- X + S + R1 
  #step3- list of parameters----
  beta <- params["beta"]
  gamma <- params["gamma"]
  sigma <- params["sigma"]
  mu <- params["mu"]
  tau1 <- params["tau1"]
  m <- params["m"]
  m1 <- params["m1"]
  c1 <- params["c1"]
  #step4- code the model equations----
  dX <-  mu*(1-m-m1-X) -beta*S*X + (tau1+gamma)*S -beta*(1-c1)*R1*X +  gamma*R1 
  dS <-  mu*(m-S) + beta*S*X - (tau1+gamma)*S - sigma*beta*(1-c1)*R1*S
  dR1 <- mu*(m1-R1)+ sigma*beta*(1-c1)*R1*S + beta*(1-c1)*R1*X - gamma*R1  
  #step5- combine results into a single vector----
  dN <- c(dX, dS, dR1)
  #step6- return as a list----
  list(dN)}
#step7- assign values to the parameters---
###
parametersa <- c(beta= 0.0007,# Estimated value
                 
                 gamma=  1/7,# gamma value from https://doi.org/10.1371/journal.pcbi.1004340
                 
                 mu=   0.0001,# Estimated value  / #0.003,# mu value from https://doi.org/10.1371/journal.pcbi.1005745
                 
                 m=   0.07,# m value from https://doi.org/10.1371/journal.pcbi.1005745
                 
                 m1=   0.001,# m1 value from https://doi.org/10.1371/journal.pcbi.1005745
                 
                 sigma= 0.25, # sigma value from https://doi.org/10.1080/17513758.2010.488300 AND https://doi.org/10.1073/pnas.0402298101
                 
                 tau1= seq(0.1, 1), # tau1 value from https://doi.org/10.1371/journal.pcbi.1005745 
                 
                 c1= 0.05) # c1 value from https://doi.org/10.1080/17513758.2010.488300 AND https://doi.org/10.1073/pnas.97.4.1938

###
parametersb <- c(beta= 0.0007,# Estimated value
                 
                 gamma=  1/7,# gamma value from https://doi.org/10.1371/journal.pcbi.1004340
                 
                 mu=   0.0005,# Estimated value  / #0.003,# mu value from https://doi.org/10.1371/journal.pcbi.1005745
                 
                 m=   0.07,# m value from https://doi.org/10.1371/journal.pcbi.1005745
                 
                 m1=   0.001,# m1 value from https://doi.org/10.1371/journal.pcbi.1005745
                 
                 sigma= 0.25, # sigma value from https://doi.org/10.1080/17513758.2010.488300 AND https://doi.org/10.1073/pnas.0402298101
                 
                 tau1= seq(0.1, 1), # tau1 value from https://doi.org/10.1371/journal.pcbi.1005745 
                 
                 c1= 0.05) # 

###
parametersc <- c(beta= 0.0007,# Estimated value
                 
                 gamma=  1/7,# gamma value from https://doi.org/10.1371/journal.pcbi.1004340
                 
                 mu=   0.001,# Estimated value  / #0.003,# mu value from https://doi.org/10.1371/journal.pcbi.1005745
                 
                 m=   0.07,# m value from https://doi.org/10.1371/journal.pcbi.1005745
                 
                 m1=   0.001,# m1 value from https://doi.org/10.1371/journal.pcbi.1005745
                 
                 sigma= 0.25, # sigma value from https://doi.org/10.1080/17513758.2010.488300 AND https://doi.org/10.1073/pnas.0402298101
                 
                 tau1= seq(0.1, 1), # tau1 value from https://doi.org/10.1371/journal.pcbi.1005745 
                 
                 c1= 0.05) # 

###
parametersd <- c(beta= 0.0007,# Estimated value
                 
                 gamma=  1/7,# gamma value from https://doi.org/10.1371/journal.pcbi.1004340
                 
                 mu=   0.0015,# Estimated value  / #0.003,# mu value from https://doi.org/10.1371/journal.pcbi.1005745
                 
                 m=   0.07,# m value from https://doi.org/10.1371/journal.pcbi.1005745
                 
                 m1=   0.001,# m1 value from https://doi.org/10.1371/journal.pcbi.1005745
                 
                 sigma= 0.25, # sigma value from https://doi.org/10.1080/17513758.2010.488300 AND https://doi.org/10.1073/pnas.0402298101
                 
                 tau1= seq(0.1, 1), # tau1 value from https://doi.org/10.1371/journal.pcbi.1005745 
                 
                 c1= 0.05) # 


###
parameterse <- c(beta= 0.0007,# Estimated value
                 
                 gamma=  1/7,# gamma value from https://doi.org/10.1371/journal.pcbi.1004340
                 
                 mu=   0.002,# Estimated value  / #0.003,# mu value from https://doi.org/10.1371/journal.pcbi.1005745
                 
                 m=   0.07,# m value from https://doi.org/10.1371/journal.pcbi.1005745
                 
                 m1=   0.001,# m1 value from https://doi.org/10.1371/journal.pcbi.1005745
                 
                 sigma= 0.25, # sigma value from https://doi.org/10.1080/17513758.2010.488300 AND https://doi.org/10.1073/pnas.0402298101
                 
                 tau1= seq(0.1, 1), # tau1 value from https://doi.org/10.1371/journal.pcbi.1005745 
                 
                 c1= 0.05) # 

###
parametersf <- c(beta= 0.0007,# Estimated value
                 
                 gamma=  1/7,# gamma value from https://doi.org/10.1371/journal.pcbi.1004340
                 
                 mu=   0.003,# Estimated value  / #0.003,# mu value from https://doi.org/10.1371/journal.pcbi.1005745
                 
                 m=   0.07,# m value from https://doi.org/10.1371/journal.pcbi.1005745
                 
                 m1=   0.001,# m1 value from https://doi.org/10.1371/journal.pcbi.1005745
                 
                 sigma= 0.25, # sigma value from https://doi.org/10.1080/17513758.2010.488300 AND https://doi.org/10.1073/pnas.0402298101
                 
                 tau1= seq(0.1, 1), # tau1 value from https://doi.org/10.1371/journal.pcbi.1005745 
                 
                 c1= 0.05) # 




#step8- state the times---- 
times <- seq(from=0, to=150, by=1)###
#step9- state the initial conditions----
xstart<- c(X= 720, S= 79, R1= 1)



###  Plot a  

#step10- compute a model trajectory----
library(tidyverse)
ode(
  func=single.modelmu,
  y=xstart,
  times=times,
  parms=parametersa
)%>%
  as.data.frame() -> data
#step11- plot the results----
plot_data <- data %>%
  gather(Legend, value,-time)%>%
  ggplot(aes(x=time,y=value,color=Legend, group=Legend))+
  geom_line(lwd=1.5) +
  # theme_classic()+#
  labs(x='Time(Days)', y='Number of individuals')
#show normal plot (with ggplot2)----
plot_data
#Step12- Interactive plot with plotly----
library(plotly)
#Show plot with plotly----
muplota <-ggplotly(plot_data)

### Plot b  

#step10- compute a model trajectory----
library(tidyverse)
ode(
  func=single.modelmu,
  y=xstart,
  times=times,
  parms=parametersb
)%>%
  as.data.frame() -> datb
#step11- plot the results----
plot_datb <- datb %>%
  gather(Legend, value,-time)%>%
  ggplot(aes(x=time,y=value,color=Legend, group=Legend))+
  geom_line(lwd=1.5) +
  # theme_classic()+#
  labs(x='Time(Days)', y='Number of individuals')
#show normal plot (with ggplot2)----
plot_datb
#Step12- Interactive plot with plotly----
library(plotly)
#Show plot with plotly----
muplotb <-ggplotly(plot_datb)

### Plot c  

#step10- compute a model trajectory----
library(tidyverse)
ode(
  func=single.modelmu,
  y=xstart,
  times=times,
  parms=parametersc
)%>%
  as.data.frame() -> datc
#step11- plot the results----
plot_datc <- datc %>%
  gather(Legend, value,-time)%>%
  ggplot(aes(x=time,y=value,color=Legend, group=Legend))+
  geom_line(lwd=1.5) +
  # theme_classic()+#
  labs(x='Time(Days)', y='Number of individuals')
#show normal plot (with ggplot2)----
plot_datc
#Step12- Interactive plot with plotly----
library(plotly)
#Show plot with plotly----
muplotc <-ggplotly(plot_datc)


### Plot d 

#step10- compute a model trajectory----
library(tidyverse)
ode(
  func=single.modelmu,
  y=xstart,
  times=times,
  parms=parametersd
)%>%
  as.data.frame() -> datd
#step11- plot the results----
plot_datd <- datd %>%
  gather(Legend, value,-time)%>%
  ggplot(aes(x=time,y=value,color=Legend, group=Legend))+
  geom_line(lwd=1.5) +
  # theme_classic()+#
  labs(x='Time(Days)', y='Number of individuals')
#show normal plot (with ggplot2)----
plot_datd
#Step12- Interactive plot with plotly----
library(plotly)
#Show plot with plotly----
muplotd <-ggplotly(plot_datd)


###   Plot e 

#step10- compute a model trajectory----
library(tidyverse)
ode(
  func=single.modelmu,
  y=xstart,
  times=times,
  parms=parameterse
)%>%
  as.data.frame() -> date
#step11- plot the results----
plot_date <- date %>%
  gather(Legend, value,-time)%>%
  ggplot(aes(x=time,y=value,color=Legend, group=Legend))+
  geom_line(lwd=1.5) +
  # theme_classic()+#
  labs(x='Time(Days)', y='Number of individuals')
#show normal plot (with ggplot2)----
plot_date
#Step12- Interactive plot with plotly----
library(plotly)
#Show plot with plotly----
muplote <-ggplotly(plot_date)



### Plot f

#step10- compute a model trajectory----
library(tidyverse)
ode(
  func=single.modelmu,
  y=xstart,
  times=times,
  parms=parametersf
)%>%
  as.data.frame() -> datf
#step11- plot the results----
plot_datf <- datf %>%
  gather(Legend, value,-time)%>%
  ggplot(aes(x=time,y=value,color=Legend, group=Legend))+
  geom_line(lwd=1.5) +
  # theme_classic()+#
  labs(x='Time(Days)', y='Number of individuals')
#show normal plot (with ggplot2)----
plot_datf
#Step12- Interactive plot with plotly----
library(plotly)
#Show plot with plotly----
muplotf <- ggplotly(plot_datf)




### Arrange Combined plots 
################################################################################
library(ggpubr)
ggarrange(plot_data,plot_datb,plot_datc,plot_datd,plot_date,plot_datf, ncol = 3, nrow = 2)





################################################################################
###########   Influence of β, µ, and tau1 on the Reproductive numbers ##########
################################################################################






