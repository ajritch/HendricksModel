#This script generates along-track plots of Dxs and 17xs
#Entirely theoretical

#code copied mostly from "big model.R"
#adapted from "model w evap balance and dxs.R" by Matt Winnick

#created Feb 4, 2016 by AJR

setwd('/Users/Annie/Documents/model/')
#plots will save here


########MODEL PARAMETERS/VARIABLES########
#change these to whatever you fancy

TempC=25 #surface temperature, degrees centigrade
TempK=TempC + 273.15 #deg K
delta18_atm_initial=-12 #initial atm. vapor d18O in VSMOW
delta18_precip_initial=-2   #initial precipitation d18O in VSMOW
rh=0.75; #relative humidity
e=seq(1,0, by = -0.2) ##evap as % of ET
t=1-e  ## transpiration as % of ET
track.length=41 #number of points along storm track

smow18=0.0020052 # smow ratio of 18O/16O
smowD=155.76*10^(-6) # smow ratio of D/H
smow17=379.9*10^(-6) # vsmow ratio of 17O/16O
theta.eq=0.529 #Barkan and Luz, 2005; equilibrium (l-v and v-l)
theta.diff=0.5185 #Barkan and Luz, 2007; diffusion
k18=1-0.01872 # kinetic fractionation factor for 18O, wet soils from Mathieu and Bariac (1996), slightly higher than lake evap (~14.3)
kD=1-0.01676 # kinetic fractionation factor for dD, wet soils from Mathieu and Bariac (1996), slightly higher than lake evap (~12.5)
k17=k18^theta.diff #Barkan and Luz, 2005, 2007


#########necessary equations###############

##run these equations and load the data above, and then each section 
#below can stand alone

##function to calculate either d18O (if d18 = True in input; default) or dD (if d18 = False)
#of ET using Craig and Gordon equations and closure assumption
craig_gordon <- function(t, h, delta_precip, ts, d18 = T){
  if(d18 == T){
    alpha18=1 - (-7.685 + (6.7123*10^3)/(ts) - (1.6664*10^6)/(ts)^2 + (0.35041*10^9)/(ts)^3)/1000 #equilibrium fractionation 18O, l-v
    r18_s <- (delta_precip/1000 + 1)*smow18  #ratio of precip sample
    Re_18 <- (1-t)*((alpha18*k18*r18_s)/(1-h+(1-t)*k18*h)) + t*r18_s*(1/(1+(1-t)*k18*(h/(1-h))))
    d18_evap <- (Re_18/smow18 - 1)*1000
    (d18_evap - delta_precip)/1000}
  else{
    alphaD=1 - ((1158.8*ts^3)/10^9 - (1620.1*ts^2)/10^6 + (794.84*ts)/10^3 - 161.04 + (2.9992*10^9)/ts^3)/1000 #equilibrium fractionation D, l-v
    rD_s <- (delta_precip/1000 + 1)*smowD 
    Re_D <- (1-t)*((alphaD*kD*rD_s)/(1-h+(1-t)*kD*h)) + t*rD_s*(1/(1+(1-t)*kD*(h/(1-h))))
    dD_evap <- (Re_D/smowD - 1)*1000
    (dD_evap - delta_precip)/1000}
}

##function to calculate d17O of ET using Craig-Gordon and closure assumption
craig_gordon17=function(t,h,delta_precip,ts) {
  alpha18=1 - (-7.685 + (6.7123*10^3)/(ts) - (1.6664*10^6)/(ts)^2 + (0.35041*10^9)/(ts)^3)/1000 #equilibrium fractionation 18O, l-v
  alpha17=alpha18^(theta.eq) #Barkan and Luz, 2005
  r17_s=(delta_precip/1000+1)*smow17
  Re_17=(1-t)*((alpha17*k17*r17_s)/(1-h+(1-t)*k17*h))+t*r17_s*(1/(1+(1-t)*k17*(h/(1-h))))
  d17_evap=(Re_17/smow17-1)*1000
  (d17_evap-delta_precip)/1000
}


##function to calculate integrated isotopic composition of ET with
#integration of Rayleigh distillation for soil moisture reservoir
isotope_e = function(resid, delta_precip, eps_e){
  delt = -resid*((1000/(1+eps_e))*resid^(eps_e)-1000)
  delta_precip + delt
}



########ATMOSPHERE COMPONENT########

## Set arrays of Nd, w, and x to iterate model over
nd <- c(seq(0.0025, 0.04, by = 0.0025),seq(0.05,1, by = 0.025), seq(1.1,10, by = .1), seq(20,100, by = 10), seq(200,500, by = 50), 1100) #quasi-logarithmic sequence of Nd's
w <- array(dim = c(track.length, length(nd)))
for(i in 1:length(nd)){
  w[,i] = seq(from = 1, by = -0.025*(1 - nd[i]/(nd[i]+1)), length.out = track.length)
}
x <- 1-w ##build array for non-dimensional distance

ts <- TempK #surface temperature
#alpha <- 1 + (-7.685 + (6.7123*10^3)/(ts) - (1.6664*10^6)/(ts)^2 + (0.35041*10^9)/(ts)^3)/1000 #equilibrium fractionation 18O, v-l
#aD <- 1 + ((1158.8*ts^3)/10^9 - (1620.1*ts^2)/10^6 + (794.84*ts)/10^3 - 161.04 + (2.9992*10^9)/ts^3)/1000 #equilibrium fractionation, D v-l
#constants for isotope value calculations:
#note that a18 != alpha18 (and same for D)
a18=1 + (-7.685 + (6.7123*10^3)/(ts) - (1.6664*10^6)/(ts)^2 + (0.35041*10^9)/(ts)^3)/1000 #equilibrium fractionation 18O, v-l
aD=1 + ((1158.8*ts^3)/10^9 - (1620.1*ts^2)/10^6 + (794.84*ts)/10^3 - 161.04 + (2.9992*10^9)/ts^3)/1000 #equilibrium fractionation, D v-l
alpha18=1 - (-7.685 + (6.7123*10^3)/(ts) - (1.6664*10^6)/(ts)^2 + (0.35041*10^9)/(ts)^3)/1000 #equilibrium fractionation 18O, l-v
alphaD=1 - ((1158.8*ts^3)/10^9 - (1620.1*ts^2)/10^6 + (794.84*ts)/10^3 - 161.04 + (2.9992*10^9)/ts^3)/1000 #equilibrium fractionation D, l-v
alpha17=alpha18^(theta.eq) #Barkan and Luz, 2005
a17=a18^theta.eq

## build empty arrays to store data for d18O and dD
delta_a <- delta_et <- delta_inf<- delta_p <- delta_D_a <- delta_et_D <- delta_D_inf <- delta_D_p <- delta_17_a <- delta_et_17 <- delta_17_inf <- delta_17_p <- array(dim=c(length(w[,1]),length(nd), length(e))) 
delta_a[1,,] <- delta18_atm_initial ## sets initial d18O atmospheric vapor
delta_p[1,,] <- delta18_precip_initial ## sets intial d18O precip value
exp <- a18+nd*a18-1 ## calculate exponent in Hendricks equation
delta_inf[1,,] <- (nd*(delta_p[1,,])-(1+nd)*(a18-1)*1000)/(a18+nd*a18-1) ## calculate d18O_infinity parameter from Hendricks equation
delta_D_a[1,,] <- delta_a[1,,]*8+10 ## sets atmospheric vapor dD
#delta_d_p[1,,] <- delta_d_a[1,,] + (aD-1)*1000 ## sets intial dD precip value
delta_D_p[1,,] <- delta_p[1,,]*8+10  #initial precip dD
exp_D <- aD+nd*aD-1 ## exponent of Hendricks for dD
delta_D_inf[1,,] <- (nd*(delta_D_p[1,,])-(1+nd)*(aD-1)*1000)/(aD+nd*aD-1) ## calculate dD_infinity parameter from Hendricks equation
delta_17_a[1,,] = (exp(0.528*log(delta_a[1,,]/1000+1)+0.000033)-1)*1000  #GMWL, Luz and Barkan, 2010
delta_17_p[1,,] = (exp(0.528*log(delta_p[1,,]/1000+1)+0.000033)-1)*1000  #GMWL, Luz and Barkan, 2010
exp_17 <- a17+nd*a17-1 ## exponent of Hendricks for d17O
delta_17_inf[1,,] <- (nd*(delta_17_p[1,,])-(1+nd)*(a17-1)*1000)/(a17+nd*a17-1) ## calculate dD_infinity parameter from Hendricks equation
rh <- rh ## this sets relative humidity

####
### RUN THE ATMOSPHERIC MODEL

## iterates over a storm track (w) for the arrays of e/t values and nd values
for(h in 1:length(e)){
  for(j in 1:length(nd)){
    for(i in 2:length(w[,1])){
      ## calculate residual soil moisture for isotope_e function
      residual <- 1-(nd[j]/(nd[j]+1))*e[h]/(1-(nd[j]/((nd[j]+1))*t[h]))
      ## if everything evaporates, set isotope value of evaporation to isotope value of precip
      if(residual == 0) {
        delta_et[i,j,h] <- delta_p[i-1,j,h]
        delta_et_D[i,j,h] <- delta_D_p[i-1,j,h]
        delta_et_17[i,j,h] <- delta_17_p[i-1,j,h]
      ## otherwise, first calculate the combined equilibrium/kinetic fractionation using craig-gordon function
      ## then calculate integrated ET value from location
      } else {
        eps_e <- craig_gordon(t[h], rh, delta_p[i-1,j,h],ts)
        eps_e_D <- craig_gordon(t[h], rh, delta_D_p[i-1,j,h],ts, d18 = F)
        eps_e_17 <- craig_gordon17(t[h], rh, delta_17_p[i-1,j,h],ts)
        delta_et[i,j,h] <- t[h]*delta_p[i-1,j,h] + e[h]*isotope_e(residual,delta_p[i-1,j,h], eps_e)
        delta_et_D[i,j,h] <- t[h]*delta_D_p[i-1,j,h] + e[h]*isotope_e(residual,delta_D_p[i-1,j,h], eps_e_D)
        delta_et_17[i,j,h] <- t[h]*delta_17_p[i-1,j,h] + e[h]*isotope_e(residual,delta_17_p[i-1,j,h], eps_e_17)
        }
      
      ## Hendricks model
      delta_inf[i,j,h] <- (nd[j]*delta_et[i,j,h]-(1+nd[j])*(a18-1)*1000)/(a18+nd[j]*a18-1)
      delta_a[i,j,h] <- (delta_a[1,j,h]-delta_inf[i,j,h])*(w[i,j]^exp[j])+delta_inf[i,j,h]
      delta_p[i,j,h] <-delta_a[i,j,h]+(a18-1)*1000
      
      delta_D_inf[i,j,h] <- (nd[j]*delta_et_D[i,j,h]-(1+nd[j])*(aD-1)*1000)/(aD+nd[j]*aD-1)
      delta_D_a[i,j,h] <- (delta_D_a[1,j,h]-delta_D_inf[i,j,h])*(w[i,j]^exp_D[j])+delta_D_inf[i,j,h]
      delta_D_p[i,j,h] <-delta_D_a[i,j,h]+(aD-1)*1000
      
      delta_17_inf[i,j,h] <- (nd[j]*delta_et_17[i,j,h]-(1+nd[j])*(a17-1)*1000)/(a17+nd[j]*a17-1)
      delta_17_a[i,j,h] <- (delta_17_a[1,j,h]-delta_17_inf[i,j,h])*(w[i,j]^exp_17[j])+delta_17_inf[i,j,h]
      delta_17_p[i,j,h] <-delta_17_a[i,j,h]+(a17-1)*1000
    }}}

#calculate deuterium, 17O excesses from dD, d17O and d18O values
Dxs <- delta_D_p - 8*delta_p
O17xs <- (log(delta_17_p/1000+1)-0.528*log(delta_p/1000+1))*10^6 #in permeg


####plots!####

library(fields) #for tim.colors


###
####plots of precip d18O , 17xs, and Dxs with varying E/ET along storm track
###

#pick Nd=1
#now use all locations along storm track:

colors = tim.colors(n = length(e)) #for varying E/ETs at each location

pdf('Precipitation d18O, Dxs, and 17Oxs with varying E-T along storm track.pdf',width=7,height=10)
par(mfrow=c(3,1))
#for d18O
plot(x[,which(nd==1)],delta_p[,which(nd==1),1],type='l',lwd=3,col=colors[1],xlab='Non-Dimensional Distance (x)',
     ylab=expression(paste(delta^18,'O of precipitation (permil)')),
     main=expression(paste(delta^18,'O, D-excess, and 17O-excess of precipitation along a storm track')))
for (i in 2:length(e)) {
  points(x[,which(nd==1)],delta_p[,which(nd==1),i],type='l',lwd=3,col=colors[i])
}
legend('bottomleft',c('E=100% ET','E=80%','E=60%','E=40%','E=20%','E=0%'),col=colors,lwd=3)

#for Dxs:
plot(x[,which(nd==1)],Dxs[,which(nd==1),1],type='l',lwd=3,col=colors[1],
     xlab='Non-Dimensional Distance (x)',ylab='Dxs of precipitation (permil)')
for (i in 2:length(e)) {
  points(x[,which(nd==1)],Dxs[,which(nd==1),i],type='l',lwd=3,col=colors[i])
}
legend('topleft',c('E=100% ET','E=80%','E=60%','E=40%','E=20%','E=0%'),col=colors,lwd=3)

#for 17xs:
plot(x[,which(nd==1)],O17xs[,which(nd==1),1],type='l',lwd=3,col=colors[1],
     xlab='Non-Dimensional Distance (x)',ylab='17Oxs of precipitation (permeg)',
     ylim=c(-45,32))
for (i in 2:length(e)) {
  points(x[,which(nd==1)],O17xs[,which(nd==1),i],type='l',lwd=3,col=colors[i])
}
legend('bottomleft',c('E=100% ET','E=80%','E=60%','E=40%','E=20%','E=0%'),col=colors,lwd=3)
dev.off()


