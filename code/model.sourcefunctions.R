#
#these are the source functions for the Hendricks model!
#load these before running the model!

#atmosphere and soil functions included

#file created 11-12-15 by AJR

#########necessary equations###############

# Calculates the fractionation factor to soil carbonate as a function
# of temperature. Default is 25C.
fract=function(Temp=25){ # From Cerling (1999)--see Figure 5
  fractionation=-(11.709-0.116*Temp+2.16e-4*Temp^2)
  fractionation
}

# Converts CO2 from ppmv to mols/cm3 for use in the soil respiration
# equation.  This conversion is a function of temperature (default is 25C).
CO2convert=function(CO2ppmv,Temp=25){ 
  CO2v=CO2ppmv/1e6
  MCO2=44.00964 # Molecular mass of bulk CO2
  SP=101325 # Standard pressure (Pa) at 1atm
  R=8.31441 # Universal Gas Constant
  TempK=Temp+273.15 # Correction to Kelvin 
  CO2m=(MCO2*SP*CO2v)/(R*TempK) # g/m3 of CO2
  CO2moles=CO2m/MCO2/1e6 # moles/cm3
  CO2moles #moles/cm3
}

# This function does the opposite of the CO2convert function. It takes
# CO2 as moles/cm3 and converts back to ppm.
CO2backconvert=function(CO2_atm,Temp=25){
  MCO2=44.00964 # Molecular mass of bulk CO2
  SP=101325 # Standard pressure (Pa) at 1atm
  R=8.31441 # Universal Gas Constant
  TempK=Temp+273.15 # Correction to Kelvin
  CO2mm3=CO2_atm*MCO2*1e6 # g/m3 of CO2
  CO2v=(CO2mm3*R*TempK)/(SP*MCO2) # absolute quantity of CO2
  CO2ppm=CO2v*1e6 # ppm CO2
  CO2ppm #ppm
}


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

#function to convert negative values to NAs (to get rid of neg. Nds)
#only used when averaging tracks
neg.toNA=function(matrix) {
  for (i in 1:dim(matrix)[1]) {
    for (j in 1:dim(matrix)[2]) {
      if (matrix[i,j]<=0) {
        matrix[i,j]=NA
      }
    }
  }
  matrix
}


