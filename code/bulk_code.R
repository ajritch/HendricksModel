#this script creates C and O isotope profiles in soils along a designated storm track
#both advection-only and eddy diffusion-only cases explored
#fixing the diffusion coefficient units, because that was a major fail before
#adapted from multiple different script components:
   #"model w evap balance and dxs.R" by Matt Winnick (atmosphere)
   #"SoilModule_d18OCodeV2.R" by D.Ibarra and Chamberlain (soil moisture)
   #"SR_preciprelationship.R" and "Soil Resp-forward.R" by J. Caves (soil respiration)
#created 4-28-15 by AJR
#last modified: 5-28-16 by AJR


setwd('/Users/Annie/Documents/model')
load("CaliforniaNevadaStormTrack_NARR_July2.RData")
##^this RData file contains 28x10 (28 points of 10 tracks) dataframes of the following variables:
#MAP.data          #mean annual precipitation rate (mm/month)
#Nd1.data          #Nd (evap/runoff) unitless DJF average, calculated with evap and NARR runoff
#Nd2.data          #Nd (ET/P-ET) unitless DJF average, calculated using lhtfl as ET
#Prw.data          #Precipitable water (w) kg/m^2 DJF average
#Temp_2m.data      #2m air temperature DFJ average (K)
#Temp_surf.data    #surface temperature JJASO average (K)
#evap.JJASO.data   #evapotranspiration (from lhtfl), JJASO average (kg/m^2)
#land.mask.data    #land mask (0 if not land, 1 if land)
#lat.data          #latitude
#lon.data          #longitude
#rhum.data         #2m relative humidity (%) DJF average
####additionally, 28x10x801 array at 801 depths (0-800cm)
#tsoil.data        #soil temperature (K) from depth 0-800cm (linearly interpolated from original NARR data)

Nd.data=Nd2.data 
evap.data=evap.JJASO.data #for coding simplification purposes

#####parameters and inputs#####

e=seq(1,0, by = -0.1) #evap as % of ET
t=1-e #transpiration as % of ET

Rpdb=0.0112372 #dimensionless PDB Belemnite Ratio
smow18=0.0020052 # smow ratio of 18O/16O
smowD=155.76*10^(-6) # smow ratio of D/H
smow17=379.9*10^(-6) # vsmow ratio of 17O/16O
theta.eq=0.529 #Barkan and Luz, 2005; equilibrium (l-v and v-l)
theta.diff=0.5185 #Barkan and Luz, 2007; diffusion
k18=1-0.01872 # kinetic fractionation factor for 18O, wet soils from Mathieu and Bariac (1996), slightly higher than lake evap (~14.3)
kD=1-0.01676 # kinetic fractionation factor for dD, wet soils from Mathieu and Bariac (1996), slightly higher than lake evap (~12.5)
k17=k18^theta.diff #Barkan and Luz, 2005, 2007

#INPUT d18O VALUES! (in permil)
initial_delta_p=-6.96 #average value from Ingraham and Taylor 1991, Traverse II
initial_delta_a=-17.476 #assuming equlibrium with precip
initial_delta_p_eddy=initial_delta_p #this is what Matt used, at least
initial_delta_a_eddy=initial_delta_a


#####soil model inputs!!!

pH=8 
delta18_atm_initial=initial_delta_p
delta18_precip_initial=initial_delta_a
maxdepth=100 #in cm, maximum soil profile depth 
zplot=seq(0,maxdepth) #for plotting profiles

#keep only subset of soil temp data above (and at) maxdepth of soil:
tsoil.data=tsoil.data[,,1:(maxdepth+1)]

#soil properties
d=2.3e-5 #cm2/sec self diffusion coefficient of soil water
Dair=0.144 #cm2/s #CO2 diffusion coefficient in air (Cerling 1991)
p=0.5 #Free-air porosity (from Cerling 1991); Cerling and Quade 1993 use 0.6
tau=0.6 #Toturosity (0.7 is uniform sand from) Barnes and Allison (1983)
Ds=Dair*p*tau #cm2/s #Bulk CO2 diffusion coefficient in the air of soils

#leaving these here as equation reminder; they're handled elsewhere
#zstar=(p*tau*d)/E
#E=(p*tau*d)/zstar 

M12=44 #Mass of 12CO2
M13=45 #Mass of 13CO2
Mair=29 #Average atomic mass of air
MCO2=44.00964 #Mass of bulk CO2
diff=sqrt(((MCO2+Mair)/(MCO2*Mair))*((M13*Mair)/(M13+Mair))) #Fractional increase in 13CO2 diffusion coefficient over bulk CO2
#diff=1.0044

Ds13=Ds/diff #cm2/s #should be Ds/1.0044 relative to 12C

delta_13_o=-27 #per mil composition of soil-respired CO2
pCO2ppm_atm=280
CO2_atm_mg.m3=543 #mg/m3 #Atmospheric CO2 #Currently set for 280ppm
CO2_atm=CO2_atm_mg.m3/1000/100^3/MCO2 #moles/cm3 CO2
#CO2_atm=1.23381e-8 #moles/cm3 This is at 280ppm
delta_13_atm=-6.5 #per mil of atmospheric CO2

z=50 #cm depth IN soil CRITICAL VALUE!! for determining depth-MAP relationship
L=100 #cm depth of soil column with linear production (for resp-MAP relationship)
#Soil parameters assume linear production with depth with a soil depth (L) of 100 cm.


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


#########basic atmosphere model for multiple storm tracks###########

#model input NARR data
w=as.matrix(Prw.data) #precipitatble water, kg/m^2
ntracks=length(w[1,]) #number of storm tracks
tracklength=length(w[,1]) #length of storm tracks
Nd=as.matrix(Nd.data) #input array from NARR data
ts=as.matrix(Temp_2m.data) #2m surface temperature (K)
rh=as.matrix(rhum.data/100) #relative humidity (fraction)


#advection only (Hendricks):
#delta_a=(delta_a_initial-delta_a_inf)*exp((alpha+alpha*Nd-1)*(-x/l))+delta_a_inf
#delta_a_inf=(Nd*delta_et-(1+Nd)*(alpha-1)*1000)/(alpha+alpha*Nd-1)
#(-x/l)=-x' ; l=-w*(dx/dw) ==> -x'=(x/w)*(dw/dx) 

#eddy diffusion only (Hendricks):
#delta_a_eddy=(delta_a_initial-delta_a_inf)*exp((sqrt(alpha+alpha*Nd)-1)*(-x/l))+delta_a_inf

#constants for isotope value calculations:
#note that a18 != alpha18 (and same for D)
a18=1 + (-7.685 + (6.7123*10^3)/(ts) - (1.6664*10^6)/(ts)^2 + (0.35041*10^9)/(ts)^3)/1000 #equilibrium fractionation 18O, v-l
aD=1 + ((1158.8*ts^3)/10^9 - (1620.1*ts^2)/10^6 + (794.84*ts)/10^3 - 161.04 + (2.9992*10^9)/ts^3)/1000 #equilibrium fractionation, D v-l
alpha18=1 - (-7.685 + (6.7123*10^3)/(ts) - (1.6664*10^6)/(ts)^2 + (0.35041*10^9)/(ts)^3)/1000 #equilibrium fractionation 18O, l-v
alphaD=1 - ((1158.8*ts^3)/10^9 - (1620.1*ts^2)/10^6 + (794.84*ts)/10^3 - 161.04 + (2.9992*10^9)/ts^3)/1000 #equilibrium fractionation D, l-v
alpha17=alpha18^(theta.eq) #Barkan and Luz, 2005
a17=a18^theta.eq

#build empty arrays to store data for d18O and dD and d17O for advection-only/eddy-only cases
delta_a <- delta_et <- delta_inf<- delta_p <- delta_a_eddy <- delta_et_eddy <- delta_inf_eddy <- delta_p_eddy <- array(dim=c(tracklength,ntracks,length(e)))
delta_D_a <- delta_D_et <- delta_D_inf <- delta_D_p <- delta_D_a_eddy <- delta_D_et_eddy <- delta_D_inf_eddy <- delta_D_p_eddy <- array(dim=c(tracklength,ntracks,length(e)))
delta_17_a <- delta_17_et <- delta_17_inf <- delta_17_p <- delta_17_a_eddy <- delta_17_et_eddy <- delta_17_inf_eddy <- delta_17_p_eddy <- array(dim=c(tracklength,ntracks,length(e)))

delta_a[1,,] = initial_delta_a #sets initial d18O atmospheric vapor for all tracks and E/ETs
delta_p[1,,] = initial_delta_p #sets intial d18O precip  #ntracksxlength(e)
delta_D_a[1,,] = delta_a[1,,]*8+10 #sets initial atmospheric vapor dD
delta_D_p[1,,] = delta_p[1,,]*8+10 #sets initial precip dD
delta_17_a[1,,] = (exp(0.528*log(delta_a[1,,]/1000+1)+0.000033)-1)*1000  #GMWL, Luz and Barkan, 2010
delta_17_p[1,,] = (exp(0.528*log(delta_p[1,,]/1000+1)+0.000033)-1)*1000  #GMWL, Luz and Barkan, 2010
#same for eddies:
delta_a_eddy[1,,] = initial_delta_a_eddy #sets initial d18O atmospheric vapor for all tracks and E/ETs
delta_p_eddy[1,,] = initial_delta_p_eddy #sets intial d18O precip  #ntracksxlength(e)
delta_D_a_eddy[1,,] = delta_a_eddy[1,,]*8+10 #sets initial atmospheric vapor dD
delta_D_p_eddy[1,,] = delta_p_eddy[1,,]*8+10 #sets initial precip dD
delta_17_a_eddy[1,,] = (exp(0.528*log(delta_a_eddy[1,,]/1000+1)+0.000033)-1)*1000  #GMWL, Luz and Barkan, 2010
delta_17_p_eddy[1,,] = (exp(0.528*log(delta_p_eddy[1,,]/1000+1)+0.000033)-1)*1000  #GMWL, Luz and Barkan, 2010


#make (mostly empty) matrix of total distance from initial point
distance.m=array(dim=c(tracklength,ntracks))
distance.m[1,]=0 #first distance zero because we haven't gone anywhere yet

###run the model!

library(geosphere) #to calculate distance between coordinates (Haversine)

#iterate over all points along each track for each E/ET ratio
for (h in 1:length(e)) {
  for(j in 1:ntracks) {
    for (i in 2:tracklength) {
      #calculate residual soil moisture for isotope_e function
      if (land.mask.data[i,j]==0) {
        residual=0   #sets 0 residual soil moisture when not over land
      } else {
        residual=1-(Nd[i,j]/(Nd[i,j]+1))*e[h]/(1-(Nd[i,j]/(Nd[i,j]+1)*t[h]))}
      #if everything evaporates, set isotope value of evaporation to isotope value of precip
      if(residual == 0){
        delta_et[i,j,h] <- delta_p[i-1,j,h]
        delta_D_et[i,j,h] <- delta_D_p[i-1,j,h]
        delta_17_et[i,j,h] = delta_17_p[i-1,j,h]
        delta_et_eddy[i,j,h]=delta_p_eddy[i-1,j,h]
        delta_D_et_eddy[i,j,h] <- delta_D_p_eddy[i-1,j,h]
        delta_17_et_eddy[i,j,h] = delta_17_p_eddy[i-1,j,h]
        
        #otherwise, first calculate the combined equilibrium/kinetic fractionation using Craig-Gordon function
      } else {
        eps_e <- craig_gordon(t[h], rh[i,j], delta_p[i-1,j,h], ts[i,j])
        eps_e_D <- craig_gordon(t[h], rh[i,j], delta_D_p[i-1,j,h], ts[i,j], d18 = F)
        eps_e_17 = craig_gordon17(t[h],rh[i,j],delta_17_p[i-1,j,h],ts[i,j])
        eps_e_eddy <- craig_gordon(t[h], rh[i,j], delta_p_eddy[i-1,j,h], ts[i,j])
        eps_e_D_eddy <- craig_gordon(t[h], rh[i,j], delta_D_p_eddy[i-1,j,h], ts[i,j], d18 = F)
        eps_e_17_eddy = craig_gordon17(t[h],rh[i,j],delta_17_p_eddy[i-1,j,h],ts[i,j])
        #then calculate integrated ET value from location
        delta_et[i,j,h] <- t[h]*delta_p[i-1,j,h] + e[h]*isotope_e(residual,delta_p[i-1,j,h], eps_e)
        delta_D_et[i,j,h] <- t[h]*delta_D_p[i-1,j,h] + e[h]*isotope_e(residual,delta_D_p[i-1,j,h], eps_e_D)
        delta_17_et[i,j,h] = t[h]*delta_17_p[i-1,j,h]+e[h]*isotope_e(residual,delta_17_p[i-1,j,h],eps_e_17) 
        delta_et_eddy[i,j,h] <- t[h]*delta_p_eddy[i-1,j,h] + e[h]*isotope_e(residual,delta_p_eddy[i-1,j,h], eps_e_eddy)
        delta_D_et_eddy[i,j,h] <- t[h]*delta_D_p_eddy[i-1,j,h] + e[h]*isotope_e(residual,delta_D_p_eddy[i-1,j,h], eps_e_D_eddy)
        delta_17_et_eddy[i,j,h] = t[h]*delta_17_p_eddy[i-1,j,h]+e[h]*isotope_e(residual,delta_17_p_eddy[i-1,j,h],eps_e_17_eddy)
      }
      
      
      ###there's an issue with some of the values after the initial point being NAs
      #this loop will make it so all of the NAs stay NAs *except* the last NA before
      #hitting land. the first ones stay NAs so that we're only calculating over land.
      #note: still necessary to keep the original initials above in the event of no NAs
      if (land.mask.data[i,j]==0) {
        delta_p[i,j,]=initial_delta_p
        delta_p[i-1,j,]=NA
        delta_a[i,j,]=initial_delta_a
        delta_a[i-1,j,]=NA
        delta_D_p[i,j,]=delta_p[i,j,]*8+10 #assume initial points are on MWL
        delta_D_p[i-1,j,]=NA
        delta_D_a[i,j,]=delta_a[i,j,]*8+10
        delta_D_a[i-1,j,]=NA
        delta_17_p[i,j,]=(exp(0.528*log(delta_p[i,j,]/1000+1)+0.000033)-1)*1000  #GMWL, Luz and Barkan, 2010
        delta_17_p[i-1,j,]=NA
        delta_17_a[i,j,]=(exp(0.528*log(delta_a[i,j,]/1000+1)+0.000033)-1)*1000  #GMWL, Luz and Barkan, 2010
        delta_17_a[i-1,j,]=NA
        #repeat for eddies (copypasta)
        delta_p_eddy[i,j,]=initial_delta_p_eddy
        delta_p_eddy[i-1,j,]=NA
        delta_a_eddy[i,j,]=initial_delta_a_eddy
        delta_a_eddy[i-1,j,]=NA
        delta_D_p_eddy[i,j,]=delta_p_eddy[i,j,]*8+10 #assume initial points are on MWL
        delta_D_p_eddy[i-1,j,]=NA
        delta_D_a_eddy[i,j,]=delta_a_eddy[i,j,]*8+10
        delta_D_a_eddy[i-1,j,]=NA
        delta_17_p_eddy[i,j,]=(exp(0.528*log(delta_p_eddy[i,j,]/1000+1)+0.000033)-1)*1000  #GMWL, Luz and Barkan, 2010
        delta_17_p_eddy[i-1,j,]=NA
        delta_17_a_eddy[i,j,]=(exp(0.528*log(delta_a_eddy[i,j,]/1000+1)+0.000033)-1)*1000  #GMWL, Luz and Barkan, 2010
        delta_17_a_eddy[i-1,j,]=NA
      } else {
        #haversine formula to calculate distance between coordinates:
        x=distVincentyEllipsoid(c(lon.data[i,j],lat.data[i,j]),c(lon.data[i-1,j],lat.data[i-1,j]))
        dx=x #meters, but the units don't matter. slightly >32km
        dw=w[i,j]-w[i-1,j]
        dw=ifelse(dw>0,0,dw) #gets rid of positive gradients (b/c we don't like explosions here)
        #Hendricks model (advection only)
        delta_inf[i,j,h] <- (Nd[i,j]*delta_et[i,j,h]-(1+Nd[i,j])*(a18[i,j]-1)*1000)/(a18[i,j]+Nd[i,j]*a18[i,j]-1)
        delta_a[i,j,h]=(delta_a[i-1,j,h]-delta_inf[i,j,h])*exp((a18[i,j]+a18[i,j]*Nd[i,j]-1)*(x/w[i,j])*(dw/dx))+delta_inf[i,j,h]
        delta_p[i,j,h] <-delta_a[i,j,h]+(a18[i,j]-1)*1000
        
        delta_D_inf[i,j,h] <- (Nd[i,j]*delta_D_et[i,j,h]-(1+Nd[i,j])*(aD[i,j]-1)*1000)/(aD[i,j]+Nd[i,j]*aD[i,j]-1)
        delta_D_a[i,j,h]=(delta_D_a[i-1,j,h]-delta_D_inf[i,j,h])*exp((aD[i,j]+aD[i,j]*Nd[i,j]-1)*(x/w[i,j])*(dw/dx))+delta_D_inf[i,j,h]
        delta_D_p[i,j,h] <-delta_D_a[i,j,h]+(aD[i,j]-1)*1000
        
        delta_17_inf[i,j,h] <- (Nd[i,j]*delta_17_et[i,j,h]-(1+Nd[i,j])*(a17[i,j]-1)*1000)/(a17[i,j]+Nd[i,j]*a17[i,j]-1)
        delta_17_a[i,j,h]=(delta_17_a[i-1,j,h]-delta_17_inf[i,j,h])*exp((a17[i,j]+a17[i,j]*Nd[i,j]-1)*(x/w[i,j])*(dw/dx))+delta_17_inf[i,j,h]
        delta_17_p[i,j,h] <-delta_17_a[i,j,h]+(a17[i,j]-1)*1000
        
        #Hendricks model (eddy diffusion only)
        delta_inf_eddy[i,j,h] <- (Nd[i,j]*delta_et_eddy[i,j,h]-(1+Nd[i,j])*(a18[i,j]-1)*1000)/(a18[i,j]+Nd[i,j]*a18[i,j]-1)
        delta_a_eddy[i,j,h]=(delta_a_eddy[i-1,j,h]-delta_inf_eddy[i,j,h])*exp((sqrt(a18[i,j]+a18[i,j]*Nd[i,j])-1)*(x/w[i,j])*(dw/dx))+delta_inf_eddy[i,j,h]
        delta_p_eddy[i,j,h] <-delta_a_eddy[i,j,h]+(a18[i,j]-1)*1000
        
        delta_D_inf_eddy[i,j,h] <- (Nd[i,j]*delta_D_et_eddy[i,j,h]-(1+Nd[i,j])*(aD[i,j]-1)*1000)/(aD[i,j]+Nd[i,j]*aD[i,j]-1)
        delta_D_a_eddy[i,j,h]=(delta_D_a_eddy[i-1,j,h]-delta_D_inf_eddy[i,j,h])*exp((sqrt(aD[i,j]+aD[i,j]*Nd[i,j])-1)*(x/w[i,j])*(dw/dx))+delta_D_inf_eddy[i,j,h]
        delta_D_p_eddy[i,j,h] <-delta_D_a_eddy[i,j,h]+(aD[i,j]-1)*1000
        
        delta_17_inf_eddy[i,j,h] <- (Nd[i,j]*delta_17_et_eddy[i,j,h]-(1+Nd[i,j])*(a17[i,j]-1)*1000)/(a17[i,j]+Nd[i,j]*a17[i,j]-1)
        delta_17_a_eddy[i,j,h]=(delta_17_a_eddy[i-1,j,h]-delta_17_inf_eddy[i,j,h])*exp((sqrt(a17[i,j]+a17[i,j]*Nd[i,j])-1)*(x/w[i,j])*(dw/dx))+delta_17_inf_eddy[i,j,h]
        delta_17_p_eddy[i,j,h] <-delta_17_a_eddy[i,j,h]+(a17[i,j]-1)*1000
      }
      
      #this is essentially the distance from the last point before landfall
      #(or first point on land, but not for our lovely NARR swath)
      if (land.mask.data[i,j]==0) {
        distance.m[i-1,j]=0
        distance.m[i,j]=0
      } else {
        distance.m[i,j]=distance.m[i-1,j]+dx
      }  
    }  
  }    
}
#rename all the resulting dataframes
delta_p.all=delta_p
delta_D_p.all=delta_D_p
delta_17_p.all=delta_17_p
delta_p_eddy.all=delta_p_eddy
delta_D_p_eddy.all=delta_D_p_eddy
delta_17_p_eddy.all=delta_17_p_eddy
distance.all=distance.m/1000
distance.fake = seq(0,by=32,length=tracklength) #approximate distances

###calculate D-excess and 17O-excess:
Dxs.all=delta_D_p.all-8*delta_p.all #in permil
O17xs.all=(log(delta_17_p.all/1000+1)-0.528*log(delta_p.all/1000+1))*10^6 #in permeg
Dxs_eddy.all=delta_D_p_eddy.all-8*delta_p_eddy.all #in permil
O17xs_eddy.all=(log(delta_17_p_eddy.all/1000+1)-0.528*log(delta_p_eddy.all/1000+1))*10^6 #in permeg

#delta_p[,,1] #quick look for E/ET=0
#^ NAs at beginning of tracks from where we're still over the ocean
#^ NaNs where Damkohler goes negative

#distance plot of delta_p

library(fields)
colors=tim.colors(ntracks)
colors_eddy=paste('gray',seq(5,86,length=10))
plot(distance.all[,1],delta_p.all[,1,1],type="o",pch=16,cex=0.5,xaxs="i",
     xlab='Distance from "Coast" (km)',ylab="delta_p",ylim=c(-30,0),col=colors[1])
for (j in 2:10){
  points(distance.all[,j],delta_p.all[,j,1],type="o",pch=16,cex=0.5,col=colors[j])
}
for(j in 1:10) {
  points(distance.all[,j],delta_p_eddy.all[,j,1],type='o',pch=15,cex=0.5,col=colors_eddy[j])
}


#plot of Nd
#for plotting distance from coast:
#convert Nd to zero when not over land
Nd.plot=Nd
for (j in 1:ntracks) {
  for (i in 2:tracklength) {
    if (land.mask.data[i,j]==0) {
      Nd.plot[i-1,j]=NA
    }}}
plot(distance.all[,1],Nd.plot[,1],type='o',pch=16,cex=0.5,col=colors[1],yaxs='i',xaxs='i',
     xlab='Distance from "Coast" (km)',ylab='Nd',ylim=c(0,20))
for (j in 2:10){
  points(distance.all[,j],Nd.plot[,j],type="o",pch=16,cex=0.5,col=colors[j])
}

##plots of Dxs and 17Oxs from "coast"
#Dxs
plot(distance.all[,1],Dxs.all[,1,1],type="o",pch=16,cex=0.5,xaxs="i",
     xlab='Distance from "Coast" (km)',ylim=c(0,132),ylab="d-excess",col=colors[1])
for (j in 2:10){
  points(distance.all[,j],Dxs.all[,j,1],type="o",pch=16,cex=0.5,col=colors[j])
}
for (j in 1:10) {
  points(distance.all[,j],Dxs_eddy.all[,j,1],type='o',pch=15,cex=0.5,col=colors_eddy[j])
}
#17xs
plot(distance.all[,1],O17xs.all[,1,1],type="o",pch=16,cex=0.5,xaxs="i",
     xlab='Distance from "Coast" (km)',ylim=c(-100,10),ylab="17O-excess (permeg)",col=colors[1])
for (j in 2:10){
  points(distance.all[,j],O17xs.all[,j,1],type="o",pch=16,cex=0.5,col=colors[j])
}
for (j in 1:10) {
  points(distance.all[,j],O17xs_eddy.all[,j,1],type='o',pch=15,cex=0.5,col=colors_eddy[j])
}


####

#######average tracks together by distance from coast#########
####

###values from each track are averaged together by the distance from 
#the coast rather than by the distance from the first data point
###note that these are approximate averages as the distance
#between points of different storm tracks are not equivalent
#ie, I'm assuming a distance of 32km between each point for each track.

#turn all of the inputs into one storm track:
#pre for pre-averaging
w.pre=as.matrix(Prw.data) #precipitatble water, kg/m^2
ntracks.pre=dim(w.pre)[2]
ntracks=1 #number of storm tracks
tracklength=dim(w.pre)[1] #length of storm track
Nd.NA=neg.toNA(as.matrix(Nd.data)) #convert negatives to NA to remove from averaging
ts.pre=as.matrix(Temp_2m.data) #2m surface temperature (K)
rh.pre=as.matrix(rhum.data/100) #relative humidity (fraction)

#if land.mask=0, cut off the previous data point and add NA to the end
#cutting off previous allows for initiation just before land, consistent with sections above
#adding NA at the end will allow for averaging of all the non-NAs for farther out on track
for (j in 1:ntracks.pre) {
  for (i in 2:tracklength) {
    if (land.mask.data[i,j]==0) {
      w.pre[,j]=c(w.pre[-(i-1),j],NA)
      Nd.NA[,j]=c(Nd.NA[-(i-1),j],NA)
      ts.pre[,j]=c(ts.pre[-(i-1),j],NA)
      rh.pre[,j]=c(rh.pre[-(i-1),j],NA)
    }
  }
}

#now average these. this average is based on distance from coast!
w=rowMeans(w.pre,na.rm=T)
Nd=rowMeans(Nd.NA,na.rm=T)
ts=rowMeans(ts.pre,na.rm=T)
rh=rowMeans(rh.pre,na.rm=T)

#now copypasta from above section.

#constants for isotope value calculations:
#note that a18 != alpha18 (and same for D)
a18=1 + (-7.685 + (6.7123*10^3)/(ts) - (1.6664*10^6)/(ts)^2 + (0.35041*10^9)/(ts)^3)/1000 #equilibrium fractionation 18O, v-l
aD=1 + ((1158.8*ts^3)/10^9 - (1620.1*ts^2)/10^6 + (794.84*ts)/10^3 - 161.04 + (2.9992*10^9)/ts^3)/1000 #equilibrium fractionation, D v-l
alpha18=1 - (-7.685 + (6.7123*10^3)/(ts) - (1.6664*10^6)/(ts)^2 + (0.35041*10^9)/(ts)^3)/1000 #equilibrium fractionation 18O, l-v
alphaD=1 - ((1158.8*ts^3)/10^9 - (1620.1*ts^2)/10^6 + (794.84*ts)/10^3 - 161.04 + (2.9992*10^9)/ts^3)/1000 #equilibrium fractionation D, l-v
alpha17=alpha18^(theta.eq) #Barkan and Luz, 2005
a17=a18^theta.eq

#build empty arrays to store data for d18O and dD and d17O for advection-only/eddy-only cases
delta_a <- delta_et <- delta_inf<- delta_p <- delta_a_eddy <- delta_et_eddy <- delta_inf_eddy <- delta_p_eddy <- array(dim=c(tracklength,length(e)))
delta_D_a <- delta_D_et <- delta_D_inf <- delta_D_p <- delta_D_a_eddy <- delta_D_et_eddy <- delta_D_inf_eddy <- delta_D_p_eddy <- array(dim=c(tracklength,length(e)))
delta_17_a <- delta_17_et <- delta_17_inf <- delta_17_p <- delta_17_a_eddy <- delta_17_et_eddy <- delta_17_inf_eddy <- delta_17_p_eddy <- array(dim=c(tracklength,length(e)))

delta_a[1,] = initial_delta_a #sets initial d18O atmospheric vapor for all tracks and E/ETs
delta_p[1,] = initial_delta_p #sets intial d18O precip  #ntracksxlength(e)
delta_D_a[1,] = delta_a[1,]*8+10 #sets initial atmospheric vapor dD
delta_D_p[1,] = delta_p[1,]*8+10 #sets initial precip dD
delta_17_a[1,] = (exp(0.528*log(delta_a[1,]/1000+1)+0.000033)-1)*1000  #GMWL, Luz and Barkan, 2010
delta_17_p[1,] = (exp(0.528*log(delta_p[1,]/1000+1)+0.000033)-1)*1000  #GMWL, Luz and Barkan, 2010
#same for eddies:
delta_a_eddy[1,] = initial_delta_a_eddy #sets initial d18O atmospheric vapor for all tracks and E/ETs
delta_p_eddy[1,] = initial_delta_p_eddy #sets intial d18O precip  #ntracksxlength(e)
delta_D_a_eddy[1,] = delta_a_eddy[1,]*8+10 #sets initial atmospheric vapor dD
delta_D_p_eddy[1,] = delta_p_eddy[1,]*8+10 #sets initial precip dD
delta_17_a_eddy[1,] = (exp(0.528*log(delta_a_eddy[1,]/1000+1)+0.000033)-1)*1000  #GMWL, Luz and Barkan, 2010
delta_17_p_eddy[1,] = (exp(0.528*log(delta_p_eddy[1,]/1000+1)+0.000033)-1)*1000  #GMWL, Luz and Barkan, 2010


###run the model!

#iterate over all points along each track for each E/ET ratio
for (h in 1:length(e)) {
  for (i in 2:tracklength) {
    #calculate residual soil moisture for isotope_e function
    residual=1-(Nd[i]/(Nd[i]+1))*e[h]/(1-(Nd[i]/(Nd[i]+1)*t[h]))
    #if everything evaporates, set isotope value of evaporation to isotope value of precip
    if(residual == 0){
      delta_et[i,h] <- delta_p[i-1,h]
      delta_D_et[i,h] <- delta_D_p[i-1,h]
      delta_17_et[i,h] = delta_17_p[i-1,h]
      delta_et_eddy[i,h]=delta_p_eddy[i-1,h]
      delta_D_et_eddy[i,h] <- delta_D_p_eddy[i-1,h]
      delta_17_et_eddy[i,h] = delta_17_p_eddy[i-1,h]
      
      #otherwise, first calculate the combined equilibrium/kinetic fractionation using Craig-Gordon function
    } else {
      eps_e <- craig_gordon(t[h], rh[i], delta_p[i-1,h], ts[i])
      eps_e_D <- craig_gordon(t[h], rh[i], delta_D_p[i-1,h], ts[i], d18 = F)
      eps_e_17 = craig_gordon17(t[h],rh[i],delta_17_p[i-1,h],ts[i])
      eps_e_eddy <- craig_gordon(t[h], rh[i], delta_p_eddy[i-1,h], ts[i])
      eps_e_D_eddy <- craig_gordon(t[h], rh[i], delta_D_p_eddy[i-1,h], ts[i], d18 = F)
      eps_e_17_eddy = craig_gordon17(t[h],rh[i],delta_17_p_eddy[i-1,h],ts[i])
      #then calculate integrated ET value from location
      delta_et[i,h] <- t[h]*delta_p[i-1,h] + e[h]*isotope_e(residual,delta_p[i-1,h], eps_e)
      delta_D_et[i,h] <- t[h]*delta_D_p[i-1,h] + e[h]*isotope_e(residual,delta_D_p[i-1,h], eps_e_D)
      delta_17_et[i,h] = t[h]*delta_17_p[i-1,h]+e[h]*isotope_e(residual,delta_17_p[i-1,h],eps_e_17) 
      delta_et_eddy[i,h] <- t[h]*delta_p_eddy[i-1,h] + e[h]*isotope_e(residual,delta_p_eddy[i-1,h], eps_e_eddy)
      delta_D_et_eddy[i,h] <- t[h]*delta_D_p_eddy[i-1,h] + e[h]*isotope_e(residual,delta_D_p_eddy[i-1,h], eps_e_D_eddy)
      delta_17_et_eddy[i,h] = t[h]*delta_17_p_eddy[i-1,h]+e[h]*isotope_e(residual,delta_17_p_eddy[i-1,h],eps_e_17_eddy)
    }
    
    
    ###there's an issue with some of the values after the initial point being NAs
    #this loop will make it so all of the NAs stay NAs *except* the last NA before
    #hitting land. the first ones stay NAs so that we're only calculating over land.
    #note: still necessary to keep the original initials above in the event of no NAs
    #only where rowMeans==0 will they all be above water==>this track will be above water
    if (rowMeans(land.mask.data)[i]==0) {
      delta_p[i,]=initial_delta_p
      delta_p[i-1,]=NA
      delta_a[i,]=initial_delta_a
      delta_a[i-1,]=NA
      delta_D_p[i,]=delta_p[i,]*8+10 #assume initial points are on MWL
      delta_D_p[i-1,]=NA
      delta_D_a[i,]=delta_a[i,]*8+10
      delta_D_a[i-1,]=NA
      delta_17_p[i,]=(exp(0.528*log(delta_p[i,]/1000+1)+0.000033)-1)*1000  #GMWL, Luz and Barkan, 2010
      delta_17_p[i-1,]=NA
      delta_17_a[i,]=(exp(0.528*log(delta_a[i,]/1000+1)+0.000033)-1)*1000  #GMWL, Luz and Barkan, 2010
      delta_17_a[i-1,]=NA
      #repeat for eddies (copypasta)
      delta_p_eddy[i,]=initial_delta_p_eddy
      delta_p_eddy[i-1,]=NA
      delta_a_eddy[i,]=initial_delta_a_eddy
      delta_a_eddy[i-1,]=NA
      delta_D_p_eddy[i,]=delta_p_eddy[i,]*8+10 #assume initial points are on MWL
      delta_D_p_eddy[i-1,]=NA
      delta_D_a_eddy[i,]=delta_a_eddy[i,]*8+10
      delta_D_a_eddy[i-1,]=NA
      delta_17_p_eddy[i,]=(exp(0.528*log(delta_p_eddy[i,]/1000+1)+0.000033)-1)*1000  #GMWL, Luz and Barkan, 2010
      delta_17_p_eddy[i-1,]=NA
      delta_17_a_eddy[i,]=(exp(0.528*log(delta_a_eddy[i,]/1000+1)+0.000033)-1)*1000  #GMWL, Luz and Barkan, 2010
      delta_17_a_eddy[i-1,]=NA
    } else {
      x=32 #km, assumed distance between points
      dx=x #meters, but the units don't matter. slightly >32km
      dw=w[i]-w[i-1]
      dw=ifelse(dw>0,0,dw) #gets rid of positive gradients (b/c we don't like explosions here)
      #Hendricks model (advection only)
      delta_inf[i,h] <- (Nd[i]*delta_et[i,h]-(1+Nd[i])*(a18[i]-1)*1000)/(a18[i]+Nd[i]*a18[i]-1)
      delta_a[i,h]=(delta_a[i-1,h]-delta_inf[i,h])*exp((a18[i]+a18[i]*Nd[i]-1)*(x/w[i])*(dw/dx))+delta_inf[i,h]
      delta_p[i,h] <-delta_a[i,h]+(a18[i]-1)*1000
      
      delta_D_inf[i,h] <- (Nd[i]*delta_D_et[i,h]-(1+Nd[i])*(aD[i]-1)*1000)/(aD[i]+Nd[i]*aD[i]-1)
      delta_D_a[i,h]=(delta_D_a[i-1,h]-delta_D_inf[i,h])*exp((aD[i]+aD[i]*Nd[i]-1)*(x/w[i])*(dw/dx))+delta_D_inf[i,h]
      delta_D_p[i,h] <-delta_D_a[i,h]+(aD[i]-1)*1000
      
      delta_17_inf[i,h] <- (Nd[i]*delta_17_et[i,h]-(1+Nd[i])*(a17[i]-1)*1000)/(a17[i]+Nd[i]*a17[i]-1)
      delta_17_a[i,h]=(delta_17_a[i-1,h]-delta_17_inf[i,h])*exp((a17[i]+a17[i]*Nd[i]-1)*(x/w[i])*(dw/dx))+delta_17_inf[i,h]
      delta_17_p[i,h] <-delta_17_a[i,h]+(a17[i]-1)*1000
      
      #Hendricks model (eddy diffusion only)
      delta_inf_eddy[i,h] <- (Nd[i]*delta_et_eddy[i,h]-(1+Nd[i])*(a18[i]-1)*1000)/(a18[i]+Nd[i]*a18[i]-1)
      delta_a_eddy[i,h]=(delta_a_eddy[i-1,h]-delta_inf_eddy[i,h])*exp((sqrt(a18[i]+a18[i]*Nd[i])-1)*(x/w[i])*(dw/dx))+delta_inf_eddy[i,h]
      delta_p_eddy[i,h] <-delta_a_eddy[i,h]+(a18[i]-1)*1000
      
      delta_D_inf_eddy[i,h] <- (Nd[i]*delta_D_et_eddy[i,h]-(1+Nd[i])*(aD[i]-1)*1000)/(aD[i]+Nd[i]*aD[i]-1)
      delta_D_a_eddy[i,h]=(delta_D_a_eddy[i-1,h]-delta_D_inf_eddy[i,h])*exp((sqrt(aD[i]+aD[i]*Nd[i])-1)*(x/w[i])*(dw/dx))+delta_D_inf_eddy[i,h]
      delta_D_p_eddy[i,h] <-delta_D_a_eddy[i,h]+(aD[i]-1)*1000
      
      delta_17_inf_eddy[i,h] <- (Nd[i]*delta_17_et_eddy[i,h]-(1+Nd[i])*(a17[i]-1)*1000)/(a17[i]+Nd[i]*a17[i]-1)
      delta_17_a_eddy[i,h]=(delta_17_a_eddy[i-1,h]-delta_17_inf_eddy[i,h])*exp((sqrt(a17[i]+a17[i]*Nd[i])-1)*(x/w[i])*(dw/dx))+delta_17_inf_eddy[i,h]
      delta_17_p_eddy[i,h] <-delta_17_a_eddy[i,h]+(a17[i]-1)*1000
    }  
  }          
}
delta_p.meanfromcoast=delta_p
delta_D_p.meanfromcoast=delta_D_p
delta_17_p.meanfromcoast=delta_17_p
delta_p_eddy.meanfromcoast=delta_p_eddy
delta_D_p_eddy.meanfromcoast=delta_D_p_eddy
delta_17_p_eddy.meanfromcoast=delta_17_p_eddy
distance.fake = seq(0,by=32,length=tracklength)

###calculate D-excess and 17O-excess:
Dxs.meanfromcoast=delta_D_p.meanfromcoast-8*delta_p.meanfromcoast #in permil
O17xs.meanfromcoast=(log(delta_17_p.meanfromcoast/1000+1)-0.528*log(delta_p.meanfromcoast/1000+1))*10^6 #in permeg
Dxs_eddy.meanfromcoast=delta_D_p_eddy.meanfromcoast-8*delta_p_eddy.meanfromcoast #in permil
O17xs_eddy.meanfromcoast=(log(delta_17_p_eddy.meanfromcoast/1000+1)-0.528*log(delta_p_eddy.meanfromcoast/1000+1))*10^6 #in permeg


##########model plots###########

library(fields)
ntracks=10
colors=tim.colors(ntracks)
colors_eddy=paste('gray',seq(5,86,length=10))

###plot Nd###
Nd=rowMeans(Nd.NA,na.rm=T)
pdf('Damkohler of each track.pdf',width=5,height=3)
plot(distance.all[,1],Nd.data[,1],pch=16,type='o',cex=0.5,xaxs="i",
     col=colors[1],xlab='Distance from "Coast" (km)',ylim=c(-5,20),
     ylab='Nd',main='Damkohler of each track')
for (j in 2:ntracks){
  points(distance.all[,j],Nd.data[,j],type="o",pch=16,cex=0.5,col=colors[j])
}
points(distance.fake,Nd,pch=16,type='o',cex=0.9,col='black')
legend('topleft',legend='average without\nnegatives',pch=16,col='black',cex=0.8)
dev.off()

###
###plot model-output d18O###

pdf('Model d18O with eddies.pdf',width=7,height=10)
par(mfrow=c(2,1))
#all evaporation
plot(distance.all[,1],delta_p.all[,1,1],type="o",pch=16,cex=0.5,xaxs="i",
     ylim=c(-35,0),col=colors[1],xlab='Distance from "Coast" (km)',
     ylab=expression(paste(delta^18,O[p],' (\u2030)')),main='All Evaporation')
for (j in 2:ntracks){
  points(distance.all[,j],delta_p.all[,j,1],type="o",pch=16,cex=0.5,col=colors[j])
}
for(j in 1:ntracks) {
  points(distance.all[,j],delta_p_eddy.all[,j,1],type='o',pch=15,cex=0.5,col=colors_eddy[j])
}
points(distance.fake,delta_p.meanfromcoast[,1],type='o',pch=16,col='black')
points(distance.fake,delta_p_eddy.meanfromcoast[,1],type='o',pch=15,col='gray50')
legend('topright',c('advection only','eddies only'),col=c('black','gray50'),pch=c(16,15))

#all tranpsiration
plot(distance.all[,1],delta_p.all[,1,11],type="o",pch=16,cex=0.5,xaxs="i",
     ylim=c(-35,0),col=colors[1],xlab='Distance from "Coast" (km)',
     ylab=expression(paste(delta^18,O[p],' (\u2030)')),main='All Transpiration')
for (j in 2:ntracks){
  points(distance.all[,j],delta_p.all[,j,11],type="o",pch=16,cex=0.5,col=colors[j])
}
for (j in 1:ntracks) {
  points(distance.all[,j],delta_p_eddy.all[,j,11],type='o',pch=15,cex=0.5,col=colors_eddy[j])
}
points(distance.fake,delta_p_eddy.meanfromcoast[,11],type='o',pch=15,col='gray50')
points(distance.fake,delta_p.meanfromcoast[,11],type='o',pch=16,col='black')
legend('topright',c('advection only','eddies only'),col=c('black','gray50'),pch=c(16,15))
dev.off()

###
###plot model-output Dxs###

pdf('Model Dxs with eddies.pdf',width=7,height=10)
par(mfrow=c(2,1))
#all evaporation
plot(distance.all[,1],Dxs.all[,1,1],type="o",pch=16,cex=0.5,xaxs="i",
     ylim=c(10,130),col=colors[1],xlab='Distance from "Coast" (km)',
     ylab='d-excess (\u2030)',main='All Evaporation')
for (j in 2:ntracks){
  points(distance.all[,j],Dxs.all[,j,1],type="o",pch=16,cex=0.5,col=colors[j])
}
for (j in 1:ntracks) {
  points(distance.all[,j],Dxs_eddy.all[,j,1],type='o',pch=15,cex=0.5,col=colors_eddy[j])
}
points(distance.fake,Dxs_eddy.meanfromcoast[,1],type='o',pch=15,col='gray50')
points(distance.fake,Dxs.meanfromcoast[,1],type='o',pch=16,col='black')
legend('topleft',c('advection only','eddies only'),col=c('black','gray50'),pch=c(16,15))

#all transpiration
plot(distance.all[,1],Dxs.all[,1,11],type="o",pch=16,cex=0.5,xaxs="i",
     ylim=c(-50,80),col=colors[1],xlab='Distance from "Coast" (km)',
     ylab='d-excess (\u2030)',main='All Transpiration')
for (j in 2:ntracks){
  points(distance.all[,j],Dxs.all[,j,11],type="o",pch=16,cex=0.5,col=colors[j])
}
for (j in 1:ntracks) {
  points(distance.all[,j],Dxs_eddy.all[,j,11],type='o',pch=15,cex=0.5,col=colors_eddy[j])
}
points(distance.fake,Dxs_eddy.meanfromcoast[,11],type='o',pch=15,col='gray50')
points(distance.fake,Dxs.meanfromcoast[,11],type='o',pch=16,col='black')
legend('topleft',c('advection only','eddies only'),col=c('black','gray50'),pch=c(16,15))
dev.off()

###
###plot model-output 17Oxs###

pdf('Model 17Oxs with eddies.pdf',width=7,height=15)
par(mfrow=c(3,1))
#all evaporation
plot(distance.all[,1],O17xs.all[,1,1],type="o",pch=16,cex=0.5,xaxs="i",
     ylim=c(-150,50),col=colors[1],xlab='Distance from "Coast" (km)',
     ylab='17O-excess (permeg)',main='All Evaporation')
for (j in 2:ntracks){
  points(distance.all[,j],O17xs.all[,j,1],type="o",pch=16,cex=0.5,col=colors[j])
}
for (j in 1:ntracks){
  points(distance.all[,j],O17xs_eddy.all[,j,1],type="o",pch=15,cex=0.5,col=colors_eddy[j])
}
points(distance.fake,O17xs_eddy.meanfromcoast[,1],type='o',pch=15,col='gray50')
points(distance.fake,O17xs.meanfromcoast[,1],type='o',pch=16,col='black')
legend('bottomleft',c('advection only','eddies only'),col=c('black','gray50'),pch=c(16,15))

#70% transpiration
plot(distance.all[,1],O17xs.all[,1,8],type="o",pch=16,cex=0.5,xaxs="i",
     ylim=c(-300,50),col=colors[1],xlab='Distance from "Coast" (km)',
     ylab='17O-excess (permeg)',main='70% Transpiration')
for (j in 2:ntracks){
  points(distance.all[,j],O17xs.all[,j,8],type="o",pch=16,cex=0.5,col=colors[j])
}
for (j in 1:ntracks){
  points(distance.all[,j],O17xs_eddy.all[,j,8],type="o",pch=15,cex=0.5,col=colors_eddy[j])
}
points(distance.fake,O17xs_eddy.meanfromcoast[,8],type='o',pch=15,col='gray50')
points(distance.fake,O17xs.meanfromcoast[,8],type='o',pch=16,col='black')
legend('bottomleft',c('advection only','eddies only'),col=c('black','gray50'),pch=c(16,15))

#all transpiration
plot(distance.all[,1],O17xs.all[,1,11],type="o",pch=16,cex=0.5,xaxs="i",
     ylim=c(-250,200),col=colors[1],xlab='Distance from "Coast" (km)',
     ylab='17O-excess (permeg)',main='All Transpiration')
for (j in 2:ntracks){
  points(distance.all[,j],O17xs.all[,j,11],type="o",pch=16,cex=0.5,col=colors[j])
}
for (j in 1:ntracks){
  points(distance.all[,j],O17xs_eddy.all[,j,11],type="o",pch=15,cex=0.5,col=colors_eddy[j])
}
points(distance.fake,O17xs_eddy.meanfromcoast[,11],type='o',pch=15,col='gray50')
points(distance.fake,O17xs.meanfromcoast[,11],type='o',pch=16,col='black')
legend('bottomleft',c('advection only','eddies only'),col=c('black','gray50'),pch=c(16,15))
dev.off()


####

########data time!########
####

moderndata=read.csv('modern water d18O dD - no USGS.csv',header=T)
head(moderndata)
library(dplyr)
library(geosphere)

#a few of the data points have NAs for lat/lon. take out those rows:
moderndata2=filter(moderndata,is.na(latitude)==F,is.na(longitude)==F)

##add a column for distance from "coast"
#where "coast" means first I&T (1991) point
#calculate three coastal distances then take the min of those as dist from coast
lati=moderndata$latitude[1] #39.43807
loni=moderndata$longitude[1] #-123.7982
lati2=40.002656
loni2=-124.024106
lati3=38.962633
loni3=-123.721038

latitudes=moderndata2$latitude
longitudes=moderndata2$longitude
distmin <- dist1 <- dist2 <- dist3 <- numeric()
for (i in 1:length(latitudes)){
  dist1[i]=distVincentyEllipsoid(c(longitudes[i],latitudes[i]),c(loni,lati))
  dist2[i]=distVincentyEllipsoid(c(longitudes[i],latitudes[i]),c(loni2,lati2))
  dist3[i]=distVincentyEllipsoid(c(longitudes[i],latitudes[i]),c(loni3,lati3))
  distmin[i]=min(dist1[i],dist2[i],dist3[i])/1000 #distance in km
}

moderndata3=cbind(moderndata2,distmin)
head(moderndata3)

###add another column for Dxs
moderndata.Dxs=mutate(moderndata3,Dxs=dD-8*d18O)
head(moderndata.Dxs)

##separate data by sample type
#rain!
Dxs_rain=filter(moderndata.Dxs,water.type=='rain')$Dxs
d18O_rain=filter(moderndata3,water.type=='rain')$d18O
dist_rain=filter(moderndata3,water.type=='rain')$distmin
lat_rain=filter(moderndata3,water.type=='rain')$lat
lon_rain=filter(moderndata3,water.type=='rain')$lon
#snow!
Dxs_snow=filter(moderndata.Dxs,water.type=='snow')$Dxs
d18O_snow=filter(moderndata3,water.type=='snow')$d18O
dist_snow=filter(moderndata3,water.type=='snow')$distmin
lat_snow=filter(moderndata3,water.type=='snow')$lat
lon_snow=filter(moderndata3,water.type=='snow')$lon
#rivers and creeks!
Dxs_river=filter(moderndata.Dxs,water.type=='river')$Dxs
d18O_river=filter(moderndata3,water.type=='river')$d18O
dist_river=filter(moderndata3,water.type=='river')$distmin
lat_river=filter(moderndata3,water.type=='river')$lat
lon_river=filter(moderndata3,water.type=='river')$lon
#groundwater!
Dxs_gw=filter(moderndata.Dxs,water.type=='GW')$Dxs
d18O_gw=filter(moderndata3,water.type=='GW')$d18O
dist_gw=filter(moderndata3,water.type=='GW')$distmin
lat_gw=filter(moderndata3,water.type=='GW')$lat
lon_gw=filter(moderndata3,water.type=='GW')$lon
#lakes!
Dxs_lake=filter(moderndata.Dxs,water.type=='lake')$Dxs
d18O_lake=filter(moderndata3,water.type=='lake')$d18O
dist_lake=filter(moderndata3,water.type=='lake')$distmin
lat_lake=filter(moderndata3,water.type=='lake')$lat
lon_lake=filter(moderndata3,water.type=='lake')$lon


####

######plots with data#########
####


###make map of where the data are; color points by isotopic value
library(maps)
library(oce)
lim=c(-22,0)
ncol=100
col=oceColorsJet(ncol)
# Draw colorbar for d18O
pdf('map of data.pdf',width=10, height=7)
par(mfrow=c(1,2))
drawPalette(zlim=lim,col=col,
            zlab=expression(paste(delta^18,"O"[precip])),
            las=1,cex.lab=1.2) # adjust palette position with mai


par(mar=c(2.5,2.5,1,1),mgp=c(2,0.7,0)) # tighten margins
map('state',xlim=c(min(moderndata2$lon)-1,max(moderndata2$lon)),
    ylim=c(min(moderndata2$lat)-1,max(moderndata2$lat)))
for (i in 1:10) {  #add storm tracks
  points(lon.data[,i],lat.data[,i],col="blue",type="o",pch=1,cex=0.5)
}
points(lon_rain,lat_rain,pch=21,cex=0.7,
       bg=col[rescale(x=d18O_rain,xlow=lim[1],
                      xhigh=lim[2],rlow=1,rhigh=ncol)])
points(lon_snow,lat_snow,pch=22,cex=0.7,
       bg=col[rescale(x=d18O_snow,xlow=lim[1],
                      xhigh=lim[2],rlow=1,rhigh=ncol)])
points(lon_river,lat_river,pch=23,cex=0.7,
       bg=col[rescale(x=d18O_river,xlow=lim[1],
                      xhigh=lim[2],rlow=1,rhigh=ncol)])
points(lon_gw,lat_gw,pch=24,cex=0.7,
       bg=col[rescale(x=d18O_gw,xlow=lim[1],
                      xhigh=lim[2],rlow=1,rhigh=ncol)])
points(lon_lake,lat_lake,pch=25,cex=0.7,
       bg=col[rescale(x=d18O_lake,xlow=lim[1],
                      xhigh=lim[2],rlow=1,rhigh=ncol)])
legend('bottomleft',pch=c(21,22,23,24,25),pt.bg='cyan',cex=0.7,bty='n',
       legend=c('rain','snow','rivers/creeks','groundwater','lakes'))
dev.off()


###
###plot model-output d18O with data###

pdf('d18O with eddies.pdf',width=7,height=15)
par(mfrow=c(3,1))
#all evaporation
plot(distance.all[,1],delta_p.all[,1,1],type="o",pch=16,cex=0.5,xaxs="i",
     ylim=c(-35,0),col=colors[1],xlab='Distance from "Coast" (km)',
     ylab=expression(paste(delta^18,O[p],' (\u2030)')),main='All Evaporation')
for (j in 2:ntracks){
  points(distance.all[,j],delta_p.all[,j,1],type="o",pch=16,cex=0.5,col=colors[j])
}
for(j in 1:ntracks) {
  points(distance.all[,j],delta_p_eddy.all[,j,1],type='o',pch=15,cex=0.5,col=colors_eddy[j])
}
points(distance.fake,delta_p.meanfromcoast[,1],type='o',pch=16,col='black')
points(distance.fake,delta_p_eddy.meanfromcoast[,1],type='o',pch=15,col='gray50')
points(dist_rain,d18O_rain,pch=21,cex=0.8,col='black',bg='blue')
points(dist_snow,d18O_snow,cex=0.8,pch=22,col='black',bg='blue')
points(dist_gw,d18O_gw,cex=0.8,pch=23,col='black',bg='blue')
points(dist_river,d18O_river,cex=0.8,pch=24,col='black',bg='blue')
points(dist_lake,d18O_lake,cex=0.8,pch=24,col='black',bg='blue')
legend('bottomleft',c('advection','eddies','rain','snow','groundwater','rivers/streams','lakes'),
       cex=0.8,pch=c(16,15,21,22,23,24,25),col=c('black','gray50',rep('black',5)),pt.bg='blue')

#70% tranpsiration
plot(distance.all[,1],delta_p.all[,1,8],type="o",pch=16,cex=0.5,xaxs="i",
     ylim=c(-35,0),col=colors[1],xlab='Distance from "Coast" (km)',
     ylab=expression(paste(delta^18,O[p],' (\u2030)')),main='70% Transpiration')
for (j in 2:ntracks){
  points(distance.all[,j],delta_p.all[,j,8],type="o",pch=16,cex=0.5,col=colors[j])
}
for (j in 1:ntracks) {
  points(distance.all[,j],delta_p_eddy.all[,j,8],type='o',pch=15,cex=0.5,col=colors_eddy[j])
}
points(distance.fake,delta_p_eddy.meanfromcoast[,8],type='o',pch=15,col='gray50')
points(distance.fake,delta_p.meanfromcoast[,8],type='o',pch=16,col='black')
points(dist_rain,d18O_rain,pch=21,cex=0.8,col='black',bg='blue')
points(dist_snow,d18O_snow,cex=0.8,pch=22,col='black',bg='blue')
points(dist_gw,d18O_gw,cex=0.8,pch=23,col='black',bg='blue')
points(dist_river,d18O_river,cex=0.8,pch=24,col='black',bg='blue')
points(dist_lake,d18O_lake,cex=0.8,pch=24,col='black',bg='blue')
legend('bottomleft',c('advection','eddies','rain','snow','groundwater','rivers/streams','lakes'),
       cex=0.8,pch=c(16,15,21,22,23,24,25),col=c('black','gray50',rep('black',5)),pt.bg='blue')

#all tranpsiration
plot(distance.all[,1],delta_p.all[,1,11],type="o",pch=16,cex=0.5,xaxs="i",
     ylim=c(-35,0),col=colors[1],xlab='Distance from "Coast" (km)',
     ylab=expression(paste(delta^18,O[p],' (\u2030)')),main='All Transpiration')
for (j in 2:ntracks){
  points(distance.all[,j],delta_p.all[,j,11],type="o",pch=16,cex=0.5,col=colors[j])
}
for (j in 1:ntracks) {
  points(distance.all[,j],delta_p_eddy.all[,j,11],type='o',pch=15,cex=0.5,col=colors_eddy[j])
}
points(distance.fake,delta_p_eddy.meanfromcoast[,11],type='o',pch=15,col='gray50')
points(distance.fake,delta_p.meanfromcoast[,11],type='o',pch=16,col='black')
points(dist_rain,d18O_rain,pch=21,cex=0.8,col='black',bg='blue')
points(dist_snow,d18O_snow,cex=0.8,pch=22,col='black',bg='blue')
points(dist_gw,d18O_gw,cex=0.8,pch=23,col='black',bg='blue')
points(dist_river,d18O_river,cex=0.8,pch=24,col='black',bg='blue')
points(dist_lake,d18O_lake,cex=0.8,pch=24,col='black',bg='blue')
legend('bottomleft',c('advection','eddies','rain','snow','groundwater','rivers/streams','lakes'),
       cex=0.8,pch=c(16,15,21,22,23,24,25),col=c('black','gray50',rep('black',5)),pt.bg='blue')
dev.off()

###
###plot model-output Dxs with data###

pdf('Dxs with eddies.pdf',width=7,height=15)
par(mfrow=c(3,1))
#all evaporation
plot(distance.all[,1],Dxs.all[,1,1],type="o",pch=16,cex=0.5,xaxs="i",
     ylim=c(-70,130),xlim=c(0,1130),col=colors[1],xlab='Distance from "Coast" (km)',
     ylab='d-excess (\u2030)',main='All Evaporation')
for (j in 2:ntracks){
  points(distance.all[,j],Dxs.all[,j,1],type="o",pch=16,cex=0.5,col=colors[j])
}
for (j in 1:ntracks) {
  points(distance.all[,j],Dxs_eddy.all[,j,1],type='o',pch=15,cex=0.5,col=colors_eddy[j])
}
points(distance.fake,Dxs_eddy.meanfromcoast[,1],type='o',pch=15,col='gray50')
points(distance.fake,Dxs.meanfromcoast[,1],type='o',pch=16,col='black')
points(dist_rain,Dxs_rain,pch=21,cex=0.8,col='black',bg='blue')
points(dist_snow,Dxs_snow,cex=0.8,pch=22,col='black',bg='blue')
points(dist_gw,Dxs_gw,cex=0.8,pch=23,col='black',bg='blue')
points(dist_river,Dxs_river,cex=0.8,pch=24,col='black',bg='blue')
points(dist_lake,Dxs_lake,cex=0.8,pch=24,col='black',bg='blue')
legend('topright',c('advection','eddies','rain','snow','groundwater','rivers/streams','lakes'),
       cex=0.8,pch=c(16,15,21,22,23,24,25),col=c('black','gray50',rep('black',5)),pt.bg='blue')

#70% transpiration
plot(distance.all[,1],Dxs.all[,1,8],type="o",pch=16,cex=0.5,xaxs="i",
     ylim=c(-70,130),xlim=c(0,1130),col=colors[1],xlab='Distance from "Coast" (km)',
     ylab='d-excess (\u2030)',main='70% Transpiration')
for (j in 2:ntracks){
  points(distance.all[,j],Dxs.all[,j,8],type="o",pch=16,cex=0.5,col=colors[j])
}
for (j in 1:ntracks) {
  points(distance.all[,j],Dxs_eddy.all[,j,8],type='o',pch=15,cex=0.5,col=colors_eddy[j])
}
points(distance.fake,Dxs_eddy.meanfromcoast[,8],type='o',pch=15,col='gray50')
points(distance.fake,Dxs.meanfromcoast[,8],type='o',pch=16,col='black')
points(dist_rain,Dxs_rain,pch=21,cex=0.8,col='black',bg='blue')
points(dist_snow,Dxs_snow,cex=0.8,pch=22,col='black',bg='blue')
points(dist_gw,Dxs_gw,cex=0.8,pch=23,col='black',bg='blue')
points(dist_river,Dxs_river,cex=0.8,pch=24,col='black',bg='blue')
points(dist_lake,Dxs_lake,cex=0.8,pch=24,col='black',bg='blue')
legend('topright',c('advection','eddies','rain','snow','groundwater','rivers/streams','lakes'),
       cex=0.8,pch=c(16,15,21,22,23,24,25),col=c('black','gray50',rep('black',5)),pt.bg='blue')

#all transpiration
plot(distance.all[,1],Dxs.all[,1,11],type="o",pch=16,cex=0.5,xaxs="i",
     ylim=c(-70,130),xlim=c(0,1130),col=colors[1],xlab='Distance from "Coast" (km)',
     ylab='d-excess (\u2030)',main='All Transpiration')
for (j in 2:ntracks){
  points(distance.all[,j],Dxs.all[,j,11],type="o",pch=16,cex=0.5,col=colors[j])
}
for (j in 1:ntracks) {
  points(distance.all[,j],Dxs_eddy.all[,j,11],type='o',pch=15,cex=0.5,col=colors_eddy[j])
}
points(distance.fake,Dxs_eddy.meanfromcoast[,11],type='o',pch=15,col='gray50')
points(distance.fake,Dxs.meanfromcoast[,11],type='o',pch=16,col='black')
points(dist_rain,Dxs_rain,pch=21,cex=0.8,col='black',bg='blue')
points(dist_snow,Dxs_snow,cex=0.8,pch=22,col='black',bg='blue')
points(dist_gw,Dxs_gw,cex=0.8,pch=23,col='black',bg='blue')
points(dist_river,Dxs_river,cex=0.8,pch=24,col='black',bg='blue')
points(dist_lake,Dxs_lake,cex=0.8,pch=24,col='black',bg='blue')
legend('topright',c('advection','eddies','rain','snow','groundwater','rivers/streams','lakes'),
       cex=0.8,pch=c(16,15,21,22,23,24,25),col=c('black','gray50',rep('black',5)),pt.bg='blue')
dev.off()


####

######Ingraham and Taylor, 1991 plots#######

library(geosphere)
library(dplyr)

IandTdata=read.csv("Ingraham and Taylor 1991 Traverse II.csv")
head(IandTdata)
#eliminate NA lat/lons
IandTdata2=filter(IandTdata,is.na(latitude)==F,is.na(longitude)==F)
###add another column for Dxs
IandTdata.Dxs=mutate(IandTdata2,Dxs=dD-8*d18O)

##separate data by sample type
#rain!
IT_Dxs_rain=filter(IandTdata.Dxs,water.type=='rain')$Dxs
IT_d18O_rain=filter(IandTdata.Dxs,water.type=='rain')$d18O
IT_dist_rain=filter(IandTdata.Dxs,water.type=='rain')$dist
IT_lat_rain=filter(IandTdata.Dxs,water.type=='rain')$lat
IT_lon_rain=filter(IandTdata.Dxs,water.type=='rain')$lon
#snow!
IT_Dxs_snow=filter(IandTdata.Dxs,water.type=='snow')$Dxs
IT_d18O_snow=filter(IandTdata.Dxs,water.type=='snow')$d18O
IT_dist_snow=filter(IandTdata.Dxs,water.type=='snow')$dist
IT_lat_snow=filter(IandTdata.Dxs,water.type=='snow')$lat
IT_lon_snow=filter(IandTdata.Dxs,water.type=='snow')$lon
#groundwater!
IT_Dxs_gw=filter(IandTdata.Dxs,water.type=='spring.GW')$Dxs
IT_d18O_gw=filter(IandTdata.Dxs,water.type=='spring.GW')$d18O
IT_dist_gw=filter(IandTdata.Dxs,water.type=='spring.GW')$dist
IT_lat_gw=filter(IandTdata.Dxs,water.type=='spring.GW')$lat
IT_lon_gw=filter(IandTdata.Dxs,water.type=='spring.GW')$lon

###
###plot model-output d18O with I&T data###

pdf('model d18O with I&T data.pdf',width=5,height=10)
par(mfrow=c(3,1))
#all evaporation
plot(distance.fake,delta_p.meanfromcoast[,1],type='o',pch=16,col='black',
     xaxs="i", ylim=c(-35,0),xlab='Distance from "Coast" (km)',
     ylab=expression(paste(delta^18,O[p],' (\u2030)')),main='All Evaporation')
points(distance.fake,delta_p_eddy.meanfromcoast[,1],type='o',pch=15,col='gray50')
points(IT_dist_rain,IT_d18O_rain,pch=21,cex=0.8,col='black',bg='blue')
points(IT_dist_snow,IT_d18O_snow,cex=0.8,pch=22,col='black',bg='light blue')
points(IT_dist_gw,IT_d18O_gw,cex=0.8,pch=23,col='black',bg='cyan')
legend('bottomleft',c('advection','eddies','rain','snow','groundwater'),
       cex=0.8,pch=c(16,15,21,22,23),col=c('black','gray50',rep('black',3)),
       pt.bg=c(NA,NA,'blue','light blue','cyan'))

#70% tranpsiration
plot(distance.fake,delta_p.meanfromcoast[,8],type='o',pch=16,col='black',
     xaxs="i", ylim=c(-35,0),xlab='Distance from "Coast" (km)',
     ylab=expression(paste(delta^18,O[p],' (\u2030)')),main='70% Transpiration')
points(distance.fake,delta_p_eddy.meanfromcoast[,8],type='o',pch=15,col='gray50')
points(IT_dist_rain,IT_d18O_rain,pch=21,cex=0.8,col='black',bg='blue')
points(IT_dist_snow,IT_d18O_snow,cex=0.8,pch=22,col='black',bg='light blue')
points(IT_dist_gw,IT_d18O_gw,cex=0.8,pch=23,col='black',bg='cyan')
legend('bottomleft',c('advection','eddies','rain','snow','groundwater'),
       cex=0.8,pch=c(16,15,21,22,23),col=c('black','gray50',rep('black',3)),
       pt.bg=c(NA,NA,'blue','light blue','cyan'))

#all tranpsiration
plot(distance.fake,delta_p.meanfromcoast[,11],type='o',pch=16,col='black',
      xaxs="i", ylim=c(-35,0),xlab='Distance from "Coast" (km)',
     ylab=expression(paste(delta^18,O[p],' (\u2030)')),main='All Transpiration')
points(distance.fake,delta_p_eddy.meanfromcoast[,11],type='o',pch=15,col='gray50')
points(distance.fake,delta_p.meanfromcoast[,11],type='o',pch=16,col='black')
points(IT_dist_rain,IT_d18O_rain,pch=21,cex=0.8,col='black',bg='blue')
points(IT_dist_snow,IT_d18O_snow,cex=0.8,pch=22,col='black',bg='light blue')
points(IT_dist_gw,IT_d18O_gw,cex=0.8,pch=23,col='black',bg='cyan')
legend('bottomleft',c('advection','eddies','rain','snow','groundwater'),
       cex=0.8,pch=c(16,15,21,22,23),col=c('black','gray50',rep('black',3)),
       pt.bg=c(NA,NA,'blue','light blue','cyan'))
dev.off()


###
###plot model-output Dxs with I&T data###

pdf('model Dxs with I&T data.pdf',width=5,height=10)
par(mfrow=c(3,1))
#all evaporation
plot(distance.fake,Dxs.meanfromcoast[,1],type='o',pch=16,col='black',
     xaxs="i", ylim=c(-5,130),xlab='Distance from "Coast" (km)',
     ylab=expression(paste(Dxs,' (\u2030)')),main='All Evaporation')
points(distance.fake,Dxs_eddy.meanfromcoast[,1],type='o',pch=15,col='gray50')
points(IT_dist_rain,IT_Dxs_rain,pch=21,cex=0.8,col='black',bg='blue')
points(IT_dist_snow,IT_Dxs_snow,cex=0.8,pch=22,col='black',bg='light blue')
points(IT_dist_gw,IT_Dxs_gw,cex=0.8,pch=23,col='black',bg='cyan')
legend('topleft',c('advection','eddies','rain','snow','groundwater'),
       cex=0.8,pch=c(16,15,21,22,23),col=c('black','gray50',rep('black',3)),
       pt.bg=c(NA,NA,'blue','light blue','cyan'))

#70% tranpsiration
plot(distance.fake,Dxs.meanfromcoast[,8],type='o',pch=16,col='black',
     xaxs="i", ylim=c(-5,130),xlab='Distance from "Coast" (km)',
     ylab=expression(paste(Dxs,' (\u2030)')),main='70% Transpiration')
points(distance.fake,Dxs_eddy.meanfromcoast[,8],type='o',pch=15,col='gray50')
points(IT_dist_rain,IT_Dxs_rain,pch=21,cex=0.8,col='black',bg='blue')
points(IT_dist_snow,IT_Dxs_snow,cex=0.8,pch=22,col='black',bg='light blue')
points(IT_dist_gw,IT_Dxs_gw,cex=0.8,pch=23,col='black',bg='cyan')
legend('topleft',c('advection','eddies','rain','snow','groundwater'),
       cex=0.8,pch=c(16,15,21,22,23),col=c('black','gray50',rep('black',3)),
       pt.bg=c(NA,NA,'blue','light blue','cyan'))

#all tranpsiration
plot(distance.fake,Dxs.meanfromcoast[,11],type='o',pch=16,col='black',
     xaxs="i", ylim=c(-5,130),xlab='Distance from "Coast" (km)',
     ylab=expression(paste(Dxs,' (\u2030)')),main='All Transpiration')
points(distance.fake,Dxs_eddy.meanfromcoast[,11],type='o',pch=15,col='gray50')
points(IT_dist_rain,IT_Dxs_rain,pch=21,cex=0.8,col='black',bg='blue')
points(IT_dist_snow,IT_Dxs_snow,cex=0.8,pch=22,col='black',bg='light blue')
points(IT_dist_gw,IT_Dxs_gw,cex=0.8,pch=23,col='black',bg='cyan')
legend('topleft',c('advection','eddies','rain','snow','groundwater'),
       cex=0.8,pch=c(16,15,21,22,23),col=c('black','gray50',rep('black',3)),
       pt.bg=c(NA,NA,'blue','light blue','cyan'))
dev.off()


####

######resp-MAP relationships#######

#here we have three different empirical, linear SR-MAP relationships


####upload data
# Cotton and Sheldon 2012 (GSA Bulletin)
Cotton12=read.csv("Cotton 2012.csv",header=T,skip=1)
head(Cotton12)

# Cotton et al., 2013 (Chemical Geology)
Cotton13=read.csv("Cotton 2013.csv",header=T,skip=1)
head(Cotton13)

# Convert to Soil Respiration Rates
#assuming soil temp=25C (in CO2convert) because that's how Jeremy does it
Cotton12=data.frame(Cotton12,sapply(Cotton12$CO2ppmv,CO2convert))
Cotton12=data.frame(Cotton12,((Ds*Cotton12[,3])/(L*z-z^2/2))*L*12*1e4*60*60*24*365) # in gC/m2/year
colnames(Cotton12)=c("SoilCO2","MAP","SoilCO2moles","SoilResp")
head(Cotton12)
#same for Cotton et al, 2013:
Cotton13=data.frame(Cotton13,sapply(Cotton13$CO2ppmv,CO2convert))
Cotton13=data.frame(Cotton13,((Ds*Cotton13[,3])/(L*z-z^2/2))*L*12*1e4*60*60*24*365) # Calculates soil respiration in gC/m2/year
colnames(Cotton13)=c("SoilCO2","MAP","SoilCO2moles","SoilResp")
head(Cotton13)

# Calculate linear regression
Cotton12lm=lm((Cotton12$SoilResp)~Cotton12$MAP)
summary(Cotton12lm)
Cot12_intercept=Cotton12lm$coef[1]
Cot12_slope=Cotton12lm$coef[2]
#same for Cotton et al, 2013
Cotton13lm=lm((Cotton13$SoilResp)~Cotton13$MAP)
summary(Cotton13lm)
Cot13_intercept=Cotton13lm$coef[1]
Cot13_slope=Cotton13lm$coef[2]

#regression line for Raich and Schlesinger, 1992:
RS_intercept=155
RS_slope=0.391

#use these values if you don't feel like running the regressions
# Cot12_slope=1.010897
# Cot12_intercept=-48.0897
# Cot13_slope=1.712771
# Cot13_intercept=-71.87511


########soil carbonate function##########


#Soil Carbonate d13C
soilcarb=function(CO2ppm_atm,MAP,Temp=298.15,MAPSR="C12",z=50,delta_13_a=-7,delta_13_o=-27){ 
  #MAP in mm/yr
  #Temp is soil temperature IN KELVIN
  #MAPSR indicates the regression to use
  
  TempC=Temp-273.15
  
  if (MAPSR=="C12") {SR=MAP*Cot12_slope+Cot12_intercept}
  if (MAPSR=="C13") {SR=MAP*Cot13_slope+Cot13_intercept}
  if (MAPSR=="RS") {SR=MAP*RS_slope+RS_intercept}
  
  prod=SR/12/365/24/60/60/1e4/L # Converts production (g/m2/yr) to (moles/cm3/s)
  
  CO2_atm=CO2convert(CO2ppm_atm,TempC) # Converts CO2 (ppm) to CO2 (moles/cm3)
  
  Sz=(prod/Ds)*(L*z-z^2/2) #S(z) Soil-respired CO2 concentration (moles/cm3)
  delta_13_s=((CO2_atm/Sz)*delta_13_a+1.0044*delta_13_o+4.4)/(1+CO2_atm/Sz) # Soil d13C value
  
  alpha=fract(TempC)
  
  dcarb=delta_13_s-alpha # Soil carbonate, assuming equilibrium
  
  # Outputs a vector of soil-respired CO2 (ppm), soil d13C, carbonate d13C, rounded to two decimal places
  result=c(round(CO2backconvert(Sz),3),round(delta_13_s,2),round(dcarb,2))
  names(result)=c("Soil CO2","Soil d13C","Carbonate d13C")
  result
}

#soilcarb(CO2ppm_atm=300,MAP=300)
#soilcarb(CO2ppm_atm=300,MAP=1000)
#soilcarb(CO2ppm_atm=1000,MAP=300)


#######soil carbon profile calculations########

###get MAP of "average from coast" track

#turn all of the inputs into one storm track:
#pre for pre-averaging
MAP.pre=as.matrix(MAP.data) #mean annual precip in mm/year
tsoil.pre=tsoil.data
ntracks.pre=dim(MAP.pre)[2]
ntracks=1 #number of storm tracks
tracklength=dim(MAP.pre)[1] #length of storm track

#if land.mask=0, cut off the previous data point and add NA to the end
#cutting off previous allows for initiation just before land, consistent with sections above
#adding NA at the end will allow for averaging of all the non-NAs for farther out on track
for (i in 2:tracklength) {
  for (j in 1:ntracks.pre) {
    if (land.mask.data[i,j]==0) {
      MAP.pre[,j]=c(MAP.pre[-(i-1),j],NA)
    } 
  }
}
#same for soil temperatures
for (i in 2:tracklength) {
  for (j in 1:ntracks.pre) {
    for (k in 1:length(zplot)) {
      if (land.mask.data[i,j]==0) {
        tsoil.pre[,j,k]=c(tsoil.pre[-(i-1),j,k],NA)
      } 
    }
  }
}

#now average these. this average is based on distance from coast!
MAP=rowMeans(MAP.pre,na.rm=T)
tsoil=apply(tsoil.pre,c(1,3),mean,na.rm=T) #28x101

###soil CO2 as a function of depth
soilCO2.1=array(dim=c(tracklength,length(zplot),3)) #for Cotton 2012
soilCO2.2=array(dim=c(tracklength,length(zplot),3)) #for Cotton 2013
soilCO2.3=array(dim=c(tracklength,length(zplot),3)) #for Raich and Schlesinger 
#currently assuming Temp=25C (function default)
for (i in 1:tracklength) {
  for (j in 1:length(zplot)) {
    soilCO2.1[i,j,]=soilcarb(CO2ppm_atm=pCO2ppm_atm,MAP=MAP[i],Temp=tsoil[i,j],z=zplot[j],delta_13_a=delta_13_atm)
    soilCO2.2[i,j,]=soilcarb(CO2ppm_atm=pCO2ppm_atm,MAP=MAP[i],Temp=tsoil[i,j],MAPSR='C13',z=zplot[j],delta_13_a=delta_13_atm)
    soilCO2.3[i,j,]=soilcarb(CO2ppm_atm=pCO2ppm_atm,MAP=MAP[i],Temp=tsoil[i,j],MAPSR='RS',z=zplot[j],delta_13_a=delta_13_atm)  
  }
}
#names(soilCO2[i,j,])=c("Soil CO2","Soil d13C","Carbonate d13C")
#^names of the three outputs of soilcarb() function


####soil CO2 concentration (combination of soil-respired and atmospheric)
#necessary for soil water model
#Cerling 1991:
CO2ppm_soil.1=soilCO2.1[,,1]+pCO2ppm_atm #for Cotton 2012
CO2ppm_soil.2=soilCO2.2[,,1]+pCO2ppm_atm #for Cotton 2013
CO2ppm_soil.3=soilCO2.3[,,1]+pCO2ppm_atm #for Raich and Schlesinger
#28 tracks, maxdepth+1 depths


######plot model soil carbonate profiles########

#Soil-Respired CO2 at beginning of track (arbitrary)
plot(soilCO2.1[1,,1],zplot,ylim=c(max(zplot),0),xlim=c(0,12000),type='l',yaxs='i',
     xlab=expression(paste('Soil-Respired ',CO[2],' (ppm)')),ylab='depth (cm)')
points(soilCO2.2[1,,1],zplot,type='l',col='red')
points(soilCO2.3[1,,1],zplot,type='l',col='blue')
legend('topright',legend=c('Cotton12','Cotton13','S&R'),col=c('black','red','blue'),lty=1)

#Soil d13C 
plot(soilCO2.1[1,,2],zplot,ylim=c(max(zplot),0),type='l',yaxs='i',
     xlab=expression(paste('Soil ',delta^13,C[VPDB],' (\u2030)')),ylab='depth (cm)')
points(soilCO2.2[1,,2],zplot,type='l',col='red')
points(soilCO2.3[1,,2],zplot,type='l',col='blue')
legend('bottomright',legend=c('Cotton12','Cotton13','S&R'),col=c('black','red','blue'),lty=1)

#Carbonate d13C
plot(soilCO2.1[1,,3],zplot,ylim=c(max(zplot),0),type='l',yaxs='i',
     xlab=expression(paste('Soil Carbonate',delta^13,C[VPDB],' (\u2030)')),ylab='depth (cm)')
points(soilCO2.2[1,,3],zplot,type='l',col='red')
points(soilCO2.3[1,,3],zplot,type='l',col='blue')
legend('bottomright',legend=c('Cotton12','Cotton13','S&R'),col=c('black','red','blue'),lty=1)



##########soil water component##########
####

# The soil profile is calculated using the equations from 
# Zimmerman et al. (1967) as outlined in Chamberlain et al. (2014)

###NOTE this is for ONE storm track (in this case the "average from coast")
#model inputs!
w.pre=as.matrix(Prw.data) #precipitatble water, kg/m^2
tracklength=length(w.pre[,1]) #length of storm tracks
ntracks.pre=dim(w.pre)[2]
ntracks=1 #number of storm tracks
Tsurf.pre=as.matrix(Temp_surf.data) #JJASO ground surface temperature (K)
tsoil.pre=tsoil.data
rh.pre=as.matrix(rhum.data/100) #relative humidity (fraction)
evap.pre=as.matrix(evap.data)  #evap rate in kg/m^2/month (mm/month over a given area)
#if land.mask=0, cut off the previous data point and add NA to the end
#cutting off previous allows for initiation just before land, consistent with sections above
#adding NA at the end will allow for averaging of all the non-NAs for farther out on track
for (j in 1:ntracks.pre) {
  for (i in 2:tracklength) {
    if (land.mask.data[i,j]==0) {
      Tsurf.pre[,j]=c(Tsurf.pre[-(i-1),j],NA)
      rh.pre[,j]=c(rh.pre[-(i-1),j],NA)
      evap.pre[,j]=c(evap.pre[-(i-1),j],NA)
    }
  }
}
#same for soil temperatures
for (i in 2:tracklength) {
  for (j in 1:ntracks.pre) {
    for (k in 1:length(zplot)) {
      if (land.mask.data[i,j]==0) {
        tsoil.pre[,j,k]=c(tsoil.pre[-(i-1),j,k],NA)
      } 
    }
  }
}
#now average these. this average is based on distance from coast!
TsurfK=rowMeans(Tsurf.pre,na.rm=T) #surface temp in Kelvin
rh=rowMeans(rh.pre,na.rm=T) #rel. humidity as fraction
tsoil=apply(tsoil.pre,c(1,3),mean,na.rm=T) #soil temp with depth in Kelvin
evap.mm.month=rowMeans(evap.pre,na.rm=T) #evap rate in mm/month
#convert evap rate to m/sec
evap.m.sec=evap.mm.month*(3.8580247e-10) #evap rate in m/sec
E=evap.m.sec
zstar=((p*tau*d)/E)/100 #in cm depth (dividing by 100 necessary for correct units)


###functions to make carb and water profiles###

# FUNCTION TO CALCULATE d18O (VSMOW) in water to d18O (VPDB) in carbonate ###
water2carb.oxygen <- function(delta18_precip = -15,TempK,pH=8,pCO2ppm=2500){
  #TempK should be soil temperature!
  #Carbon section Equil gas Zhang et al. (1995); CaCO3 Bottinga (1968)
  #Equilibrium for carbon isotopes
  TempC=TempK-273.15
  DeltaC1 = -(0.1141*TempC) + 10.78    #Bicarb-CO2g
  DeltaC2 = -(0.052*TempC) + 7.22    #Carb-CO2g
  DeltaC3 = (0.0049*TempC) - 1.31    #CO2aq-CO2g
  DeltaC4 = -2.4612 + (7.6663*10^3/TempK) - (2.988*10^6/TempK^2) #CO2g-CaCO3s
  #Kinetic for carbon isotopes
  DeltaCK = -0.81 #Kinectic fractionation CO2aq-CO2gas @ 21C for 5C use -0.95
  #Oxygen isotope section Equil gas Beck et al. (2005); CaCO3-H2O Kim and O'Neil (1997)
  DeltaO1 = (2.59*10^6/TempK^2) + 1.89    #Bicarb-H2O
  DeltaO2 = (2.39*10^6/TempK^2) - 2.7    #Carb-H2O
  DeltaO3 = (2.52*10^6/TempK^2) + 12.12    #CO2aq-H2O
  DeltaO4 = (18.03*10^3/TempK) - 32.42    #CaCO3-H2O
  DeltaO5 = DeltaO1 - DeltaO2             #Bicarb-Carb
  # Constants needed for calculation of carbonate species
  Ko = 0.032  #Henrys Law constant for CO2gas to CO2aq units mol/(kg*atm)
  pK1 = 3404.71/TempK + 0.032786*TempK - 14.8435 # Harned and Davis (1943) 
  pK2 = 2902.39/TempK + 0.02379*TempK - 6.4980 # Harned and Davis (1943) 
  K1 = 10^(-pK1)
  K2 = 10^(-pK2)
  # Calculate CO2aq from Henry's Law
  pCO2 = pCO2ppm/(1*10^6) #Conversion to atm
  CO2 = Ko * pCO2  # CO2 is in mol/kg PCO2 is in atm
  # Calculate carbonate species give CO2 conc and pH
  H = 10^(-pH)
  DIC = CO2 * (1 + K1/H + K1*K2/H^2)  
  BICARB = DIC / (1 + H/K1 + K2/H)
  CARB = DIC / (1 + H/K2 + H^2/K1*K2)
  # Calculate the oxygen isotope compostion of calcite with speciation
  XBICARB = BICARB/(BICARB + CARB)
  Delta18Ocalcite = XBICARB * DeltaO5 + delta18_precip + DeltaO2
  # function outputs a vector of Delta18Ocalcite based on each temperature value used
  #Delta18Ocalcite #VSMOW!
  #convert to VPDB:
  d18O_calciteVPDB=(Delta18Ocalcite-30.91)/1.03091
  d18O_calciteVPDB
}
# water2carb.oxygen(-5,284)

# FUNCTION TO CALCULATE SOIL WATER PROFILE ###
# spits out profile of soil water d18O given delta18_precip, zstar, maxdepth, and tempK_surface
soil.profile <- function(delta18_precip = -15.0,TsurfK,rh,zstar = 50,maxdepth=200){
  #defaults to zstar=50, maxdepth=200cm if none given *in function call*
  #temperature should be ground surface temp, not air 2m temp!
  # kinetic and equilibrium fracionation factors
  alpha_eq = ((-7.685+6712.3/(TsurfK)-1666400/((TsurfK)^2)+350410000/((TsurfK)^3)) + 1000)/1000
  alpha_kin = ((14.2*(1-rh)) + 1000)/1000
  
  R_precip = (delta18_precip/1000)*0.0020052 + 0.0020052
  R_atm = R_precip/alpha_eq
  
  R_sur = alpha_eq*((1-rh)*alpha_kin*R_precip + rh*R_atm)
  #z = seq(1,maxdepth,1) # assuming soil depth of 2 m max (z is depth in cm)
  z=seq(0,maxdepth)
  # zstar = (p*t*d)/E # this is z*
  profile <- (((R_sur - R_precip)*exp(-z/zstar) + R_precip)/0.0020052 - 1)*1000
  profile #vector of soil water d18O with depth
}
#soil.profile(-7,284,0.85,50,100)

####calculate soil profiles

#calculate soil water profile at every point along track for each E/T ratio
soilwater18=array(dim=c(tracklength,length(zplot),length(e)))
soilcarb18=array(dim=c(tracklength,length(zplot),length(e)))
#and also for eddies!
soilwater18_eddy=array(dim=c(tracklength,length(zplot),length(e)))
soilcarb18_eddy=array(dim=c(tracklength,length(zplot),length(e)))
for (i in 1:tracklength) {
  for (j in 1:length(zplot)) {
    for (h in 1:length(e)) {
      soilwater18[i,,h]=soil.profile(delta_p.meanfromcoast[i,h],TsurfK[i],rh[i],zstar[i],maxdepth)
      soilwater18_eddy[i,,h]=soil.profile(delta_p_eddy.meanfromcoast[i,h],TsurfK[i],rh[i],zstar[i],maxdepth)
      soilcarb18[i,j,h]=water2carb.oxygen(soilwater18[i,j,h],tsoil[i,j],pH)
      soilcarb18_eddy[i,j,h]=water2carb.oxygen(soilwater18_eddy[i,j,h],tsoil[i,j],pH)
    }
  }
}
#these have dim 28x219x11 #28 points, 219 depth values, 11 E/T ratios

##FOR COTTON 2012
#calculate soil water profile at every point along track for each E/T ratio
soilwater18.1=array(dim=c(tracklength,length(zplot),length(e)))
soilcarb18.1=array(dim=c(tracklength,length(zplot),length(e)))
#and also for eddies!
soilwater18_eddy.1=array(dim=c(tracklength,length(zplot),length(e)))
soilcarb18_eddy.1=array(dim=c(tracklength,length(zplot),length(e)))
for (i in 1:tracklength) {
  for (j in 1:length(zplot)) {
    for (h in 1:length(e)) {
      soilwater18.1[i,,h]=soil.profile(delta_p.meanfromcoast[i,h],TsurfK[i],rh[i],zstar[i],maxdepth)
      soilwater18_eddy.1[i,,h]=soil.profile(delta_p_eddy.meanfromcoast[i,h],TsurfK[i],rh[i],zstar[i],maxdepth)
      soilcarb18.1[i,j,h]=water2carb.oxygen(soilwater18.1[i,j,h],tsoil[i,j],pH,CO2ppm_soil.1[i,j])
      soilcarb18_eddy.1[i,j,h]=water2carb.oxygen(soilwater18_eddy.1[i,j,h],tsoil[i,j],pH,CO2ppm_soil.1[i,j])
    }
  }
}

##COTTON 2013
soilwater18.2=array(dim=c(tracklength,length(zplot),length(e)))
soilcarb18.2=array(dim=c(tracklength,length(zplot),length(e)))
#and also for eddies!
soilwater18_eddy.2=array(dim=c(tracklength,length(zplot),length(e)))
soilcarb18_eddy.2=array(dim=c(tracklength,length(zplot),length(e)))
for (i in 1:tracklength) {
  for (j in 1:length(zplot)) {
    for (h in 1:length(e)) {
      soilwater18.2[i,,h]=soil.profile(delta_p.meanfromcoast[i,h],TsurfK[i],rh[i],zstar[i],maxdepth)
      soilwater18_eddy.2[i,,h]=soil.profile(delta_p_eddy.meanfromcoast[i,h],TsurfK[i],rh[i],zstar[i],maxdepth)
      soilcarb18.2[i,j,h]=water2carb.oxygen(soilwater18.2[i,j,h],tsoil[i,j],pH,CO2ppm_soil.2[i,j])
      soilcarb18_eddy.2[i,j,h]=water2carb.oxygen(soilwater18_eddy.2[i,j,h],tsoil[i,j],pH,CO2ppm_soil.2[i,j])
    }
  }
}

###RAICH and SCHLESINGER 1992
soilwater18.3=array(dim=c(tracklength,length(zplot),length(e)))
soilcarb18.3=array(dim=c(tracklength,length(zplot),length(e)))
#and also for eddies!
soilwater18_eddy.3=array(dim=c(tracklength,length(zplot),length(e)))
soilcarb18_eddy.3=array(dim=c(tracklength,length(zplot),length(e)))
for (i in 1:tracklength) {
  for (j in 1:length(zplot)) {
    for (h in 1:length(e)) {
      soilwater18.3[i,,h]=soil.profile(delta_p.meanfromcoast[i,h],TsurfK[i],rh[i],zstar[i],maxdepth)
      soilwater18_eddy.3[i,,h]=soil.profile(delta_p_eddy.meanfromcoast[i,h],TsurfK[i],rh[i],zstar[i],maxdepth)
      soilcarb18.3[i,j,h]=water2carb.oxygen(soilwater18.3[i,j,h],tsoil[i,j],pH,CO2ppm_soil.3[i,j])
      soilcarb18_eddy.3[i,j,h]=water2carb.oxygen(soilwater18_eddy.3[i,j,h],tsoil[i,j],pH,CO2ppm_soil.3[i,j])
    }
  }
}

#######plot model soil water profiles########


#we have 11 E/ET ratios, 28 points along track, 3 MAP-resp relationships, and eddies vs advection.
#arbitrarily choosing distance of 576km from coast (approximately Hay data from below)

point=which(distance.fake==576) #arbitrary point on storm track (19)
#dim(soilcarb18) #28 points, 219 depths, 11 E/ETs

#start with advection only, constant soil CO2
plot(soilcarb18[point,,1],zplot,yaxs='i',ylim=c(max(zplot),0),lwd=3,type='l',lty=1,col='red',
     xlim=c(-25,0),xlab=expression(paste(delta^18,'O'[VPDB],' (\u2030)')),ylab='depth (cm)')
points(soilcarb18[point,,8],zplot,yaxs='i',lwd=3,type='l',lty=1,col='green')
points(soilcarb18[point,,11],zplot,yaxs='i',lwd=3,type='l',lty=1,col='blue')
#Cotton 2012, advection only
points(soilcarb18.1[point,,1],zplot,yaxs='i',lwd=3,type='l',lty=2,col='red')
points(soilcarb18.1[point,,8],zplot,yaxs='i',lwd=3,type='l',lty=2,col='green')
points(soilcarb18.1[point,,11],zplot,yaxs='i',lwd=3,type='l',lty=2,col='blue')
#Cotton 2013, advection only
points(soilcarb18.2[point,,1],zplot,yaxs='i',lwd=3,type='l',lty=3,col='red')
points(soilcarb18.2[point,,8],zplot,yaxs='i',lwd=3,type='l',lty=3,col='green')
points(soilcarb18.2[point,,11],zplot,yaxs='i',lwd=3,type='l',lty=3,col='blue')
#Raich and Schlesinger, advection only
points(soilcarb18.3[point,,1],zplot,yaxs='i',lwd=3,type='l',lty=4,col='red')
points(soilcarb18.3[point,,8],zplot,yaxs='i',lwd=3,type='l',lty=4,col='green')
points(soilcarb18.3[point,,11],zplot,yaxs='i',lwd=3,type='l',lty=4,col='blue')
#constant soil CO2, eddies only
points(soilcarb18_eddy[point,,1],zplot,yaxs='i',lwd=1,type='l',lty=1,col='red')
points(soilcarb18_eddy[point,,8],zplot,yaxs='i',lwd=1,type='l',lty=1,col='green')
points(soilcarb18_eddy[point,,11],zplot,yaxs='i',lwd=1,type='l',lty=1,col='blue')
#Cotton 2012, advection only
points(soilcarb18_eddy.1[point,,1],zplot,yaxs='i',lwd=1,type='l',lty=2,col='red')
points(soilcarb18_eddy.1[point,,8],zplot,yaxs='i',lwd=1,type='l',lty=2,col='green')
points(soilcarb18_eddy.1[point,,11],zplot,yaxs='i',lwd=1,type='l',lty=2,col='blue')
#Cotton 2013, advection only
points(soilcarb18_eddy.2[point,,1],zplot,yaxs='i',lwd=1,type='l',lty=3,col='red')
points(soilcarb18_eddy.2[point,,8],zplot,yaxs='i',lwd=1,type='l',lty=3,col='green')
points(soilcarb18_eddy.2[point,,11],zplot,yaxs='i',lwd=1,type='l',lty=3,col='blue')
#Raich and Schlesinger, advection only
points(soilcarb18_eddy.3[point,,1],zplot,yaxs='i',lwd=1,type='l',lty=4,col='red')
points(soilcarb18_eddy.3[point,,8],zplot,yaxs='i',lwd=1,type='l',lty=4,col='green')
points(soilcarb18_eddy.3[point,,11],zplot,yaxs='i',lwd=1,type='l',lty=4,col='blue')
#####THERE seems to be NO sensitivity to soil pCO2 in the d18O values



#######modern soil profiles########

soildata=read.csv('modern soils.csv')
head(soildata)

#add column for distance from "coast"
library(geosphere)
library(dplyr)
lati4=37.7
loni4=-122.5
lati5=36.27
loni5=-121.8
latitudes=soildata$latitude
longitudes=soildata$longitude
distmin <- dist1 <- dist2 <- dist3 <- dist4 <- dist5 <- numeric()
for (i in 1:length(latitudes)){
  dist1[i]=distVincentyEllipsoid(c(longitudes[i],latitudes[i]),c(loni,lati))
  dist2[i]=distVincentyEllipsoid(c(longitudes[i],latitudes[i]),c(loni2,lati2))
  dist3[i]=distVincentyEllipsoid(c(longitudes[i],latitudes[i]),c(loni3,lati3))
  dist4[i]=distVincentyEllipsoid(c(longitudes[i],latitudes[i]),c(loni4,lati4))
  dist5[i]=distVincentyEllipsoid(c(longitudes[i],latitudes[i]),c(loni5,lati5))
  distmin[i]=min(dist1[i],dist2[i],dist3[i],dist4[i],dist5[i])/1000 #distance in km
}

soildata2=cbind(soildata,distmin) #added column for distance from coast
head(soildata2)

#our comparable distances for the modern soil data (left) and model output (right) are:
#640.1107     #640  #Newark Valley
#791.6856     #800  #South of Wendover
#586.1828     #576  #Hay  
#367.1008     #352  #Pendall et al 1994 
#551.1161     #544  #Quade et al 1989
#447.7820     #448  #Quade et al 1989

#I'm going to plot these six for data-model comparison
point1=which(distance.fake==352)
point2=which(distance.fake==448)
point3=which(distance.fake==544)
point4=which(distance.fake==576)
point5=which(distance.fake==640)
point6=which(distance.fake==800)


###
#####d18O plots!
###


#Profile1: distance 367.1008/352
profile1.1=soilcarb18[point1,,1] #all evap
profile1.2=soilcarb18[point1,,8] #70% transpiration
profile1.3=soilcarb18[point1,,11] #all transpiration
profile1.1_eddy=soilcarb18_eddy[point1,,1] #all evap, eddies only
profile1.2_eddy=soilcarb18_eddy[point1,,8] #70% transpiration, eddies only
profile1.3_eddy=soilcarb18_eddy[point1,,11] #all transpiration, eddies only
profile1.d18O=filter(soildata2,latitude==37.6947)$d18O
profile1.depth=filter(soildata2,latitude==37.6947)$depth

#Profile2: distance 448
profile2.1=soilcarb18[point2,,1] #all evap
profile2.2=soilcarb18[point2,,8] #70% transpiration
profile2.3=soilcarb18[point2,,11] #all transpiration
profile2.1_eddy=soilcarb18_eddy[point2,,1] #all evap, eddies only
profile2.2_eddy=soilcarb18_eddy[point2,,8] #70% transpiration, eddies only
profile2.3_eddy=soilcarb18_eddy[point2,,11] #all transpiration, eddies only
profile2.d18O=filter(soildata2,latitude==36.2469)$d18O
profile2.depth=filter(soildata2,latitude==36.2469)$depth

#Profile3: distance 544
profile3.1=soilcarb18[point3,,1] #all evap
profile3.2=soilcarb18[point3,,8] #70% transpiration
profile3.3=soilcarb18[point3,,11] #all transpiration
profile3.1_eddy=soilcarb18_eddy[point3,,1] #all evap, eddies only
profile3.2_eddy=soilcarb18_eddy[point3,,8] #70% transpiration, eddies only
profile3.3_eddy=soilcarb18_eddy[point3,,11] #all transpiration, eddies only
profile3.d18O=filter(soildata2,latitude==36.2981)$d18O
profile3.depth=filter(soildata2,latitude==36.2981)$depth

#Profile4: distance 576
profile4.1=soilcarb18[point4,,1] #all evap
profile4.2=soilcarb18[point4,,8] #70% transpiration
profile4.3=soilcarb18[point4,,11] #all transpiration
profile4.1_eddy=soilcarb18_eddy[point4,,1] #all evap, eddies only
profile4.2_eddy=soilcarb18_eddy[point4,,8] #70% transpiration, eddies only
profile4.3_eddy=soilcarb18_eddy[point4,,11] #all transpiration, eddies only
profile4.d18O=filter(soildata2,latitude==39.572481)$d18O
profile4.depth=filter(soildata2,latitude==39.572481)$depth

#Profile5: distance 640
profile5.1=soilcarb18[point5,,1] #all evap
profile5.2=soilcarb18[point5,,8] #70% transpiration
profile5.3=soilcarb18[point5,,11] #all transpiration
profile5.1_eddy=soilcarb18_eddy[point5,,1] #all evap, eddies only
profile5.2_eddy=soilcarb18_eddy[point5,,8] #70% transpiration, eddies only
profile5.3_eddy=soilcarb18_eddy[point5,,11] #all transpiration, eddies only
profile5.d18O=filter(soildata2,latitude==39.876944)$d18O
profile5.depth=filter(soildata2,latitude==39.876944)$depth

#Profile6: distance 800
profile6.1=soilcarb18[point6,,1] #all evap
profile6.2=soilcarb18[point6,,8] #70% transpiration
profile6.3=soilcarb18[point6,,11] #all transpiration
profile6.1_eddy=soilcarb18_eddy[point6,,1] #all evap, eddies only
profile6.2_eddy=soilcarb18_eddy[point6,,8] #70% transpiration, eddies only
profile6.3_eddy=soilcarb18_eddy[point6,,11] #all transpiration, eddies only
profile6.d18O=filter(soildata2,latitude==40.594903)$d18O
profile6.depth=filter(soildata2,latitude==40.594903)$depth


### d18O profile plots!

pdf('Modern d18O Profiles with eddies and Data.pdf',width=10,height=8)
par(mfrow=c(2,3))
plot(profile1.1,zplot,ylim=c(maxdepth,0),xlim=c(-32,3),type='l',lwd=3,col='red',
     xlab=expression(paste(delta^18,'O'[VPDB],' (\u2030)')),ylab='depth (cm)',
     main='Pendall et al, 1994 (~350km)')
points(profile1.2,zplot,type='l',lwd=3,col='green')
points(profile1.3,zplot,type='l',lwd=3,col='blue')
points(profile1.1_eddy,zplot,type='l',col='red')
points(profile1.2_eddy,zplot,type='l',col='green')
points(profile1.3_eddy,zplot,type='l',col='blue')
points(profile1.d18O,profile1.depth,pch=16)
legend('bottomleft',pch=c(16,rep(NA,5)),lty=c(NA,rep(1,5)),lwd=c(rep(1,4),3),
       col=c('black','red','green','blue','black','black'),cex=0.8,
       legend=c('data','all E','70% T','all T','advection','eddies'))

plot(profile2.1,zplot,ylim=c(maxdepth,0),xlim=c(-32,3),type='l',lwd=3,col='red',
     xlab=expression(paste(delta^18,'O'[VPDB],' (\u2030)')),ylab='depth (cm)',
     main='Quade et al, 1989 (~450km)')
points(profile2.2,zplot,type='l',lwd=3,col='green')
points(profile2.3,zplot,type='l',lwd=3,col='blue')
points(profile2.1_eddy,zplot,type='l',col='red')
points(profile2.2_eddy,zplot,type='l',col='green')
points(profile2.3_eddy,zplot,type='l',col='blue')
points(profile2.d18O,profile2.depth,pch=16)
legend('bottomleft',pch=c(16,rep(NA,5)),lty=c(NA,rep(1,5)),lwd=c(rep(1,4),3),
       col=c('black','red','green','blue','black','black'),cex=0.8,
       legend=c('data','all E','70% T','all T','advection','eddies'))

plot(profile3.1,zplot,ylim=c(maxdepth,0),xlim=c(-32,3),type='l',lwd=3,col='red',
     xlab=expression(paste(delta^18,'O'[VPDB],' (\u2030)')),ylab='depth (cm)',
     main='Quade et al, 1989 (~550km)')
points(profile3.2,zplot,type='l',lwd=3,col='green')
points(profile3.3,zplot,type='l',lwd=3,col='blue')
points(profile3.1_eddy,zplot,type='l',col='red')
points(profile3.2_eddy,zplot,type='l',col='green')
points(profile3.3_eddy,zplot,type='l',col='blue')
points(profile3.d18O,profile3.depth,pch=16)
legend('bottomright',pch=c(16,rep(NA,5)),lty=c(NA,rep(1,5)),lwd=c(rep(1,4),3),
       col=c('black','red','green','blue','black','black'),cex=0.8,
       legend=c('data','all E','70% T','all T','advection','eddies'))

plot(profile4.1,zplot,ylim=c(maxdepth,0),xlim=c(-32,3),type='l',lwd=3,col='red',
     xlab=expression(paste(delta^18,'O'[VPDB],' (\u2030)')),ylab='depth (cm)',
     main='Hay (~580km)')
points(profile4.2,zplot,type='l',lwd=3,col='green')
points(profile4.3,zplot,type='l',lwd=3,col='blue')
points(profile4.1_eddy,zplot,type='l',col='red')
points(profile4.2_eddy,zplot,type='l',col='green')
points(profile4.3_eddy,zplot,type='l',col='blue')
points(profile4.d18O,profile4.depth,pch=16)
legend('bottomright',pch=c(16,rep(NA,5)),lty=c(NA,rep(1,5)),lwd=c(rep(1,4),3),
       col=c('black','red','green','blue','black','black'),cex=0.8,
       legend=c('data','all E','70% T','all T','advection','eddies'))

plot(profile5.1,zplot,ylim=c(maxdepth,0),xlim=c(-32,3),type='l',lwd=3,col='red',
     xlab=expression(paste(delta^18,'O'[VPDB],' (\u2030)')),ylab='depth (cm)',
     main='Newark Valley (~640km)')
points(profile5.2,zplot,type='l',lwd=3,col='green')
points(profile5.3,zplot,type='l',lwd=3,col='blue')
points(profile5.1_eddy,zplot,type='l',col='red')
points(profile5.2_eddy,zplot,type='l',col='green')
points(profile5.3_eddy,zplot,type='l',col='blue')
points(profile5.d18O,profile5.depth,pch=16)
legend('bottomright',pch=c(16,rep(NA,5)),lty=c(NA,rep(1,5)),lwd=c(rep(1,4),3),
       col=c('black','red','green','blue','black','black'),cex=0.8,
       legend=c('data','all E','70% T','all T','advection','eddies'))

plot(profile6.1,zplot,ylim=c(maxdepth,0),xlim=c(-32,3),type='l',lwd=3,col='red',
     xlab=expression(paste(delta^18,'O'[VPDB],' (\u2030)')),ylab='depth (cm)',
     main='South of Wendover (~800km)')
points(profile6.2,zplot,type='l',lwd=3,col='green')
points(profile6.3,zplot,type='l',lwd=3,col='blue')
points(profile6.1_eddy,zplot,type='l',col='red')
points(profile6.2_eddy,zplot,type='l',col='green')
points(profile6.3_eddy,zplot,type='l',col='blue')
points(profile6.d18O,profile6.depth,pch=16)
legend('bottomright',pch=c(16,rep(NA,5)),lty=c(NA,rep(1,5)),lwd=c(rep(1,4),3),
       col=c('black','red','green','blue','black','black'), cex=0.8,
       legend=c('data','all E','70% T','all T','advection','eddies'))
dev.off()


###
######d13C plots!
####


#extract d13C data:

#Profile1: distance 352
profile1.d13C=filter(soildata2,latitude==37.6947)$d13C
profile1.depth=filter(soildata2,latitude==37.6947)$depth
#Profile2: distance 448
profile2.d13C=filter(soildata2,latitude==36.2469)$d13C
profile2.depth=filter(soildata2,latitude==36.2469)$depth
#Profile3: distance 544
profile3.d13C=filter(soildata2,latitude==36.2981)$d13C
profile3.depth=filter(soildata2,latitude==36.2981)$depth
#Profile4: distance 576
profile4.d13C=filter(soildata2,latitude==39.572481)$d13C
profile4.depth=filter(soildata2,latitude==39.572481)$depth
#Profile5: distance 640
profile5.d13C=filter(soildata2,latitude==39.876944)$d13C
profile5.depth=filter(soildata2,latitude==39.876944)$depth
#Profile6: distance 800
profile6.d13C=filter(soildata2,latitude==40.594903)$d13C
profile6.depth=filter(soildata2,latitude==40.594903)$depth

#plot for multiple points and 3 MAP-resp relationships vs depth
pdf('Modern d13C Profiles.pdf',width=10,height=8)
par(mfrow=c(2,3))
plot(soilCO2.1[point1,,3],zplot,ylim=c(max(zplot),0),yaxs='i',type='l',col='red',
     xlab=expression(paste('Soil-Respired ',CO[2],' (ppm)')),ylab='depth (cm)',
     main='Pendall et al, 1994 (~350km)')
points(soilCO2.2[point1,,3],zplot,type='l',col='blue')
points(soilCO2.3[point1,,3],zplot,type='l',col='green')
points(profile1.d13C,profile1.depth,pch=16)
legend('bottomright',legend=c('Cotton12','Cotton13','S&R'),col=c('red','blue','green'),lty=1)

plot(soilCO2.1[point2,,3],zplot,ylim=c(max(zplot),0),yaxs='i',type='l',col='red',xlim=c(-13,2),
     xlab=expression(paste('Soil-Respired ',CO[2],' (ppm)')),ylab='depth (cm)',
     main='Quade et al, 1989 (~450km)')
points(soilCO2.2[point2,,3],zplot,type='l',col='blue')
points(soilCO2.3[point2,,3],zplot,type='l',col='green')
points(profile2.d13C,profile2.depth,pch=16)
legend('bottomright',legend=c('Cotton12','Cotton13','S&R'),col=c('red','blue','green'),lty=1)

plot(soilCO2.1[point3,,3],zplot,ylim=c(max(zplot),0),yaxs='i',type='l',col='red',xlim=c(-13,2),
     xlab=expression(paste('Soil-Respired ',CO[2],' (ppm)')),ylab='depth (cm)',
     main='Quade et al, 1989 (~550km)')
points(soilCO2.2[point3,,3],zplot,type='l',col='blue')
points(soilCO2.3[point3,,3],zplot,type='l',col='green')
points(profile3.d13C,profile3.depth,pch=16)
legend('bottomright',legend=c('Cotton12','Cotton13','S&R'),col=c('red','blue','green'),lty=1)

plot(soilCO2.1[point4,,3],zplot,ylim=c(max(zplot),0),yaxs='i',type='l',col='red',xlim=c(-13,2),
     xlab=expression(paste('Soil-Respired ',CO[2],' (ppm)')),ylab='depth (cm)',
     main='Hay (~580km)')
points(soilCO2.2[point4,,3],zplot,type='l',col='blue')
points(soilCO2.3[point4,,3],zplot,type='l',col='green')
points(profile4.d13C,profile4.depth,pch=16)
legend('bottomright',legend=c('Cotton12','Cotton13','S&R'),col=c('red','blue','green'),lty=1)

plot(soilCO2.1[point5,,3],zplot,ylim=c(max(zplot),0),yaxs='i',type='l',col='red',xlim=c(-13,2),
     xlab=expression(paste('Soil-Respired ',CO[2],' (ppm)')),ylab='depth (cm)',
     main='Newark Valley (~640km)')
points(soilCO2.2[point5,,3],zplot,type='l',col='blue')
points(soilCO2.3[point5,,3],zplot,type='l',col='green')
points(profile5.d13C,profile5.depth,pch=16)
legend('bottomright',legend=c('Cotton12','Cotton13','S&R'),col=c('red','blue','green'),lty=1)

plot(soilCO2.1[point6,,3],zplot,ylim=c(max(zplot),0),yaxs='i',type='l',col='red',xlim=c(-13,2),
     xlab=expression(paste('Soil-Respired ',CO[2],' (ppm)')),ylab='depth (cm)',
     main='South of Wendover (~800km)')
points(soilCO2.2[point6,,3],zplot,type='l',col='blue')
points(soilCO2.3[point6,,3],zplot,type='l',col='green')
points(profile6.d13C,profile6.depth,pch=16)
legend('bottomright',legend=c('Cotton12','Cotton13','S&R'),col=c('red','blue','green'),lty=1)
dev.off()

###OOOF THESE ARE STILL BAD

###########remaking map########
#this time including profile locations and initial lat/lon points
library(maps)
library(oce)
library(dplyr)

latis=c(lati,lati2,lati3,lati4,lati5)
lonis=c(loni,loni2,loni3,loni4,loni5)
proflats=summarise(group_by(soildata2,latitude),longitude=mean(longitude))$latitude
proflons=summarise(group_by(soildata2,latitude),longitude=mean(longitude))$longitude

lim=c(-22,0)
ncol=100
col=oceColorsJet(ncol)
# Draw colorbar for d18O
pdf('map of data with profile locations.pdf',width=10, height=7)
par(mfrow=c(1,2))
drawPalette(zlim=lim,col=col,
            zlab=expression(paste(delta^18,"O"[precip])),
            las=1,cex.lab=1.2) # adjust palette position with mai


par(mar=c(2.5,2.5,1,1),mgp=c(2,0.7,0)) # tighten margins
map('state',xlim=c(min(moderndata2$lon)-1,max(moderndata2$lon)),
    ylim=c(min(moderndata2$lat)-1,max(moderndata2$lat)))
for (i in 1:10) {  #add storm tracks
  points(lon.data[,i],lat.data[,i],col="blue",type="o",pch=1,cex=0.5)
}
points(lon_rain,lat_rain,pch=21,cex=0.7,
       bg=col[rescale(x=d18O_rain,xlow=lim[1],
                      xhigh=lim[2],rlow=1,rhigh=ncol)])
points(lon_snow,lat_snow,pch=22,cex=0.7,
       bg=col[rescale(x=d18O_snow,xlow=lim[1],
                      xhigh=lim[2],rlow=1,rhigh=ncol)])
points(lon_river,lat_river,pch=23,cex=0.7,
       bg=col[rescale(x=d18O_river,xlow=lim[1],
                      xhigh=lim[2],rlow=1,rhigh=ncol)])
points(lon_gw,lat_gw,pch=24,cex=0.7,
       bg=col[rescale(x=d18O_gw,xlow=lim[1],
                      xhigh=lim[2],rlow=1,rhigh=ncol)])
points(lon_lake,lat_lake,pch=25,cex=0.7,
       bg=col[rescale(x=d18O_lake,xlow=lim[1],
                      xhigh=lim[2],rlow=1,rhigh=ncol)])
points(lonis,latis,pch=16,col='black',cex=0.8)
points(proflons,proflats,pch=15,col='black',cex=0.8)
legend('bottomleft',pch=c(21,22,23,24,25,15,16),pt.bg='cyan',cex=0.7,bty='n',
       legend=c('rain','snow','rivers/creeks','groundwater','lakes','soil profiles','coast points'))
dev.off()


###
#####porosity sensitivity analysis#####
####

#testing sensitivity of d18O, d13C of soil carbonate to soil porosity

point=which(distance.fake==576) #arbitrary point on storm track (19)

###OXYGEN
#run the code multiple times for multiple values of porosity
#oxygen not sensitive to MAP-Resp, so we're omitting here
#save the different values here

##BE CAREFUL not to overwrite these once already saved
#p=0.4
O18_p40=soilcarb18
O18_p40_eddy=soilcarb18_eddy
#p=0.5
O18_p50=soilcarb18
O18_p50_eddy=soilcarb18_eddy
#p=0.6
O18_p60=soilcarb18
O18_p60_eddy=soilcarb18_eddy

#pdf('O18 porosity sensitivity.pdf',width=15,height=7)
par(mfrow=c(1,3))
#all evaporation
plot(O18_p40[point,,1],zplot,yaxs='i',ylim=c(max(zplot),0),lwd=2,type='l',lty=1,col='red',
     xlim=c(-25,0),xlab=expression(paste(delta^18,'O'[VPDB],' (\u2030)')),ylab='depth (cm)',
     main='All Evaporation')
points(O18_p50[point,,1],zplot,lwd=2,type='l',lty=1,col='green')
points(O18_p60[point,,1],zplot,lwd=2,type='l',lty=1,col='blue')
legend('bottomright',legend=c('p=0.4','p=0.5','p=0.6'),col=c('red','green','blue'),lty=1)
#70% transpiration
plot(O18_p40[point,,8],zplot,yaxs='i',ylim=c(max(zplot),0),lwd=2,type='l',lty=1,col='red',
     xlim=c(-25,0),xlab=expression(paste(delta^18,'O'[VPDB],' (\u2030)')),ylab='depth (cm)',
     main='70% Transpiration')
points(O18_p50[point,,8],zplot,lwd=2,type='l',lty=1,col='green')
points(O18_p60[point,,8],zplot,lwd=2,type='l',lty=1,col='blue')
#all transpiration
plot(O18_p40[point,,11],zplot,yaxs='i',ylim=c(max(zplot),0),lwd=2,type='l',lty=1,col='red',
     xlim=c(-25,0),xlab=expression(paste(delta^18,'O'[VPDB],' (\u2030)')),ylab='depth (cm)',
     main='All Transpiration')
points(O18_p50[point,,11],zplot,lwd=2,type='l',lty=1,col='green')
points(O18_p60[point,,11],zplot,lwd=2,type='l',lty=1,col='blue')
dev.off()


###CARBON
#like oxygen run the code multiple times for multiple values of porosity
#also save the different MAP-Resp relationships
#save the different values here

##BE CAREFUL not to overwrite these once already saved
#p=0.4
C13_p40.Cot12=soilCO2.1[,,3]
C13_p40.Cot13=soilCO2.2[,,3]
C13_p40.RS=soilCO2.3[,,3]
#p=0.5
C13_p50.Cot12=soilCO2.1[,,3]
C13_p50.Cot13=soilCO2.2[,,3]
C13_p50.RS=soilCO2.3[,,3]
#p=0.6
C13_p60.Cot12=soilCO2.1[,,3]
C13_p60.Cot13=soilCO2.2[,,3]
C13_p60.RS=soilCO2.3[,,3]

pdf('C13 porosity sensitivity.pdf',width=15,height=7)
par(mfrow=c(1,3))
#Cotton 2012
plot(C13_p40.Cot12[point,],zplot,yaxs='i',ylim=c(max(zplot),0),lwd=2,type='l',lty=1,col='red',
     xlim=c(-13,2),xlab=expression(paste(delta^13,'C'[VPDB],' (\u2030)')),ylab='depth (cm)',
     main='Cotton 2012')
points(C13_p50.Cot12[point,],zplot,type='l',lwd=2,lty=1,col='green')
points(C13_p60.Cot12[point,],zplot,type='l',lwd=2,lty=1,col='blue')
legend('bottomright',legend=c('Cotton12','Cotton13','R&S'),col=c('red','green','blue'),lty=1)
#Cotton 2013
plot(C13_p40.Cot13[point,],zplot,yaxs='i',ylim=c(max(zplot),0),lwd=2,type='l',lty=1,col='red',
     xlim=c(-13,2),xlab=expression(paste(delta^13,'C'[VPDB],' (\u2030)')),ylab='depth (cm)',
     main='Cotton 2013')
points(C13_p50.Cot13[point,],zplot,type='l',lwd=2,lty=1,col='green')
points(C13_p60.Cot13[point,],zplot,type='l',lwd=2,lty=1,col='blue')
#Raich and Schlesinger 1992
plot(C13_p40.RS[point,],zplot,yaxs='i',ylim=c(max(zplot),0),lwd=2,type='l',lty=1,col='red',
     xlim=c(-13,2),xlab=expression(paste(delta^13,'C'[VPDB],' (\u2030)')),ylab='depth (cm)',
     main='Raich and Scheslinger, 1992')
points(C13_p50.RS[point,],zplot,type='l',lwd=2,lty=1,col='green')
points(C13_p60.RS[point,],zplot,type='l',lwd=2,lty=1,col='blue')
dev.off()
