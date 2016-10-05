#function for the Hendricks model
#calcuates d18O, dD, d17O of vapor and precip along a storm track
#calculates advection-only and eddy diffusion-only scenarios

#this code is adpated from "model w evap balance and dxs.R" by Matt Winnick
#created: 11/12/15 by AJR
#last modified: 12/6/15 by AJR

#NOTE!!!! initial delta values must be defined before running this function!
#other model input parameters listed below because there's less likely to need changing

library(geosphere)

###model parameters
Rpdb=0.0112372 #dimensionless PDB Belemnite Ratio
smow18=0.0020052 # smow ratio of 18O/16O
smowD=155.76*10^(-6) # smow ratio of D/H
smow17=379.9*10^(-6) # vsmow ratio of 17O/16O
theta.eq=0.529 #Barkan and Luz, 2005; 17O equilibrium (l-v and v-l)
theta.diff=0.5185 #Barkan and Luz, 2007; 17O diffusion
k18=1-0.01872 # kinetic fractionation factor for 18O, wet soils from Mathieu and Bariac (1996), slightly higher than lake evap (~14.3)
kD=1-0.01676 # kinetic fractionation factor for dD, wet soils from Mathieu and Bariac (1996), slightly higher than lake evap (~12.5)
k17=k18^theta.diff #Barkan and Luz, 2005, 2007


Hendricksmodel.onetraj=function(revdf) {
  #takes in a dataframe that's optimized for model and has E/T specified per point
  #NA rows removed to match NARR's 3-hourly resolution
  #rows already in temporal order, with HYSPLIT initialization point last (at bottom)
  
  t=revdf$t
  e=1-t #evaporation as a fraction of total ET
  
  tracklength=nrow(revdf) #length of storm track
  
  w=revdf$Prw #precipitatble water, kg/m^2
  ntracks=1 #number of storm tracks
  Nd=revdf$Nd #calculated as E/(P-E) from NARR. LOTS OF NEGATIVES (like, 90%) :/
  ts=revdf$Temp_2m #2m surface temperature (K)
  rh=(revdf$rhum)/100 #relative humidity (fraction)
  
  #constants for isotope value calculations:
  #note that a18 != alpha18 (and same for D)
  a18=1 + (-7.685 + (6.7123*10^3)/(ts) - (1.6664*10^6)/(ts)^2 + (0.35041*10^9)/(ts)^3)/1000 #equilibrium fractionation 18O, v-l
  aD=1 + ((1158.8*ts^3)/10^9 - (1620.1*ts^2)/10^6 + (794.84*ts)/10^3 - 161.04 + (2.9992*10^9)/ts^3)/1000 #equilibrium fractionation, D v-l
  alpha18=1 - (-7.685 + (6.7123*10^3)/(ts) - (1.6664*10^6)/(ts)^2 + (0.35041*10^9)/(ts)^3)/1000 #equilibrium fractionation 18O, l-v
  alphaD=1 - ((1158.8*ts^3)/10^9 - (1620.1*ts^2)/10^6 + (794.84*ts)/10^3 - 161.04 + (2.9992*10^9)/ts^3)/1000 #equilibrium fractionation D, l-v
  alpha17=alpha18^(theta.eq) #Barkan and Luz, 2005
  a17=a18^theta.eq
  
  #build empty vectors to store data for d18O and dD and d17O for advection-only/eddy-only cases
  delta_a <- delta_et <- delta_inf<- delta_p <- delta_a_eddy <- delta_et_eddy <- delta_inf_eddy <- delta_p_eddy <- numeric(length=tracklength)
  delta_D_a <- delta_D_et <- delta_D_inf <- delta_D_p <- delta_D_a_eddy <- delta_D_et_eddy <- delta_D_inf_eddy <- delta_D_p_eddy <- numeric(length=tracklength)
  delta_17_a <- delta_17_et <- delta_17_inf <- delta_17_p <- delta_17_a_eddy <- delta_17_et_eddy <- delta_17_inf_eddy <- delta_17_p_eddy <- numeric(length=tracklength)
  
  delta_a[1] = initial_delta_a #sets initial d18O atmospheric vapor
  delta_p[1] = initial_delta_p #sets intial d18O precip
  delta_D_a[1] = delta_a[1]*8+10 #sets initial atmospheric vapor dD
  delta_D_p[1] = delta_p[1]*8+10 #sets initial precip dD
  delta_17_a[1] = (exp(0.528*log(delta_a[1]/1000+1)+0.000033)-1)*1000  #GMWL, Luz and Barkan, 2010
  delta_17_p[1] = (exp(0.528*log(delta_p[1]/1000+1)+0.000033)-1)*1000  #GMWL, Luz and Barkan, 2010
  #same for eddies-only case:
  delta_a_eddy[1] = initial_delta_a_eddy #sets initial d18O atmospheric vapor
  delta_p_eddy[1] = initial_delta_p_eddy #sets intial d18O precip
  delta_D_a_eddy[1] = delta_a_eddy[1]*8+10 #sets initial atmospheric vapor dD
  delta_D_p_eddy[1] = delta_p_eddy[1]*8+10 #sets initial precip dD
  delta_17_a_eddy[1] = (exp(0.528*log(delta_a_eddy[1]/1000+1)+0.000033)-1)*1000  #GMWL, Luz and Barkan, 2010
  delta_17_p_eddy[1] = (exp(0.528*log(delta_p_eddy[1]/1000+1)+0.000033)-1)*1000  #GMWL, Luz and Barkan, 2010
  
  
  ###run the model!
  
  for (i in 2:tracklength) {
    #calculate residual soil moisture for isotope_e function
    residual=1-(Nd[i]/(Nd[i]+1))*e[i]/(1-(Nd[i]/(Nd[i]+1)*t[i]))
    #if everything evaporates, set isotope value of evaporation to isotope value of precip
    if(residual == 0){
      delta_et[i] <- delta_p[i-1]
      delta_D_et[i] <- delta_D_p[i-1]
      delta_17_et[i] = delta_17_p[i-1]
      delta_et_eddy[i]=delta_p_eddy[i-1]
      delta_D_et_eddy[i] <- delta_D_p_eddy[i-1]
      delta_17_et_eddy[i] = delta_17_p_eddy[i-1]
      
      #otherwise, first calculate the combined equilibrium/kinetic fractionation using Craig-Gordon function
    } else {
      eps_e <- craig_gordon(t[i], rh[i], delta_p[i-1], ts[i])
      eps_e_D <- craig_gordon(t[i], rh[i], delta_D_p[i-1], ts[i], d18 = F)
      eps_e_17 = craig_gordon17(t[i],rh[i],delta_17_p[i-1],ts[i])
      eps_e_eddy <- craig_gordon(t[i], rh[i], delta_p_eddy[i-1], ts[i])
      eps_e_D_eddy <- craig_gordon(t[i], rh[i], delta_D_p_eddy[i-1], ts[i], d18 = F)
      eps_e_17_eddy = craig_gordon17(t[i],rh[i],delta_17_p_eddy[i-1],ts[i])
      #then calculate integrated ET value from location
      delta_et[i] <- t[i]*delta_p[i-1] + e[i]*isotope_e(residual,delta_p[i-1], eps_e)
      delta_D_et[i] <- t[i]*delta_D_p[i-1] + e[i]*isotope_e(residual,delta_D_p[i-1], eps_e_D)
      delta_17_et[i] = t[i]*delta_17_p[i-1]+e[i]*isotope_e(residual,delta_17_p[i-1],eps_e_17) 
      delta_et_eddy[i] <- t[i]*delta_p_eddy[i-1] + e[i]*isotope_e(residual,delta_p_eddy[i-1], eps_e_eddy)
      delta_D_et_eddy[i] <- t[i]*delta_D_p_eddy[i-1] + e[i]*isotope_e(residual,delta_D_p_eddy[i-1], eps_e_D_eddy)
      delta_17_et_eddy[i] = t[i]*delta_17_p_eddy[i-1]+e[i]*isotope_e(residual,delta_17_p_eddy[i-1],eps_e_17_eddy)
    }
    
    
    ###there's an issue with some of the values after the initial point being NAs
    #this loop will make it so all of the NAs stay NAs *except* the last NA before
    #hitting land. the first ones stay NAs so that we're only calculating over land.
    #note: still necessary to keep the original initials above in the event of no NAs
    #only where rowMeans==0 will they all be above water==>this track will be above water
    if (revdf$Land[i]==0) {
      delta_p[i]=initial_delta_p
      delta_p[i-1]=NA
      delta_a[i]=initial_delta_a
      delta_a[i-1]=NA
      delta_D_p[i]=delta_p[i]*8+10 #assume initial points are on MWL
      delta_D_p[i-1]=NA
      delta_D_a[i]=delta_a[i]*8+10
      delta_D_a[i-1]=NA
      delta_17_p[i]=(exp(0.528*log(delta_p[i]/1000+1)+0.000033)-1)*1000  #GMWL, Luz and Barkan, 2010
      delta_17_p[i-1]=NA
      delta_17_a[i]=(exp(0.528*log(delta_a[i]/1000+1)+0.000033)-1)*1000  #GMWL, Luz and Barkan, 2010
      delta_17_a[i-1]=NA
      #repeat for eddies (copypasta)
      delta_p_eddy[i]=initial_delta_p_eddy
      delta_p_eddy[i-1]=NA
      delta_a_eddy[i]=initial_delta_a_eddy
      delta_a_eddy[i-1]=NA
      delta_D_p_eddy[i]=delta_p_eddy[i]*8+10 #assume initial points are on MWL
      delta_D_p_eddy[i-1]=NA
      delta_D_a_eddy[i]=delta_a_eddy[i]*8+10
      delta_D_a_eddy[i-1]=NA
      delta_17_p_eddy[i]=(exp(0.528*log(delta_p_eddy[i]/1000+1)+0.000033)-1)*1000  #GMWL, Luz and Barkan, 2010
      delta_17_p_eddy[i-1]=NA
      delta_17_a_eddy[i]=(exp(0.528*log(delta_a_eddy[i]/1000+1)+0.000033)-1)*1000  #GMWL, Luz and Barkan, 2010
      delta_17_a_eddy[i-1]=NA
    } else {
      #haversine formula to calculate distance between coordinates:
      x=distVincentyEllipsoid(c(revdf$Lon[i],revdf$Lat[i]),c(revdf$Lon[i-1],revdf$Lat[i-1]))
      dx=x #meters, but the units don't matter. slightly >32km
      dw=w[i]-w[i-1]
      dw=ifelse(dw>0,0,dw) #gets rid of positive gradients (b/c we don't like explosions here)
      
      #Hendricks model (advection only)
      delta_inf[i] <- (Nd[i]*delta_et[i]-(1+Nd[i])*(a18[i]-1)*1000)/(a18[i]+Nd[i]*a18[i]-1)
      delta_a[i]=(delta_a[i-1]-delta_inf[i])*exp((a18[i]+a18[i]*Nd[i]-1)*(x/w[i])*(dw/dx))+delta_inf[i]
      delta_p[i] <-delta_a[i]+(a18[i]-1)*1000
      
      delta_D_inf[i] <- (Nd[i]*delta_D_et[i]-(1+Nd[i])*(aD[i]-1)*1000)/(aD[i]+Nd[i]*aD[i]-1)
      delta_D_a[i]=(delta_D_a[i-1]-delta_D_inf[i])*exp((aD[i]+aD[i]*Nd[i]-1)*(x/w[i])*(dw/dx))+delta_D_inf[i]
      delta_D_p[i] <-delta_D_a[i]+(aD[i]-1)*1000
      
      delta_17_inf[i] <- (Nd[i]*delta_17_et[i]-(1+Nd[i])*(a17[i]-1)*1000)/(a17[i]+Nd[i]*a17[i]-1)
      delta_17_a[i]=(delta_17_a[i-1]-delta_17_inf[i])*exp((a17[i]+a17[i]*Nd[i]-1)*(x/w[i])*(dw/dx))+delta_17_inf[i]
      delta_17_p[i] <-delta_17_a[i]+(a17[i]-1)*1000
      
      #Hendricks model (eddy diffusion only)
      delta_inf_eddy[i] <- (Nd[i]*delta_et_eddy[i]-(1+Nd[i])*(a18[i]-1)*1000)/(a18[i]+Nd[i]*a18[i]-1)
      delta_a_eddy[i]=(delta_a_eddy[i-1]-delta_inf_eddy[i])*exp((sqrt(a18[i]+a18[i]*Nd[i])-1)*(x/w[i])*(dw/dx))+delta_inf_eddy[i]
      delta_p_eddy[i] <-delta_a_eddy[i]+(a18[i]-1)*1000
      
      delta_D_inf_eddy[i] <- (Nd[i]*delta_D_et_eddy[i]-(1+Nd[i])*(aD[i]-1)*1000)/(aD[i]+Nd[i]*aD[i]-1)
      delta_D_a_eddy[i]=(delta_D_a_eddy[i-1]-delta_D_inf_eddy[i])*exp((sqrt(aD[i]+aD[i]*Nd[i])-1)*(x/w[i])*(dw/dx))+delta_D_inf_eddy[i]
      delta_D_p_eddy[i] <-delta_D_a_eddy[i]+(aD[i]-1)*1000
      
      delta_17_inf_eddy[i] <- (Nd[i]*delta_17_et_eddy[i]-(1+Nd[i])*(a17[i]-1)*1000)/(a17[i]+Nd[i]*a17[i]-1)
      delta_17_a_eddy[i]=(delta_17_a_eddy[i-1]-delta_17_inf_eddy[i])*exp((sqrt(a17[i]+a17[i]*Nd[i])-1)*(x/w[i])*(dw/dx))+delta_17_inf_eddy[i]
      delta_17_p_eddy[i] <-delta_17_a_eddy[i]+(a17[i]-1)*1000
    }  
  }
  
 
  #for now only returning d18O of precip (for advection only and eddies only)
  outputdf=cbind(revdf,d18O.p=delta_p,d18O.p_eddy=delta_p_eddy)
  #uncomment to add all the isotopic values!
  #   outputdf=cbind(revdf,d18O.atm=delta_a,d18O.p=delta_p,dD.atm=delta_D_a,dD.p=delta_D_p,
  #                  d17O.atm=delta_17_a,d17O.p=delta_17_p,d18O.atm_eddy=delta_a_eddy,d18O.p_eddy=delta_p_eddy)
  
  return(outputdf)
}

