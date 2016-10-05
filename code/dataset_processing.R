#this prepares NARR data for use with HYSPLIT-generated trajectories
#   in the vapor transport model

#loops over all the desired years
#outputs .RData files for NARR data required for that year

#created 10-20-15 by AJR
#last modified 11-27-15 by AJR


library(ncdf)

#switch to directory with the NARR netcdf files
setwd("/Users/Annie/Desktop/NARR_for_processing")

###Import datasets!###

#land mask
ncdf = open.ncdf("land.nc") # 
land.data=get.var.ncdf(nc=ncdf,varid="land") # 

#latitude/longitude
lat.data=get.var.ncdf(nc=ncdf,varid="lat")
lon.data=get.var.ncdf(nc=ncdf,varid='lon')

save(land.data,lat.data,lon.data,file='NARRgriddata.RData')
rm(land.data,lat.data,lon.data) 
#clear these out so I can just save the image later!


years=seq(2005,2014,1) #ADJUST FOR YEARS DESIRED

#loop and save over each year!!!!

#for each year, save an RData file containing
 #model input climatological variables
for (year in years) {

  #Air Temperature at 2m (Kelvin)
  filename=paste('air.2m.',year,'.nc',sep='')
  dataname=paste('Temp_2m.',year,sep='')
  ncdf=open.ncdf(filename)
  assign(dataname,get.var.ncdf(nc=ncdf,varid='air'))  #Kelvin
  
  
  #Precipitable water (w) (kg/m^2)
  filename=paste('pr_wtr.',year,'.nc',sep='')
  dataname=paste('Prw.',year,sep='')
  ncdf=open.ncdf(filename)
  assign(dataname,get.var.ncdf(nc=ncdf,varid='pr_wtr'))  #kg/m2
  
  
  # Relative humidity (rhum) #% <-PERCENT!
  filename=paste('rhum.2m.',year,'.nc',sep='')
  dataname=paste('rhum.',year,sep='')
  ncdf=open.ncdf(filename)
  assign(dataname,get.var.ncdf(nc=ncdf,varid='rhum'))  #PERCENT!!!
  
  
  # 3-hourly accumulated evaporation total (evap) # kg/m^2
  filename=paste('evap.',year,'.nc',sep='')
  #dataname=paste('evap.',year,sep='')
  ncdf=open.ncdf(filename)
  #assign(dataname,get.var.ncdf(nc=ncdf,varid='evap')) #kg/m2/(3hrs) which is mm/(3hrs)
  evap=get.var.ncdf(nc=ncdf,varid='evap') #kg/m2/(3hrs) which is mm/(3hrs)
  #don't need to save evap.year file; only need for Damkohler
  
  
  # Precipitation rate (prate) #kg/m2/(3hrs)
  filename=paste('prate.',year,'.nc',sep='')
  #dataname=paste('Pratemm.',year,sep='')
  ncdf=open.ncdf(filename)
  #assign(dataname,(get.var.ncdf(nc=ncdf,varid='prate'))*(3600*3)) #kg/m2/(3hrs) which is mm/(3hrs)
  precip=(get.var.ncdf(nc=ncdf,varid='prate'))*(3600*3) #kg/m2/(3hrs) which is mm/(3hrs)
  #we don't actually have to save a year prate file; only needed for Damkohler
  
  
  # Time (hours since 1800-1-1 00:00:0.0)
  #filename=paste('prate.',year,'.nc',sep='') #commented out because continuation from prate above
  dataname=paste('time.',year,sep='')
  #ncdf=open.ncdf(filename) #commented out because continuation from prate above
  assign(dataname,get.var.ncdf(nc=ncdf,varid='time')) #hours since 1800-1-1 00:00:0.0
  
  
  #Damkohler number (evap/(p-evap)) #optimal format determined in cees RStudio
  dataname=paste('Nd.',year,sep='')
  assign(dataname,(evap/(precip-evap)))
  
  
  ###save the .RData file for each year!
  #create name of file we're saving to
  savename=paste('NARR',year,'.RData',sep='')
  #remove the things we don't want to save in the year file
  rm(dataname,evap,filename,ncdf,precip)
  
  #save everything that contains the year (variable.year)
  save(list=ls()[grepl(year,ls())],file=savename)
  
  #remove all the files from that year so that they're not double-saved!
  rm(list=ls()[grepl(year,ls())])
  
}
