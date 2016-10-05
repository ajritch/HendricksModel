#latlonNARR function
#function to convert HYSPLIT lat/lons to NARR grid

#created 10-22-15 by AJR
#last modified 12-5-15 by AJR

#this function will take in a df with lat/lon and spit out a df with NARR-friendly lat/lon
#where "NARR-friendly" means its regridded to the NARR grid

library(RANN) #for nearest neighbor search
library(ncdf) #to get lat/lons from file
library(dplyr) #for data structure manipulation

#load NARR coordinates 
#setwd("/Users/Annie/Documents/300/HYSPLIT/")
#load('NARRgriddata.RData')

#necessary for function!
NARRlat=as.numeric(lat.data)
NARRlon=as.numeric(lon.data)
NARRlatlon=cbind(NARRlat,NARRlon) #has negatives for lon  #96673x2


#this function needs to run AFTER timechanger!!!!
#^because timechanger adds the extra first row for initialization
  #and we'll want to NARR-ify the initialization coords, too!
latlonNARR=function(df) {
  raw.latlon=cbind(df$Lat,df$Long)
  
  #nearest neighbor search between two coordinate sets
  nearest=nn2(NARRlatlon,raw.latlon,k=1) #k=1 to get just one nearest neighbor
  #^the output length is the same as the second variable in the nn2 function
  nn.idx=nearest[[1]] #one-column vector of the indices (row) of the nearest neighbor in NARRlatlon
  #nn.dists=nearest[[2]] #same, but for distance from traj.latlon to nearest neighbor in NARRlatlon

  newcoords=NARRlatlon[nn.idx,] #extract NARR coordinates
  df$Lat=newcoords[,1] #replace coordinates in the dataframe
  df$Long=newcoords[,2]
  names(df)[4]='Lon'  #switching name to my preferable version
  
  #spit out new df
  return(df)
}