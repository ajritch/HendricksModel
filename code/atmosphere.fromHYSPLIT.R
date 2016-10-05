#This code runs the atmosphere component of the model along individual HYSPLIT trajectories
#the bulk of the model and data manipulation handled in loaded files.

#Created: 11-11-15 by AJR
#Last Modified: 12-4-15 by AJR


library(RANN) #for nearest neighbor search in latlonNARR function
library(ncdf) #to get lat/lons from file
library(dplyr) #for data structure manipulation
library(geosphere) #to calculate distance between coordinates (Haversine)

#load functions for syncing with NARR data
#NARR data must already by saved as .RData files

setwd("/Users/Annie/documents/model/") #CHANGE to wherever function files stored
load('NARRgriddata.RData') #these necessary for latlonNARR() function
#^this loads land.data, lat.data, lon.data
source("timechanger.R") #function to help merge time on the datasets
source("latlonNARR.R") #function to regrid HYSPLIT output onto NARR grid
source("NARRgrabber.R") #function to load and grab NARR data along trajectories
#^^^IMPORTANT! before RUNNING this function user MUST switch to directory with NARR .RData files!!!
source('model.sourcefunctions.R') #functions used within the model
source('Hendricksmodel.onetraj.R') #function that runs the model itself


#beginning only with atmosphere portion
#first, establish a dataframe for testing!
setwd('/Users/Annie/Documents/model/HYSPLIT output/RData files')
load('Pyramid.SON.precip.RData') #list of 367 precipitation trajectories
rawlist=Pyramid.SON #CHANGE based on name of list loaded above

####NOTE uncomment below if list files are not saved elsewhere
  ##these steps are very computationally expensive, so it's best
  ##to use the pre-generated files for each location.
# #prepare list for use in model
# list.timechanged=lapply(rawlist,timechanger) #change time format to sync datasets
# list.fixedlatlon=lapply(list.timechanged,latlonNARR) #convert HYSPLIT output to NARR coordinates
# rm(list.timechanged) #clear space
# setwd('/Users/Annie/Desktop/NARR_for_processing') #CHANGE to directory with NARR .RData files
# #list.withNARR=lapply(list.fixedlatlon,NARRgrabber) #grab NARR data associated with each trajectory point
# #^files aren't closing properly due to use of lapply() 
# 
# #NOTE! best if ~10GB free on hard drive
# list.withNARR=list()
# for (i in 1:length(list.fixedlatlon)) { #CHANGE TO PROPER LENGTH
#   list.withNARR[i][[1]]=NARRgrabber(list.fixedlatlon[i][[1]])
# }
# setwd('/Users/Annie/Desktop/HYSPLIT output/RData files')
# save(list.withNARR,file='list.Pyramid.SON.withNARR.RData') #CHANGE to appropriate name
# 
# #clear out remaining files to save space
# rm(Nd.2013,Nd.2014)
# rm(Prw.2013,Prw.2014)
# rm(rhum.2013,rhum.2014)
# rm(Temp_2m.2013,Temp_2m.2014)
# rm(time.2013,time.2014)
# rm(list.fixedlatlon) #clear space

#now load saved lists to run the model
rm(list.withNARR)
#setwd('/Users/Annie/Desktop/HYSPLIT output/RData files')
load('list.Pyramid.SON.withNARR.RData') #CHANGE based on location/season


list.noNA=lapply(list.withNARR,filter,is.na(Prw)==F) #get rid of all NA'd points
#^^(NARR only has data for every third hour but HYSPLIT outputs more)

#reorder so temporally (switch backhours to temporal order)
reorder=function(dataframe) {
  reversed=dataframe[rev(rownames(dataframe)),]
  return(reversed)
}
list.ready=lapply(list.noNA,reorder) #model-ready list!!!!


###TESTING
# #picking out a trajectory with only 2004/2005
# test=jj.coord[[100]] #NOTE these don't include initialization point yet!
# test.timechanged=timechanger(test)
# test.fixedlatlon=latlonNARR(test.timechanged)
# setwd('/Users/Annie/Documents/300/HYSPLIT/NARR 2004') #directory with NARR .RData files
# test.withNARR=NARRgrabber(test.fixedlatlon)
# 
# #this leaves NA's for all but every third point! delete all of these NA rows
# test.noNA=filter(test.withNARR,is.na(Prw)==F)
# 
# #now, reverse the row order because the last point of trajectory is currently at top
# reversed=test.noNA[rev(rownames(test.noNA)),]
# 
# #Now test with some multiple trajectories!
# test=jj.coord[c(100,102)]
# test2=lapply(test,timechanger)
# test3=lapply(test2,latlonNARR)
# setwd('/Users/Annie/Documents/300/HYSPLIT/NARR 2004')
# test4=lapply(test3,NARRgrabber)
# test.noNA=lapply(test4,filter,is.na(Prw)==F)


#########model time!###########

#inputs to the model

initial_delta_p=-6.96 #these numbers essentially from Ingraham and Taylor transect
initial_delta_a=-17.476 
initial_delta_p_eddy=initial_delta_p #this is what Matt used, at least
initial_delta_a_eddy=initial_delta_a

#df=reversed
#tracklength=nrow(df) #length of storm track
#t=rep(0.7,times=tracklength) #T as fraction of total ET

#write a function to bind the T/ET to the prepared df so that lapply() can be used efficiently
#not fully sure yet that this will be useful
tbinder=function(df,t) {
  #t is a T/ET fraction value
  vector=rep(t,times=nrow(df))
  outputdf=cbind(df,t=vector)
  return(outputdf)
}

#input=tbinder(df,0.7)
#results=Hendricksmodel.onetraj(input) #IT WORKS!

#test model on list
input.list=lapply(list.ready,tbinder,0.7)
results=lapply(input.list,Hendricksmodel.onetraj) #IT WORKS but tons of NAs, NaNs....

Pyramid.SON.modelresults=results
save(Pyramid.SON.modelresults,file='Pyramid.SON.modelresults.RData')

plot(results[[1]]$Lon,results[[1]]$d18O.p,type='l')
points(results[[1]]$Lon,results[[1]]$d18O.p_eddy,type='l',col='red')
for (i in 2:length(results)) {
  points(results[[i]]$Lon,results[[i]]$d18O.p,type='l')
  points(results[[i]]$Lon,results[[i]]$d18O.p_eddy,type='l',col='red')
}
#meh this is kinda nasty