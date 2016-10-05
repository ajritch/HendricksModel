# Script to filter trajectories by precipitation and plot them on a map.
# Created on:  5/20/2015
# Last Modified:  3/15/16 by AJR
# BASED ON: "Precip plotter.R" by Jeremy Caves
# modified from: "Trajectory Plotter-NARR.R" (last mod. 11-8-15) by Jeremy Caves 

# This script takes HYSPLIT output, loads it into a list, which can then be 
# analyzed for precipitation.  This code then plots this data as weighted by
# precipitation amount over the past hours. Because this script does not alter 
# the files, it works in the primary HYSPLIT output directory.  The code also
# allows one to plot all of the non-precipitating trajectories.  Currently,
# it is configured for plotting in the Western US.

# Notes:  Currently, this code cannot analyze the 31st day 
#         of each month.

# Currently set-up for: Coso Basin
# Months: Jan-Dec
# Years: 1980-2014
# Reanalysis Dataset: NARR

library(maps)
library(gplots)
library(abind)
library(fields)
library(maptools)
library(ncdf)

#name the directories where HYSPLIT output files are located
wd1000='/Users/Annie/Desktop/HYSPLIT output/1000Verdi_1980.2014'
wd1500='/Users/Annie/Desktop/HYSPLIT output/1500Verdi_1980.2014'
wd2000='/Users/Annie/Desktop/HYSPLIT output/2000Verdi_1980.2014'

#name the filenames, excluding altitude prefix and date postscript
filename='verdi' 

#number of altitudes being combined
naltitudes=3



#### Set up Time Vector ####
first.time=1980010100 # Starts January 1, 1980
years=35
year.interval=1000000 #Do not change
months=12
month.interval=10000 #Do not change
hour.interval=6
day.interval=100 # Do not change
year=numeric(length=years)
month=numeric(length=months)
hour=numeric(length=24/hour.interval)
days=array(dim=c(31,length(hour),length(month),years))
dim(days)

for(l in 1:years){
  year[l]=first.time+year.interval*(l-1)
  for (k in 1:length(month)){
    month[k]=year[l]+month.interval*(k-1)
    for (j in 1:length(hour)){
      hour[j]=month[k]+hour.interval*(j-1)
      for (i in 1:dim(days)[1]){
        days[i,j,k,l]=hour[j]+day.interval*(i-1)
      }
    }
  }
}

# days
days=days[-31,,,] #Removes 31st of each month
days=as.numeric(days)
length(days)
days=sort(days) #Sorts the days

# Eliminate 29 and 30 of all Februarys
days=days[-c((233:240)                    # Year 1
               ,((233+1440):(240+1440))     # 2 
               ,((233+1440*2):(240+1440*2)) # 3
               ,((233+1440*3):(240+1440*3)) # 4
               ,((233+1440*4):(240+1440*4)) # 5
               ,((233+1440*5):(240+1440*5)) # 6
               ,((233+1440*6):(240+1440*6)) # 7
               ,((233+1440*7):(240+1440*7)) # 8
               ,((233+1440*8):(240+1440*8)) # 9
               ,((233+1440*9):(240+1440*9)) # 10
               ,((233+1440*10):(240+1440*10)) # 11
               ,((233+1440*11):(240+1440*11)) # 12
               ,((233+1440*12):(240+1440*12)) # 13
               ,((233+1440*13):(240+1440*13)) # 14
               ,((233+1440*14):(240+1440*14)) # 15
               ,((233+1440*15):(240+1440*15)) # 16
               ,((233+1440*16):(240+1440*16)) # 17
               ,((233+1440*17):(240+1440*17)) # 18
               ,((233+1440*18):(240+1440*18)) # 19
               ,((233+1440*19):(240+1440*19)) # 20
               ,((233+1440*20):(240+1440*20)) # 21
               ,((233+1440*21):(240+1440*21)) # 22
               ,((233+1440*22):(240+1440*22)) # 23
               ,((233+1440*23):(240+1440*23)) # 24
               ,((233+1440*24):(240+1440*24)) # 25
               ,((233+1440*25):(240+1440*25)) # 26
               ,((233+1440*26):(240+1440*26)) # 27
               ,((233+1440*27):(240+1440*27)) # 28
               ,((233+1440*28):(240+1440*28)) # 29
               ,((233+1440*29):(240+1440*29)) # 30
               ,((233+1440*30):(240+1440*30)) # 31
               ,((233+1440*31):(240+1440*31)) # 32
               ,((233+1440*32):(240+1440*32)) # 33
               ,((233+1440*33):(240+1440*33)) # 34
               ,((233+1440*34):(240+1440*34)) # 35
) # Need to add more rows if more years are added
]
length(days)


#### Load HYSPLIT trajectories into arrays ####
back.plot=24*7 #Number of back-hours to plot for the trajectory
met.variables=20 #Number of meteorological variables computed by HYSPLIT
jj1000=list() # Create list to store all trajectories at 1000m
jj1500=list() #list to store 1500m trajectories
jj2000=list() #list to store 2000m trajectories
runnumber=length(days) # days.run*length(hour)*length(month)*length(year)
for (i in 1:runnumber){
#   jj[[i]]=as.matrix(read.table(file=paste("1000wendover",days[i],sep=""),sep="",skip=6,nrows=back.plot,colClasses="numeric")) #"1000westwendover for GDAS files"
#   jj[[i]]=as.matrix(read.table(file=paste("1000hay",days[i],sep=""),sep="",skip=6,nrows=back.plot,colClasses="numeric")) # Hay
#   jj[[i]]=as.matrix(read.table(file=paste("1000swendover",days[i],sep=""),sep="",skip=6,nrows=back.plot,colClasses="numeric")) # South of Wendover
#   jj[[i]]=as.matrix(read.table(file=paste("1000dubois",days[i],sep=""),sep="",skip=6,nrows=back.plot,colClasses="numeric")) # Dubois
#   jj[[i]]=as.matrix(read.table(file=paste("1000oupico",days[i],sep=""),sep="",skip=6,nrows=back.plot,colClasses="numeric")) # Oupico
#   jj[[i]]=as.matrix(read.table(file=paste("1000newark",days[i],sep=""),sep="",skip=6,nrows=back.plot)) # Newark Valley
#   jj[[i]]=as.matrix(read.table(file=paste("1000hay",days[i],sep=""),sep="",skip=6,nrows=back.plot)) # Hay
#   jj[[i]]=as.matrix(read.table(file=paste("1000swendover",days[i],sep=""),sep="",skip=6,nrows=back.plot)) # South of Wendover
#   jj[[i]]=as.matrix(read.table(file=paste("1000oupico",days[i],sep=""),sep="",skip=6,nrows=back.plot)) # Oupico
#   jj[[i]]=as.matrix(read.table(file=paste("1000pyramid",days[i],sep=""),sep="",skip=6,nrows=back.plot)) # Pyramid Lake
  
  setwd(wd1000)
  jj1000[[i]]=as.matrix(read.table(file=paste('1000',filename,days[i],sep=""),sep="",skip=6,nrows=back.plot)) 
  colnames(jj1000[[i]])=c("U","U","Year","Month","Day","Hour","U","U",
                       "Back.Hour","Lat","Long","AGL","Pressure","Temp",
                       "Rainfall","MixDepth","RH","SH","H2OMix","Terrain")
  setwd(wd1500)
  jj1500[[i]]=as.matrix(read.table(file=paste('1500',filename,days[i],sep=""),sep="",skip=6,nrows=back.plot)) 
  colnames(jj1500[[i]])=c("U","U","Year","Month","Day","Hour","U","U",
                          "Back.Hour","Lat","Long","AGL","Pressure","Temp",
                          "Rainfall","MixDepth","RH","SH","H2OMix","Terrain")
  setwd(wd2000)
  jj2000[[i]]=as.matrix(read.table(file=paste('2000',filename,days[i],sep=""),sep="",skip=6,nrows=back.plot)) 
  colnames(jj2000[[i]])=c("U","U","Year","Month","Day","Hour","U","U",
                         "Back.Hour","Lat","Long","AGL","Pressure","Temp",
                         "Rainfall","MixDepth","RH","SH","H2OMix","Terrain")
}
# for (i in 1:runnumber){
#   colnames(jj[[i]])=c("U","U","Year","Month","Day","Hour","U","U",
#                       "Back.Hour","Lat","Long","AGL","Pressure","Temp",
#                       "Rainfall","MixDepth","RH","SH","H2OMix","Terrain")
# }

##sometimes the file doesn't read if there's nothing there
 #skip those days and just keep adding on to jj list:
###ONLY DO THIS IF FIRST RUN ABOVE DIDN'T COMPLETE
newstart=i  #CHANGE FOR NEW START NUMBER
for (i in (newstart+1):runnumber) {
  setwd(wd1000)
  jj1000[[i]]=as.matrix(read.table(file=paste('1000',filename,days[i],sep=""),sep="",skip=6,nrows=back.plot)) 
  colnames(jj1000[[i]])=c("U","U","Year","Month","Day","Hour","U","U",
                         "Back.Hour","Lat","Long","AGL","Pressure","Temp",
                         "Rainfall","MixDepth","RH","SH","H2OMix","Terrain")
  setwd(wd1500)
  jj1500[[i]]=as.matrix(read.table(file=paste('1500',filename,days[i],sep=""),sep="",skip=6,nrows=back.plot)) 
  colnames(jj1500[[i]])=c("U","U","Year","Month","Day","Hour","U","U",
                          "Back.Hour","Lat","Long","AGL","Pressure","Temp",
                          "Rainfall","MixDepth","RH","SH","H2OMix","Terrain")
  setwd(wd2000)
  jj2000[[i]]=as.matrix(read.table(file=paste('2000',filename,days[i],sep=""),sep="",skip=6,nrows=back.plot)) 
  colnames(jj2000[[i]])=c("U","U","Year","Month","Day","Hour","U","U",
                         "Back.Hour","Lat","Long","AGL","Pressure","Temp",
                         "Rainfall","MixDepth","RH","SH","H2OMix","Terrain")
}

maxhours1000=numeric(length=runnumber)
maxhours1500=numeric(length=runnumber)
maxhours2000=numeric(length=runnumber)
for (i in 1:runnumber){
  #maxhours[i]=dim(jj[[i]])[1]
  maxhours1000[i]=ifelse(length(dim(jj1000[[i]]))==0,0,dim(jj1000[[i]])[1])
  maxhours1500[i]=ifelse(length(dim(jj1000[[i]]))==0,0,dim(jj1000[[i]])[1])
  maxhours2000[i]=ifelse(length(dim(jj2000[[i]]))==0,0,dim(jj2000[[i]])[1])
  #^this makes it zero hours if there's no HYSPLIT output (like for Newark 32097)
}
max(maxhours1000)
max(maxhours1500)
max(maxhours2000)

###COMBINE SEPARATE ALTITUDES INTO ONE LIST
jj=c(jj1000,jj1500,jj2000)


#### Create index vectors for precip and non-precip trajectories ####
back.hour=6 #This is the time over which to calculate precipitation
precip=numeric(length=runnumber*naltitudes)
for (i in 1:length(precip)){
  if (sum(as.numeric(jj[[i]][1:back.hour,15])>0)) {precip[i]=sum(as.numeric(jj[[i]][1:back.hour,15]))}
  else {precip[i]=0}
}
# precip
precip.day=which(precip>0)
noprecip=which(precip==0)
length(precip.day)
length(noprecip)
length(precip.day)+length(noprecip)

#### Separate into precipitating vs. non-precipitating trajectories ####
# Precipitating
jj.coord=list()
for (i in 1:length(precip.day)){
  jj.coord[[i]]=jj[[precip.day[i]]][,c(3:6,9:20)]
}
jj.coord[[1]][1:10,]

# Non-precipitating
jj.noprecip=list()
for (i in 1:length(noprecip)){
  jj.noprecip[[i]]=jj[[noprecip[i]]][,c(3:6,9:20)]
}

#### Eliminate positive longitudes (ie, past the dateline) ####
# Only does this for precipitating trajectories
for (i in 1:length(precip.day)){
  temptraj=as.data.frame(jj.coord[[i]])
  temptraj.Long=numeric(length=nrow(temptraj))
  for (j in 1:nrow(temptraj)){
    temptraj$Long[j]=ifelse(temptraj$Long[j]>0,NA,temptraj$Long[j])
    temptraj2=cbind(temptraj,temptraj.Long)
    temptraj2=na.omit(temptraj)
    temptraj2=temptraj2[,-21]
  }
  jj.coord[[i]]=temptraj2
}

#### Partition by Season ####
Month.vector=numeric(length=length(precip.day))
for (i in 1:length(precip.day)){
  Month.vector[i]=jj.coord[[i]]$Month[1]
}
length(Month.vector)
MAM.index=which(Month.vector>=3 & Month.vector<=5)
length(MAM.index)
JJA.index=which(Month.vector>=6 & Month.vector<=8)
length(JJA.index)
SON.index=which(Month.vector>=9 & Month.vector<=11)
length(SON.index)
D.index=which(Month.vector>=12)
JF.index=which(Month.vector<=2)
DJF.index=sort(c(D.index,JF.index))
length(DJF.index)
#NOTE some months are of value NA so the length totals don't necessarily match

MAM.traj=list()
for (i in 1:length(MAM.index)){
  MAM.traj[[i]]=jj.coord[[MAM.index[i]]]
}

JJA.traj=list()
for (i in 1:length(JJA.index)){
  JJA.traj[[i]]=jj.coord[[JJA.index[i]]]
}

SON.traj=list()
for (i in 1:length(SON.index)){
  SON.traj[[i]]=jj.coord[[SON.index[i]]]
}

DJF.traj=list()
for (i in 1:length(DJF.index)){
  DJF.traj[[i]]=jj.coord[[DJF.index[i]]]
}


jj.coord=jj.coord[which(is.na(Month.vector)==F)]
#^this just removes the ones with NA month values....

#### Plot all trajectories that produce precipitation over the site ####

# All Precipitating Trajectories
par(mfrow=c(1,1))
map(database="state",xlim=c(-130,-90),ylim=c(30,50))
for (i in 1:length(precip.day)){
  lines(jj.coord[[i]][,7],jj.coord[[i]][,6],lwd=0.1)      
}
#points(locations$Long,locations$Lat,col="red",pch=16)
points(jj.coord[[1]][1,7],jj.coord[[1]][1,6],col="blue",pch=16)

# DJF Trajectories
par(mfrow=c(1,1))
map(database="state",xlim=c(-130,-90),ylim=c(30,50))
for (i in 1:length(DJF.index)){
  lines(DJF.traj[[i]][,7],DJF.traj[[i]][,6],lwd=0.1)      
}
points(DJF.traj[[1]][1,7],DJF.traj[[1]][1,6],col="blue",pch=16)

# MAM Trajectories
par(mfrow=c(1,1))
map(database="state",xlim=c(-130,-90),ylim=c(30,50))
for (i in 1:length(MAM.index)){
  lines(MAM.traj[[i]][,7],MAM.traj[[i]][,6],lwd=0.1)      
}
points(MAM.traj[[1]][1,7],MAM.traj[[1]][1,6],col="blue",pch=16)

# JJA Trajectories
par(mfrow=c(1,1))
map(database="state",xlim=c(-130,-90),ylim=c(30,50))
for (i in 1:length(JJA.index)){
  lines(JJA.traj[[i]][,7],JJA.traj[[i]][,6],lwd=0.1)      
}
points(JJA.traj[[1]][1,7],JJA.traj[[1]][1,6],col="blue",pch=16)

# SON Trajectories
par(mfrow=c(1,1))
map(database="state",xlim=c(-130,-90),ylim=c(30,50))
for (i in 1:length(SON.index)){
  lines(SON.traj[[i]][,7],SON.traj[[i]][,6],lwd=0.1)      
}
points(SON.traj[[1]][1,7],SON.traj[[1]][1,6],col="blue",pch=16)

#### Plot all trajectories that don't produce precipitation ####
# jj.noprecip=list()
# for (i in 1:length(noprecip)){
#   jj.noprecip[[i]]=jj[[noprecip[i]]][,c(3:6,9:20)]
# }

map(database="state",xlim=c(-130,-90),ylim=c(30,50))
for (i in 1:length(noprecip)){
  lines(jj.noprecip[[i]][,7],jj.noprecip[[i]][,6],lwd=0.1)      
}
points(jj.coord[[1]][1,7],jj.coord[[1]][1,6],col="blue",pch=16)

#### Write precipitation and non-precipitation lists to file ####
#setwd("/Users/jeremycaves/Desktop/Box Sync/EESS 300 - Spring 2015 Terrestrial Paleoclimate/HYSPLIT/R Output")
#setwd("/Users/jeremycaves/Desktop/Box Sync/SoCo 2015/Guidebook/Dubois HYSPLIT/R and NCL Files")
setwd('/Users/Annie/Desktop/HYSPLIT output/RData files/Hari')
#CHANGE NAMES AS APPROPRIATE!!!!
Verdi.precip=jj.coord
Verdi.MAM=MAM.traj
Verdi.JJA=JJA.traj
Verdi.SON=SON.traj
Verdi.DJF=DJF.traj
Verdi.noprecip=jj.noprecip
save(Verdi.precip,file="verdi.precip_annual.RData")
save(Verdi.MAM,file="verdi.MAM.precip.RData")
save(Verdi.JJA,file="verdi.JJA.precip.RData")
save(Verdi.SON,file="verdi.SON.precip.RData")
save(Verdi.DJF,file="verdi.DJF.precip.RData")
save(Verdi.noprecip,file="verdi.noprecip.RData")

# Load any presaved files
#load("Coso.SON.precip.RData")

#### Write trajectories to a file ####
# NEED TO FIX THESE NEXT LINES OF CODE
#setwd("/Users/jeremycaves/Desktop/Assignments/Mongolia/HYSPLIT/Results/Contour Figure")

#Adds storm track ID to each storm track
# Currently set to only keep lat and long and add storm track ID
traj.analyze=SON.traj #jj.coord #SON.traj
jj.coord.st=list()
for (i in 1:length(traj.analyze)){ 
  ST=matrix(rep(i,length(traj.analyze[[i]][,5])),ncol=1,nrow=length(traj.analyze[[i]][,5]))
  jj.coord.st[[i]]=cbind((traj.analyze[[i]][,6:7]),ST)
}

st=as.matrix(jj.coord.st[[1]])
for (i in 2:length(jj.coord.st)){
  st=rbind(st,jj.coord.st[[i]])
}
dim(st)

# write.table(jj.coord.st[[1]][,],file="swendover.traj.txt",
#             row.names=FALSE,col.names=FALSE,sep=",")
# for (i in 2:length(precip.day)){ 
#   write.table(jj.coord.st[[i]][,],file="swendover.traj.txt",append=TRUE,
#               col.names=FALSE,row.names=FALSE,sep=",")
# }

#Should fix this step sometime in the future
# big <- read.table('swendover.traj.txt',sep = ',', header = FALSE)
big=st
# big=as.data.frame(big)
colnames(big) = c("Lat","Long","ST")

##CREATE AN EMPTY GRID TO STORE HISTOGRAM DATA
grid_p <- array(0, dim = c(180*2-1,90*2-1))
grid_n <- array(0, dim = c(180*2-1,90*2-1))

for(i in 1:max(big$ST)){
  # MAKE ARRAY OF COORDINATES
  coords <- big[which(big$ST == i),1:2]
  # BIN COORDINATES INTO 0.5X0.5 DEGREE GRID
  coords <- floor(coords*2)/2
  # REMOVE DUPLICATES SO THAT A GRID SPACE CAN ONLY BE COUNTED ONCE PER STORM TRACK
  coords<-unique(coords)
  # RECORD HISTOGRAM COUNT IN THE GRID, SEPERATE MATRICIES FOR NEGATIVE AND POSITIVE LON'S
  for(k in 1:dim(coords)[1]){
    if(coords[k,2]<0){
      grid_n[abs(coords[k,2]*2)-1,coords[k,1]*2-1] <- grid_n[abs(coords[k,2]*2)-1,coords[k,1]*2-1] + 1}
    else{
      grid_p[coords[k,2]*2-1,coords[k,1]*2-1] <- grid_p[coords[k,2]*2-1,coords[k,1]*2-1] + 1}
  }
}    

##CONVERT GRID TO PERCENTAGES
grid_p <- grid_p/max(big$ST)
grid_n <- grid_n/max(big$ST)
grid_all <- rbind(grid_n[359:1,], grid_p)

##SET CONTOUR LEVELS
levels = c(1,.7,.5,.4,.3,.2,.15, 0.1, .05,.025)
#levels=c(625,300,150,75,40,20,10,5)
colors = c('black', 'violet', 'purple', 'blue','cyan', 'green', 'yellow', 'orange', 'red',"darkred")
#colors=rainbow(length(levels))

##PLOT BACKGROUND MAP
map('state',interior=TRUE)#,xlim =c(75,120),ylim=c(20,60),col='gray')

##PLOT CONTOUR LINES
contour(c(seq(-180, -1, by = 0.5),seq(1, 180, by = 0.5)), seq(1, 90, by = 0.5), 
        grid_all, add = TRUE, 
        levels = levels, 
        col = colors, 
        drawlabels = FALSE)

##MAKE NetCDF FILE OF HISTOGRAM GRID FOR NCL MAPPING
x <- dim.def.ncdf("lon", "degrees_east", 
                  c(seq(-180, -1, by = 0.5),seq(1, 180, by = 0.5)))
y <- dim.def.ncdf("lat", "degrees_north", seq(1,90, by = 0.5))

storm <- var.def.ncdf("storm", "percent", list(x,y), -1)

histogram <- create.ncdf("histverdi.SON.nc", storm) #CHANGE
put.var.ncdf(histogram, storm, grid_all)
close.ncdf(histogram)

#### Plot Specific Humidity ####
plot(jj.coord[[1]][,10]~jj.coord[[1]][,1],type="l",xlab="Back Hour",
     ylab="Specific Humidity",bty="n",ylim=c(0,15)
     #,xlim=c(40,120)
)
for (i in 2:length(precip.day)){
  lines(jj.coord[[i]][,10]~jj.coord[[i]][,1])
}
abline(v=jj.coord[[1]][1,3],lty="dashed")

meanSH=numeric(length=back.plot)
meanSHarray=array(dim=c(back.plot,length(precip.day)))
for (i in 1:back.plot){
  for (j in 1:length(precip.day)){
    meanSHarray[i,j]=jj.coord[[j]][i,10]
  }
}
for (i in 1:back.plot){
  meanSH[i]=mean(meanSHarray[i,])
}
length(meanSH)
meanSH=cbind(jj.coord[[1]][,1],meanSH)
colnames(meanSH)=c("Back.Hour","SH")
lines(x=meanSH[,1],y=meanSH[,2],lwd=5,col="blue")

#Writes all non-precipitation producing trajectories to a csv file (for Tseterleg)
#Adds ID to each trajectory
# Currently set to only keep lat and long and add trajectory ID
jj.coord.st=array(dim=c(back.plot,3,length(noprecip)))
for (i in 1:length(noprecip)){
  ST=matrix(rep(i,back.plot),ncol=1,nrow=back.plot)
  jj.coord.st[,,i]=cbind(jj.noprecip[,2:3,i],ST)
}

#This for-loop writes all trajectories to a single file
write.table(jj.coord.st[1:back.plot,,1],file="tlgalltrajnoprecip.txt",
            row.names=FALSE,col.names=FALSE,sep=",")
for (i in 2:length(noprecip)){
  write.table(jj.coord.st[1:back.plot,,i],file="tlgalltrajnoprecip.txt",append=TRUE,
              col.names=FALSE,row.names=FALSE,sep=",")
}

#This for-loop writes each trajectory to a different csv file
for (i in 1:length(precip.day)){
  write.table(jj.coord[1:back.plot,,i],file=paste("tg",i,sep=""),sep=",",
              row.names=FALSE,quote=FALSE)
}

#Turns all of the trajectories into a single matrix
alltraj=matrix(ncol=2,nrow=(back.plot*dim(jj.coord)[3]))
s=seq(from=1,to=(back.plot*dim(jj.coord)[3]),by=back.plot)
e=seq(from=back.plot,to=(back.plot*dim(jj.coord)[3]),by=back.plot)
for (i in 1:dim(jj.coord)[3]){
  alltraj[s[i]:e[i],]=jj.coord[,2:3,i]
}
colnames(alltraj)=c("Lat","Long")
#write.table(alltraj,file="dzeregalltraj.txt",sep=",",row.names=FALSE,col.names=FALSE) #Writes "alltraj" to a .txt file