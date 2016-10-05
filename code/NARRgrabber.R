#grab NARR variables

#this function takes (modified) HYSPLIT output, regridded to NARR grid, and
  #tacks on the associated NARR meteorological variables to the dataframe
#the following meteorological variables are added:
    #Air Temperature at 2m (Temp_2m.year) (Kelvin)
    #Relative Humidity at 2m (rhum.year)(kg/kg)
    #Precipitable Water (Prw.year) (kg/m^2)
    #Time (time.year) (hours since 1800-1-1 00:00:0.0)
    #Damkohler Number (evap/(prate-evap))

#created: 11-4-15 by AJR
#last modified: 11-10-15 by AJR



#CHANGE WORKING DIRECTORY to the one with the NARR .RData files!
#setwd('/Users/Annie/Desktop/NARR_for_processing')



#begin with function for one row at a time; this function is called later in NARRgrabber()
#this function uses the lat/lon/time of a row of a df and picks out appropriate NARR data
NARRgrabber.onerow=function(vector) {
  #year.pre=vector$Year 
  year.pre=as.numeric(vector['Year'])
  year=ifelse(year.pre<16,year.pre+2000,year.pre+1900) #turn into an actual year format
  #^these lines repeat in NARRgrabber() as well. Here we need for individual year
     #in other function we need just to download appropriate files
  
  #select the indices we'll be using to dig out NARR data using lat/lons
  #dlati=which(lat.data==vector$Lat,arr.ind=T) #desired lat index (dlati)
  dlati=which(lat.data==as.numeric(vector['Lat']),arr.ind=T) #desired lat index (dlati)
  #dloni=which(lon.data==vector$Lon,arr.ind=T) #desired lon index (dloni)
  dloni=which(lon.data==as.numeric(vector['Lon']),arr.ind=T) #desired lon index (dloni)
  #^these give vectors where indices correspond to indices in lat.data and lon.data
  #ie lat.data[dlati[1],dlati[2]] returns desired lat VALUE
  
  #there are non-unique lat/lons!
  #length(as.numeric(lat.data)) #96673
  #length(unique(as.numeric(lat.data))) #96130
  #length(unique(as.numeric(lon.data))) #95662
  
  #desired coordinate indices because there are non-unique lats and lons
  #keep only the indices present in lat and lon (keep only duplicates)
  coordis=rbind(dlati,dloni)
  dcoordis=coordis[duplicated(coordis),,drop=F] #desired coordinate indices (dcoordis)
  
  
  yearfilelist=ls(envir=.GlobalEnv)[grepl(year,ls(envir=.GlobalEnv))] #list of all the files for given year
  yeartime=get(yearfilelist[grepl("time",yearfilelist)],envir=.GlobalEnv) #list of times for given year
  
  #get desired time index
  #dtimei=which(yeartime==vector$Time)
  dtimei=which(yeartime==as.numeric(vector['Time']))
  #dtimei=which(yeartime==1796997)
  
  #NARR is every three hours, so 2/3 of our HYSPLIT hours aren't extractable
  #for those that are extractable, pick out met values using associated lat/lon/time indices
  if (length(dtimei)==1) { #if NOT compatible, dtimei is length zero
    #extract from NARR files using desired lat/lon/time indices
    vector['Nd']=get(yearfilelist[grepl('Nd',yearfilelist)])[dcoordis[1],dcoordis[2],dtimei]
    vector['Prw']=get(yearfilelist[grepl('Prw',yearfilelist)])[dcoordis[1],dcoordis[2],dtimei]
    vector['rhum']=get(yearfilelist[grepl('rhum',yearfilelist)])[dcoordis[1],dcoordis[2],dtimei]
    vector['Temp_2m']=get(yearfilelist[grepl('Temp_2m',yearfilelist)])[dcoordis[1],dcoordis[2],dtimei]
  } #otherwise it all remains as NAs!
  
  #pick out the row's land mask from desired coord indices:
  vector['Land']=land.data[dcoordis[1],dcoordis[2]]
  
  return(vector)
  
}





#take in dataframe (with fixed coords and time), append appropriate NARR output variables
NARRgrabber=function(df) {
  #add five extra columns to the df for all of the synced info!!!!
  #also, make it a matrix so that it's compatible with apply() 
  dfnew=cbind(df,Land=NA,Nd=NA,Prw=NA,rhum=NA,Temp_2m=NA)
  #^leaving as NAs; the non-three-hourly times will remain NAs
  
  ###bookkeeping of necessary NARR data files
  #THIS IS CRITICAL so as to not exceed RAM space
  
  year.pre=df[1,]$Year #get year of initialization point
  year=ifelse(year.pre<16,year.pre+2000,year.pre+1900) #turn into an actual year format
  openname=paste('NARR',year,'.RData',sep='')
  #^this is for the year of the initialization point
  #the trajectory may extend into previous year, so load that up too!
  year.prev=year-1
  openname.prev=paste('NARR',year.prev,'.RData',sep='')
  
  #Now, we don't want to have too much gunk floating around in our environment!
  #get rid of variables from years past that we no longer need
  year.twoago=year-2
  rm(list=ls(envir=.GlobalEnv)[grepl(year.twoago,ls(envir=.GlobalEnv))],envir=.GlobalEnv) 
  #^removes anything from more than a year ago
  
  #if the files (variables) for a year have already been loaded, do nothing!
  #otherwise load the appropriate .RData file into environment
  dummyname=paste('Nd.',year,sep='') #arbitrarily choosing Nd
  if (exists(dummyname)==F) {
    load(openname,envir=.GlobalEnv)
  }
  #do the same for the previous year's file (previous year=year prior to initialization year)
  dummyname.prev=paste('Nd.',year.prev,sep='')
  if (exists(dummyname.prev)==F) {
    load(openname.prev,envir=.GlobalEnv)
  }

  

  ###grab appropriate NARR met variables using NARRgrabber.onerow() (defined above)
  outputdf=as.data.frame(t(apply(as.matrix(dfnew),1,NARRgrabber.onerow)))
  
  return(outputdf)
  
  
}


