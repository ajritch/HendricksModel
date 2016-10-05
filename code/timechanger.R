#timechanger function
#converts Year/Month/Day/Hour to Hours since 1800-1-1 00:00:00.0

#created 10-22-15 by AJR
#last modified 11-10-15



#takes in dataframe of precipitating trajectories and their met. variables
#dataframe has 16 columns; add a 17th with the new time!
#also adds a new top row of empties for the initialization point
# ^this is currently commented out as it is unnecessary
timechanger=function(df) {
  vec=df[1,1:4]
  backhours=df$Back.Hour
  #vec=c(Year,Month,Day,Hour)
  Year=as.numeric(vec[1])
  Month=as.numeric(vec[2])
  Day=as.numeric(vec[3])
  Hour=as.numeric(vec[4])
  #Year is in years after 1900/2000 (ie, 2004 has Year=4, 1980 has Year=80)
  Y2K=0 #year 2000
  Y79=79 #year 1979
  
  
  #calculate number of hours since a known reference
  #must execute as either before or after year 2000
  if (Year<16) { #Execute this part for years after (and including) 2000
    #establish whether or not it's a leap year
    yeardiff=Year-Y2K #difference between Year and 2000
    remainder=yeardiff%%4
    Leap=F #default not a leap year
    if (remainder==0) Leap=T  #leap year True if ....if it's a leap year... :P
    
    
    #how many leap years have passed between Year and 2000? (2000 counts as one)
    yeardiv=(yeardiff-1)/4
    nleapyears=1+floor(yeardiv) #add one to include 2000 in the count
    
    #add on hours from 2000-01-01 to YEAR-01-01 (at midnight)
    Y2K.hours=1753152
    #number of hours at Jan 1, midnight, of year at hand
    Yearhours.Jan1=Y2K.hours+(yeardiff*365*24)+(nleapyears*24)
    
  } else { #Execute this part for years before 2000
    #establish whether or not it's a leap year
    yeardiff=Year-Y79 #number of years between Year and 1979
    remainder=yeardiff%%4
    Leap=F
    if (remainder==0) Leap=T #indicates leap year TRUE
    
    #how many leap years have passed between Year and 1979??
    yeardiv=(yeardiff-1)/4
    yeardiv=ifelse(yeardiv>0,yeardiv,0)
    nleapyears=ceiling(yeardiv) #it's magic. it works. don't question.

    #add on hours from 1979-01-01 to YEAR-01-01 (at midnight)
    Y79.hours=1569072
    #number of hours at Jan 1, midnight, of year at hand
    Yearhours.Jan1=Y79.hours+(yeardiff*365*24)+(nleapyears*24)
  }
  
  
    
  #add on hours until midnight of first day of each month
  #months are all different, so just do it by month, I guess
    
  JanHours=31*24
  FebHours=28*24
  FebHours.leap=29*24
  MarHours=31*24
  AprHours=30*24
  MayHours=31*24
  JunHours=30*24
  JulHours=31*24
  AugHours=31*24
  SepHours=30*24
  OctHours=31*24
  NovHours=30*24
  
  MonthHoursVec=c(0,JanHours,FebHours,MarHours,AprHours,MayHours,JunHours,
                  JulHours,AugHours,SepHours,OctHours,NovHours)
  MonthHoursVec.leap=c(0,JanHours,FebHours.leap,MarHours,AprHours,MayHours,JunHours,
                       JulHours,AugHours,SepHours,OctHours,NovHours)
  #^zeros are placeholders because there are no hours before start of Jan
  
  #hours until beginning of first day of the month, depending on leap year
  MonthHours=sum(MonthHoursVec[1:Month])
  if (Leap==T) { #if current year is a leap year
    MonthHours=sum(MonthHoursVec.leap[1:Month])
  }
  
  #add on hours until start of the day
  DayHours=(Day-1)*24
  
  #add on hours until start of the hour
  HourHours=Hour 
  
  #add them all together!
  #this is hour since 1800-01-01 00:00:00.0 at our given point
  PointHour=sum(Yearhours.Jan1,MonthHours,DayHours,HourHours)


  #okay, we've got the hour for our first point along trajectory.
  #now add an additional column to the d.f. for the PointHour (call it "Time")
  Time=PointHour+backhours   #backhours MUST start with a zero!!!! (at initialization point)
  #newdf=cbind(df,Time) #17 met variables total now, with Time being #17
  #keeping everything is making gigantic lists. trim unnecessary junk:
  newdf=cbind(df[c('Year','Back.Hour','Lat','Long')],Time)
  
  #commenting this out unless initialization point isn't included for some reason
  # #while we're manipulating the dataframe, add the initialization point
  # #add a first column for Back.Hour=0, with appropriate Time
  # #everything else (Lat/Long, most importantly) can be filled in later
  # toprow=c(rep(NA,4),0,rep(NA,11),(PointHour+1))
  # newdf2=rbind(toprow,newdf)
  
  #awesome!
  #return(newdf2) #only if initialization point isn't included from TrajectoryPlotter
  return(newdf)

}