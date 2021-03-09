# This file is to test the functioning of the activity_metric_functions20201203.r (or older).
# I will later transfer bits of this script to the example markdown and to a 
# file where we produce results for the paper. 

#	4. Functions require the following information. 
#		ClockTime: time of the event as recorded by local time (15:30 becomes 15.5)
#		day: day of the month at which the event ocurred
#		month: number of the month of the year at which the event ocurred
#		lat: latitude in degrees, positive (north) and negative (south) of the 
#		equator
#		long: longitude in degrees, positive (west), and negative (east) of the
#	 	equator
#		timeZone


library(activity)
library(insol)
library(pbapply)
library(circular)
source("activity_metric_functions20201203.R")
# Load data ----

# need to add the errorbarplot# source("01_errorbarplot.r")
# source("activity_metric_functions20201203.r") # This file has changed names 2015. 

# Aguti - diurnal
# Armadillo - nocturnal
# peccary looks crepuscular. 
# Brocket - cathemeral

data(BCItime)

table(BCItime$species)
#Fit with confidence limits (limited reps to speed up)

par(mfrow=c(3,3))
species.classifications<-list()
for (i in 1:length(unique(BCItime$species))){
species <- subset(BCItime, species==unique(BCItime$species)[i])
time.rad<-2*pi*species$time
mod <- fitact(time.rad, sample="data", reps=10)
plot(mod, "radians", "density", main=paste(species$species[1]), las=1)
datetime <- as.POSIXct(species$date, format = "%d/%m/%Y %H:%M")
species.classifications[[i]]<-ActType(datetime, 9.15, -80, -5, 18)
}
names(species.classifications)<- unique(BCItime$species)

# Activity type ----

# First we need to have the date of the event, including latitude
# For activity type, we need to know how many observations took place when. 
time.clock <- species$time*24
datetime <- as.POSIXct(species$date, format = "%d/%m/%Y %H:%M")
day      <- format(datetime, "%d")
month    <-format(datetime, "%m")


#ActClass: returns a dataframe with number of observations per part of the day. 
#Act type: runs a chi square test on the distribution of the observations across
#  		   the day. Then classifies the animals as: 
# NOCTURNAL: More than 90% Observations during the night 
# MOSTLY NOCTURNAL: 70-90% observations during the night
# DIURNAL: less than 10 % Observations during the night (why not 90% during day?)
# MOSTLY DIURNAL: less than 10-30% Observations during the night
# crepsUSCULAR: 50% of records during twilight
# CATHERMERAL: all others or any with a non significant p value. 
#datetime, lat, long, tmz, creps=18
# input: datetime = a datetime object
#        lat = latitude - a unique value, or one as long as the datetime object. 
#        long = longitude - a unique value, or one as long as the datetime object.
#        tmz = time zone - a unique value, or one as long as the datetime object.
# 		   creps= a value for changing the definition of twilight representing the 
#		     number of degrees the sun travels below the horizon (civil, nautical or 
#		     astronomical. Default is astronomical twilight (18º below the horizon)). 
# Requires input from SunInfo. and the insol package. 
# Change order of ourcome for the table, as is comes out alphabetical

ActClass(datetime=datetime, lat=9.15, long=-80, tmz=-5)
#         Time of day Total observations % observations
# 1               day              12681           71.2
# 2    day crepuscule               1314            7.4
# 3             night               2647           14.9
# 4 night crepsuscule               1178            6.6

ActType(datetime, 9.15, -80, -5, 18)

# Chi-squared test for given probabilities
# 
# data:  chivec
# X-squared = 3390.4, df = 2, p-value < 2.2e-16
# 
# [1] "mostly diurnal"


# Activity peaks ----
# tp.est
# requires tp.est, tp.def, match.tp, tp.boot
# There is also tp.plot. Although I think you just need to add the points from 
# the table to the plot. 
#Estimates positions, heights and their SEs for 1st derivative turning points of a von Mises kernel density
#Output table with columns:
#	1: turning point position
#	2: turning point height
#	3: position SE
#	4: height SE
#	5: pairwise significance of differences between given row and the following
#	6: groupings based on significance of pairwise comparisons
#	7: flags "true" turning points (max/min within groups exluding those intermediate groups) 
Sig.turning.points<-tp.est(time.rad) # Works (as in gives a table, but many warnings)

warnings()
# 1: In min(x, na.rm = T) : no non-missing arguments to min; returning Inf

tp.plot(mod, Sig.turning.points, main=paste(species$species[1]), las=1) #xunit doesn't work

# We need to make sure people can get this in hours and freq! 


onset.end<-onset(mod=mod, tp=Sig.turning.points)
onset.end
abline(v=onset.end[1,1], col="red")
abline(v=onset.end[2,1], col="red")
abline(v=onset.end[3,1], col="red")
abline(v=onset.end[4,1], col="red")


# Polygons for change in sunrise/sunset # 


info <- SunInfo(datetime, 9.15, -80, -5, ) #Difference in sunrise-sunset in the study time

maxHrise <- max(info$Hrise)*2*pi/24
minHrise <- min(info$Hrise)*2*pi/24
maxHset <- max(info$Hset)*2*pi/24
minHset <- min(info$Hset)*2*pi/24


#Visual representation
#Adds polygons to better distinguish Sunrise and sunset


polygon(x=c(0,maxHrise,maxHrise,0), y=c(0,0,10,10), col=rgb(0.19, 0.19, 0.19, 0.3), border=NA)
polygon(x=c(minHset,(2*pi),(2*pi),minHset),y=c(0,0,10,10), col=rgb(0.19, 0.19, 0.19, 0.3), border=NA)
polygon(x=c(0,minHrise, minHrise,0), y=c(0,0,10,10), col=rgb(0.19, 0.19, 0.19, 0.3), border=NA)
polygon(x=c(maxHset,(2*pi),(2*pi),maxHset),y=c(0,0,10,10), col=rgb(0.19, 0.19, 0.19, 0.3), border=NA)
box()

