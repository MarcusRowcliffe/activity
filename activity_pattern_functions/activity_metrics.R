# Functions 
# Lets check

###########################################################################
## Sun Information: We use package insol to generate data on sunrise and ##
## sunset, we have tested this and worked best in out area of study.     ##
## given time, date, latitude and longitude. 					                   ##
## Input:                            						                         ##
## 	lat=latitude, long=longitude, tmz=timezone used in the study,        ##
##    clocktime in decimals,  datetime= object of class (POSIXct)        ##
##  creps: value of twilight you wish to use (civil, astronomical, etc)  ##
##         default=18													                           ##
## Output:                  											                       ##
##	dataframe ClockTime,Hrise,Hset,Daylength,NightTwil, Midday, Midnight ##
###########################################################################

SunInfo<- function(datetime, lat, long, tmz, creps=NULL){
  if (is.null(creps)) creps <- 18 else creps<-creps
  
  jd<-JD(datetime) #string date 8/31/56, 8-31-1956, 31 8 56, 083156,31Aug56, or August 31 1956.
  insol<-daylength(lat=lat, long=long, jd=jd, tmz=tmz) 
  Hrise<-insol[,1]  
  Hset<-insol[,2]
  Midday<-Hrise+((Hset-Hrise)/2)	
  Daylength<-(Hset-Hrise)
  Night.with.Twil<-(24-Daylength)
  Twilightlength<- 2*Night.with.Twil*creps/100 
  Nightlength<-Night.with.Twil-Twilightlength
  mn <- Hrise- (Night.with.Twil/2)
  mn24<-mn+24
  Midnight<-list()
  for (i in 1:length(mn))
    if (mn[i]<0) Midnight[i]<- mn24[i] else Midnight[i]<- mn[i] 
  Midnight<-unlist(Midnight)
  data.frame(cbind(Hrise,Hset,Daylength,Twilightlength, Nightlength, Night.with.Twil, Midday, Midnight))
}


###########################################################################
## Activity Classiffication: Using Information about the sun, classifies ##
## events as diurnal (ocurred between sunrise and sunset), crepsuscular  ##
## (between the time the sun reaches the horizon and travels 18 degrees) ##
## Input:
##	 timedate, lat, long, tmz,  clocktime, creps: civil (6dg),      ##
##   nautical (12dg), or astronomical (18dg) twilight                    ##
## output: table with number of observations per type of day		         ##
########################################################################### 


ActClass <- function (datetime, lat, long, tmz, creps=18){
  if (creps>18) message("Warning: creps > 18. When the sun is 18 degrees below the surface, astronomical twilight has not yet started")
  y<-SunInfo(datetime, lat, long, tmz, creps)
  Hset<-y$Hset
  Hrise<-y$Hrise
  twilight<-y$Twilightlength
  hour<-as.numeric(format(as.POSIXct(datetime, format = "%d/%m/%Y %H:%M"), "%H"))
  minute<-as.numeric(format(as.POSIXct(datetime, format = "%d/%m/%Y %H:%M"), "%M"))
  clocktime<-hour+minute/60
  class<-NULL
  
  for (i in 1:length(clocktime)){
    if (Hrise[i] < clocktime[i] && Hset[i] >= clocktime[i]) class[i]<- "day" else
      if (clocktime[i]>Hset[i] && (Hset[i] + (twilight[i]/2))>= clocktime[i]) class[i]<- "night crepsuscule" else 
        if (Hrise[i] - (twilight[i]/2) < clocktime[i] && Hrise[i] >= clocktime[i]) class [i]<-"day crepuscule" else
          class[i]<- "night"}
  act.table<-data.frame(table((class)))
  res<-data.frame(matrix(nrow=4, ncol=3))
  names(res)<-c("Time of day", "Total observations", "% observations")
  res[ ,1] <- c("day", "day crepuscule", "night", "night crepsuscule")
  res[1,2] <- sum(class==res[1,1])
  res[2,2] <- sum(class==res[2,1])
  res[3,2] <- sum(class==res[3,1])
  res[4,2] <- sum(class==res[4,1])
  res[ ,3] <- round(res[,2]*100/sum(res[,2]),1)
  res
}

################################################################################
# Activity type: Chisquare test to determine the uniformity of the distribution #
# In the case that the pattern is not uniform, we follow Gomez et al 2005 to   #
# classify the animals as nocturnal, diurnal or crepsuscular. If the pattern is #
# uniform, we consider the species cathermeral. 							   #
# NOCTURNAL: More than 90% Observations during the night (they say dark Does this mean when the sun is below the horizon?)
# MOSTLY NOCTURNAL: 70-90% observations during the night
# DIURNAL: less than 10 % Observations during the night (why not 90% during day?)
# MOSTLY DIURNAL: less than 10-30% Observations during the night
# crepsUSCULAR: 50% of records during twilight
# CATHERMERAL: all others
# However, they defined crepuscule as anything that occured between an hour before 
# and an hour after sunset. and the same for sunrise. Our def does not include 
# any time during the "day' 
#################################################################################


ActType<-function(datetime, lat, long, tmz, creps=18){
  y<-ActClass(datetime, lat, long, tmz, creps)
  w<-SunInfo(datetime, lat, long, tmz , creps)
  avblDay<- mean(w$Daylength)/24
  avblNight<-mean(w$Nightlength)/24
  avblTwilight<-mean(w$Twilightlength)/24
  expectedprop<-c(avblDay, avblNight, avblTwilight)
  
  chivec<-c(y[1,2], y[3,2], y[2,2]+y[4,2])	
  

  res<- list()
  res[[1]]<-y
  res[[2]]<-chisq.test(chivec, p=expectedprop)
  res[[3]]<-(ifelse(res[[2]]$p.value>= 0.05, "cathemeral", 
                    ifelse(chivec[2]/sum(chivec)<=0.1, "diurnal",
                          ifelse(chivec[3]/sum(chivec)>=0.5, "crepuscular",
                                ifelse(chivec[1]/sum(chivec)<=0.1, "nocturnal",
                                      ifelse(chivec[2]/sum(chivec)>0.1 && chivec[2]/sum(chivec)<=0.3, "mostly diurnal", 
                                            ifelse(chivec[1]/sum(chivec)>0.1 && chivec[1]/sum(chivec)<=0.3, "mostly nocturnal", "cathemeral" )))))))
names(res)<- c("Classification of activity", "chisq test", "activity type")
res
}


########################################
# Definition of Turning points
#######################################

#Estimates positions, heights and their SEs for 1st derivative turning points of a von Mises kernel density
#Output table with columns:
#	1: turning point position
#	2: turning point height
#	3: position SE
#	4: height SE
#	5: pairwise significance of differences between given row and the following
#	6: groupings based on significance of pairwise comparisons
#	7: flags "true" turning points (max/min within groups excluding those intermediate groups) 

tp.def <- function(mod, nx=1000){
  numderiv <- function(y) (y[-1] - y[-length(y)]) 
  x<-seq(0,2*pi,pi/nx) 
  hts <- dvmkern(x, dat=mod@data, adj=mod@adj, bw=mod@bw) 
  nd1 <- numderiv(hts)*nx/pi
  nd2 <- numderiv(c(nd1[2*nx],nd1))*length(x)/pi
  i <- which(sign(nd1) != sign(c(nd1[-1],nd1[1])))
  tp1 <- x[i+1]
  ht1 <- hts[i+1]
  j <- which(sign(nd2) != sign(c(nd2[-1],nd2[1])))
  tp2 <- x[j]+pi/(2*nx)
  ht2 <- (hts[j]+hts[j+1])/2
  ord <- c(rep(1,length(tp1)), rep(2,length(tp2)))[order(c(tp1,tp2))]
  tp <- sort(c(tp1,tp2))
  ht <- c(ht1,ht2)[order(c(tp1,tp2))]
  k <- (2*1:(length(tp)/2))
  i1 <- which(tp %in% tp1)[1] 
  if(i1/2-floor(i1/2)==0.5) k <- k-1
  tp <- tp[k]
  ht <- ht[k]
  n2 <- which(ord[k]==2)
  truetp <- rep(c(TRUE,FALSE), c(length(tp),length(n2)))
  truetp[n2] <- FALSE
  tp <- c(tp, tp[n2])
  ht <- c(ht, ht[n2])
  truetp <- truetp[order(tp)]
  ht <- ht[order(tp)]
  tp <- sort(tp)
  cbind(tp,ht,truetp)
}

match.tp <- function(btp,tp){	
  bn <- length(btp)
  n <- length(tp)
  difmat <- matrix(abs(rep(btp,n)-rep(tp,each=bn)), nrow=bn)
  difmat[difmat>pi] <- 2*pi-difmat[difmat>pi]	
  q <- apply(difmat, 1, function(x) which(x==min(x)))
  if(is.list(q)){
    qi <- which(lapply(q,length)>1)
    q[qi] <- lapply(q[qi], function(x) sample(x,1))
    q <- unlist(q)}
  
  tm <- difmat/t(apply(difmat, 1, function(x) x==min(x)))
  tm[tm==Inf] <- NA
  w <- unlist(apply(tm,2,function(x) which(x!=min(x,na.rm=T))))  # There is not min or max. 
  pkres <- tp[q]
  q[w] <- NA
  q
}

tp.boot <- function(i,btps,mod, nx=1000){	
  x <- mod@data
  smod <- fitact(sample(x,length(x),replace=T)) # Think this has changed! 
  tps <- tp.def(smod,nx)
  bm <- 2*(1:(length(btps[,1])/2))
  m <- 2*(1:(length(tps[,1])/2))
  bpk1st <- ifelse(min(btps[bm,2]-btps[bm-1,2])>=0,0,1)
  pk1st <- ifelse(min(tps[m,2]-tps[m-1,2])>=0,0,1)
  
  qpk <- match.tp(btps[bm-bpk1st,1], tps[m-pk1st,1]) ##generates NAS
  qtr <- match.tp(btps[bm+bpk1st-1,1], tps[m+pk1st-1,1]) ##generates NAS
  res <- htres <- matrix(nrow=length(btps[,1]), ncol=2)
  res[bm-bpk1st,1] <- tps[m-pk1st,1][qpk]
  res[bm+bpk1st-1,1] <- tps[m+pk1st-1,1][qtr]
  res[bm-bpk1st,2] <- tps[m-pk1st,2][qpk]
  res[bm+bpk1st-1,2] <- tps[m+pk1st-1,2][qtr]
  res
}

tp.est <- function(x, nx=100, its=100){	
  ptest <- function(mns,ses){	
    tval <- abs(mns[-1]-mns[-length(mns)])/sqrt(ses[-1]^2+ses[-length(ses)]^2)
    pt(tval,Inf,lower.tail=F)}
  funky <- function(i,tp) {	
    q <- tapply(tp[,2], tp[,6], mean)
    j <- i+(-1:1)
    if(i==1) j[1] <- length(q)
    if(i==length(q)) j[3] <- 1
    df <- diff(q[j])
    if(sign(df[1])==sign(df[2])) NA else
      if(df[1]>0) which(tp[,2]==max(tp[tp[,6]==i,2])) else
        which(tp[,2]==min(tp[tp[,6]==i,2]))}
  ######Find all turning points and their SEs
  mod<-fitact(x)
  btps <- tp.def(mod,nx)
  n <- length(btps[,1])
  bootres <- pbsapply(1:its, tp.boot, btps, mod, nx) ##nx
  # There are warnings here! 
  tpSE <- apply(bootres[1:n,],1,sd.circular,na.rm=T) # sd.circular
  htSE <- apply(bootres[(n+1):(2*n),],1,sd,na.rm=T)
  i <- btps[,3]==1
  res <- cbind(btps[i,1:2], tpSE=tpSE[i], htSE=htSE[i])
  ht.p=ptest(c(res[,2],res[1,2]), c(res[,4],res[1,4]))
  names(ht.p) <- NULL
  sig <- ht.p<0.05
  gp <- cumsum(sig)+1-sig
  if(!sig[length(sig)]) gp[gp==max(gp)] <- 1
  tptab <- cbind(res, ht.p, gp)
  ######Add column defining whether a "true" turning point:
  ######	excludes groups between a trough and a peak, 
  ######	and finds max/minima within groups (peaks/troughs respectively)
  j <- sapply(1:max(tptab[,6]), funky, tptab)
  incl <- vector(length=nrow(tptab)) 
  incl[j] <- 1
  cbind(tptab,incl)
  
}

###############################################################################
# Plot Turning point function: function of a model class actmod, a tpest result 
# with the possibility of plotting the error of the turning points (default)
# and plotting only the significant turning points (plot.sig=TRUE)
###############################################################################

tp.plot<-function(mod, tpest, plot.error=TRUE, plot.sig=FALSE, ...){
  if (plot.sig==TRUE) tpest<-tpest[tpest[,7]==1,] else tpest<-tpest
  sigtpest<-tpest[tpest[,7]==1,]
  errorx<-data.frame(ul=(tpest[,1]+tpest[,3]), ll=(tpest[,1]-tpest[,3]))
  errory<-data.frame(ul=(tpest[,2]+tpest[,4]), ll=(tpest[,2]-tpest[,4]))
  plot(mod, xunit = c("radians"),yunit = c("density")) 
  if (plot.error==FALSE) {
    points(x=tpest[,1], y=tpest[,2])
    points(x=sigtpest[,1], y=sigtpest[,2], pch=16)
  } else {
    points(x=tpest[,1], y=tpest[,2])
    points(x=sigtpest[,1], y=sigtpest[,2], pch=16)
    segments(x0=tpest[,1], y0=errory[,1], x1 = tpest[,1], y1 = errory[,2])
    segments(x0=errorx[,1], y0=tpest[,2], x1 = errorx[,2], y1 =tpest[,2])
  }
}


# Make degree of bimodality and amplitude a result of tp.est. 


##Degree of bimodality
##Either use a set of times, then get turning points OR make it a function 
## in which people just introduce a set of values (tp locations), and dob is 
## calculated from that point. This way, we give the user the choice to either
## let the program decide what is and is not bimodal, or ignore minor peaks
## and get DOB anyway.
## x is a dataframe containing a column for the location and height of the 
## turning points 
DOB <- function(x, tp.cal=TRUE){
  if (tp.cal==TRUE) x <- tp.est(x) else x <- x
  if (length(x[x[,7]==1,]) > 4 || length(x[x[,7]==1,])< 4) print("Species is not bimodal") else {
    #set the dataframe in order: highest ht to lowest ht
    # x[,1] becomes peak 1
    # x[,2] becomes peak 2
    # x[,3] becomes the lull
    # x[,4] becomes the minimum value
    order.ht <- order(x[,1], decreasing=TRUE)
    y <- x[order.ht,]
    res <- (((y[1,1]-y[1,3])+(y[1,2]-y[1,3]))/2)/y[1,1]
    res
  }}


# Need to add bootstrapping
# Add ratio and activity time and resting time. 
onset <- function(time.rad=NULL, mod=NULL, tp=NULL) {	
  if (is.null(time.rad) & is.null(mod)) warning ("time.rad and mod are both null. Input either time in radians or an acmod object")
  if (!is.null(time.rad) & sum(time.rad>2*pi)>0) ("time.rad should be expressed in radians. Some values are larger than 2*pi")
  is.actmod <- function(x) inherits(x, "actmod")
  if (!is.null(mod) &   !is.actmod(mod)) ("mod if not an actmod object. Please see actfit function")
  if (is.null(mod)) mod<-fitact(time.rad)
  if (is.null(tp)) tp <- tp.est(time.rad)
  numderiv <- function(y) (y[-1] - y[-length(y)]) 
  # True turning points
  n <- sum(tp[,7])							#number of true turning points
  ttp <- tp[tp[,7]==1,1]
  htp <- tp[tp[,7]==1,2]
  ttp <- c(ttp,ttp[1])		
  htp <- c(htp, htp[1])				#true turning points, wrapped
res <- data.frame(matrix(ncol=2 , nrow=n))						#empty results vector
for(i in 1:n) {
  if (i<n) { period1 <- mod@pdf[,1] >=ttp[i] & mod@pdf[,1] < ttp[i+1]} else 
  { period1 <- mod@pdf[,1] >=ttp[i] | mod@pdf[,1] < ttp[1]}
  rate.change.x<-numderiv(mod@pdf[,1][period1])
  rate.change.y<-numderiv(mod@pdf[,2][period1])
  abs.slope <- abs(rate.change.y/rate.change.x)
  sel<-abs.slope==max(abs.slope)
  res[i, 1] <- mod@pdf[,1][period1][sel]
  res[i, 2] <- mod@pdf[,2][period1][sel]
  if (dim(mod@pdf)[2]==5){
    res[i, 3] <- mod@pdf[,3][period1][sel]
    res[i, 4] <- mod@pdf[,4][period1][sel]
    res[i, 5] <- mod@pdf[,5][period1][sel]}}
if(dim(res)[2]==2) names(res)<-c("Radian time", "pdf") else 
     names(res)<-c("Radian time", "pdf", "std. error", "lower 95% CL", "upper 95% CL")
 
message("This function returns the points of highest change of rate between true turning points")
message("These represent onset and end of a population's activity only after a period of low/no activity")
message("Cathermeral species have no onset or end of activity")
res
}
 


# i 
# btps - is a table with the definition of the turning points 
# mod is the model. From there there is a sampling procedure. 
# tp.boot <- function(i,btps,mod, nx=1000){
#   x <- mod@data
#   smod <- fitact(sample(x,length(x),replace=T)) # Think this has changed! 
#   tps <- tp.def(smod,nx)
#   bm <- 2*(1:(length(btps[,1])/2))
#   m <- 2*(1:(length(tps[,1])/2))
#   bpk1st <- ifelse(min(btps[bm,2]-btps[bm-1,2])>=0,0,1)
#   pk1st <- ifelse(min(tps[m,2]-tps[m-1,2])>=0,0,1)
#   
#   qpk <- match.tp(btps[bm-bpk1st,1], tps[m-pk1st,1]) ##generates NAS
#   qtr <- match.tp(btps[bm+bpk1st-1,1], tps[m+pk1st-1,1]) ##generates NAS
#   res <- htres <- matrix(nrow=length(btps[,1]), ncol=2)
#   res[bm-bpk1st,1] <- tps[m-pk1st,1][qpk]
#   res[bm+bpk1st-1,1] <- tps[m+pk1st-1,1][qtr]
#   res[bm-bpk1st,2] <- tps[m-pk1st,2][qpk]
#   res[bm+bpk1st-1,2] <- tps[m+pk1st-1,2][qtr]
#   res
# }

