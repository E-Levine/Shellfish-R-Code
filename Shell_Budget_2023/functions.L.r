standard.header<-function(){
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#                                                                              #
#                     STANDARD    HEADER                                       #
#______________________________________________________________________________#
#Remove all objects that currently exist in memory
rm(list=ls(all=TRUE))
#Load required packages
require(graphics)
require(lattice)
#source standard functions
dir.func<- "/home/dhennen/DRH_general_R_code/"
source(paste(dir.func,"functions.L.r",sep=""))
#or better yet - use the new package!!
library(devtools)
load_all("/home/dhennen/DRH_general_R_code/DanFunc")
library(DanFunc)

options(stringsAsFactors = F) #turn off factor reading!

#if a multi-plot is required here is a good template
plot.new()
resize.win(Width=15, Height=10, xpinch=30, ypinch=30)
par(mar=c(4,5,4,2)+0.1,cex.axis=1.5,cex.lab=2,cex.main=2) #bottom, left, top, right
#or for equal sized plots with wide outer margins...
#par(mar=c(0,0,0,0)+0.5,oma=c(8,8,8,8),cex.axis=1.5,cex.lab=2,cex.main=2) #bottom, left, top, right 
l.out<-matrix(c(1:6),3,2, byrow=T)
#l.out<-matrix(c(1,1,1,1,2,3),3,2, byrow=T)
l1=layout(l.out)  
layout.show(l1)
  
dir.graphics<-"/home/dhennen/"
filename<-paste(dir.graphics,"out.put1",sep="")
savePlot(filename=paste(filename,"png",sep="."),type="png")
#savePlot(filename=paste(filename,"tif",sep="."),type="tif")

#to get R to see a mounted drive on a linux box
setwd("~/.gvfs")
dir()
#or
test1<-read.csv("~/.gvfs/home2 on net/dhennen/test_2012_clam_db/REGION_MEANS.CSV")
list.files("~/.gvfs/home2 on net/dhennen/")
#Useful legend template
p1<-par()$usr
legend(
  x=0.75*(p1[2]-p1[1])+p1[1]
  ,y=0.9*(p1[4]-p1[3])+p1[3]
  ,legend=c("My Legend","My Legend2")
  ,col=c("col1","col2")
  ,lty=c(1,2)
  ,lwd=c(2,2)
  ,cex=2
  ,bty="n"
)  
#This is how to get math and evaluated variables into a label - bquote() and .() and ~ to separate the strings is key
l1<-bquote("Fish > 40 cm (" ~ bar(n[y]) == .(mean(M40n)) ~ ")")
#Note this doesn't really work in a legend - you have to monkey with text()!

#Function to clear memory after reading a big file... 
save(list = ls(all = TRUE), file = "MyData.RData")
clearMemory() #get rid of the memory pig you just made!
load("MyData.RData")
#This may work better - example
interim_object <- data.frame(rep(1:100, 10),
                             rep(101:200, 10),
                             rep(201:300, 10))
object.size(interim_object) # Reports the memory size allocated to the object
rm("interim_object") # Removes only the object itself and not necessarily the memory allotted to it
gc() # Force R to release memory it is no longer using
ls() # Lists all the objects in your current workspace
rm(list = ls()) # If you want to delete all the objects in the workspace and start with a clean slate

#little timer
ptm <- proc.time()
#process you want to time
proc.time() - ptm 

#Example of an oracle querry
#require(ROBDC)
#channel1 <- odbcConnect("nova","dhennen","$$dhe2###QWE")
#sel1<-paste("select a.year",sep="")
#from1<-paste("from oblnh",sep="")
#where1<-paste("where a.link3(+) = b.link3 and b.link3 = c.link3"
#  ,"and b.catdisp = '1' ;",sep="")
#put the query together to make an sql query that will run with channel 1
#qry<-paste(sel1,from1,where1,sep=" ")
#qry  
#now run it and store the result
#working<-sqlQuery(channel1,qry,max=1000000)  
  
#see how much memory objects are using...
tmp=as.data.frame(rev(sort( sapply(mget(ls()),object.size)))[1:20])
names(tmp)=c("Size");tmp

#This converts a list to data.frame pretty nicely 
MyNewDataFrame=data.frame(matrix(unlist(Mylist), nrow=Myrows, byrow=T),stringsAsFactors=FALSE)

#use in a function to get a package that may not be loaded (or even downloaded) before the call
if(require("lme4")){
  print("lme4 is loaded correctly")
} else {
  print("trying to install lme4")
  install.packages("lme4")
  if(require(lme4)){
    print("lme4 installed and loaded")
  } else {
    stop("could not install lme4")
  }
}


} #end function header
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#adjust window size
resize.win <- function(Width=10, Height=10, xpinch=20, ypinch=30)
{
  if(length(dev.list())>0) {dev.off()}
  platform <- sessionInfo()$platform 
  if (grepl("linux",platform)) { 
	x11(width=Width, height=Height) 
  } else if (grepl("pc",platform)) { 
	windows(width=Width, height=Height) 
  } else if (grepl("apple", platform)) { 
      quartz(width=Width, height=Height) 
  } 
}
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#useful for getting packages in functions 
check.package=function(package1){
  if(suppressWarnings(require(package1,character.only=TRUE))){
    print(paste0(package1," is loaded correctly"))
  } else {
    print(paste0("trying to install ",package1))
    install.packages(package1,character.only=TRUE)
    if(require(package1,character.only=TRUE)){
      print(paste0(package1," installed and loaded"))
    } else {
      stop(paste0("could not install ", package1))
    }
  }}
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearMemory=function(){
  rm(list = ls(envir = globalenv()),envir = globalenv()) #clear Vars from global enviroment
  gc()  #grabage colector
  cat("\014") #clc
  .rs.restartR() #clear session
}
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Filter<-function(x) return(x[,colSums(is.na(x))<nrow(x)]) #drop columns that are all NA
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
frame.class<-function(x){
#prints the names and class of each column in a dataframe (x)
for(i in 1:dim(x)[2]){
  print(c(names(x[,i]),class(x[,i])) )
  }
}
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#Check the currently installed packeges and their version
check.packages=function(){
  ip <- as.data.frame(installed.packages()[,c(1,3:4)])
  rownames(ip) <- NULL
  ip <- ip[is.na(ip$Priority),1:2,drop=FALSE]
  print(ip, row.names=FALSE)
}
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# logit and inverse logit
logitx<-function(x){
  if(x>=1 | x<0) {
    return("MUST BE BETWEEN 0 & 1 FOOL!")
  } else {return(log((x/(1-x))))}
}
elogitx<-function(x){
  return((exp(x)/(1+exp(x))))
}
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
co.var<-function(x)(
       100*sd(x)/mean(x)
    )
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cv<-function(x)(
  sd(x)/mean(x)
)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
se<-function(x) {
      sd(x)/sqrt(length(x))
    }
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
na.replace<-function(x,rep1=0){
  #This function will replace NA's in
  #vector or matrix x with rep1 (0 is the default)
  if(!is.null(ncol(x))) {
    ret=apply(x,2, rep1,FUN=function(x,...) {ifelse(is.finite(x),x,rep1)})
  } else ret=ifelse(is.finite(x),x,rep1)
  
  return(ret)
}
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
closest<-function(vec,x){
  #find the index of the closest value to x in vec
  which(abs(vec-x)==min(abs(vec-x)))
}
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# a function to make bins
makebins<-function(x,binsize=10,first.bin=0) return(trunc(x[x>=first.bin]/binsize)*binsize+binsize/2)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mid<-function(x){return(((max(x)-min(x))/2+min(x)))}
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
std.lnorm<-function(x1,cv.x){
  #give a mean and cv in real space and this will give back the appropriate
  #parameters for using the R dlnorm functions
  sx1<-log(x1)
  ss1<-sqrt(log(1+cv.x^2))
  return(data.frame("meanlog"=sx1,"sdlog"=ss1))
}
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
RandLN<-function(n,x1,cv.x){
  #give a mean and cv in real space and this will give back the n random log 
  #normal deviates
  sx1<-log(x1)
  ss1<-sqrt(log(1+cv.x^2))
  rez<-rlnorm(n,sx1,ss1)
  return(rez)
}
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
logify<-function(a){
  #transform normal a into a log normal variable with bias correction
  CF<-(se(a)^2)/2
  a2<-exp(a-CF)
  return(a2)
}
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
make.bin.fixed=function(minx,maxx,brks,x) {
  # minx is the minimum value for the first bin
  # maxx is the maximum value for the last bin
  # brks is the bin interval (size of bin)
  # x is the variable you want to bin
  # this routine will return a vector of bin
  # assignments for each value in x, but the
  # levels will be fixed by minx, maxx and brks.
  # DRH 1/13/12
  hist.brks<-c(-Inf,seq(from=minx,to=maxx,by=brks),Inf)
  cuts <- cut(x, breaks=hist.brks, labels=seq(from=minx,to=(maxx+brks),by=brks))
  return(cuts)
}
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
find.row.num=function(to.match,array1,col1) {
  #find the row number corresponding to the first occurence of an item of interest
  #to.match = the item to match
  #array1 is the array that you want to match to
  #col1 is the column you expect the match to occur in...
  #example:  (using the iris data set)
  #   my.row<-find.row.num(3.6,iris,2)
  #   my.row
  row.n<-c()
  for(j in 1:dim(array1)[1]) {
     if(array1[j,col1]==to.match) {
        row.n<-j
        break
     }
  }
  return(row.n)
}
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
find.row.num.vec=function(to.match,vec1) {
  #find the row number corresponding to the first occurence of an item of interest
  #to.match = the item to match
  #vec1 is the vector that you want to match to
  row.n<-c()
  for(j in 1:length(vec1)) {
     if(vec1[j]==to.match) {
        row.n<-j
        break
     }
  }
  return(row.n)
}
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
find.row.num.all=function(to.match,array1,col1) {
  #find the row number corresponding to all occurences of an item of interest
  #to.match = the item to match
  #array1 is the array that you want to match to
  #col1 is the column you expect the match to occur in...
  #example:  (using the iris data set)
  #   my.rows<-find.row.num.all(3.6,iris,2)
  #   my.rows
  row.n<-c()
  k<-1
  for(j in 1:dim(array1)[1]) {
     if(array1[j,col1]==to.match) {
        row.n[k]<-j
        k<-k+1
     }
  }
  return(row.n)
}
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
edist<-function(slat,slon,xlat,xlon)
{
#/*-------------
# *
# * Name:        e_dist
# * Input:       two lat, lon pairs
# *              (double) slat, slon, xlat, xlon
# * Return:      distance in meters
# * Result:      distance calculated
# *
# * Note:        Adapted from FORTRAN code written by B. Unger from
# *              Andoyer-Lambert distance formula, Journal of the
# *              Washington Academy of Sciences, May 15, 1942
# *
# *              Calculates distance between two points on the globe
# *              with corrections for the GRS 80 ellipsoid.
# *
#/*-------------

  FLAT<-1.0/298.257      # flattening at poles
  R<-6378137.0           # radius of earth (at equator)
  DTOR<-0.0174532925     # decimal degrees to radians

# account for negative longitude - suppress the warning for a vector
suppressWarnings(if((slon * xlon) < 0.0)
  dlo <- (180.0 - (abs(slon))) + (180.0 - abs(xlon))
 else dlo <- abs(slon-xlon)
)
# account for negative latitude - suppress the warning for a vector
suppressWarnings( if((slat * xlat) < 0.0)
 {
    if(slat < 0.0)
    {
      slat <- slat + 90.0
    }
    else
    {
      xlat <- (-xlat + 90.0) * (1.0)
      slat <- 90.0 - slat
    }
  }
)
# translate to radians for use in trig functions
 slatr <- slat * DTOR
 slonr <- slon * DTOR
 xlatr <- xlat * DTOR
 xlonr <- xlon * DTOR
 dlor  <- dlo  * DTOR

 s <- sin(slatr) * sin(xlatr)
 c <- cos(slatr) * cos(xlatr)
 arc  <- acos(s + (c * cos(dlor)))

# corrections for ellipsoid
 t1 <- (FLAT * R)/4.0
 t2 <- ((3.0*(sin(arc))-arc)*(1.0+s+c)*(1.0+s-c))/(1.0+cos(arc))
 t3 <- ((3.0*(sin(arc))+arc)*(1.0-s+c)*(1.0-s-c))/(1.0-cos(arc))

 sig <- t1 * (t2 - t3)
 range <- (R * arc) + sig
 return(range)
}
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#*******************************************************************
meters=function(XLO,YLO,XHI,YHI){
  # This subroutine computes the difference in two positions in meters.
  # This method utilizes ellipsoidal rather than spherical
  # parameters.  I believe that the original approach and code
  # for this came from Ed McKay.
  # The reference used by Ed McKay for this was:
  #       A Course in Higher Geodesy' by P.W. Zakatov, Israel Program
  #       for Scientific Translations, Jerusalem, 1962
  #       Warren T. Dewhurst
  #       11/1/89
  # Note that this subroutine is set up for +west longitude
  AXIS1 = 6378137.0
  E2 = 0.0066943800229
  RHOSEC = 206264.806247
  #REAL*8 W, LM, LP, AVLAT
  #REAL*8 LAT1S, LAT2S, LONG1S, LONG2S
  #REAL*8 DLAT, DLONG
  LAT1=YLO
  LONG1=abs(XLO)
  LAT2=YHI
  LONG2=abs(XHI)      
  # Change into sec.ddd and convert to +west longitude
  LAT1S =    LAT1*60.*60./RHOSEC
  LONG1S = -LONG1*60.*60./RHOSEC
  LAT2S =    LAT2*60.*60./RHOSEC
  LONG2S = -LONG2*60.*60./RHOSEC
  DLAT  = ( LAT2S -  LAT1S)*RHOSEC
  DLONG = (LONG2S - LONG1S)*RHOSEC
  AVLAT = (LAT1S + LAT2S)/2.0
  W  = sqrt(1.00 - E2*sin(AVLAT)^2)
  LM = AXIS1*(1.00 - E2)/(W^3*RHOSEC)
  LP = AXIS1*cos(AVLAT)/(W*RHOSEC)
  LATMTR = LM*DLAT
  LONMTR = LP*DLONG
  DX=LONMTR
  DY=LATMTR
  return(data.frame(DX,DY))
}
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pad.zero=function(num){
  #this function just adds some zeros in front of a number <1000 to get to
  #three digits
  num<-as.numeric(as.character(num))
  if(num<10) {
    pad.num<-paste("00",num,sep="")
    } else if(num>9 & num<100) {
    pad.num<-paste("0",num,sep="")
    } else  {pad.num<-paste(num) }
return(pad.num)
}
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
decimal.degrees=function(pos) {
  #convert from degree minutes to decimal degrees in latitude or longitude
  d.pos<-c()
  for(i in 1:length(pos)) {
      #clip the minutes first
      if(pos[i]<0) {
        pos[i]<--pos[i]
        neg<-T
      } else neg<-F
      c.pos<-substr(paste(pos[i]),3,8)
      m.pos<-(as.numeric(c.pos)/60.0 )
      #clip the degrees last
      c.pos<-substr(paste(pos[i]),1,2)
      h.pos<-as.numeric(c.pos)
      #now put it all together
      d.pos[i]<-as.numeric(h.pos)+as.numeric(m.pos)
      if(neg) d.pos[i]<--d.pos[i]
    }
  return(d.pos)
}
#  lon<-decimal.degrees(d.lon[1:50])
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
degree.min.sec=function(pos){
#convert from decimal dgrees to degree minutes 
  out=c();neg=F
  for(i in 1:length(pos)){
    #get the decimal part
    if(pos[i]<0) {pos=abs(pos);neg=T}
    deg=floor(pos[i])
    min=floor(60*(pos[i]-deg))
    sec=3600*(pos[i]-deg-(min/60))
    #remove decimal point
    sec=gsub(".","",sec,fixed = T)
    #clip this to 4 places
    sec=substr(sec,1,min(4,nchar(sec)))
    if(neg) deg=deg*-1
    out[i]=paste(deg,min,sec,sep=".")
  }
  return(out)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dev.table=function(glm1,pathname,filename) {
library(xtable)
#_______________________________________________________________________________
# Here is a function to produce a deviance table for any glm object.  The
#arguments are a glm object, a pathname where you want the table to be
#written and a filename (optional).  The default pathname is "deviance_table.csv".
#Remember to add / at the end of your pathname...
#                       DRH             4/13/11
#_______________________________________________________________________________
  anova.dev1<-anova(glm1)
  names(anova.dev1)<-c("Df","Deviance","Resid.Df","Resid.Dev")
  delta_dev<-NA
  perc_dev<-NA
  sig<-""
  for(i in 2:dim(anova.dev1)[1]) {
    delta_dev[i]<-anova.dev1$Resid.Dev[i-1]-anova.dev1$Resid.Dev[i]
    perc_dev[i]<-delta_dev[i]/(anova.dev1$Resid.Dev[1]-anova.dev1$Resid.Dev[dim(anova.dev1)[1]])*100.
    if(perc_dev[i]>=10.) sig[i]="**"
    else if(perc_dev[i]>=5.0 & perc_dev[i]<10.0) sig[i]="*"    #significant
    else if(perc_dev[i]>=4.5 & perc_dev[i]<5.0) sig[i]="."   #near significant
    else    sig[i]=""     #not significant
  }
  dev1<-data.frame(rownames(anova.dev1),anova.dev1$Resid.Dev,delta_dev,perc_dev,sig)
  names(dev1)<-c("Model","Residual deviance","Delta deviance","% of deviance explained","significance")
  if(filename=="") filename<-"deviance_table.csv"
  write.csv(dev1,file=paste(pathname,filename,sep=""),quote=F,na="",row.names=F )
  #also produce an x table version of this for use in latex
  write(print(xtable(dev1,caption=paste("Deviance table for ...."))),
    file=paste(pathname,substr(filename,1,(nchar(filename)-3)),"tex",sep=""))
  }
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
read.sas.data<-function(type1){
#_______________________________________________________________________________
# type1= "reg" or "tow" or "age"
# "reg" is the regional file which reports the stats of interest by a user
#defined (in the sas code) region.  "tow" gives a tow by tow break down.
#_______________________________________________________________________________
#Now read send this to the R compiler from windows to load your data frame.
  load("B:/ClamDataBase/SASdataFiles/tempfile.Rdata")
  objects()



#There are 2 types of files one might get: regional data or tow by tow...
  if(type1=="reg") {
      load("B:/ClamDataBase/SASdataFiles/tempfile.Rdata")
      ############################################################################
      wrk.t.data<-reg.data$region_data_last_run_403
      #wrk.t.data<-reg.data$region_data_last_run_409

      names(wrk.t.data)
      #select particular columns as follows:
      table.out<-data.frame(wrk.t.data[,c(4,1,2,3,5,6,9,10,11,12,13,14,15,16
              ,21,22,26,27,28)])
      #[1] "REGION"   "CRUISE6"  "LENGRP"   "SVSPP"    "REGNAM"   "SURVAREA"
      #[7] "LWSARCYR" "REV_DATE"
      #[9] "KGPERTOW" "SE_KGPER" "NPERTOW"  "SE_NPER"  "STOCKNPE" "SE_STK_N"
      #[15] "STOCKKGP" "SE_STK_K"
      #[17] "FISHABLE" "SE_FIS_N" "FISHABL2" "SE_FIS_K" "NTOWS"    "POSTOW"
      #[23] "BORROWED" "BORROW_S"
      #[25] "STRATWT"  "STRATSUR" "CVN"      "CVKG"     "CVSTKN"   "CVSTKKG"
      #[31] "CVFISHN"  "CVFISHKG"
      #[33] "REGAREA"  "NREGSTRA" "TOTALN"   "TOTALKG"  "TOTALSTK" "TOTALST2"
      #[39] "TOTALFIS" "TOTALFI2"
      #[41] "PERC_N"   "PERC_KG"  "PERCSTOC" "PERCSTO2" "PERCFISH" "PERCFIS2"
      names(table.out)
      dim(table.out)
      table.out[1:50,]
      year<-round(table.out$CRUISE6/100,0)
      x.data<-data.frame(table.out,year)
  }
  if(type1=="tow") {
      load("B:/ClamDataBase/SASdataFiles/tempfile.Rdata")
      ##########################################################################
      wrk.t.data<-tow.data$tow_data_last_run_403
      names(wrk.t.data)
      #select particular columns as follows:
      table.out<-data.frame(wrk.t.data[,c(1,2,3,4,5,7,11,12,13,14,23,24,25,26,27)])
      # [1] "REGION"   "CRUISE6"  "STATION"  "LENGRP"   "STRATUM"
      # [6] "TOW"      "REGNAM"   "SVSPP"    "DOPDISTB" "DIST2USE"
      #[11] "LAT"      "LON"      "DEPTH"    "TNMS"     "REV_DATE"
      #[16] "JWSTCODE" "RANDLIKE" "STATYPE"  "HAUL"     "GEARCOND"
      #[21] "LWSARCYR" "DISTANCE" "NUMLEN"   "NPERTOW"  "KGPERTOW"
      #[26] "STOCKNPE" "STOCKKGP" "FISHABLE" "FISHABL2" "TOTALN"
      #[31] "TOTALKG"  "TOTASTOC" "TOTALSTO" "TOTALFIS" "TOTALFI2"
      #[36] "PERC_N"   "PERC_KG"  "PERCSTOC" "PERCSTO2" "PERCFISH"
      #[41] "PERCFIS2" "POSTOW"   "BORROWED" "ORIG_CRU" "ORIG_TOW"
      #[46] "ORIG_STA"
      names(table.out)
      dim(table.out)
      table.out[1:50,]
      year<-round(table.out$CRUISE6/100,0)
      x.data<-data.frame(table.out,year)
    }
    if(type1=="age") {
      ############################################################################
      load("B:/ClamDataBase/Ages/tempfile.Rdata")
      wrk.t.data<-age.data1
      #print(names(wrk.t.data))
      #select particular columns as follows:
      table.out<-data.frame( wrk.t.data[,c(7,3,4,5,12,16)]  )
      #[1] "CRUISE"   "STATION"  "YEAR"     "LENGTH"   "AGE"      "CRUISE6"
      #[7] "STRATUM"  "TOW"      "SVGEAR"   "MONTH"    "DAY"      "AVGDEPTH"
      #[13] "SVSPP"    "REGION"   "REV_DATE" "REGNAM"
      print(names(table.out))
      x.data<-data.frame(table.out)
    }

return(x.data)
}
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
check.lim<-function(x,min){
  #check to see that there are at least min finite values in x.
  #if yes return TRUE, otherwise FALSE
  checksum<-ifelse(is.finite(x),1,0)
  if(sum(checksum)>min) {
    return(TRUE)
    } else {return(FALSE)}
}
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
blank.plot<-function(string){
  #plot nothing but string in the middle of the plotting region - a useful
  #placeholder for automated plotting loops
  plot(x=c(0,1),y=c(0,1)
    ,main=""
    ,axes=F
    ,ylab=""
    ,xlab=""
    ,col="white"
  )
  text(x=0.5,y=0.5
    ,labels=string
    ,cex=2.
  )
}
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
m.to.psi<-function(x) {return((x)*1.46)}
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
increment=function(x,inc=1) {return(x+inc)}
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
last <- function(x) { tail(x, n = 1) }
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
m.to.nm<-function(x) {return((x)* 0.00053995680346)}
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nm.to.m<-function(x) {return((x)* (1/0.00053995680346))}
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
m2.to.nm2<-function(x) {return((x)* 0.00053995680346^2)} #2.915533496e-7
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nm.to.ft<-function(x) {return((x)* 6076.1154856)}
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#1 square mile [nautical] = 3 429 904 square meter
nm2.to.m2=function(x) {return((x)* 3429904)}
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
get.region<-function(strata1){
  #A function to assign a region name to a stata number, using lookup tables
  #which can be adjusted for changes to the region/strata relationship over time
  #assign strata to region
  dir.region.lookup <- paste("/home/dhennen/"
                  ,"DRH_general_R_code/Region Lookup/",sep="")
  lookup<-read.csv(paste(dir.region.lookup,"lookup_strata_region.csv",sep="")
        ,stringsAsFactors=F
        ,as.is=T
      )
  names(lookup)
  if(mean(strata1)>1000) REG.CODE<-lookup$REGION[match(((strata1-6000)/10),lookup$STRATA)]
  if(mean(strata1)<100) REG.CODE<-lookup$REGION[match(strata1,lookup$STRATA)]
  rm(lookup)
  #now get the region name
  lookup<-read.csv(paste(dir.region.lookup,"lookup_region_name.csv",sep="")
        ,stringsAsFactors=F
        ,as.is=T
  )
  names(lookup)
  Region<-lookup$REGNAM[match(REG.CODE,lookup$REGION)]
  rm(lookup)
  return(Region)
}
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# to make a gray scale
#*******************************************************************************
#        col=gray(seq(0.1,0.9,length=10)))
#*******************************************************************************
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
error.bar <- function(x, y, upper, lower=upper, length=0.1,col1="black",...){
    if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
    stop("vectors must be same length")
    arrows(x
      ,y+upper
      ,x
      ,y-lower
      ,angle=90
      ,code=3
      ,length=length
      ,col=col1
      , ...)
}
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
asm.ci2<-function(x,cv.x){
  #generate approximate 95% confidence intervals based on log normal variable x
  #and it's cv. cv.x
  s<-sqrt(log(1+cv.x^2))
  s<-ifelse(is.finite(s),s,0)
  lci<-x*exp(-1.96*(s))
  uci<-x*exp(1.96*(s))
  return(data.frame("lci"=lci,"uci"=uci)) 
}
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#Both of these confidence intervals are set to 95% - change 'em if you want to!
asm.ci<-function(x,cv.x,bounds=95){
  #generate approximate 95% confidence intervals based on log normal variable x
  #and it's cv. cv.x
  s<-sqrt(log(1+cv.x^2))
  s<-ifelse(is.finite(s),s,0)
  p<-(1-(bounds/100))/2
  Z<-qnorm(p)
  lci<-x*exp(Z*(s))
  uci<-x*exp(-Z*(s))
  return(data.frame("lci"=lci,"uci"=uci)) 
}
gamma.ci<-function(ts,cv,bounds=95){
  #ts is your data and cv its cv
  #lp and up are the lower and upper probabilities for your gamma CI's
  cv2<-(1/cv)^2
  scale=ts/cv2
  shape=ts/scale
  lp<-(1-(bounds/100))/2
  up<-1-lp
  LB=qgamma(p=lp,shape=shape,scale=scale)
  UB=qgamma(p=up,shape=shape,scale=scale)
  return(data.frame("lci"=LB,"uci"=UB))
}
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fix.NJ<-function(x) {
  #NJ has a region code of 8, when it should be 3 or 4 to make the regions line up
  #S to N
  x$REGION<-with(x,ifelse((REGNAM=="NJ" | REGION==8),3,REGION))
  return(x)
}  
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot.region.trend.ci<-function(
#plot a panel of regional trends - does not plot "other"
#This should be useful in the build a bridge phase of assessment development.
#This application is specific to Larry's .sas database output circa 2012
#requires these variables:
  plot.var=NULL               #what variable should be plotted
  ,plot.var.lci=NULL           #confidence region of that plotting variable
  ,plot.var.uci=NULL
  ,plot.var.name=NULL         #name for y axis label
  ,reg.data=NULL             #data frame from Larry's .sas database
  ,path.name=NULL             #where should the data be saved
  ,split.yr=NULL             #Year to split the survey time series (end year of the first series)
  ) {
    plot.new()
    resize.win(Width=15, Height=10, xpinch=30, ypinch=30)    
    par(oma=c(1,2,0,0))
    #exclude "other"
    reg.data<-fix.NJ(reg.data)    
    reg.data2<-data.frame(reg.data,plot.var,plot.var.lci,plot.var.uci)
    reg.data2<-reg.data2[reg.data2$REGNAM!="OTH",,drop=T]
    if(length(unique(reg.data2$REGNAM))>2) reg.data2<-OrderRegionsSC(reg.data2)
    reg.lev<-levels(as.factor(reg.data2$REGNAM))[levels(as.factor(reg.data2$REGNAM))!="OTH",drop=T]
    if(length(reg.lev)>2) {
      par(mar=c(4,5,4,2)+0.1,cex.axis=1.5,cex.lab=2,cex.main=2) #bottom, left, top, right      
      l.out<-matrix(c(1:6),3,2, byrow=T)
    } 
    if(length(reg.lev)==2) {
      par(mar=c(4,5,4,2)+0.1,cex.axis=1.5,cex.lab=2,cex.main=2) #bottom, left, top, right      
      l.out<-matrix(c(1:2),2,1, byrow=T)
    } 
    if(length(reg.lev)==1) {
      par(mar=c(5,5,4,2)+0.1,cex.axis=1.5,cex.lab=2,cex.main=2) #bottom, left, top, right      
      l.out<-matrix(c(1:1),1,1, byrow=T)
    }
      
    nf <- layout(l.out)    
    
  for(i in 1:length(reg.lev)){
      t.dat<-c()
      t.dat<-subset(reg.data2,REGNAM==reg.lev[i])
      #find a scale
      yl<-c(0,(max(na.omit(t.dat$plot.var.uci))+.1*max(na.omit(t.dat$plot.var.uci))))
      m1<-paste(t.dat$REGNAM[1])
      #See if a survey split is required
      if(!is.null(split.yr)) {
        pnt=which(t.dat$YR==split.yr)
        YR1=t.dat$YR[1:pnt]; YR2=t.dat$YR[(pnt+1):length(t.dat$YR)];
        plot.var1=t.dat$plot.var[1:pnt]; plot.var2=t.dat$plot.var[(pnt+1):length(t.dat$plot.var)];
        plot.var1.lci=t.dat$plot.var.lci[1:pnt]; plot.var2.lci=t.dat$plot.var.lci[(pnt+1):length(t.dat$plot.var.lci)];
        plot.var1.uci=t.dat$plot.var.uci[1:pnt]; plot.var2.uci=t.dat$plot.var.uci[(pnt+1):length(t.dat$plot.var.uci)];
      } else {
        YR1=t.dat$YR; YR2=NULL;
        plot.var1=t.dat$plot.var; plot.var2=NULL;
        plot.var1.lci=t.dat$plot.var.lci; plot.var2.lci=NULL;
        plot.var1.uci=t.dat$plot.var.uci; plot.var2.uci=NULL;       
      }
      
      plot(plot.var1~YR1
        ,type="b"
        ,main=m1
        ,ylab=""
        ,xlab=""
        ,pch=20
        ,ylim=yl
        ,xlim=c(min(pretty(c(min(t.dat$YR),max(t.dat$YR)))),max(pretty(c(min(t.dat$YR),max(t.dat$YR)))))
      )
      lines(plot.var1.lci~YR1
        ,lty=2
      )
      lines(plot.var1.uci~YR1
        ,lty=2
      )      
      if(!is.null(split.yr)){
        points(plot.var2~YR2,pch=20,cex=1.5,col="red")
        if(length(YR2)==1) {
          error.bar(x=YR2,y=plot.var2,upper=(plot.var2.uci-plot.var2),lower=(plot.var2-plot.var2.lci),col="red")
        }  else {
          lines(plot.var2~YR2,lty=1,col="red")
          lines(plot.var2.lci~YR2,lty=2,col="red") 
          lines(plot.var2.uci~YR2,lty=2,col="red")           
        }
      }
    }  
    
    mtext(text="Year"
      ,side=1
      ,outer = T
      ,line=-1
      ,cex=2
    )
    mtext(text=plot.var.name
      ,side=2
      ,outer = T
      ,line=-2
      ,cex=2
    )
    filename<-path.name
    savePlot(filename=paste(filename,"png",sep="."),type="png")
}
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
compare.region.trend<-function(
#plot a panel of regional trends, comparing 2 versions - does not plot "other"
#This should be useful in the build a bridge phase of assessment development.
#This application is specific to Larry's .sas database output circa 2012
#requires these variables:
  plot.var               #what variable should be plotted
  ,plot.var2
  ,plot.var.se           #std. error of that plotting variable
  ,plot.var.se2
  ,plot.var.name         #name for y axis label
  ,reg.data
  ,reg.data.2               #data frames from Larry's .sas database
  ,path.name             #where should the data be saved
  ) {
    plot.new()
    resize.win(Width=20, Height=15, xpinch=30, ypinch=30)
    par(mar=c(4,5,4,2)+0.1,cex.axis=1.5,cex.lab=2,cex.main=2) #bottom, left, top, right
    l.out<-matrix(c(1:6),3,2, byrow=T)
    nf <- layout(l.out)

    #find a scale
    y1<-max(plot.var)+.1*max(plot.var)
    y2<-max(plot.var2)+.1*max(plot.var)
    yl<-c(0,max(y1,y2))


    reg.data2<-data.frame(reg.data,plot.var,plot.var.se)
    reg.data2.2<-data.frame(reg.data.2,plot.var2,plot.var.se2)

    #exclude "other"
    reg.lev<-levels(as.factor(reg.data$REGION))

    for(i in 2:length(reg.lev)){
      t.dat<-c()
      t.dat<-subset(reg.data2,REGION==reg.lev[i])
      t.dat2<-c()
      t.dat2<-subset(reg.data2.2,REGION==reg.lev[i])
      m1<-paste(t.dat$REGNAM[1])
      plot(plot.var~YR
        ,data=t.dat
        ,type="b"
        ,main=m1
        ,ylab=""
        ,xlab=""
        ,ylim=yl
      )
      with(t.dat,error.bar(x=YR
          ,y=plot.var
          ,upper=plot.var.se
          ,length=0.05
        )
      )
      points(plot.var2~YR
        ,data=t.dat2
        ,type="b"
        ,lty="dashed"
        )
      with(t.dat2,error.bar(x=YR
          ,y=plot.var2
          ,upper=plot.var.se2
          ,length=0.05
        )
      )
      legend( x  = "topright"
        ,legend=c("2009","2012")
        ,lty=c(1,2)
      )
    }

    mtext(text="Year"
      ,side=1
      ,outer = T
      ,line=-2
      ,cex=2
    )
    mtext(text=plot.var.name
      ,side=2
      ,outer = T
      ,line=-2
      ,cex=2
    )
    filename<-path.name
    savePlot(filename=paste(filename,"pdf",sep="."),type="pdf")
    savePlot(filename=paste(filename,"tif",sep="."),type="tif")
}
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot.clam.map.slow<-function(
  shell.strata=T
  ,state.boundary=T
) {
#a function that plots a background map onto which the user can put whatever.
#The map is of the Northeast coast of the US, and potionally, the 3 mile
#state boundary and the shellfish strata
#Load required packages
  library(gpclib)
  library(maptools)
  gpclibPermit()
  #get some shape files
  dir.maps<-"/home/dhennen/Clam/ARC_data/"
  fname<-"shstr"
  strata1<-readShapeLines(paste(dir.maps,fname,sep=""))
  fname<-"dcw"
  coast<-readShapePoly(paste(dir.maps,fname,sep=""))
  fname<-"3mileME-SC"
  state.h2o<-readShapeLines(paste(dir.maps,fname,sep=""))
  #first plot the coast line
  plot(coast
    ,col="grey"
    ,bg="white"
    ,ylim=c(35,43)
    ,xlim=c(-76.5,-65.5)
    ,axes=T
  )
  #lines(tnms.lines)
  if(state.boundary) lines(state.h2o,col="red")
  if(shell.strata) {
    lines(strata1
      ,col="dark grey"
    )
  }
}
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot.clam.map<-function(map=1,axes=F,SC=T,...){
  #Much faster way to do this
  require(RColorBrewer)
  options(stringsAsFactors = F)
  if(map==1){  
    #only get the map data if you need to
    if(!exists("clamap")){
      library(rworldmap)
      assign("clamap",getMap(resolution = "low"), envir = .GlobalEnv)
    }
    plot(clamap, xlim = c(-76, -66), ylim = c(35.5, 42.5),col="gray80",axes=axes,...)
    #only read the shape file for strata in if you haven't been here before!
    if(!exists("strata")) {
      library(gpclib)
      library(maptools)
      gpclibPermit()
      #dir.maps<-"/home/dhennen/Clam/ARC_data/"
      #fname<-"shstr"
      #assign("strata1",readShapeLines(paste(dir.maps,fname,sep=""))
      #  , envir = .GlobalEnv) 
      CurrentStrata.dir="/home/dhennen/Clam/Assessment Stuff/Survey Analysis/StrataMadeIn2018/"
      if(SC) load(paste0(CurrentStrata.dir,"SCstrata.RData")) #objects is called "tosave"
      if(!SC) load(paste0(CurrentStrata.dir,"OQstrata.RData")) #objects is called "tosave"
      strata=tosave #transfer the object you just loaded into a usefully named memory location
      #assign("strata1",strata, envir = .GlobalEnv) 
    }
    
    col=brewer.pal(12,"Set3")
    col2=data.frame("Strata"=strata@data$Group.1,"col"=
                      unlist(lapply(col[1:length(strata@data$Group.1)],FUN=function(x) transparent(x,20))))
    col3=data.frame("Strata"=strata@data$Group.1,"col"=col)
    #                  unlist(lapply(col1[1:length(strata@data$Group.1)],FUN=function(x) transparent(x,90))))
    plot(strata,add=T,col=col2$col) #plot the strata polygons

  }
  
  if(map==2) {
    name1<-"map_base"
    require(png)
    dir.maps<-"/home/dhennen/Clam/ARC_data/"
    map1 <- readPNG(paste(dir.maps,name1,".png",sep=""))
    par1 = (read.table(file=paste(dir.maps,name1,".par",sep="")))
    #Get the plot information so the image will fill the plot box, and draw it  
    #Set up the plot area 
    plot(y=c(par1[3,1],par1[4,1]),x=c(par1[1,1],par1[2,1])
      , type='n', main="", xlab="", ylab="",axes=F)  
    #Get the plot information so the image will fill the plot box, and draw it  
    rasterImage(map1, par1[1,1], par1[3,1], par1[2,1], par1[4,1],interpolate=F)
  }
  if(map==3) {
    library (ggmap)
    mapImageData <- get_googlemap(center= c(lon=-72.0, lat=40.0)
      ,zoom=7
      ,size=c(640,640)  
      ,scale=2  
      ,maptype="satellite"  
    )
    ggmap(mapImageData)      
  }
  
}
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
generate.plot.clam.map<-function(){
  #Make sure to change the name if you are making a new background!  Then add it to the 
  #function above with a new code!
  name1<-"map_base"
  library(gpclib)
  library(maptools)
  gpclibPermit()
  #get some shape files
  par(mar=c(0,0,0,0),cex.axis=1.5,cex.lab=2,cex.main=2) #bottom, left, top, right

  dir.maps<-"/home/dhennen/Clam/ARC_data/"
  fname<-"shstr"
  strata1<-readShapeLines(paste(dir.maps,fname,sep=""))
  fname<-"dcw"
  coast<-readShapePoly(paste(dir.maps,fname,sep=""))
  fname<-"3mileME-SC"
  state.h2o<-readShapeLines(paste(dir.maps,fname,sep=""))
  #first plot the coast line
  plot(coast
    ,col="grey"
    ,bg="white"
    ,ylim=c(35,43)
    ,xlim=c(-76.5,-65.5)
    ,axes=F
  )
  if(T) {
    lines(strata1
      ,col="dark grey"
    )
  }
  savePlot(filename=paste(filename,name1,"png",sep="."),type="png")
  p1<-par()
  write.table(p1$usr,file=paste(dir.maps,name1,"par",sep="."))
  
  #This is a way to save this map, but it is fast so no real need to do so!
  mapImageData <- get_googlemap(center= c(lon=-72.0, lat=40.0)
    ,zoom=7
    ,size=c(640,640)  
    ,scale=2  
    ,maptype="satellite"  
  )
  ggmap(mapImageData)
  save(mapImageData, file=paste(dir.maps,"savedMap.rda",sep=""))

  
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#/**
# \name{weighted.median}
# \alias{weighted.median}
#
# \title{Weighted Median Value}
#
# \usage{weighted.median(x, w, na.rm=TRUE, ties=NULL)}
#
# \description{
#   Compute a weighted median of a numeric vector.
# }
#
# \arguments{
#   \item{x}{a numeric vector containing the values whose weighted median is
#            to be computed.}
#   \item{w}{a vector of weights the same length as \code{x} giving the
#weights
#            to use for each element of \code{x}. Default value is equal
#weight
#            to all values.}
#   \item{na.rm}{a logical value indicating whether \code{NA} values in
#            \code{x} should be stripped before the computation proceeds.}
#   \item{ties}{a character string specifying how to solve ties between two
#            \code{x}'s that are satisfying the weighted median criteria.
#            Note that at most two values can satisfy the criteria.
#            When \code{ties} is \code{"min"}, the smaller value of the two
#            is returned and when it is \code{"max"}, the larger value is
#            returned.
#            If \code{ties} is \code{"mean"}, the mean of the two values is
#            returned and if it is \code{"both"}, both values are returned.
#            Finally, if \code{ties} is \code{"weighted"} (or \code{NULL}) a
#            weighted average of the two are returned, where the weights are
#            weights of all values \code{x[i] <= x[k]} and \code{x[i] >= x[k]},
#            respectively. Default value is \code{NULL}.}
# }
#
# \value{
#   Returns the weighted median.
# }
#
# \details{
#  For the \code{n} elements \code{x = c(x[1], x[2], ..., x[n])} with positive
#  weights \code{w = c(w[1], w[2], ..., w[n])} such that \code{sum(w) = S},
#  the \emph{weighted median} is defined as the element \code{x[k]} for which
#initial  the total weight of all elements \code{x[i] < x[k]} is less or equal to
#  \code{S/2} and for which the total weight of all elements
#  \code{x[i] > x[k]} is less or equal to \code{S/2} (c.f. [1]).
#
#  If \code{w} is missing then all elements of \code{x} are given the same
#  positive weight. If all weights are zero, \code{NA} is returned.
#
#  When all the weights are the same and \code{ties} is \code{"weighted"}
#(or \code{NULL}) \code{weighted.median} gives the same result as \code{median}.
#
#  If one or more weights are \code{Inf}, it is the same as these weights
#  have the same weight and the others has zero. This makes things easier
#  for cases where the weights are result of a division with zero.
#
#  The weighted median solves the following optimization problem:
#
#  \deqn{\alpha^* = \arg_\alpha \min \sum_{k=1}{K} w_k |x_k-\alpha|}
#  where \eqn{x=(x_1,x_2,\ldots,x_K)} are scalars and
#  \eqn{w=(w_1,w_2,\ldots,w_K)} are the corresponding "weights" for
#  each individual \eqn{x} value.
# }
#
# \examples{
#   x <- 1:10
#   n <- length(x)
#   median(x)                            # 5.5
#   weighted.median(x)                   # 5.5
#   w <- rep(1, n)
#   weighted.median(x, w)                # 5.5 (default)
#   weighted.median(x, ties="weighted")  # 5.5 (default)
#   weighted.median(x, ties="min")       # 5
#   weighted.median(x, ties="max")       # 6
#
#   # Pull the median towards zero
#   w[1] <- 5
#   weighted.median(x, w)                # 3.5
#   y <- c(rep(0,w[1]), x[-1])           # Only possible for integer weights
#   median(y)                            # 3.5
#
#   # Put even more weight on the zero
#   w[1] <- 8.5
#   weighted.median(x, w)                # 2
#
#   # All weight on the first value
#   w[1] <- Inf
#   weighted.median(x, w)                # 1
#
#   # All weight on the first value
#   w[1] <- 1
#   w[n] <- Inf
#   weighted.median(x, w)                # 10
#
#   # All weights set to zero
#   w <- rep(0, n)
#   weighted.median(x, w)                # NA
# }
#
# \seealso{
#   \code{\link[base]{median}}, \code{\link[base]{mean}} and
#   \code{\link[base]{weighted.mean}}
# }
#
# \reference{
#   [1]  T.H. Cormen, C.E. Leiserson, R.L. Rivest, Introduction to Algorithms,
#        The MIT Press, Massachusetts Institute of Technology, 1989.
# }
#
# \author{Henrik Bengtsson, hb at maths.lth.se with help from
#         Roger Koenker, roger at ysidro.econ.uiuc.edu}
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
weighted.median <- function(x, w, na.rm=TRUE, ties=NULL) {
  if (missing(w))
    w <- rep(1, length(x));

  # Remove values that are NA's
  if (na.rm == TRUE) {
    keep <- !(is.na(x) | is.na(w));
    x <- x[keep];
    w <- w[keep];
  } else if (any(is.na(x)))
    return(NA);

  # Assert that the weights are all non-negative.
  if (any(w < 0))
    stop("Some of the weights are negative; one can only have positive
weights.");

  # Remove values with weight zero. This will:
  #  1) take care of the case when all weights are zero,
  #  2) make sure that possible tied values are next to each others, and
  #  3) it will most likely speed up the sorting.
  n <- length(w);
  keep <- (w > 0);
  nkeep <- sum(keep);
  if (nkeep < n) {
    x <- x[keep];
    w <- w[keep];
    n <- nkeep;
  }

  # Are any weights Inf? Then treat them with equal weight and all others
  # with weight zero.
  wInfs <- is.infinite(w);
  if (any(wInfs)) {
    x <- x[wInfs];
    n <- length(x);
    w <- rep(1, n);
  }

  # Are there any values left to calculate the weighted median of?
  if (n == 0)
    return(NA);

  # Order the values and order the weights accordingly
  ord <- order(x);
  x <- x[ord];
  w <- w[ord];

  wcum <- cumsum(w);
  wsum <- wcum[n];
  wmid <- wsum / 2;

  # Find the position where the sum of the weights of the elements such that
  # x[i] < x[k] is less or equal than half the sum of all weights.
  # (these two lines could probably be optimized for speed).
  lows <- (wcum <= wmid);
  k  <- sum(lows);

  # Two special cases where all the weight are at the first or the
  # last value:
  if (k == 0) return(x[1]);
  if (k == n) return(x[n]);

  # At this point we know that:
  #  1) at most half the total weight is in the set x[1:k],
  #  2) that the set x[(k+2):n] contains less than half the total weight
  # The question is whether x[(k+1):n] contains *more* than
  # half the total weight (try x=c(1,2,3), w=c(1,1,1)). If it is then
  # we can be sure that x[k+1] is the weighted median we are looking
  # for, otherwise it is any function of x[k:(k+1)].

  wlow  <- wcum[k];    # the weight of x[1:k]
  whigh <- wsum - wlow;  # the weight of x[(k+1):n]
  if (whigh > wmid)
    return(x[k+1]);

  if (is.null(ties) || ties == "weighted") {  # Default!
    (wlow*x[k] + whigh*x[k+1]) / wsum;
  } else if (ties == "max") {
    x[k+1];
  } else if (ties == "min") {
    x[k];
  } else if (ties == "mean") {
    (x[k]+x[k+1])/2;
  } else if (ties == "both") {
    c(x[k], x[k+1]);
  }
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
###
# HISTORY:
# 2002-02-07
# * Optimized the code for speed.
# * Added support for zero and Inf weights.
# * Added the 'ties' argument.
# * Created!
# * Thanks to the following people for helping me out with this one:
#   - David Brahm, brahm at alum.mit.edu
#   - David Eppstein, eppstein at ics.uci.edu
#   - Frank E Harrell Jr, fharrell at virginia.edu
#   - Markus Jantti, markus.jantti at iki.fi
#   - Roger Koenker, roger at ysidro.econ.uiuc.edu
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
max2=function(x){max(x[!is.na(x)])}
min2=function(x){min(x[!is.na(x)])}
mean2<-function(x){return(mean(na.omit(x)))}
median2<-function(x){return(median(na.omit(x)))}
sd2<-function(x){return(sd(na.omit(x)))}
sum2=function(x){sum(x[!is.na(x)])}
standardize<-function(x){return((na.omit(x)-mean(na.omit(x)))/sd(na.omit(x)))}
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
get.set.up.tows<-function() {
  load("X:/Depletion/SCSetupTowList.Rdata")
  SC.setups<-SCSetupTowList
  return(SC.setups)
}
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
geomean<-function(x){exp(mean(log(x)))  }
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
geomean2 <- function(data) {
    log_data <- log(data)
    gm <- exp(mean(log_data[is.finite(log_data)]))  #to deal with 0's
    return(gm)
}
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Hmean=function(x){1/mean(1/x)}
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Hmean2=function(x){1/mean(1/x[is.finite(1/x)])} #deal with 0's
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
var.ratio<-function(x1,x2){
  #Calculate the variance of ratio x1/x2; where x1 and x2 are vectors of random
  #variables (x2>0)
  v.x1<-var(x1)
  v.x2<-var(x2)
  E.x1<-mean(x1)
  E.x2<-mean(x2)
  E2.x1<-mean(x1)^2
  E2.x2<-mean(x2)^2
  cov.x1x2<-cov(x1,x2)
  var.x1x2<-((E.x1/E.x2)^2)*((v.x1/E2.x1)+(v.x2/E2.x2)-((2*cov.x1x2)/(E.x1*E.x2)))
  return(var.x1x2)
}
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
wt.var.ratio<-function(x1,x2,w){
  #Calculate the variance of ratio x1/x2; where x1 and x2 are vectors of random
  #variables (x2>0) with weights: w
  v.x1<-weighted.var(x1,w)
  v.x2<-weighted.var(x2,w)
  E.x1<-weighted.mean(x1,w)
  E.x2<-weighted.mean(x2,w)
  E2.x1<-mean(x1)^2
  E2.x2<-mean(x2)^2
  #print(data.frame(x1,x2))
  #print(w)
  cov.x1x2<-cov.wt(data.frame(standardize(x1),standardize(x2)),w)
  #print(cov.x1x2)
  var.x1x2<-((E.x1/E.x2)^2)*((v.x1/E2.x1)+(v.x2/E2.x2)-((2*cov.x1x2$cov[1,2])/(E.x1*E.x2)))
  return(var.x1x2)
}
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#Calculate the variance of a ratio given the means, cvs and correlation between 
var.ratio.stats<-function(x1,x2,cv1,cv2,cor){
  #Calculate the variance of ratio x1/x2; where x1 and x2 are vectors of random
  #variables (x2>0)
  v.x1<-(cv1*x1)^2
  v.x2<-(cv2*x2)^2
  E.x1<-x1
  E.x2<-x2
  E2.x1<-x1^2
  E2.x2<-x2^2
  cov.x1x2<-cor*(sqrt(v.x1)*sqrt(v.x2)) #covariance is equal to the correlation times the product of the std devs
  var.x1x2<-((E.x1/E.x2)^2)*((v.x1/E2.x1)+(v.x2/E2.x2)-((2*cov.x1x2)/(E.x1*E.x2)))
  return(var.x1x2)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
weighted.var <- function(x, w, na.rm = TRUE) {
    if (na.rm) {
        w <- w[i <- !is.na(x)]
        x <- x[i]
    }
    sum.w <- sum(w)
    sum.w2 <- sum(w^2)
    mean.w <- sum(x * w) / sum(w)
    (sum.w / (sum.w^2 - sum.w2)) * sum(w * (x - mean.w)^2, na.rm =na.rm)
}
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
stack.vec<-function(y.range,x,data1){
#transpose and stack a long dataframe of time series like data for analysis
#y.range is the column numbers that constitute the vector of interest
#x is the dependent predictor e.g. year in a time series
#data1 is the dataframe from which y.range is to be drawn
for(i in 1:dim(data1)[1]){
    #transpose the selected portion of data1
    dt.t<-t(data1[i,y.range])
    #attach dependent variable
    if (i==1) {
        dtt<-data.frame("Y"=dt.t[,1],"X"=x)
      } else{
        dtt<-rbind(dtt,data.frame("Y"=dt.t[,1],"X"=x))
      }
  }
  return(dtt)
}
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#Remove any leading zeros from alphanumeric vector d
rm.lead0=function(d) gsub("(^|[^0-9])0+", "\\1", d, perl = TRUE)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot.convolve <- function(
  #_____________________________________________________________________________
  ##############################################################################
  #Plot and convolve two overlapping distributions  
  a    #these must be vectors of random deviates from a distribution
  ,b   #using the r functions like rnorm() works well...
  ,aname=NULL  #name for the a distribution - will go in legend
  ,bname=NULL  #name for the b distribution - will go in legend
  ,colors=c("blue","red","purple")  #colors for each distribution and the 
  #overalp between them
  ,b1=1000 #the number of breaks in the plot of the historgram of a and b
  ,xlim=NULL #predefined limits (if you want to restrict your plot)
  ,ylim=NULL
  ,xlab=NULL #other plotting info
  ,ylab=NULL
  ,main=NULL
  ,put.legend="topright" #where should the legend go?
  ,plot.overlap=T #should the probability of overlap be plotted?
  ) {                       #DRH 2/5/13
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # 
  #             #EXAMPLE# 
  #             n<-1000000
  #             xA<-rgamma(n,shape=500)
  #             xB<-rlnorm(n,5.9189,0.1668)
  #             plot.convolve(a=xA,b=xB,aname="Biomass in 2011",bname="Threshold"
  #             ,xlab="Biomass",ylab="Probability density")
  #_____________________________________________________________________________
  #clear arrays and make hist objects
  ahist=NULL
  bhist=NULL
  ahist=hist(a,breaks=b1,plot=F)
  bhist=hist(b,breaks=b1,plot=F)
  #some jerking around to ge the histograms on the same scale...
  dist1 = ahist$breaks[2]-ahist$breaks[1] 
  br1 = seq((min(ahist$breaks,bhist$breaks)-dist1)
      ,(max(ahist$breaks,bhist$breaks)+dist1),dist1)
  bhist=hist(b,breaks=br1,plot=F)
  ahist=hist(a,breaks=br1,plot=F)

  #rescaling the density portion of the histogram to force sum(density)=1
  ahist$density<-ahist$density/sum(ahist$density)
  bhist$density<-bhist$density/sum(bhist$density)
    
  if(is.null(xlim)){
    xlim = c(min(ahist$breaks,bhist$breaks),max(ahist$breaks,bhist$breaks))
  }
  if(is.null(ylim)){
    ylim = c(0,max(ahist$density,bhist$density))
  }
  overlap = ahist
  #the area of overlap is easily defined by taking advantage of the hist properties
  #Note: the approximation to true overlap (based on integrals) degrades as b1 gets small
  #200 seems more or less adequate for most functions, but there is no reason to skimp! 
  for(i in 1:length(overlap$density)){
    if(ahist$density[i] > 0 & bhist$density[i] > 0){
      overlap$density[i] = min(ahist$density[i],bhist$density[i])
    } else {
      overlap$density[i] = 0
    }
  }
  if(is.null(main) & plot.overlap) p1<-paste("P[overlap]~",round(sum(overlap$density),3),sep="")
  if(!(is.null(main))) p1<-main
  print(paste("P[overlap]~",round(sum(overlap$density),3),sep=""))
  par(cex.main=3)
  plot(ahist
    ,freq=F
    ,xlim=xlim
    ,ylim=ylim
    ,col=colors[1]
    ,xlab=xlab
    ,ylab=ylab
    ,main=p1
    ,border="transparent"
  )
  plot(bhist
    ,freq=F
    ,xlim=xlim
    ,ylim=ylim
    ,col=colors[2]
    ,add=T
    ,border="transparent"
  )
  plot(overlap
    ,freq=F
    ,xlim=xlim
    ,ylim=ylim
    ,col=colors[3]
    ,add=T
    ,border="transparent"
  )
  if(!(is.null(aname)&is.null(bname))){
      leg1<-legend(x=put.legend
      ,legend=(c(aname,bname,"Overlap") )
      ,fill=colors
      ,bty="n"
      ,cex=2
    )  
  }
}
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#Warning: this pretty slow...
Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot.status.check <- function(
  #_____________________________________________________________________________
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  #Plot and convolve two overlapping distributions
  a    #these must be vectors of random deviates from a distribution
  ,b   #using the r functions like rnorm() works well...
  ,aname=NULL  #name for the a distribution - will go in legend
  ,bname=NULL  #name for the b distribution - will go in legend
  ,colors=c("blue","red","purple")  #colors for each distribution and the
  #overalp between them
  ,b1=1000 #the number of breaks in the plot of the historgram of a and b
  ,xlim=NULL #predefined limits (if you want to restrict your plot)
  ,ylim=NULL
  ,xlab=NULL #other plotting info
  ,ylab=NULL
  ,main=NULL
  ,sub1=NULL  #subtext for plot
  ,put.legend="topright" #where should the legend go?
  ,status.test="overfishing"  #could be "overfished"
  ) {                       #DRH 2/5/13
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  #
  #             #EXAMPLE#
  #             n<-1000000
  #             xA<-rgamma(n,shape=500)
  #             xB<-rlnorm(n,5.9189,0.1668)
  #             plot.convolve(a=xA,b=xB,aname="Biomass in 2011",bname="Threshold"
  #             ,xlab="Biomass",ylab="Probability density")
  #_____________________________________________________________________________
  #clear arrays and make hist objects
  ahist=NULL
  bhist=NULL
  ahist=hist(a,breaks=b1,plot=F)
  bhist=hist(b,breaks=b1,plot=F)
  #some jerking around to ge the histograms on the same scale...
  dist1 = ahist$breaks[2]-ahist$breaks[1]
  br1 = seq((min(ahist$breaks,bhist$breaks)-dist1)
      ,(max(ahist$breaks,bhist$breaks)+dist1),dist1)
  bhist=hist(b,breaks=br1,plot=F)
  ahist=hist(a,breaks=br1,plot=F)

  #rescaling the density portion of the histogram to force sum(density)=1
  ahist$density<-ahist$density/sum(ahist$density)
  bhist$density<-bhist$density/sum(bhist$density)
  
  p<-0
  a.star<-a[a<max(b)]
  b.star<-b[b>min(a)]
  print(length(a.star))
  print(length(b.star))
  for (i in 1:length(b.star)) {
      p<-p+length(a.star[a.star<=b.star[i]]) #instances of overfishing
    }
  p2<-p/(length(a)^2)
  
  if(is.null(xlim)){
    xlim = c(min(ahist$breaks,bhist$breaks),max(ahist$breaks,bhist$breaks))
  }
  if(is.null(ylim)){
    ylim = c(0,max(ahist$density,bhist$density))
  }
  overlap = ahist
  #the area of overlap is easily defined by taking advantage of the hist properties
  #Note: the approximation to true overlap (based on integrals) degrades as b1 gets small
  #200 seems more or less adequate for most functions, but there is no reason to skimp!
  for(i in 1:length(overlap$density)){
    if(ahist$density[i] > 0 & bhist$density[i] > 0){
      overlap$density[i] = min(ahist$density[i],bhist$density[i])
    } else {
      overlap$density[i] = 0
    }
  }
  if(is.null(main)) p1<-paste("P[",status.test,"]~",round(p2,3),sep="")
  if(!(is.null(main))) p1<-main
  print(paste("P[overlap]~",round(sum(overlap$density),3),sep=""))
  print(paste("P[overfishing]~",round(p2,3),sep=""))
  par(cex.main=3)
  plot(ahist
    ,freq=F
    ,xlim=xlim
    ,ylim=ylim
    ,col=colors[1]
    ,xlab=xlab
    ,ylab=ylab
    ,main=p1
    ,sub=sub1
    ,border="transparent"
  )
  plot(bhist
    ,freq=F
    ,xlim=xlim
    ,ylim=ylim
    ,col=colors[2]
    ,add=T
    ,border="transparent"
  )
  plot(overlap
    ,freq=F
    ,xlim=xlim
    ,ylim=ylim
    ,col=colors[3]
    ,add=T
    ,border="transparent"
  )
  if(!(is.null(aname)&is.null(bname))){
      leg1<-legend(x=put.legend
      ,legend=(c(aname,bname,"Overlap") )
      ,fill=colors
      ,bty="n"
      ,cex=2
    )
  }
  return(data.frame("Overlap"=round(sum(overlap$density),3),"Probability"=round(p2,5)))
}
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot.status.check.F <- function(
  #_____________________________________________________________________________
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  #Plot and convolve two overlapping distributions
  a    #these must be vectors of random deviates from a distribution
        #using the r functions like rnorm() works well...
  ,b.val  #single threshold value 
  ,aname=NULL  #name for the a distribution - will go in legend
  ,bname=NULL  #name for the b distribution - will go in legend
  ,colors=c("black","red")  #colors for each distribution and the
  #overalp between them
  ,b1=1000 #the number of breaks in the plot of the historgram of a and b
  ,xlim=NULL #predefined limits (if you want to restrict your plot)
  ,ylim=NULL
  ,xlab=NULL #other plotting info
  ,ylab=NULL
  ,main=NULL
  #,sub1=NULL  #subtext for plot
  ,put.legend="topright" #where should the legend go?
  ,status.test="overfishing"  #could be "overfished"
  ,leg=T #put a legend?
  ) {                       #DRH 2/5/13
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  #
  #             #EXAMPLE#
  #             n<-1000000
  #             xA<-rgamma(n,shape=500)
  #             xB<-rlnorm(n,5.9189,0.1668)
  #             plot.convolve(a=xA,b=xB,aname="Biomass in 2011",bname="Threshold"
  #             ,xlab="Biomass",ylab="Probability density")
  #_____________________________________________________________________________
  #clear arrays and make hist objects
  ahist=NULL
  ahist=hist(a,breaks=b1,plot=F)
  #rescaling the density portion of the histogram to force sum(density)=1
  ahist$density<-ahist$density/sum(ahist$density)
  
  p<-0

  p<-length(a[a>=b.val]) #instances of overfishing
  p2<-p/(length(a))

  if(is.null(main)) p1<-paste("P[",status.test,"]~",round(p2,3),sep="")
  if(!(is.null(main))) p1<-main
  print(paste("P[overfishing]~",round(p2,3),sep=""))
  if(is.null(xlim)){
    xlim = c(0,max(c(max(ahist$breaks),(b.val+0.5*b.val))))
  }
  if(is.null(ylim)){
    ylim = c(0,max(ahist$density))
  }
  
  par(cex.main=3)
  plot(ahist
    ,freq=F
    ,xlim=xlim
    ,ylim=ylim
    ,col=colors[1]
    ,xlab=xlab
    ,ylab=ylab
    ,main=p1
 #   ,sub=sub1
    ,border="transparent"
  )
  abline(v=b.val,lwd=3,lty=1,col=colors[2])
  if(!(is.null(aname)&is.null(bname))&leg){
      leg1<-legend(x=put.legend
      ,legend=(c(aname,bname) )
      ,fill=colors
      ,bty="n"
      ,cex=2
    )
  }
  return(data.frame("Probability"=round(p2,5)))
}
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot.correlated.status.check <- function(
  #_____________________________________________________________________________
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  #Plot and conrast two correlated overlapping lognormal distributions, will make
  #a plot and return the probability of exceeding the threshold ie overfishing 
  #or overfished status 
  mu=NULL #vector of means of the two distributions one for the estimate, one
  #for the reference point (given in real space, not logged)
  ,cv1=NULL #vector of cv's
  ,cor1=NULL #correlation of the two above
  ,aname=NULL  #name for the a marginal distribution A[,1] - will go in legend
  ,bname=NULL  #name for the b marginal distribution A[,2] - will go in legend
  ,colors=c("blue","red","purple")  #colors for each marginal distribution and the
  #overalp between them
  ,b1=1000 #the number of breaks in the plot of the historgram of a and b
  ,xlim=NULL #predefined limits (if you want to restrict your plot)
  ,ylim=NULL
  ,xlab=NULL #other plotting info
  ,ylab=NULL
  ,main=NULL
  ,sub1=NULL  #subtext for plot
  ,put.legend="topright" #where should the legend go?
  ,status.test="overfishing"  #could be "overfished"
  ,leg=T #put a legend?
  ) {                       #DRH 2/5/13
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  #
  #                           #EXAMPLE#
  #                plot.correlated.status.check(mu=c(1000,500)
  #                   ,cv1=c(0.22,0.18)
  #                   ,cor1=0.9
  #                   ,aname="B2011"
  #                   ,bname="Threshold"
  #                   ,status.test="overfished"
  #                  )
  #_____________________________________________________________________________
  #uses some routines from the functions file by DRH...
  library(MASS)
  #covariance = cor1*(var1*var2)
  dx1<-std.lnorm(mu[1],cv1[1])
  dx2<-std.lnorm(mu[2],cv1[2])
  var1<-dx1[,2]^2
  var2<-dx2[,2]^2
  Sigma<-matrix(c(var1,rep(cor1*(dx1[,2]*dx2[,2]),2),var2),2,2)
  n<-10000000  #one million iterations ought to be enough!
  #check
  #Sigma = var(mvrnorm(n,mu,Sigma))
  if(cv1[1]>0 & cv1[2]>0 ) {A<-mvrnorm(n,c(dx1[,1],dx2[,1]),Sigma)}
  a<-logify(A[,1])
  b<-logify(A[,2])

  #clear arrays and make hist objects
  ahist=NULL
  bhist=NULL
  ahist=hist(a,breaks=b1,plot=F)
  bhist=hist(b,breaks=b1,plot=F)
  #some jerking around to ge the histograms on the same scale...
  dist1 = ahist$breaks[2]-ahist$breaks[1]
  br1 = seq((min(ahist$breaks,bhist$breaks)-dist1)
      ,(max(ahist$breaks,bhist$breaks)+dist1),dist1)
  bhist=hist(b,breaks=br1,plot=F)
  ahist=hist(a,breaks=br1,plot=F)

  #rescaling the density portion of the histogram to force sum(density)=1
  ahist$density<-ahist$density/sum(ahist$density)
  bhist$density<-bhist$density/sum(bhist$density)

  if(status.test=="overfished") p<-length(a[a<=b]) #instances of overfishing only need to compare accross
  if(status.test=="overfishing") p<-length(a[a>=b])
  p2<-p/length(a)    #rows because each row is one draw from the correlated dists

  if(is.null(xlim)){
    xlim = c(min(ahist$breaks,bhist$breaks),max(ahist$breaks,bhist$breaks))
  }
  if(is.null(ylim)){
    ylim = c(0,max(c(ahist$density,bhist$density)))
  }
  overlap = ahist
  #the area of overlap is easily defined by taking advantage of the hist properties
  #Note: the approximation to true overlap (based on integrals) degrades as b1 gets small
  #200 seems more or less adequate for most functions, but there is no reason to skimp!
  for(i in 1:length(overlap$density)){
    if(ahist$density[i] > 0 & bhist$density[i] > 0){
      overlap$density[i] = min(ahist$density[i],bhist$density[i])
    } else {
      overlap$density[i] = 0
    }
  }
  if(is.null(main)) p1<-paste("P[",status.test,"]~",round(p2,3),sep="")
  if(!(is.null(main))) p1<-main
  print(paste("P[overlap]~",round(sum(overlap$density),3),sep=""))
  print(paste("P[",status.test,"]~",round(p2,5),sep=""))
  par(cex.main=3)
  plot(ahist
    ,freq=F
    ,xlim=xlim
    ,ylim=ylim
    ,col=colors[1]
    ,xlab=xlab
    ,ylab=ylab
    ,main=p1
    ,sub=sub1
    ,border="transparent"
  )
  plot(bhist
    ,freq=F
    ,xlim=xlim
    ,ylim=ylim
    ,col=colors[2]
    ,add=T
    ,border="transparent"
  )
  plot(overlap
    ,freq=F
    ,xlim=xlim
    ,ylim=ylim
    ,col=colors[3]
    ,add=T
    ,border="transparent"
  )
  if(!(is.null(aname)&is.null(bname))&leg){
      leg1<-legend(x=put.legend
      ,legend=(c(aname,bname,"Overlap") )
      ,fill=colors
      ,bty="n"
      ,cex=2
    )
  } 
  return(data.frame("Overlap"=round(sum(overlap$density),3),"Probability"=round(p2,5)))
}
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Status checks without pictures - for each year in the projection
correlated.status.check<-function(mu,cv1,cor1){
  library(MASS)
  #generate correalted random lognormals and check for overlap 
  #_____________________________________________________________________________
  #mu=vector of means of the two distributions one for the estimate, one
  #for the reference point (given in real space, not logged)
  #cv1=vector of cv's
  #cor1=correlation of the two above  
  #_____________________________________________________________________________
  dx1<-std.lnorm(mu[1],cv1[1])
  dx2<-std.lnorm(mu[2],cv1[2])
  var1<-dx1[,2]^2
  var2<-dx2[,2]^2
  Sigma<-matrix(c(var1,rep(cor1*(dx1[,2]*dx2[,2]),2),var2),2,2)
  n<-1000000  #one million iterations ought to be enough!
  #check
  #Sigma = var(mvrnorm(n,mu,Sigma))
  A<-mvrnorm(n,c(dx1[,1],dx2[,1]),Sigma)

  a<-logify(A[,1])
  b<-logify(A[,2])

  p<-length(a[a<=b]) #instances of overfished
  p2<-p/length(a)
  return(p2)
}
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# F Status checks without pictures - for each year in the projection
F.status.check<-function(mu,cv1,thresh,n=1000000){
  #generate correalted random lognormals and check for overlap 
  #_____________________________________________________________________________
  #mu=vector of means of the two distributions one for the estimate, one
  #for the reference point (given in real space, not logged)
  #cv1=vector of cv's 
  #threshold F for comparison
  #n = # of comparisons
  #_____________________________________________________________________________
  dx1<-std.lnorm(mu,cv1)
  var1<-dx1[,2]
  A<-rnorm(n,dx1[,1],var1)

  a<-logify(A)

  p<-length(a[a>=thresh]) #instances of overfishing
  p2<-p/length(a)
  return(p2)
}
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
get.seed<-function(){
  #This function generates a randomish seed for use in simulation settings
  library(microbenchmark)
  f<- function() NULL
  res <- microbenchmark(NULL,f(), times=1)
  res.c<-paste(res$time[1])
  l<-nchar(res.c)
  res.2<-as.numeric(substr(res.c,l,l))*as.numeric(substr(res.c,1,2))^3
  if(res.2>0) {rez<-res.2} else {rez<-res$time[1]}
  #print(rez)
  return(rez)
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
test.arg<-function(arg1){
  if(!(length(arg1)>0)){print(str(arg1))}
  return((length(arg1)>0)) #return false if null argument
}
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
reversibly.standardize<-function(x){
    m1<-mean(na.omit(x))
    s1<-sd(na.omit(x))
    vec<-(na.omit(x)-m1)/s1
    return(data.frame("st.x"=vec,"mean"=m1,"sd"=s1))
}
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
performance<-function(dat1,advice){
  #track some useful NSTOW performance metrics
  status<-dat1$SSB[length(dat1$SSB)]/dat1$Unfished.SSB
  if(length(status)==0) status=0.0
  rez<-data.frame("Stability"=abs(sum(diff(dat1$Index)))/sum(dat1$Index)  
    ,"Status"=status
    ,"Advice"=advice$advice,"cor"=advice$correlation,"Extinct"=0)
  return(rez)
}
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fill.maturity<-function(n){
  #fill a maturity vector
  set.seed(get.seed()) 
  mat<-rbeta(n,2,5)
  return(cumsum(sort(mat))/sum(mat))
}
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fill.weight<-function(n){
  #fill a maturity vector
  set.seed(get.seed()) 
  wt<-rbeta(n,2,5)
  return(cumsum(sort(wt)))
}
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fill.selx<-function(n){
  #fill a maturity vector
  set.seed(get.seed())  
  slx<-rbeta(n,2,5)
  return(sort(slx))
}
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pool.cv<-function(x){return(sqrt(mean2(na.omit(x)^2)))}
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pool.sd<-function(v,n){
#v is the variance of x : n is the sample size of each x
vx<-na.omit(v)
if(length(vx)!=length(n)) return("OMIT NA's!!")
ps<-sum((n-1)*vx)/(sum(n)-length(n))
return(sqrt(ps))
}
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#Loads some sensor data into memory
get.2012.sensor.data<-function(){
dir.data<-paste("/home/dhennen//"
  ,"Clam/Assessment Stuff/Survey Analysis/Sensors/data/",sep="")

s.data<-read.csv(paste(dir.data,"sensor.data.2012.csv",sep=""))
#print(dim(s.data))
names(s.data)<-c("Cruise","Station","Time","Depth","Lat","Lon"
    ,"SOG","Ambient.pressure","Manifold.pressure","Star.temp"
    ,"Star.tilt.x","Star.tilt.y","Star.tilt.z")
names(s.data)
selx<-c(33,53,108,117,127,136,150,162,170,178,182,184,172) #selectivity stations don't have sensors
s.d2<-s.data[!s.data$Station%in%selx,]
#print(dim(s.d2))
#unique(s.d2$Station)
return(s.d2)
}
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## return the names of the objects (from a vector of list of
## names of objects) that are functions and have debug flag set
isdebugged_safe <- function(x,ns=NULL)  {
    g <- if (is.null(ns)) get(x) else getFromNamespace(x,ns)
    is.function(g) && isdebugged(g)
}

which_debugged <- function(objnames,ns=NULL) {
    if (!length(objnames)) return(character(0))
    objnames[sapply(objnames,isdebugged_safe,ns=ns)]
}

all_debugged <- function(where=search(), show_empty=FALSE) {
    ss <- setNames(lapply(where,function(x) {
        which_debugged(ls(x,all.names=TRUE))
        }),gsub("package:","",where))
    ## find attached namespaces
    ## (is there a better way to test whether a namespace exists with a given name??)
    ns <- unlist(sapply(gsub("package:","",where),
                 function(x) {
                     if (inherits({n <- try(getNamespace(x),silent=TRUE)},
                         "try-error")) NULL else x
                 }))
    ss_ns <- setNames(lapply(ns,function(x) {
        objects <- ls(getNamespace(x),all.names=TRUE)
        which_debugged(objects,ns=x)
        }),ns)
    if (!show_empty) {
        ss <- ss[sapply(ss,length)>0]
        ss_ns <- ss_ns[sapply(ss_ns,length)>0]
    }
    ## drop overlaps
    for (i in names(ss))
        ss_ns[[i]] <- setdiff(ss_ns[[i]],ss[[i]])
    list(env=ss,ns=ss_ns)
}

undebug_all <- function(where=search()) {
    aa <- all_debugged(where)
    lapply(aa$env,undebug)
    ## now debug namespaces
    invisible(mapply(function(ns,fun) {
        undebug(getFromNamespace(fun,ns))
    },names(aa$ns),aa$ns))
}
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
make.tex.table<-function(x,pathname,filename,caption=NULL
  ,display=NULL,tab.env="longtable",align=c(rep("l",2),rep("c",(dim(x)[2]-1)))
  ,hline=NULL,add.to.row=NULL,caption.placement="top",...) {
#Make a latex table out of dataframe x.  User supplies path and file name
#for the file to be created.  Also requires an argument called
#caption which is a string that is the caption you want to go with your table.
#There is also an argument called display which is a Character vector of length 
#equal to the number of columns of the resulting table indicating the format for 
#the corresponding columns. Since the row names are printed in the first column, 
#the length of display is one greater than ncol(x) if x is a data.frame. These 
#values are passed to the formatC function. Use "d" (for integers), 
#"f", "e", "E", "g", "G", "fg" (for reals), or "s" (for strings). 
#"f" gives numbers in the usual xxx.xxx format; "e" and "E" give n.ddde+nn or 
#n.dddE+nn (scientific format); "g" and "G" put x[i] into scientific format only 
#if it saves space to do so. "fg" uses fixed format as "f", but digits as number 
#of significant digits. Note that this can lead to quite long result strings. 
#Default depends on the class of x.
#                       DRH             8/13/13
#_______________________________________________________________________________ 
library(xtable)
write.csv(x,file=paste(pathname,filename,sep=""),quote=F,na="",row.names=F )
#also produce an x table version of this for use in latex
xtab1<-xtable(x,caption=caption,display=display,align=align
	      ,label=substr(filename,1,(nchar(filename)-4)),...)
#if(!is.null(hline)) hline<-c(-1,0,hline,dim(x)[1])
print.xtable(xtab1
	      ,file=paste(pathname,substr(filename,1,(nchar(filename)-3)),"tex",sep="") 
	      ,include.rownames=FALSE,hline.after=hline,add.to.row=add.to.row
        ,tabular.environment=tab.env,caption.placement=caption.placement,... 
  )	# ,floating=F    
}
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## This function will plot a polygon representing the desired confidence interval from
## the user supplied gam fit to the current plotting region
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot.ci.gam<-function(
  m1=NULL # a fit from a gam object
  ,x.var=NULL #the independent variable from fit1
  ,ci=0.95  #the desired confidence interval (e.g. 95%)
  ,col1="#1F78B4" #color of the confidence region
  ,type="response"
  )
{
  ci.lev <- qnorm((1-(1-ci)/2))  #get the zscore from ci
  fit1<-predict(m1,type=type,se=T)$fit
  se <- predict(m1,type=type,se=T)$se.fit
  lcl <- fit1 - ci.lev * se
  ucl <- fit1 + ci.lev * se
  x.to<-order(x.var)
  x.from<-order(x.var, decreasing=T)
  x.polygon <- c( x.var[x.to] , x.var[x.from] )
  y.polygon <- c( ucl[x.to] , lcl[x.from] )
  #plot polygon to current plot window
  polygon( x.polygon , y.polygon , col = col1 , border = NA )
}
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#Used to control tick marks on log scale axes
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
minor.ticks.axis <- function(ax,n,t.ratio=0.5,mn,mx,...){
#ax is the axis you want to manipulate (1=bottom, 2=left, 3=above and 4=right)  
#n is the number of minor ticks you want
#t.ratio is the ratio between the major and the minor ticks,
#mn and mx are the nim and max (optional) in log space (eg mn=0 and mx=8)  
#and with ... you can pass extra parameters to axis  
  
  lims <- par("usr")
  if(ax %in% c(1,3)) lims <- lims[1:2] else lims <- lims[3:4]
  
  major.ticks <- pretty(lims,n=5)
  
  if(missing(mn)) mn <- min(major.ticks)
  if(missing(mx)) mx <- max(major.ticks)
  
  major.ticks <- major.ticks[major.ticks >= mn & major.ticks <= mx]
  
  labels <- sapply(major.ticks,function(i)
    as.expression(bquote(10^ .(i)))
  )
  axis(ax,at=major.ticks,labels=labels,...)
  
  n <- n+2
  minors <- log10(pretty(10^major.ticks[1:2],n))-major.ticks[1]
  minors <- minors[-c(1,n)]
  
  minor.ticks = c(outer(minors,major.ticks,`+`))
  minor.ticks <- minor.ticks[minor.ticks > mn & minor.ticks < mx]
   
  axis(ax,at=minor.ticks,tcl=par("tcl")*t.ratio,labels=FALSE)
} 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#Convert a named color to a transparent equivalent
transparent<-function(col,deg){
  #col is the color, deg is the amount of transparency desired from 1 (almost invisible) 
  #to 100 (opaque).
  
  #rescale deg to rgb (1:255)
  deg<-(deg/100)*255
  return(rgb(red=col2rgb(col)[1,1],green=col2rgb(col)[2,1],blue=col2rgb(col)[3,1],alpha=deg
      ,maxColorValue=255))  
}
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
col2hex<-function(col){
  #convert a color from names "col" to a hexadecimal representation
  return(rgb(t(col2rgb(col))/255))
}  
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
set.up.OQ.run<-function(){  
  ## First step is to draw population parameters from a random distribution 
  ## or set them up on your own if you prefer...
  nyears<-100    #round(runif(1,20,50),0) #years of simulated data
  nages<-300               #round(runif(1,200,300),0) #ages
  len<-94.2*(1-exp(-0.04*(seq(1,nages)+8.7))) #(see Kilada etal 2007)
  len[1:10]<-len[11]/10*seq(1,10)+0.1  #making growth linear below 10 yrs.  
  Mat<-1/(1+exp(-5.92+0.0927*len)) #vector of maturity at age
  Mat<-1-Mat #values above are from Thorarinsdottir and Jacobson, 2005.
  Mat[150:length(Mat)]<-0.9999 #make animals older than 150 fully mature - this is my addition -DRH
  Wt<-(exp(-9.21)*len^2.82)/1e6 #fill.weight(nages) #vector of weight at age (NEFSC 2013) - convert g to 1000 mt
  med.recr<-NULL #median recruitment for constant or purely random recruitment
  rvar<-runif(1,0.05,0.5) #standard deviation of the recruitment error parameter
  Rscale<-5   #100000000/1.5
  type<-3 #oq type recruitment - see documentation for this project or comments in Spawn code   
  steepness<-0.5 #this should be a random variable to test....
  #round(runif(1,1,2),0) #1=Bev-Holdt, 2=Ricker
  #if(type==1) steepness<-runif(1,0.2,0.99)  #stock recruit parameter for Bev-Holdt type recruitment
  #if(type==2) steepness<-runif(1,0.2,2.0)  #stock recruit parameter for Ricker type recruitment
  #Bzero<-runif(1,100000,10000000)  #stock recruit parameter - mainly affects the asymptote of the B-H
  M<-0.02         #runif(1,0.1,0.5) #natural mortality
  F.dat<-rep(runif(1,0.01,0.1),nyears) #vector of F in each year
  catch<-NULL #vector of catch (weight) in each year
  #st.pop<-runif(nages,10000,Bzero) #population in year one
  selectivity.I<-rep(1,nages) #selx for index 
  #I.error<-runif(1,0.0,0.5) #standard deviation of the index error parameter
  I.error<-0.0
  selectivity.F<-1-(1/(1+exp(-7.63+(0.105*len)))) #selx for fishery from Thorarinsdottir and Jacobson, 2005.     
  selectivity.F[1]<-0  #Need to do this or large recruitment classes will swamp last years total population 
  #and force the extinction trigger. 
  #F.error<-runif(1,0.0,0.25) #standard deviation of the fishery error parameter
  F.error<-0.
  catchab<-1
  Bzero<-100 #fixed for simulation testing at first!
  h_bound<-0.5 #upper bound on steepness parameter for Bev-Holt
  F_bound<-0.15 #upper bound on F
  sigk<-0.0075 #stdev in growth parameter k : sigk=0.0075,sigL=0.0075,sigt0=1.62
  sigL<-0.0075 #stdev in growth parameter Linfinity
  sigt0<-1.62
  burnin<-1000
  ##### set up the recovery time specific variables
  Bstar_vec<-seq(0.25,0.5,0.05)  
  Bdelta_vec<-seq(0.00,0.5,0.1)
  
  interval<-30 #interval for big recruit classes 
  #catchab<-runif(1,0.1,0.8) #index catchability
  return(list("nages"=nages,"nyears"=nyears,"len"=len,"Mat"=Mat,"Wt"=Wt,"rvar"=rvar
  ,"med.recr"=med.recr,"steepness"=steepness,"Bzero"=Bzero,"Rscale"=Rscale
  ,"type"=type,"M"=M,"F.dat"=F.dat,"catch"=catch,"selectivity.I"=selectivity.I
  ,"I.error"=I.error,"selectivity.F"=selectivity.F,"F.error"=F.error,"catchab"=catchab
  ,"h_bound"=h_bound,"F_bound"=F_bound,"Bstar_vec"=Bstar_vec,"Bdelta_vec"=Bdelta_vec
  ,"burnin"=burnin,"sigk"=sigk,"sigL"=sigL,"sigt0"=sigt0,"interval"=interval)  
  )    
}
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#Quahog specific functions
  change.growth<-function(nages,sigk=0.0075,sigL=0.0075,sigt0=1.62){
  #function to change the growth curve a bit given the number of ages and parameters for the 
  #stdev of the change in L infinity and k
    len<-(97.425+rnorm(1,0,sigL))*(1-exp((-0.045+rnorm(1,0,sigk))*(seq(1,nages)+(11.72+rnorm(1,0,sigt0))))) #(see Kilada etal 2007)
    len[1:10]<-len[11]/10*seq(1,10)+0.1  #making growth linear below 10 yrs.
    return(len)
  }
  change.selx.F<-function(len){
  #function to calculate the selectivity at age given the len at age vector  
    selectivity.F<-1-(1/(1+exp(-7.63+(0.105*len)))) #selx for fishery from Thorarinsdottir and Jacobson, 2005.     
    selectivity.F[1]<-0  #Need to do this or large recruitment classes will swamp last years total population 
    return(selectivity.F)      
  }
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## Function to calculate generation times (mean age of parents for a given cohort)
### NOTE: use a very large max age to reduce bias in numerator
# Liz Brooks
# Version 1.0
# Created 22 November 2013
# Last Modified: 
#------- Gen Time ----------------------------------
gen.time<-function(nages,fec.age,mat.age,M.age, F.mult, sel.age, spawn.time ) {
gt=0.0
top=0.0
cum.survive=1.0
z=0.0
for (i in 1:(nages-1)  ) {
  z=M.age[i] + F.mult*sel.age[i]
  z.ts=(M.age[i]+F.mult*sel.age[i])*spawn.time
  top=top+cum.survive*i*fec.age[i]*mat.age[i]*exp(-z.ts)
  cum.survive=cum.survive*exp(-z )
  }
z= M.age[nages] + F.mult*sel.age[nages]
z.ts=(M.age[nages]+F.mult*sel.age[nages])*spawn.time
top=top + fec.age[nages]*mat.age[nages]*cum.survive*exp(-z.ts)/( 1- exp(-z ) )
bottom=0.0
cum.survive=1.0
z=0.0
for (i in 1:(nages-1)  ) {
  z=M.age[i] + F.mult*sel.age[i]
  z.ts=(M.age[i]+F.mult*sel.age[i])*spawn.time
  bottom=bottom+cum.survive*fec.age[i]*mat.age[i]*exp(-z.ts)
  cum.survive=cum.survive*exp(-z )
  }
  z= M.age[nages] + F.mult*sel.age[nages]
  z.ts=(M.age[nages]+F.mult*sel.age[nages])*spawn.time
  bottom=bottom + fec.age[nages]*mat.age[nages]*cum.survive*exp(-z.ts)/( 1- exp(-z ) )
  gt=top/bottom
  return(gt)
  }
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
is.odd<-function(x){ #check if a number is odd
  val<-T 
  if((abs(x) %% 2) < sqrt(.Machine$double.eps)){
     val<-F  
   }
  return(val)
}
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
threeD<-function(x,y,z,bins=9,x.lab,y.lab,z.lab){
  #Divide x into bins and plot y vs z for each bin  
  #dir.func<- "/home/dhennen/DRH_general_R_code/"
  #source(paste(dir.func,"functions.L.r",sep=""))
  library(mgcv)  
  plot.new()
  resize.win(Width=8, Height=10, xpinch=30, ypinch=30)
  par(mar=c(4,5,4,2)+0.1,cex.axis=1.5,cex.lab=2,cex.main=2) #bottom, left, top, right
  plots<-bins
  if(is.odd(bins)) {plots<-bins+1}
  l.out<-matrix(c(1:plots),(plots/2),2, byrow=T)
  nf <- layout(l.out)
  
  xbin<-makebins(x,binsize=(range(x)[2]-range(x)[1])/bins,first.bin=range(x)[1])  
  xtest<-sort(unique(xbin))
  for(i in 1:bins){
    plot(y[xbin==xtest[i]]~z[xbin==xtest[i]]
      ,main=paste(x.lab,"=",round(xtest[i],3))
      ,ylab=y.lab
      ,xlab=z.lab
      )
  #add a spline to help observe trend
  z2<-z[xbin==xtest[i]]
  m1<-gam(y[xbin==xtest[i]]~s(z2,k=5))
  x2<-seq(min(z2),max(z2),(max(z2)-min(z2))/1000)
  lines(x2,predict(m1,data.frame("z2"=x2),type="response")
    ,col="red"
    ,lwd=2
  )
  }     
}
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
reverse<-function(x){return(x[seq(length(x),1)])} #can just use rev()!
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
distancePointLine <- function(x, y, slope, intercept) {
   ## x, y is the point to test.
   ## slope, intercept is the line to check distance.
   ##
   ## Returns distance from the line.
   ##
   ## Returns 9999 on 0 denominator conditions.
   x1 <- x-10
   x2 <- x+10
   y1 <- x1*slope+intercept
   y2 <- x2*slope+intercept
   distancePointSegment(x,y, x1,y1, x2,y2)
  }
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
distancePointSegment <- function(px, py, x1, y1, x2, y2) {
 ## px,py is the point to test.
 ## x1,y1,x2,y2 is the line to check distance.
 ##
 ## Returns distance from the line, or if the intersecting point on the line nearest
 ## the point tested is outside the endpoints of the line, the distance to the
 ## nearest endpoint.
 ##
 ## Returns 9999 on 0 denominator conditions.
 lineMagnitude <- function(x1, y1, x2, y2) sqrt((x2-x1)^2+(y2-y1)^2)
 ans <- NULL
 ix <- iy <- 0   # intersecting point
 lineMag <- lineMagnitude(x1, y1, x2, y2)
 if( lineMag < 0.00000001) {
   warning("short segment")
   return(9999)
 }
 u <- (((px - x1) * (x2 - x1)) + ((py - y1) * (y2 - y1)))
 u <- u / (lineMag * lineMag)
 if((u < 0.00001) || (u > 1)) {
   ## closest point does not fall within the line segment, take the shorter distance
   ## to an endpoint
   ix <- lineMagnitude(px, py, x1, y1)
   iy <- lineMagnitude(px, py, x2, y2)
   if(ix > iy)  ans <- iy
   else ans <- ix
 } else {
   ## Intersecting point is on the line, use the formula
   ix <- x1 + u * (x2 - x1)
   iy <- y1 + u * (y2 - y1)
   ans <- lineMagnitude(px, py, ix, iy)
 }
 ans
}
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bin.axis<-function(x=NULL,side=1){
  #____________________________________________________________________________________
  #function to draw major and minor ticks for an axis of x an object made by make.bins()
  big.tics<-2
  #lab1<-pretty(x,n=big.tics)
  lab1<-round(c(min(x),((max(x)-min(x))/2+min(x)),max(x)),2)
  #print(lab1)
  #if(lab1[length(lab1)]>max(x)) lab1<-lab1[1:(length(lab1)-1)] #clip the end of you go over the max value
  at1<-c()
  for(i in 1:length(lab1)){at1[i]<-closest(lab1[i],sort(unique(x)))} #match up the tics to the labels
  #print(sort(unique(x)))
  axis(side=side 
       ,at=at1
       ,labels=lab1
       ,lwd=2 # width of the long axis line is zero, makes invisible
       ,lwd.ticks=1 # width of the etick lines also zero, makes them invisible
  )        
  axis(side=side 
       ,at=seq(1,length(unique(x)))
       ,labels=rep("",length(seq(1,length(unique(x)))   )) 
       ,lwd=0 # width of the long axis line is zero, makes invisible
       ,col.ticks="black"
       ,lwd.ticks=0.5 # width of the etick lines also zero, makes them invisible
       ,tck=-0.02  
  )
} 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
scaled.axis<-function(x=null,side=2){
  #_______________________________________________________________________________________________
  #function to draw a reasonable axis for a box plot of vector x against an object from make.bins()
  big.tics<-2 #5 is default
  axis(side=side 
       ,at=pretty(x,n=big.tics)
       ,labels=pretty(x,n=big.tics) # labels will be for the years
       ,lwd=2 # width of the long axis line is zero, makes invisible
       ,lwd.ticks=1 # width of the etick lines also zero, makes them invisible
  )          
  inte<-(pretty(x,n=big.tics)[2]-pretty(x,n=big.tics)[1])/4
  at1<-seq(0,pretty(x,n=big.tics)[length(pretty(x,n=big.tics))]
          ,inte)
  lab1<-rep("",length(seq(0,pretty(x,n=big.tics)[length(pretty(x,n=big.tics))]
          ,inte)))
  axis(side=side 
       ,at=at1
       ,labels=lab1        
       ,lwd=0 # width of the long axis line is zero, makes invisible
       ,col.ticks="black"
       ,lwd.ticks=0.5 # width of the etick lines also zero, makes them invisible
       ,tck=-0.02
  )   
}
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
filled.contour3 <-
  function (x = seq(0, 1, length.out = nrow(z)),
            y = seq(0, 1, length.out = ncol(z)), z, xlim = range(x, finite = TRUE), 
            ylim = range(y, finite = TRUE), zlim = range(z, finite = TRUE), 
            levels = pretty(zlim, nlevels), nlevels = 20, color.palette = cm.colors, 
            col = color.palette(length(levels) - 1), plot.title, plot.axes, 
            key.title, key.axes, asp = NA, xaxs = "i", yaxs = "i", las = 1, 
            axes = TRUE, frame.plot = axes,mar,rev=F, ...) 
{

  # http://wiki.cbr.washington.edu/qerm/index.php/R/Contour_Plots
  
  # modification by Ian Taylor of the filled.contour function
  # to remove the key and facilitate overplotting with contour()
  # further modified by Carey McGilliard and Bridget Ferris
  # to allow multiple plots on one page

  if (missing(z)) {
    if (!missing(x)) {
      if (is.list(x)) {
        z <- x$z
        y <- x$y
        x <- x$x
      }
      else {
        z <- x
        x <- seq.int(0, 1, length.out = nrow(z))
      }
    }
    else stop("no 'z' matrix specified")
  }
  else if (is.list(x)) {
    y <- x$y
    x <- x$x
  }
  if (any(diff(x) <= 0) || any(diff(y) <= 0)) 
    stop("increasing 'x' and 'y' values expected")
 # mar.orig <- (par.orig <- par(c("mar", "las", "mfrow")))$mar
 # on.exit(par(par.orig))
 # w <- (3 + mar.orig[2]) * par("csi") * 2.54
 # par(las = las)
 # mar <- mar.orig
 if(rev) col=rev(col)
    
 plot.new()
 # par(mar=mar)
  plot.window(xlim, ylim, "", xaxs = xaxs, yaxs = yaxs, asp = asp)
  if (!is.matrix(z) || nrow(z) <= 1 || ncol(z) <= 1) 
    stop("no proper 'z' matrix specified")
  if (!is.double(z)) 
    storage.mode(z) <- "double"
# RV - 10-03-2012
# note replacement of .Internal(filledcontour(as.double(x),...)
# with .filled.contour() as of R-2.15.0
  .filled.contour(as.double(x), as.double(y), z, as.double(levels), 
                          col = col)

                          
                          
  if (missing(plot.axes)) {
    if (axes) {
      title(main = "", xlab = "", ylab = "")
      Axis(x, side = 1)
      Axis(y, side = 2)
    }
  }
  else plot.axes
  if (frame.plot) 
    box()
  if (missing(plot.title)) 
    title(...)
  else plot.title
  invisible()
}
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
filled.legend <-
  function (x = seq(0, 1, length.out = nrow(z)), y = seq(0, 1, 
    length.out = ncol(z)), z, xlim = range(x, finite = TRUE), 
    ylim = range(y, finite = TRUE), zlim = range(z, finite = TRUE), 
    levels = pretty(zlim, n=nlevels), nlevels = 20, color.palette = cm.colors, 
    col = color.palette(length(levels) - 1), plot.title, plot.axes, 
    key.title, key.axes, asp = NA, xaxs = "i", yaxs = "i", las = 1, rev=F,
    axes = TRUE, frame.plot = axes,side=side, ...) 
{
  
    # x = seq(0, 1, length.out = nrow(z)); y = seq(0, 1, 
    # length.out = ncol(z)); z; xlim = range(x, finite = TRUE); 
    # ylim = range(y, finite = TRUE); zlim = range(z, finite = TRUE); 
    # levels = pretty(zlim, n=nlevels); nlevels = 20; color.palette = cm.colors; 
    # col = color.palette(length(levels) - 1); plot.title=NULL; plot.axes=NULL; 
    # key.title=NULL; key.axes=NULL; asp = NA; xaxs = "i"; yaxs = "i"; las = 1; rev=F;
    # axes = TRUE; frame.plot = axes;side=side;  
    
  # modification of filled.contour by Carey McGilliard and Bridget Ferris
  # designed to just plot the legend
    #modified by DRH to make axis appear on left or right
    #and to reverse the order of labels if rev=T
    if (missing(z)) {
        if (!missing(x)) {
            if (is.list(x)) {
                z <- x$z
                y <- x$y
                x <- x$x
            }
            else {
                z <- x
                x <- seq.int(0, 1, length.out = nrow(z))
            }
        }
        else stop("no 'z' matrix specified")
    }
    if (is.list(x)) {
        y <- x$y
        x <- x$x
    }
    if (any(diff(x) <= 0) || any(diff(y) <= 0)) stop("increasing 'x' and 'y' values expected")
  #  mar.orig <- (par.orig <- par(c("mar", "las", "mfrow")))$mar
  #  on.exit(par(par.orig))
  #  w <- (3 + mar.orig[2L]) * par("csi") * 2.54
    #layout(matrix(c(2, 1), ncol = 2L), widths = c(1, lcm(w)))
  #  par(las = las)
  #  mar <- mar.orig
  #  mar[4L] <- mar[2L]
  #  mar[2L] <- 1
  #  par(mar = mar)
   # plot.new()
    plot.window(xlim = c(0, 1), ylim = range(levels), xaxs = "i", 
        yaxs = "i")
    rect(0, levels[-length(levels)], 1, levels[-1L], col = col)
    if (missing(key.axes)) {
        if (axes) 
            if(!rev) axis(side) else axis(side,labels = rev(levels), at=levels)
    } else key.axes
    box()
}
    #
#    if (!missing(key.title)) 
#        key.title
#    mar <- mar.orig
#    mar[4L] <- 1
#    par(mar = mar)
#    plot.new()
#    plot.window(xlim, ylim, "", xaxs = xaxs, yaxs = yaxs, asp = asp)
#    if (!is.matrix(z) || nrow(z) <= 1L || ncol(z) <= 1L) 
#        stop("no proper 'z' matrix specified")
#    if (!is.double(z)) 
#        storage.mode(z) <- "double"
#    .Internal(filledcontour(as.double(x), as.double(y), z, as.double(levels), 
#        col = col))
#    if (missing(plot.axes)) {
#        if (axes) {
#            title(main = "", xlab = "", ylab = "")
#            Axis(x, side = 1)
#            Axis(y, side = 2)
#        }
#    }
#    else plot.axes
#    if (frame.plot) 
#        box()
#    if (missing(plot.title)) 
#        title(...)
#    else plot.title
#    invisible()
#}
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MakeLetter <- function(a="A", where="topleft", cex=2)
     legend(where, pt.cex=0, bty="n", title=a, cex=cex, legend=NA)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PlotRegression <- function (fit,labs=NULL,CI=T,...) {
# put the results of a simple linear regression into 1 plot.
#fit is the result of an lm fit (e.g. fit<-lm(your data model here))
#labs are optional quoted labels for x and y (labs[1] and labs[2], respectively )  
if(!is.null(labs)){x<-labs[1];y<-labs[2];} else{x <- names(fit$model)[2];y = names(fit$model)[1];}
plot(fit$model[[1]]~fit$model[[2]]
    ,main = paste("Adj R2 = ",signif(summary(fit)$adj.r.squared, 5),
                     "; Intercept =",signif(fit$coef[[1]],5 ),
                     "; Slope =",signif(fit$coef[[2]], 5),
                     "; P =",signif(summary(fit)$coef[2,4], 5)) 
    ,xlab=x 
    ,ylab=y
    ,...
    )
abline(fit,col="red",lwd=2)
if(CI) {plot.ci.gam(m1=fit,x.var=fit$model[[2]] #the independent variable from fit1
  ,ci=0.95  #the desired confidence interval (e.g. 95%)
  ,col1=transparent("gray",50)) #color of the confidence region
}
}
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
convert.factors.to.strings.in.dataframe <- function(dataframe)
    {
        class.data  <- sapply(dataframe, class)
        factor.vars <- class.data[class.data == "factor"]
        for (colname in names(factor.vars))
        {
            dataframe[,colname] <- as.character(dataframe[,colname])
        }
        return (dataframe)
    }
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#  Code to determine if a point is in a complex polygon of any shape and
# with any number of vertices
#
# Original C code and nice explanation from:
#   www.codeproject.com/Tips/84226/Is-a-Point-inside-a-Polygon
# and
# http://alienryderflex.com/polygon/
#  
#   Adapted for R as an include file by Larry Jacobson and Jiashen Tang
# on March 30, 2014
#   
#   Input variables:
#
# polyX    =  horizontal coordinates of corners (vector set in calling code)
# polyY    =  vertical coordinates of corners  (vector set in calling code)
# x,y  =   x and y coordinates of points set to be tested
#
# Output variables:
#   The code will return a boolean value (TRUE) if the point x,y is inside the polygon, or
# 0 (FALSE) if it is not.  If the point is exactly on the edge of the polygon,
# then inpoly may be set to TRUE or FALSE.
#
pointInPolygon<-function(polyX=0,polyY=0,x=0,y=0)
{
  if (length(polyX)!=length(polyY)) stop("Coordinates of polygon do not match!")
  if (any(is.na(c(polyX,polyY,x,y)))) stop("Missing value in input data!")
  polySides <-length(polyX)
  j <- polySides
  oddNodes <- vector(mode='logical',length=length(x))
  for (i in 1:polySides)
  {
    Ind_1 <- (polyY[i]<y & polyY[j]>=y) | (polyY[j]<y & polyY[i]>=y)
    Ind_2 <- polyX[i]+(y-polyY[i])/(polyY[j]-polyY[i])*(polyX[j]-polyX[i])<x
    Ind_3 <- Ind_1 & Ind_2
    oddNodes[Ind_3] <- !oddNodes[Ind_3]
    j <- i
  }
  return(oddNodes)    
}
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#Assign regions based on stat area...
assign.region.by.stat.area<-function(area) {
    if (area %in% c(511, 512, 513, 514, 515)) {
      Region = 'MNE'      
    } else if (area %in% c(522, 525, 551, 552, 561, 562)) {
      Region = 'GBK'      
    } else if (area %in% c(521, 526, 533, 534, 537, 538, 539, 541)) {
      Region = 'SNE'
    } else if (area %in% c(612, 613, 616)) {
      Region = 'LI'      
    } else if (area %in% c(614, 615, 622, 623, 624)) {
      Region = 'NJ'
    } else if (area %in% c(621, 625, 626, 627, 628, 629)) {
      Region = 'DMV'
    } else if (area >=630 & area <= 640) {Region = 'SVA'}
  return(Region)     
}      
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#Assign regions based on stat area... vertorized 
assign.region.by.stat.area2<-function(area) {
    Region<-ifelse(area %in% c(511, 512, 513, 514, 515),"MNE"
      ,ifelse(area %in% c(522, 525, 551, 552, 561, 562),"GBK" 
         ,ifelse(area %in% c(521, 526, 533, 534, 537, 538, 539, 541),"SNE"
            ,ifelse(area %in% c(612, 613, 616),"LI"
               ,ifelse(area %in% c(614, 615, 622, 623, 624),"NJ"
                  ,ifelse(area %in% c(621, 625, 626, 627, 628, 629),"DMV"
                    ,ifelse((area >=630 & area < 640),"SVA","OTH")))))))
  return(Region)     
}
#area.test<-c(612,615,552,643,999,621)
#assign.region.by.stat.area2(areatest)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gammaCI<-function(ts,cv,lp=0.1,up=0.9){
  #ts is your data and cv its cv
  #lp and up are the lower and upper probabilities for your gamma CI's
  cv2<-(1/cv)^2
  scale=ts/cv2
  shape=ts/scale
  LB=qgamma(p=lp,shape=shape,scale=scale)
  UB=qgamma(p=up,shape=shape,scale=scale)
  return(data.frame("lci"=LB,"uci"=UB))
}
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# color.column - color the columns of a table for LaTeX output.
#
# Takes an xtable object and a column name as input, and returns the
# xtable with the chosen column stringified and prepended with color
# values using the \cellcolor command from the xcolor package. For an
# easy copy-paste into a LaTeX document, use xtable's print function
# with the following parameter to preserve the \cellcolor command:
# "sanitize.text.function=function(x){x}". I also like to add
# "include.rownames=F". Here's a worked example:
#
# Within your LaTeX preamble, be sure to add the xcolor package with a
# 'table' argument and you're good to go:
#
# \usepackage[table]{xcolor}
#
color.column <- function(xt, column, limits=range(xt[, column], na.rm=T), low='white', high='blue') {
## Make sure the colorspace library is loaded to convert rgb to hex.
library(colorspace)
color.f <- colorRamp(c(low, high))
vals <- xt[, column]
## Rescale the values from (limits[1], limits[2]) to (0, 1)
vals <- (vals - limits[1])/diff(limits) * diff(c(0,1)) + 0
## Clip values outside the desired range.
vals <- ifelse(!is.finite(vals) | vals %inside% c(0,1), vals, NA)
 
clrs <- color.f(vals)
clr.string <- hex(RGB(clrs/255))
clr.string <- substring(clr.string, 2) # Remove the hash prefix.
 
indx <- which(colnames(xt) == column)
# Grab the digits and display values set by xtable
digt <- attr(xt, 'digits')[indx+1]
disply <- attr(xt, 'display')[indx+1]
formatted.vals <- formatC(xt[, column], digits=digt, format=disply)
 
formatted.column <- paste('\\cellcolor[HTML]{', clr.string, '}', formatted.vals, sep='')
xt[, column] <- formatted.column
xt
}
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
alter.cell <- function(xt, nrow, ncol, change=NULL, col=NULL) {
  # alter.cell - alter the text of a cell of a table for LaTeX output.
  # xt is an xtable object, row and col are the row and column number you want to alter
  # change is the change you want to make: "ital", "bold" or "boldital" (both) are options so far.
  # col is a color e.g. "red"
  # NOTE: Transparent color doesn't work  
  # NOTE: Requires \usepackage[table]{xcolor} in LaTex preamble!!
  clr.string <- (col2hex(col))
  clr.string <- substring(clr.string, 2) # Remove the hash prefix.
  if(is.numeric(ncol)) { indx.c<-ncol
  } else if(is.finite(which(colnames(xt) == ncol))) {indx.c <- which(colnames(xt) == ncol)
  } else(return("ERROR unknown column reference in alter.cell"))
  if(is.numeric(nrow)) {indx.r<-nrow
  } else if(is.finite(which(rownames(xt) == nrow))) {indx.r <- which(rownames(xt) == nrow)
  } else(return("ERROR unknown row reference in alter.cell"))  
  # Grab the digits and display values set by xtable
  digt <- attr(xt, 'digits')[indx.c+1]
  disply <- attr(xt, 'display')[indx.c+1]
  
  formatted.vals <- formatC(paste(xt[nrow, ncol]), digits=digt, format=disply)
   
  if(!is.null(change)){ 
    if(change=="bold") fmt.string <- paste('\\textbf{', formatted.vals, '}', sep='')
    if(change=="ital" | change=="italics") fmt.string <- paste('\\textit{', formatted.vals, '}', sep='')  
    if(change=="boldital" | change=="italbold") fmt.string <- paste('\\textit{\\textbf{', formatted.vals, '}}', sep='')  
  } else {fmt.string <- formatted.vals}
  if(!is.null(col)) { formatted.cell <- paste('\\cellcolor[HTML]{',clr.string, '}',fmt.string, sep='')
  } else {formatted.cell <- fmt.string }
  if(!is.null(attr(xt[,ncol],"levels"))) {
    attr(xt[,ncol],"levels")[length(levels(xt[,ncol]))+1]<-formatted.cell} #need to do this to allow a new value to a cell!
  xt[nrow, ncol] <- formatted.cell
  #The digits argument is getting messed up by this change for some reason affects the column of change (but not the row??)
  if(class(try(paste(round(suppressWarnings(as.numeric(xt[!(is.na(as.numeric(xt[,ncol]))), ncol]))
                     ,attr(xt, 'digits')[!(is.na(suppressWarnings(as.numeric(xt[,ncol])))), ncol]))
               ,silent=T))!="try-error"){
  xt[!(is.na(suppressWarnings(as.numeric(xt[,ncol])))), ncol]=
    paste(round(suppressWarnings(as.numeric(xt[!(is.na(as.numeric(xt[,ncol]))), ncol]))
              ,attr(xt, 'digits')[!(is.na(suppressWarnings(as.numeric(xt[,ncol])))), ncol]))
  } #test to make sure this is possible first...
  return(xt)
}
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
alter.cell.old <- function(xt, nrow, ncol, change=NULL, col=NULL) {
  # alter.cell - alter the text of a cell of a table for LaTeX output.
  # xt is an xtable object, row and col are the row and column number you want to alter
  # change is the change you want to make: "ital", "bold" or "boldital" (both) are options so far.
  # col is a color e.g. "blue"
  # NOTE: Transparent color doesn't work  
  clr.string <- (col2hex(col))
  clr.string <- substring(clr.string, 2) # Remove the hash prefix.
  if(is.numeric(ncol)) { indx.c<-ncol
  } else if(is.finite(which(colnames(xt) == ncol))) {indx.c <- which(colnames(xt) == ncol)
  } else(return("ERROR unknown column reference in alter.cell"))
  if(is.numeric(nrow)) {indx.r<-nrow
  } else if(is.finite(which(rownames(xt) == nrow))) {indx.r <- which(rownames(xt) == nrow)
  } else(return("ERROR unknown row reference in alter.cell"))  
  # Grab the digits and display values set by xtable
  digt <- attr(xt, 'digits')[indx.c+1]
  disply <- attr(xt, 'display')[indx.c+1]
  
  formatted.vals <- formatC(paste(xt[nrow, ncol]), digits=digt, format=disply)
  if(!is.null(change)){ 
    if(change=="bold") fmt.string <- paste('\\textbf{', formatted.vals, '}', sep='')
    if(change=="ital") fmt.string <- paste('\\textit{', formatted.vals, '}', sep='')  
    if(change=="boldital" | change=="italbold") fmt.string <- paste('\\textit{\\textbf{', formatted.vals, '}}', sep='')  
  } else {fmt.string <- formatted.vals}
  formatted.cell <- paste('\\cellcolor[HTML]{',clr.string, '}',fmt.string, sep='')
  
  attr(xt[,ncol],"levels")[length(levels(xt[,ncol]))+1]<-formatted.cell #need to do this to allow a new value to a cell!
  xt[nrow, ncol] <- formatted.cell
  return(xt)
}
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
###Fit and plot a logistic curve for selectivity
log.fit <- function(
  x
  ,y
  ,xlab="x"
  ,ylab="y"
  ,main="Logistic selectivity fit s50 with shown with C.I."
  ,strt=list(phi1 = 0.2, phi2 = -.99, phi3 = 3) #you may have to adjust these to get a fit
  ,ci=0.95  #the desired confidence interval (e.g. 95%)
  ,col1="red" #color of the line
  ){
  log.ss <- nls(y ~ phi1/(1 + exp(-(phi2 + phi3*x))),
    start=strt,
    trace=TRUE)
  #this is another way to do it
  #log.ss <- nls(y ~ SSlogis(x, phi1, phi2, phi3))
  C <- summary(log.ss)$coef[1]
  A <- exp((summary(log.ss)$coef[2]) * (1/summary(log.ss)$coef[3]))
  K <- (1 / summary(log.ss)$coef[3])
  plot(y ~ x, main = main, xlab=xlab, ylab=ylab,pch=20,cex=.25)
  xdat<-data.frame(x=min(x):max(x))
  fit1<-predict(log.ss,xdat)  
  lines((x=min(x):max(x)),fit1 , col=col1)
  library(car)
  s50<-deltaMethod(log.ss, "-phi2/phi3")
  rownames(s50)<-"s50"
  #print(s50)  
  abline(v= -coef(log.ss)[2]/coef(log.ss)[3], lty=1)
  #plot the se around the s50 estimate
  ci.lev <- qnorm((1-(1-ci)/2))  #get the zscore from ci
  lcl <- s50[1] - ci.lev * s50[2]
  ucl <- s50[1] + ci.lev * s50[2]
  abline(v=lcl,lty=2)
  abline(v=ucl,lty=2) 
  
  r1 <- sum((x - mean(x))^2)
  r2 <- sum(residuals(log.ss)^2)
  r_sq <- (r1 - r2) / r1
  out <- data.frame(s50, "R.sq"=r_sq)
  dat<-data.frame("data"=(x=min(x):max(x)),"fit"=fit1)
  return(list(out,dat))
}
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
###Fit and plot a logistic curve for selectivity
log.fit2 <- function(
  x
  ,y
  ,xlab="x"
  ,ylab="y"
  ,main="Logistic selectivity fit s50 with shown with C.I."
  ,strt=list(phi1 = 0.2, phi2 = -.99, phi3 = 3) #you may have to adjust these to get a fit
  ,ci=0.95  #the desired confidence interval (e.g. 95%)
  ,col1="red" #color of the line
  ){
  #log.ss <- nls(y ~ phi1/(1 + exp(-(phi2 + phi3*x))),
  #  start=strt,
  #  trace=TRUE)
  #this is another way to do it
  log.ss <- nls(y ~ SSlogis(x, phi1, phi2, phi3))
  plot(y ~ x, main = main, xlab=xlab, ylab=ylab,pch=20,cex=.25)
  xdat<-data.frame(x=min(x):max(x))
  fit1<-predict(log.ss,xdat)  
  lines((x=min(x):max(x)),fit1 , col=col1)
  library(car)
  s50<-summary(log.ss)$coef[2,1]
  #print(s50)  
  abline(v=s50, lty=1)
  #plot the se around the s50 estimate
  ci.lev <- qnorm((1-(1-ci)/2))  #get the zscore from ci
  lcl <- s50 - ci.lev * summary(log.ss)$coef[2,2] 
  ucl <- s50 + ci.lev * summary(log.ss)$coef[2,2]
  abline(v=lcl,lty=2)
  abline(v=ucl,lty=2) 
  
  r1 <- sum((x - mean(x))^2)
  r2 <- sum(residuals(log.ss)^2)
  r_sq <- (r1 - r2) / r1
  out <- data.frame(s50, "R.sq"=r_sq)
  dat<-data.frame("data"=(x=min(x):max(x)),"fit"=fit1)
  return(list(out,dat))
}
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SSlogis2 <- function (input, Asym, xmid, scal) 
  {
      .expr1 <- xmid - input
      .expr3 <- exp(.e2 <- .expr1/scal)
      .expr4 <- 1 + .expr3
      .value <- Asym/.expr4
      .actualArgs <- as.list(match.call()[c("Asym", "xmid", "scal")])
      if (all(unlist(lapply(.actualArgs, is.name)))) {
          .expr10 <- .expr4^2
          .grad <- array(0, c(length(.value), 3L), list(NULL, c("Asym", 
              "xmid", "scal")))
          .grad[, "Asym"] <- 1/.expr4
          .grad[, "xmid"] <- -(xm <- Asym * .expr3/scal/.expr10)
          .grad[, "scal"] <- xm * .e2
          dimnames(.grad) <- list(NULL, .actualArgs)
          attr(.value, "gradient") <- .grad
      }
      .value
  }
  function (mCall, data, LHS) 
  {
      xy <- data.frame(sortedXyData(mCall[["input"]], LHS, data))
      if (nrow(xy) < 4) {
          stop("too few distinct input values to fit a logistic model")
      }
      z <- xy[["y"]]
      if (min(z) <= 0) {
          z <- z - 1.05 * min(z)
      }
      z <- z/(1.05 * max(z))
      xy[["z"]] <- log(z/(1 - z))
      aux <- coef(lm(x ~ z, xy))
      pars <- as.vector(coef(nls(y ~ 1/(1 + exp((xmid - x)/scal)), 
          data = xy, start = list(xmid = aux[1L], scal = aux[2L]), 
          algorithm = "plinear")))
      value <- c(pars[3L], pars[1L], pars[2L])
      names(value) <- mCall[c("Asym", "xmid", "scal")]
      value
  }
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
###Fit and plot a logistic curve for selectivity
log.fit3 <- function(
  x
  ,y
  ,xlab="x"
  ,ylab="y"
  ,main="Logistic selectivity fit s50 with shown with C.I."
  ,strt=list(phi1 = 0.2, phi2 = -.99, phi3 = 3) #you may have to adjust these to get a fit
  ,ci=0.95  #the desired confidence interval (e.g. 95%)
  ,col1="red" #color of the confidence region
  ){
  #log.ss <- nls(y ~ phi1/(1 + exp(-(phi2 + phi3*x))),
  #  start=strt,
  #  trace=TRUE)
  #this is another way to do it
  log.ss <- nls(y ~ SSlogis2(x, phi1, phi2, phi3))
  plot(y ~ x, main = main, xlab=xlab, ylab=ylab,pch=20,cex=.25)
  xdat<-data.frame(x=min(x):max(x))
  fit1<-predict(log.ss,xdat)  
  lines((x=min(x):max(x)),fit1 , col=col1)
  library(car)
  s50<-summary(log.ss)$coef[2,1]
  #print(s50)  
  abline(v=s50, lty=1)
  #plot the se around the s50 estimate
  ci.lev <- qnorm((1-(1-ci)/2))  #get the zscore from ci
  lcl <- s50 - ci.lev * summary(log.ss)$coef[2,2] 
  ucl <- s50 + ci.lev * summary(log.ss)$coef[2,2]
  abline(v=lcl,lty=2)
  abline(v=ucl,lty=2) 
  
  r1 <- sum((x - mean(x))^2)
  r2 <- sum(residuals(log.ss)^2)
  r_sq <- (r1 - r2) / r1
  out <- data.frame(s50, "R.sq"=r_sq)
  dat<-data.frame("data"=(x=min(x):max(x)),"fit"=fit1)
  return(list(out,dat))
}
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#Replace missing column(s) with 0's
fill.col<-function(x,y,fill.val=0){
  y2<-x[1:dim(y)[1],]
  y2[,1:dim(y)[2]]<-y
  y2[,(dim(y)[2]+1):dim(x)[2]]<-rep(rep(fill.val,dim(y)[1]),dim(x)[2]-dim(y)[2])
  return(y2)
}
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#Replace missing column(s) in a matrix with 0's
fill.col.mat<-function(x,y,fill.val=0){
  y2<-matrix(data=x[1,],nrow=dim(y)[1],ncol=dim(x)[2]) 
  y2[,1:dim(y)[2]]<-y
  y2[,(dim(y)[2]+1):dim(x)[2]]<-rep(rep(fill.val,dim(y)[1]),dim(x)[2]-dim(y)[2])
  rownames(y2)<-rownames(y)
  return(y2)
}
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#This is based on the xtable code - should be useful for text cleaning before porting to LaTex!
#Escapes disallowed characters.  
sanitize <- function(str) { 
  result <- str
  #these just need to be escaped
  result <- gsub("$", "\\$", result, fixed = TRUE)
  result <- gsub(">", "$>$", result, fixed = TRUE)
  result <- gsub("<", "$<$", result, fixed = TRUE)
  result <- gsub("|", "$|$", result, fixed = TRUE)
  result <- gsub("%", "\\%", result, fixed = TRUE)
  result <- gsub("\\%", "\\\\%", result, fixed = TRUE)  #need 4 escapes for %!
  result <- gsub("&", "\\&", result, fixed = TRUE)
  result <- gsub("-", "$-$", result, fixed = TRUE)
  result <- gsub(";", "\\;", result, fixed = TRUE)
  result <- gsub("*", "\\*", result, fixed = TRUE)  
  result <- gsub("#", "\\#", result, fixed = TRUE)
  result <- gsub("(", "$($", result, fixed = TRUE) 
  result <- gsub(")", "$)$", result, fixed = TRUE)   
  result <- gsub("~", "\\~{}", result, fixed = TRUE)
  result <- gsub("/", "\\/", result, fixed = TRUE)
  result <- gsub("'", "\\'", result, fixed = TRUE)  
  #close up leading and trailing spaces in parentheses and leading spaces to semicolons
  result <- gsub("\\( ", "\\(", result, fixed = TRUE)  
  result <- gsub(" \\)", "\\)", result, fixed = TRUE)
  result <- gsub(" \\;", "\\;", result, fixed = TRUE)  
  #replace XXth, etc with Latex compatible superscripting 
  result <- gsub("([0-9])(th)", "\\1$^{th}$", result,perl=T)  
  result <- gsub("([0-9])(rd)", "\\1$^{rd}$", result,perl=T)  
  result <- gsub("([0-9])(nd)", "\\1$^{nd}$", result,perl=T)  
  result <- gsub("([0-9])(st)", "\\1$^{st}$", result,perl=T)  
  
  #these are problematic because they are meaningful in common LaTex math scripting 
  #result <- gsub("^", "\\verb|^|", result, fixed = TRUE)
  #result <- gsub("_", "\\_", result, fixed = TRUE)  
  #result <- gsub("{", "\\{", result, fixed = TRUE)
  #result <- gsub("}", "\\}", result, fixed = TRUE)
  #result <- gsub("\\\\", "SANITIZE.BACKSLASH", result)
  #result <- gsub("SANITIZE.BACKSLASH", "$\\backslash$",
  #   result, fixed = TRUE)
  return(result)
}
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#Same as above with a differences accounting for the different way certain items
#are treated in xtable
sanitize.cap <- function(str) { 
  result <- str
  #these just need to be escaped
  result <- gsub("% ", "\\% ", result, fixed = TRUE)
  result <- gsub("\\\\", "\\", result, fixed = TRUE)
  
  #replace XXth, etc with Latex compatible superscripting 
  result <- gsub("([0-9])(th)", "\\1$^{th}$", result,perl=T)  
  result <- gsub("([0-9])(rd)", "\\1$^{rd}$", result,perl=T)  
  result <- gsub("([0-9])(nd)", "\\1$^{nd}$", result,perl=T)  
  result <- gsub("([0-9])(st)", "\\1$^{st}$", result,perl=T)  
  
  return(result)
}
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fac_to_num<-function(x){return(as.numeric(paste(x)))}
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fac_to_str<-function(x){return((paste(x)))}
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
OrderRegionsSC<-function(dt1){
  #dt1 is clam db output
  if(exists("dt1$Region")){
    lev1<-c("SVA","DMV","NJ","LI","SNE","GBK","All")
    lev1<-lev1[lev1%in%unique(dt1$Region)]
    dt1$Region<-factor(dt1$Region,levels=lev1)
    dt1<-dt1[order(dt1$Region),]
  } else{
    lev1<-c("SVA","DMV","NJ","LI","SNE","GBK","All")
    lev1<-lev1[lev1%in%unique(dt1$REGNAM)]
    dt1$REGNAM<-factor(dt1$REGNAM,levels=lev1)
    dt1<-dt1[order(dt1$REGNAM),]
  }    
}
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
OrderRegionsOQ<-function(dt1){
  #dt1 is clam db output
  if(exists("dt1$Region")){
    lev1<-c("SVA","DMV","NJ","LI","SNE","GBK","All")
    lev1<-lev1[lev1%in%unique(dt1$Region)]
    dt1$Region<-factor(dt1$Region,levels=lev1)
    dt1<-dt1[order(dt1$Region),]
  } else{
    lev1<-c("SVA","DMV","NJ","LI","SNE","GBK","All")
    lev1<-lev1[lev1%in%unique(dt1$REGNAM)]
    dt1$REGNAM<-factor(dt1$REGNAM,levels=lev1)
    dt1<-dt1[order(dt1$REGNAM),]
  }  
}
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
strata.med.depth<-function(x){
  #returns the median depth for each shellfish strata (m)
  depth<-0
  if(x%in%seq(80,96)) depth<-(27-9)/2+9
  if(x%in%c(1,5,9,13,17,21,25,29,33,37,41,45,49)) depth<-(46-27)/2+27
  if(x%in%c(2,6,10,14,18,22,26,30,34,38,42,46,50)) depth<-(55-46)/2+46
  if(x%in%c(3,7,11,15,19,23,27,31,35,39,43,47,51)) depth<-(73-55)/2+55
  if(x%in%c(4,8,12,16,20,24,28,32,36,40,44,48,52)) depth<-(110-73)/2+73
  if(x%in%c(68,69,72,73)) depth<-(46)/2
  return(depth)
}
strata.max.depth<-function(x){
  #returns the max depth for each shellfish strata (m)
  depth<-0
  if(x%in%seq(80,96)) depth<-(27)
  if(x%in%c(1,5,9,13,17,21,25,29,33,37,41,45,49)) depth<-(46)
  if(x%in%c(2,6,10,14,18,22,26,30,34,38,42,46,50)) depth<-(55)
  if(x%in%c(3,7,11,15,19,23,27,31,35,39,43,47,51)) depth<-(73)
  if(x%in%c(4,8,12,16,20,24,28,32,36,40,44,48,52)) depth<-(110)
  if(x%in%c(68,69,72,73)) depth<-46
  return(depth)
}
strata.min.depth<-function(x){
  #returns the min depth for each shellfish strata (m)
  depth<-0
  if(x%in%seq(80,96)) depth<-9
  if(x%in%c(1,5,9,13,17,21,25,29,33,37,41,45,49)) depth<-27
  if(x%in%c(2,6,10,14,18,22,26,30,34,38,42,46,50)) depth<-46
  if(x%in%c(3,7,11,15,19,23,27,31,35,39,43,47,51)) depth<-55
  if(x%in%c(4,8,12,16,20,24,28,32,36,40,44,48,52)) depth<-73
  if(x%in%c(68,69,72,73)) depth<-(46)
  return(depth)
}
strata.area=function(x){
  #Return the area (nm) of each strata in x
  areaCheck=c(1137,169,122,141,690,  22,  6,  15,  1894,  190,  246,  270,  1149,  205,  387,  202,  703,  240,  266,  104,  1693
  ,305,  724,  454,  647,  190,  442,  152,  1078,  667,  932,  631,  361,  207,  614,  704,  660,  268,  946,  516,  580
  ,356,  442,  400,  407,  205,  873,  1063,  222,  156,  114,  345,  269,  295,  386,  214,  176,  303,  512,  801,  588
  ,731,  655,  933,  184,  272,  196,  380,  902,  544,  168,  472,  526,  443,  739,  362,  123,  216,  336,  348,  199
  ,356,  484,  343,  117,  342,  165,  97,  215,  278,  490)
  strata=c( 1, 2, 3, 4, 5, 6, 7, 8, 9,10, 11, 12, 13, 14, 15, 16, 17, 18,19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31
 ,32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62
 ,63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96  )
  area=areaCheck[match(x,strata)]
  #should preserve order 
  return(area)
}
#test x=c(58,72,15); strata.area(x);

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cascade<- function(x) return(exp(-cumsum(x)))
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Halibut_Catch_Ramp<-function(maxCatch=5725,noise=500,fishingyrs=200){
  #noise<-500 #noise around the made up catch ramp prior to start of catch reporting
  #fishingyrs<-200
  #maxCatch<-1e4    #5725 is the highest observed value - informs the asymptote of the curve
  A<-1e6;k=.07;
  catchramp<-maxCatch/(1+A*exp(-k*seq(1,fishingyrs)))+noise*rnorm(fishingyrs)
  catchramp<-ifelse(catchramp<=0,abs(catchramp),catchramp)
  catchramp[51:100]<-abs(100*rnorm(50))
  catchramp[1:50]<-abs(10*rnorm(50))
  return(catchramp)
}
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
renderPlot.tex<-function(filename,texname,dir.upload){
 #Make a standalone latex file and compile it, naming it texname and putting it in upload 
 texname2<-paste(dir.upload,"latex/",filename,sep="")
 write.table("\\documentclass[]{report}",file=texname2,quote=F,sep="",row.names=F,col.names=F)
 write.table("\\usepackage[pdftex]{graphicx}",file=texname2,append=T,quote=F,sep="",row.names=F,col.names=F)
 write.table("\\usepackage{multicol}",file=texname2,append=T,quote=F,sep="",row.names=F,col.names=F)
 write.table("\\usepackage{tabu}",file=texname2,append=T,quote=F,sep="",row.names=F,col.names=F)  
 write.table("\\usepackage{setspace}",file=texname2,append=T,quote=F,sep="",row.names=F,col.names=F)  
 write.table("\\usepackage{longtable}",file=texname2,append=T,quote=F,sep="",row.names=F,col.names=F)  
 write.table("\\usepackage{pdflscape}",file=texname2,append=T,quote=F,sep="",row.names=F,col.names=F)  
 write.table("\\usepackage[width=\\textwidth]{caption}",file=texname2,append=T,quote=F,sep="",row.names=F,col.names=F)  
 write.table("\\usepackage{threeparttablex}",file=texname2,append=T,quote=F,sep="",row.names=F,col.names=F)  
 write.table("\\usepackage{tabularx}",file=texname2,append=T,quote=F,sep="",row.names=F,col.names=F)    
 write.table("\\usepackage{amsmath}",file=texname2,append=T,quote=F,sep="",row.names=F,col.names=F)  
 write.table("\\usepackage{rotating}",file=texname2,append=T,quote=F,sep="",row.names=F,col.names=F)  
 write.table("\\usepackage[table]{xcolor}",file=texname2,append=T,quote=F,sep="",row.names=F,col.names=F)
 write.table("\\usepackage[showframe=false]{geometry}",file=texname2,append=T,quote=F,sep="",row.names=F,col.names=F)
 write.table("\\usepackage{adjustbox}",file=texname2,append=T,quote=F,sep="",row.names=F,col.names=F)
 write.table("\\usepackage{float}",file=texname2,append=T,quote=F,sep="",row.names=F,col.names=F)
 write.table("\\usepackage{etoolbox}",file=texname2,append=T,quote=F,sep="",row.names=F,col.names=F)  
 write.table("\\begin{document}",file=texname2,append=T,quote=F,sep="",row.names=F,col.names=F)    
 write.table("\\newcommand*{\\graphwidth}{5in}",file=texname2,append=T,quote=F,sep="",row.names=F,col.names=F) 
 write.table(paste("\\input{",texname,"}",sep=""),file=texname2,append=T,quote=F,sep="",row.names=F,col.names=F)  
 write.table("\\end{document}",file=texname2,append=T,quote=F,sep="",row.names=F,col.names=F) 
 #now compile and clean up...  
 call1<-paste(
   "cd '",dir.upload,"figures/","';"
   ," ls -l;"
   ," pdflatex -halt-on-error ",texname2,";","mv *.log *.aux ../latex"
   ,sep="")
 system(call1);
} 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
renderTable.tex<-function(filename,texname,dir.upload){
  #Make a standalone latex file and compile it, naming it texname and putting it in upload 
  texname2<-paste(dir.upload,"latex/",filename,sep="")
  write.table("\\documentclass[]{report}",file=texname2,quote=F,sep="",row.names=F,col.names=F)
  write.table("\\usepackage[pdftex]{graphicx}",file=texname2,append=T,quote=F,sep="",row.names=F,col.names=F)
  write.table("\\usepackage{multicol}",file=texname2,append=T,quote=F,sep="",row.names=F,col.names=F)
  write.table("\\usepackage{tabu}",file=texname2,append=T,quote=F,sep="",row.names=F,col.names=F)  
  write.table("\\usepackage{setspace}",file=texname2,append=T,quote=F,sep="",row.names=F,col.names=F)  
  write.table("\\usepackage{longtable}",file=texname2,append=T,quote=F,sep="",row.names=F,col.names=F)  
  write.table("\\usepackage{pdflscape}",file=texname2,append=T,quote=F,sep="",row.names=F,col.names=F)  
  write.table("\\usepackage[width=\\textwidth]{caption}",file=texname2,append=T,quote=F,sep="",row.names=F,col.names=F)  
  write.table("\\usepackage{threeparttablex}",file=texname2,append=T,quote=F,sep="",row.names=F,col.names=F)  
  write.table("\\usepackage{tabularx}",file=texname2,append=T,quote=F,sep="",row.names=F,col.names=F)    
  write.table("\\usepackage{amsmath}",file=texname2,append=T,quote=F,sep="",row.names=F,col.names=F)  
  write.table("\\usepackage{rotating}",file=texname2,append=T,quote=F,sep="",row.names=F,col.names=F)  
  write.table("\\usepackage[table]{xcolor}",file=texname2,append=T,quote=F,sep="",row.names=F,col.names=F)
  write.table("\\usepackage[showframe=false]{geometry}",file=texname2,append=T,quote=F,sep="",row.names=F,col.names=F)
  write.table("\\usepackage{adjustbox}",file=texname2,append=T,quote=F,sep="",row.names=F,col.names=F)
  write.table("\\usepackage{float}",file=texname2,append=T,quote=F,sep="",row.names=F,col.names=F)
  write.table("\\usepackage{etoolbox}",file=texname2,append=T,quote=F,sep="",row.names=F,col.names=F)  
  write.table("\\begin{document}",file=texname2,append=T,quote=F,sep="",row.names=F,col.names=F)    
  write.table("\\newcommand*{\\graphwidth}{5in}",file=texname2,append=T,quote=F,sep="",row.names=F,col.names=F) 
  write.table(paste("\\input{",texname,"}",sep=""),file=texname2,append=T,quote=F,sep="",row.names=F,col.names=F)  
  write.table("\\end{document}",file=texname2,append=T,quote=F,sep="",row.names=F,col.names=F) 
  #now compile and clean up...  
  call1<-paste(
    "cd '",dir.upload,"tables/","';"
    ," ls -l;"
    ," pdflatex -halt-on-error ",texname2,";","mv *.log *.aux ../latex"
    ,sep="")
  system(call1);
} 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Solve.F<-function(F1=0.1,selectivity=NULL,M=NULL,Wt=NULL,N.age=NULL,catch1=NULL)
{
  for (n in 1:20) 
  {
    delta=F1*0.0001
    fofFa=sum(((F1*selectivity)/(F1*selectivity+M))*N.age*Wt*(1-exp(-1*(F1*selectivity+M))))
    fprimeFhi=sum((((F1+delta)*selectivity)/((F1+delta)*selectivity+M))*N.age*Wt*(1-exp(-1*((F1+delta)*selectivity+M))))
    fprimeFlo=sum((((F1-delta)*selectivity)/((F1-delta)*selectivity+M))*N.age*Wt*(1-exp(-1*((F1-delta)*selectivity+M))))          
    fofF=fofFa-catch1      
    if(!is.finite(fofF)) {print(paste("ERROR fofFa = ",fofFa," catch1 = ",catch1,sep="")); break;}
    fprimeF=(fprimeFhi-fprimeFlo)/(2.0*delta)
    F1=F1-(fofF/fprimeF)
    if(fofF<abs(0.0001)) {break}
    if (F1<0) F1=0.5*(F1+(fofF/fprimeF))
    if (F1>3) F1=3
  } #end newton raphson for loop
  return(F1)
}
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
getTimeFlag<-function(){tm<-paste(Sys.time());tm<-gsub("-","_",tm);tm<-gsub(":","",tm);tm<-gsub(" ","_",tm);return(tm)}
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot.lognormal=function(mn,sd,ylab="",xlab="",main="",col.curve="red",ylim=NULL,xlim=NULL,axes=T
                        ,r=NULL #this is the range of x values over which to plot
  ){
  if(is.null(r)) r <- seq(0,1,0.01) else r=seq(r[1],r[2],0.01)
  d <- dlnorm(r, meanlog = mn, sdlog = sd)
  #print(log(quantile(d,probs=c(0.1,0.25,0.5,0.75,0.9))))
  #print(ecdf(d))
  if(is.null(ylim)) ylim=range(pretty(d))
  if(is.null(xlim)) xlim=range(pretty(r))
  plot(d~r
       ,col=col.curve
       ,type="l"
       ,lwd=3
       ,lty=1
       ,ylab=ylab
       ,xlab=xlab
       ,main=main
       ,ylim=ylim
       ,xlim=xlim
       ,axes=axes
  )
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
estBetaParams <- function(mu, var) {
  alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
  beta <- alpha * (1 / mu - 1)
  if(alpha<0 & beta<0) return(params = list(alpha = abs(alpha), beta = abs(beta)))
  return(params = list(alpha = (alpha), beta = (beta)))
}
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot.beta=function(mn,sd,alpha=NULL,beta=NULL,ylab="",xlab="",main="",col.curve="red",ylim=NULL,xlim=NULL,axes=T
                        ,r=NULL #this is the range of x values over which to plot
){
  if(is.null(r)) r <- seq(0,1,0.01) else r=seq(r[1],r[2],0.01)
  if(is.null(alpha) | is.null(beta) ) { pars=estBetaParams(mn,sd^2);alpha=pars$alpha;beta=pars$beta}
  d <- dbeta(r, shape1 = alpha, shape2 = beta)
  #print(log(quantile(d,probs=c(0.1,0.25,0.5,0.75,0.9))))
  #print(ecdf(d))
  if(is.null(ylim)) ylim=range(pretty(d))
  if(is.null(xlim)) xlim=range(pretty(r))
  plot(d~r
       ,col=col.curve
       ,type="l"
       ,lwd=3
       ,lty=1
       ,ylab=ylab
       ,xlab=xlab
       ,main=main
       ,ylim=ylim
       ,xlim=xlim
       ,axes=axes
  )
}
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mid=function(x){((x[2]-x[1])/2)+x[1]} #return the midpoint of two values in a vector (e.g. use with range(my.vector))
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#Both of these confidence intervals are set to 95% - change 'em if you want to!
asm.ci<-function(x,cv.x,bounds=95){
  #generate approximate 95% confidence intervals based on log normal variable x
  #and it's cv. cv.x
  s<-sqrt(log(1+cv.x^2))
  s<-ifelse(is.finite(s),s,0)
  p<-(1-(bounds/100))/2
  Z<-qnorm(p)
  lci<-x*exp(Z*(s))
  uci<-x*exp(-Z*(s))
  return(data.frame("lci"=lci,"uci"=uci)) 
}
gamma.ci<-function(ts,cv,bounds=95){
  #ts is your data and cv its cv
  #lp and up are the lower and upper probabilities for your gamma CI's
  cv2<-(1/cv)^2
  scale=ts/cv2
  shape=ts/scale
  lp<-(1-(bounds/100))/2
  up<-1-lp
  LB=qgamma(p=lp,shape=shape,scale=scale)
  UB=qgamma(p=up,shape=shape,scale=scale)
  return(data.frame("lci"=LB,"uci"=UB))
}
normal.ci<-function(x,cv.x,bounds=95){
  #generate approximate (bounds)% confidence intervals based on normal variable x
  #and it's cv. cv.x
  s<-cv.x*x
  s<-ifelse(is.finite(s),s,0)
  error <- qnorm((1-(bounds/100))/2)*s
  lci<-x+error
  uci<-x-error
  return(data.frame("lci"=lci,"uci"=uci)) 
}
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
LineReplace=function(dnm,Find1,Replace1,ext=NULL,recursive=F,backup=T){
  #Find (Find1) and replace with (Replace1) in all files (with extension - or containing - ext) 
  #in directory (dnm) - this can be done recursively or not. You will by default also create backups
  #of any files that you alter, but you can turn that off with backup=F  
  #               NOTE: Find1, Replace1 and ext all must be quoted strings! 
  filenames=list.files(dnm,recursive=recursive)
  for( f in filenames ){
    if(!is.null(ext)){
      if(grepl(ext,x=f,fixed=T)){
        #Create a backup (if needed)
        if(backup){
          setwd(dnm)
          Sys.chmod(list.dirs("."), "774") #permissions needed to make the changes
          tmpdir=gsub(basename(paste0(dnm,"/",f)),"",paste0(dnm,"/",f),fixed=T)
          dir.create(paste0(tmpdir,"BackUp")) #create backup directory
          file.copy(f,paste0(tmpdir,"BackUp/",basename(f)))
        }
        x <- readLines(f)
        y <- gsub(Find1,Replace1, x ,fixed=T)
        cat(y, file=f, sep="\n")
      }  
    }
  }
}  

#Example
#dnm="/home/dhennen/Clam/Assessment Stuff/SS3/2016/runs/North/Explorations/TestReplaceFunc/"
#Find1=" 3 31 14.1093 10.3 -1 3 1 # SR_LN(R0)" 
#Replace1=" 3 31 14.1093 10.3 -1 3 -1 # SR_LN(R0)" 
#ext=".ctl"
#lineReplace(dnm,Find1,Replace1,ext)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

blockReplace=function(dnm,Find1,Replace1,ext=NULL,recursive=F,backup=T){
  #Find a block of text (Find1) and replace with another block (Replace1) in all files (with extension - ext or containing - ext) 
  #in directory (dnm) - this can be done recursively or not. You will by default also create backups
  #of any files that you alter, but you can turn that off with backup=F  
  #               NOTE: Find1, Replace1 and ext all must be quoted strings! 
  #******* WARNING: Watch out for automatic indentation - diable or use shift+tab to outdent!!!!!
  filenames=list.files(dnm,recursive=recursive)
  for( f in filenames ){
    if(!is.null(ext)){
      if(grepl(ext,x=f,fixed=T)){
        #Create a backup (if needed)
        if(backup){
          setwd(dnm)
          Sys.chmod(list.dirs("."), "774") #permissions needed to make the changes
          tmpdir=gsub(basename(paste0(dnm,"/",f)),"",paste0(dnm,"/",f),fixed=T)
          if(!dir.exists(paste0(tmpdir,"BackUp"))) dir.create(paste0(tmpdir,"BackUp")) #create backup directory
          file.copy(f,paste0(tmpdir,"BackUp/",basename(f)))
        }
        x=readChar(f, file.info(f)$size)
        #x <- readLines(f)
        y <- gsub(Find1,Replace1, x ,fixed=T)
        cat(y, file=f, sep="\n")
      }  
    }
  }
} 




#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#double normal selectivity from the sel.line function in
# the r2ss package, orginally written by Tommy Garrison.
#_______________________________________________________________________________
# The first argument (x) is a vector of sorted lengths or ages
# (midpoints).
#The second argument (sp) is a vector of six parameters in the
# same order and the same scale as in the SELEX24.XLSW and
# SS3 control files.  The vector contains six parameters:
#     PEAK=length/age at peak of dome
#     TOP=width of peak (length or age units)
#     ASC_WIDTH=width of ascending limb (log length or age units)
#     DSC_WIDTH=width of descending limb (log length or age units)
#     INIT=vertical offset at smallest size/age, units are
#          logits of proportions so 0 means offset to 0.5,
#          use a value <= -999 to omit from calculations
#     FINAL=vertical offset at largest size/age, details same as INIT
# NOTE-the order of the parameters in the sp vector is important and
#  the names for parameters mentioned above are not required.
sel <- function(x,sp) {
  sel <- rep(NA, length(x))
  startbin <- 1
  peak <- sp[1]
  upselex <- exp(sp[3])
  downselex <- exp(sp[4])
  final <- sp[6]
  if (sp[5] < -1000) {
    j1 <- -1001 - round(sp[5])
    sel[1:j1] <- 1e-06
  }
  if (sp[5] >= -1000) {
    j1 <- startbin - 1
    if (sp[5] > -999) {
      point1 <- 1/(1 + exp(-sp[5]))
      t1min <- exp(-(x[startbin] - peak)^2/upselex)
    }
  }
  if (sp[6] < -1000)  j2 <- -1000 - round(sp[6])
  if (sp[6] >= -1000) j2 <- length(x)
  peak2 <- peak + 2 + (0.99 * x[j2] - peak - 2)/(1 +exp(-sp[2]))
  if (sp[6] > -999) {
    point2 <- 1/(1 + exp(-final))
    t2min <- exp(-(x[j2] - peak2)^2/downselex)
  }
  t1 <- x - peak
  t2 <- x - peak2
  join1 <- 1/(1 + exp(-(20/(1 + abs(t1))) * t1))
  join2 <- 1/(1 + exp(-(20/(1 + abs(t2))) * t2))
  if (sp[5] > -999) asc <- point1 + (1 - point1) * (exp(-t1^2/upselex) - t1min)/(1 - t1min)
  if (sp[5] <= -999) asc <- exp(-t1^2/upselex)
  if (sp[6] > -999) dsc <- 1 + (point2 - 1) * (exp(-t2^2/downselex) - 1)/(t2min - 1)
  if (sp[6] <= -999) dsc <- exp(-(t2)^2/downselex)
  sel[(j1 + 1):j2] <- asc * (1 - join1) + join1 * (1 - join2 + dsc * join2)
  if (startbin > 1 && sp[5] >= -1000) sel[1:startbin] <- (x[1:startbin]/x[startbin])^2 *sel[startbin]
  if (j2 < length(x)) sel[(j2 + 1):length(x)] <- sel[j2]
  #normalize to max = 1
  sel=sel/max(sel)
  return(sel)
}
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MakeTexSymbols=function(){  #define some common latex terms
  assign("sosHead.tx"," \\\\textbf{State of Stock: }{}", envir = .GlobalEnv)
  assign("fig1.2.ref"," (Figures \\\\ref{SSB_plot1}-\\\\ref{F_plot1}){}", envir = .GlobalEnv)
  assign("fig1.ref"," Figure \\\\ref{SSB_plot1}{}", envir = .GlobalEnv)
  assign("fig2.ref"," Figure \\\\ref{F_plot1}{}", envir = .GlobalEnv)
  assign("fig3.ref"," Figure \\\\ref{Recr_plot1}{}", envir = .GlobalEnv)  
  assign("fig4.ref"," Figure \\\\ref{Fish_plot1}{}", envir = .GlobalEnv)
  assign("fig5.ref"," Figure \\\\ref{Surv_plot1}{}", envir = .GlobalEnv)
  assign("tab1.ref"," Table \\\\ref{Catch_Status_Table}{}", envir = .GlobalEnv)  
  assign("tab2.ref"," Table \\\\ref{BRP_Table}{}", envir = .GlobalEnv)
  assign("tab3.ref"," Table \\\\ref{Proj_Table}{}", envir = .GlobalEnv)  
  assign("SSBMSY.tx"," $SSB_{MSY}${}", envir = .GlobalEnv)
  assign("FMSY.tx", " $F_{MSY}${}", envir = .GlobalEnv)
  assign("FMSYproxy.tx", " $F_{MSY}${} \\\\textit{proxy}{}", envir = .GlobalEnv)
  assign("EMSY.tx", " $E_{MSY}${}", envir = .GlobalEnv)
  assign("EMSYproxy.tx", " $E_{MSY}${} \\\\textit{proxy}{}", envir = .GlobalEnv) 
  assign("RMSY.tx", " $R_{MSY}${}", envir = .GlobalEnv)
  assign("RMSYproxy.tx", " $R_{MSY}${} \\\\textit{proxy}{}", envir = .GlobalEnv)   
  assign("F30.tx", " $F_{30\\%}${}", envir = .GlobalEnv)
  assign("F35.tx", " $F_{35\\%}${}", envir = .GlobalEnv) 
  assign("F40.tx", " $F_{40\\%}${}", envir = .GlobalEnv) 
  assign("F45.tx", " $F_{45\\%}${}", envir = .GlobalEnv) 
  assign("F50.tx", " $F_{50\\%}${}", envir = .GlobalEnv)
  assign("E30.tx", " $E_{30\\%}${}", envir = .GlobalEnv)
  assign("E35.tx", " $E_{35\\%}${}", envir = .GlobalEnv) 
  assign("E40.tx", " $E_{40\\%}${}", envir = .GlobalEnv) 
  assign("E45.tx", " $E_{45\\%}${}", envir = .GlobalEnv) 
  assign("E50.tx", " $E_{50\\%}${}", envir = .GlobalEnv)  
  assign("FFull.tx", " $F_{Full}${}", envir = .GlobalEnv)
  assign("EFull.tx", " $E_{Full}${}", envir = .GlobalEnv) #fully selected exploitation rate
  assign("Ftarg.tx", " $F_{Target}${}", envir = .GlobalEnv) 
  assign("Fthresh.tx", " $F_{Threshold}${}", envir = .GlobalEnv)
  assign("Frebuild.tx", " $F_{Rebuild}${}", envir = .GlobalEnv)
  assign("F0.1.tx", " $F_{0.1}${}", envir = .GlobalEnv)
  assign("Fbar.tx", " $\\bar{F}${}", envir = .GlobalEnv)  
  assign("Etarg.tx", " $E_{Target}${}", envir = .GlobalEnv) 
  assign("Ethresh.tx", " $E_{Threshold}${}", envir = .GlobalEnv)
  assign("Erebuild.tx", " $E_{Rebuild}${}", envir = .GlobalEnv)
  assign("E0.1.tx", " $E_{0.1}${}", envir = .GlobalEnv)
  assign("Ebar.tx", " $\\bar{E}${}", envir = .GlobalEnv)    
  assign("SSBbar.tx", " $\\bar{SSB}${}", envir = .GlobalEnv)  
  assign("SSBMSY.tx", " $SSB_{MSY}${}", envir = .GlobalEnv)
  assign("SSB0.tx", " $SSB_{0}${}", envir = .GlobalEnv)
  assign("SSBMSYproxy.tx", " $SSB_{MSY}${} \\\\textit{proxy}{}", envir = .GlobalEnv)  
  assign("SSBtarg.tx", " $SSB_{Target}${}", envir = .GlobalEnv) 
  assign("SSBthresh.tx", " $SSB_{Threshold}${}", envir = .GlobalEnv)
  assign("BMSY.tx", " $B_{MSY}${}", envir = .GlobalEnv)
  assign("BMSYproxy.tx", " $B_{MSY}${} \\\\textit{proxy}{}", envir = .GlobalEnv)  
  assign("Btarg.tx", " $B_{Target}${}", envir = .GlobalEnv) 
  assign("Bthresh.tx", " $B_{Threshold}${}", envir = .GlobalEnv)
  assign("F.Status.tx", " \\\\textit{Overfishing}{}", envir = .GlobalEnv)  
  assign("B.Status.tx", " \\\\textit{Overfished}{}", envir = .GlobalEnv)  
  assign("half.tx", " $\\\\dfrac{1}{2}${}", envir = .GlobalEnv)   
  assign("Fratio.tx", " $\\\\frac{F}{F_{Threshold}}${}", envir = .GlobalEnv)  
  assign("Bratio.tx", " $\\\\frac{B}{B_{Threshold}}${}", envir = .GlobalEnv)  
  assign("SSBratio.tx", " $\\\\frac{SSB}{SSB_{Threshold}}${}", envir = .GlobalEnv)  
  assign("ProjHead.tx"," \\\\textbf{Projections: }{}", envir = .GlobalEnv)  
  assign("SpecComHead.tx"," \\\\textbf{Special Comments: } \\\\begin{itemize}{}", envir = .GlobalEnv)  
  assign("item.tx"," \\\\item{}", envir = .GlobalEnv)
  assign("beginitem.tx"," \\\\begin{itemize}{}", envir = .GlobalEnv)    
  assign("enditem.tx"," \\\\end{itemize}{}", envir = .GlobalEnv)  
  assign("RefHead.tx"," \\\\textbf{References: }{}", envir = .GlobalEnv)
  assign("lbreak.tx",  " \\\\linebreak{}", envir = .GlobalEnv)
  assign("indent.tx",  " \\\\hspace*{0.5cm}", envir = .GlobalEnv)
  assign("rho.tx",  " \\\\textrho{}", envir = .GlobalEnv)
  assign("alpha.tx",  " \\\\textalpha{}", envir = .GlobalEnv)
  assign("beta.tx",  " \\\\textbeta{}", envir = .GlobalEnv)  
  assign("gamma.tx",  " \\\\textgamma{}", envir = .GlobalEnv)
  assign("delta.tx",  " \\\\textdelta{}", envir = .GlobalEnv)
  assign("epsilon.tx",  " \\\\textepsilon{}", envir = .GlobalEnv)
  assign("mu.tx",  " \\\\textmu{}", envir = .GlobalEnv)
  assign("lambda.tx",  " \\\\textlambda{}", envir = .GlobalEnv)
  assign("sigma.tx",  " \\\\textsigma{}", envir = .GlobalEnv) 
  assign("lt.tx",  " \\\\textless{}", envir = .GlobalEnv) 
  assign("gt.tx",  " \\\\textgreater{}", envir = .GlobalEnv)  
  assign("le.tx",  " $\\\\leq${}", envir = .GlobalEnv) 
  assign("ge.tx",  " $\\\\geq${}", envir = .GlobalEnv)
  assign("rhoSSB.tx",  " $SSB_{\\\\rho}${}", envir = .GlobalEnv) 
  assign("rhoB.tx",  " $B_{\\\\rho}${}", envir = .GlobalEnv)  
  assign("rhoF.tx",  " $F_{\\\\rho}${}", envir = .GlobalEnv)
  assign("PortalLink.tx", "\\\\href{http://www.nefsc.noaa.gov/saw/sasi/sasi_report_options.php}{SASINF}{}" ,envir = .GlobalEnv)
  assign("RhoDecisionTab.ref"," Table \\\\ref{RhoDecision_tab}{}", envir = .GlobalEnv)
}
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
color.column2 <- function(xt, column, col='gray80') {
  # color.column - color the columns of a table for LaTeX output.
  #xt is an xtable object, column is the column name you want to color, col is a color
  # Make sure the colorspace library is loaded to convert rgb to hex.
  clr.string <- (col2hex(col))
  clr.string <- substring(clr.string, 2) # Remove the hash prefix.
  
  indx <- which(colnames(xt) == column)
  # Grab the digits and display values set by xtable
  digt <- attr(xt, 'digits')[indx+1]
  disply <- attr(xt, 'display')[indx+1]
  
  formatted.vals <- formatC(paste(xt[, column]), digits=digt, format=disply)
  
  formatted.column <- paste('\\cellcolor[HTML]{', clr.string, '}', formatted.vals, sep='')
  xt[, column] <- formatted.column
  xt
}
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
RandomWalk=function(n=1000,min=-1,max=1,slope=0,seed=NA){
  #Generate a random walk with steps between min and max. Add an optional offset with slope if you like.
  if(!is.na(seed)) set.seed(seed)
  return(cumsum(sample(c(min,max), n, TRUE)+slope))
}
#plot(RandomWalk(slope=0.02),type="l",lwd=2)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
RandomWalkIndex=function(n=1000,drift=1,sd=1){
  #Generate a random walk with mean = drift and variance = sd^2
  cumsum(rnorm(n, mean=drift, sd=sd))
}
#plot(RandomWalkIndex(drift=0.04),type="l",lwd=2)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ProjectionStrings=function(){
#add a projection string! These are from Burton
#CRS("+proj=utm +zone=18 +ellps=GRS80 +datum=NAD83 +units=m +no_defs") #this is the midatlantic
#CRS("+proj=utm +zone=19 +ellps=GRS80 +datum=NAD83 +units=m +no_defs") #this is GBK
#CRS("+proj=longlat +ellps=GRS80 +datum=NAD83 +no_defs")
}
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#Convert from old strata to new strata - requires some additional arguments - see bleow:
StratumConvert=function(dt1,svspp=403){
  #dt1 is a data.frame that has at least the columns "Strata2" (Strata with the 6 and 0 removed), "Lat", "Lon", "Depth", "CRUISE6"
  #svspp is the species code for surfclams (403) or quahogs (409)
  stratum2=dt1$Strata2 #Need to make sure we don't change the strata %in% 2018;
  stratum=paste(stratum2)
  lat=dt1$Lat
  lon=dt1$Lon
  depth=na.replace(dt1$Depth)
  cruise6=dt1$CRUISE6
  
  # CONVERT LATITUDE & LONGITUDE TO DECIMAL DEGREES;
  #-- latitude first; 
  degrees=floor(lat/100);
  minutes=(lat-degrees*100)/60;
  lat=degrees+minutes;
  #-- longitude too; 
  degrees=floor(lon/100);
  minutes=(lon-degrees*100)/60;
  lon=degrees+minutes;
  lon=-1*lon	
  
  #The formula for determ%in%%in%g if a po%in%t is one side or the other of a line is: ;
  # d=(lon-x[1])*(y[2]-y[1])-(lat-y[1])*(x[2]-x[1]) ;
  for(i in 1:nrow(dt1)) {
    if(stratum[i] %in% c('47') ) { #Need to determine the split of stratum[i] 47;
      if( ( (lon[i]-69.23)*(41-40)-(lat[i]-40)*(69.03-69.23) ) > 0 ) { stratum[i] = '471';		
      } else {  stratum[i] = '472'}
    }
    
    if(stratum[i] %in% c('73') ) { #Need to determine the split of stratum[i] 73;
      if( ( (lon[i]-66.8)*(41.9-41.35)-(lat[i]-41.35)*(67.5-66.8) ) > 0) { stratum[i] = '73';
      } else {  stratum[i] = '74'}		
    }
    
    #The Hudson canyon quahog split is complicated!;
    #This is a line along the bottom of Hudson Canyon;
    if(svspp %in% c(409) & stratum[i] %in% c('25', '26') & lat[i] >= 39.3 & lat[i]<= 40.2 ) { #Need to determine the split of strata 26 27 28;
      if( ( (lon[i]-72)*(40.2-39.3)-(lat[i]-39.3)*(73.75-72) ) < 0) { #above the line;
        if (stratum[i] %in% c('26')) { stratum[i] = '30' 
        } else { if(stratum[i] %in% c('25')) { stratum[i] = '29'}	}
      }
    }
    #This is the same line along the bottom of Hudson Canyon; 
    if(svspp %in% c(409) & stratum[i] %in% c( '31', '32') ) { #Need to determine the split of strata 31 32;
      if( ( (lon[i]-72)*(40.2-39.3)-(lat[i]-39.3)*(73.75-72) ) < 0) { #below the line;
        if(stratum[i] %in% c('31')) { stratum[i] = '27'
        } else { if(stratum[i] %in% c('32')) { stratum[i] = '28'} }
      }		
    }
    # This is a t%in%y line segment that is completely conta%in%ed %in% stratum[i] 25 & 26;
    if(svspp %in% c(409) & stratum[i] %in% c('25', '26') & lat[i] >= 40.2 & lat[i]<= 40.25) { #Need to determine the split of stratum[i] 26;
      if( ( (lon[i]-73.75)*(40.25-40.2)-(lat[i]-40.25)*(73.775-73.75) ) < 0) { #right of the line;
        if(stratum[i] %in% c('25')) { stratum[i] = '29'
        } else { if(stratum[i] %in% c('26')) { stratum[i] = '30'} }
      }	
    }
    # This is a short line segment that is completely conta%in%ed %in% stratum[i] 25 & 26;
    if(svspp %in% c(409) & stratum[i] %in% c('25', '26') & lat[i] >= 40.25 & lat[i]<= 40.5) { #Need to determine the split of stratum[i] 26;
      if( ( (lon[i]-73.775)*(40.5-40.25)-(lat[i]-40.25)*(73.825-73.775) ) < 0) { #right of the line;
        if(stratum[i] %in% c('25')) { stratum[i] = '29'
        } else { if(stratum[i] %in% c('26')) { stratum[i] = '30'} }
      }		
    }
    
    #Finally we have to split the southernmost big strata for oq only;
    if(svspp %in% c(409) & stratum[i] %in% c('17') ) { #Need to determine the split of strata 17;
      if( ( (lon[i]-74.29)*(38.6-38.94)-(lat[i]-38.94)*(74.57-74.29) ) < 0) { stratum[i] = '0'} #left of the line;	
    }
    if(svspp %in% c(409) & stratum[i] %in% c('13') & lat[i] >= 38.41) { #Need to determine the split of stratum[i] 13;
      if( ( (lon[i]-74.57)*(38.41-38.6)-(lat[i]-38.6)*(74.64-74.57) ) < 0) { stratum[i] = '0'} #left of the line;	
    }
    if(svspp %in% c(409) & stratum[i] %in% c('13') & lat[i] >= 38.15 & lat[i]<= 38.41) { #Need to determine the split of stratum[i] 13;
      if( ( (lon[i]-74.64)*(38.15-38.41)-(lat[i]-38.41)*(74.67-74.64) ) < 0) { stratum[i] = '0'} #left of the line;	
    }
    if(svspp %in% c(409) & stratum[i] %in% c('13') & lat[i]<= 38.15) { #Need to determine the split of stratum[i] 13;
      if( ( (lon[i]-74.67)*(37.83-38.15)-(lat[i]-38.15)*(74.87-74.67) ) < 0) { stratum[i] = '0'} #left of the line;	
    }
    
    #Now we can assign the new strata by the old strata;	
    if (svspp %in% c(403) & (cruise6[i]/100)<2018) {		 #Surfclam		
      if (stratum[i] %in% c('05','09','81') ) { stratum2[i] = "1S"} 
      if (stratum[i] %in% c('84', '85', '86', '87')) { stratum2[i] = "2S"}  #NJ DMV %in%shore
      if (stratum[i] %in% c('13', '17', '21', '25', '29')) { stratum2[i] = "3S"}   #DMV, NJ & LI mid depth
      if (stratum[i] %in% c('10', '14', '18', '22')) { stratum2[i] = "4S"}   #NJ, DMV & LI off shore
      if (stratum[i] %in% c('88', '89', '90', '91', '92', '93')) { stratum2[i] = "5S"}   #LI & NJ %in%shore
      if (stratum[i] %in% c('45', '46', '95', '96')) { stratum2[i] = "6S"}  #SNE
      #GBK-surfs
      if (stratum[i] %in% c('53', '54')) { stratum2[i] = "7S" }  #NW deep
      if (stratum[i] %in% c('67', '69', '70')) { stratum2[i] = "8S" }  #West shallow
      if (stratum[i] %in% c('57', '58', '59', '60')) { stratum2[i] = "9S" }  #S deep
      if (stratum[i] %in% c('65', '66')) { stratum2[i] = "10S" }  #NE deep
      if (stratum[i] %in% c('68', '72', '73')) { stratum2[i] = "11S" }  #middle shallow
      if (stratum[i] %in% c('71', '74')) { stratum2[i] = "12S" }  #middle east shallow
      if (!stratum[i] %in% c('05','09','81','84', '85', '86', '87','13', '17', '21', '25', '29'
                             ,'10', '14', '18', '22','88', '89', '90', '91', '92', '93','45', '46', '95', '96'
                             ,'53', '54','67', '69', '70','57', '58', '59', '60','65', '66','68', '72', '73'
                             ,'71', '74')) {stratum2[i] = '0'}
    }
    
    if (svspp %in% c(409) & (cruise6[i]/100)<2018) { #Quahog
      if (stratum[i] %in% c('10', '11', '12', '14', '15', '16', '18', '19', '20')) { stratum2[i] = "1Q" }  #DMV deep
      if (stratum[i] %in% c('13', '17', '21', '22', '25', '26')) { stratum2[i] = "2Q"}   #NJ mid depth
      if (stratum[i] %in% c('23', '24', '27', '28', '31', '32', '35', '36'))	{ stratum2[i] = "3Q" }  #LI & NJ deep
      if (stratum[i] %in% c('29', '30', '33', '34')) { stratum2[i] = "4Q" }  #NJ & LI mid
      if (stratum[i] %in% c('92', '93', '94', '95', '37', '41')) { stratum2[i] = "5Q"  } #DMV/SVA
      if (stratum[i] %in% c('38', '39', '40', '46', '471', '48')) { stratum2[i] = "6Q" }  #NJ DMV %in%shore
      
      #GBK-Quahogs
      if (stratum[i] %in% c('53','54','55','56','472','56')) { stratum2[i] = "7Q"  } #West side deep
      if (stratum[i] %in% c('70')) { stratum2[i] = "8Q" }  #middle shallow
      if (stratum[i] %in% c('57','58','59','60')) { stratum2[i] = "9Q" }  #S deep
      if (stratum[i] %in% c('65','66')) { stratum2[i] = "10Q" }  #NE deep
      if (stratum[i] %in% c('74')) { stratum2[i] = "11Q" } #East side mid-depth
      if (stratum[i] %in% c('61','62')) { stratum2[i] = "12Q"  } #SE deep		
      if(!stratum[i] %in% c('10', '11', '12', '14', '15', '16', '18', '19', '20','13', '17', '21', '22', '25', '26'
                            ,'23', '24', '27', '28', '31', '32', '35', '36', '29', '30', '33', '34','92', '93'
                            , '94', '95', '37', '41','38', '39', '40', '46', '471', '48', '53','54','55','56','472','56'
                            ,'70', '57','58','59','60','65','66','74','61','62') ) { stratum2[i] = '0'}
    }
    if(depth[i]>80) { stratum2[i] = '0'}
  }
  return(stratum2)
  
}




