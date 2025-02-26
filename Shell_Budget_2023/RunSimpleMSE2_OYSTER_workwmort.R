################################################################################
#                                                                              #
#   SPAWN project - simulate a population with specific dynamics               #
#  DRH and JJD 4/24/13                                                    #
#______________________________________________________________________________#
#Remove all objects that currently exist in memory
rm(list=ls(all=TRUE))
#source functions
require(microbenchmark)
library(data.table)
library(reshape2)
library(dplyr)
library(ggplot2)
library(plotly)
library(metR)
library(latticeExtra)
library(RColorBrewer)
library(lattice)
library(extrafont)
library(cowplot)
#setwd("C:/Users/Admin/Documents/Presentation/OysterSeminarPresentationBase/Oyster_MSE")
wd <- "C:/Users/laura.solinger/Documents/Publications/OysterManuscript2_CultchSeed/Oyster_MSE/"
#______________________________________________________________________________
source(paste0(wd,"functions.L.r"))
source(paste0(wd,"SPAWNSC2_OYSTER_workwmortdermomat.r"))
source(paste0(wd,"SRSC1_OYSTER.r"))
load(paste0(wd,"ShellRockDEBay.RData"))
load(paste0(wd,"HighMortDEBay.RData"))
load(paste0(wd,"MidMortDEBay.RData"))

#source("./SRSC1_OYSTER_workwmortdermomatWORKS.r")
dir.graphics<-wd
dir.data<-wd
#______________________________________________________________________________
# Rscript ./Surfclam/MSE/Rcode/RunSimpleMSE2.r &
#______________________________________________________________________________

### specify which oyster bed you want and the name of the output file

  st <- ShellRock ## Specify which bed you want to model
  outputname <- "ShellRockSims" ### Note, if this file already exists, will add additional simulations onto existing CSV
  
  fold <- TRUE
  
  if(fold){
####
  nmonths<-600    #round(runif(1,20,50),0) #months of simulated data
  burnin<-st$burnin
  nages<-st$nages               #round(runif(1,200,300),0) #ages
  len<-st$len           #
  Mat<-st$Mat #vector of maturity at age
  Wt<-st$Wt #fill.weight(nages) #vector of weight at age (NEFSC 2013) - convert g to 1000 mt
  med.recr<-st$med.recr #median recruitment for constant or purely random recruitment
  rvar<-st$rvar #standard deviation of the recruitment error parameter
  Rscale<-st$Rscale   #reduce overall recruitment scale to 1/4 of what it was!
  type<-st$type #oq type recruitment - see documentation for this project or comments in Spawn code   
  steepness<-st$steepness #this should be a random variable to test....
  hlo<-st$hlo #min steepness
  M<-st$M         #runif(1,0.1,0.5) #natural mortality
  D.B<-st$D.B #disarticulation of boxes
  D.C<-st$D.C #decay of cultch
  C.C<-st$C.C #correction factor of cultch
  C.L<-st$C.L #correction factor for live
  C.B<-st$C.B #correction factor for boxes
  shell.plant<-st$shell.plant # shell plant
  a.wt<-st$a.wt # a in length width conversion
  b.wt<-st$b.wt # b in length width conversion
  wt.1<-st$wt.1 # width at age by length
  saatage<-st$saatage # surface area at age in hectares
  rec.a.sa<-st$rec.a.sa
  rec.b.sa<-st$rec.b.sa
  shell.correct<-st$shell.correct
  Mlo<-st$Mlo
  Mhi=st$Mhi
  F.dat<-st$F.dat #vector of F in each month
  F.dat<-st$D.dat
  F.dat<-st$J.dat #seed fishery (J)
  catch<-st$catch #vector of catch (weight) in each month
  selectivity.I<-st$selectivity.I #selx for index 
  I.error<-st$I.error
  selectivity.F<-st$selectivity.F #selx for fishery from Thorarinsdottir and Jacobson, 2005.     
  Selectivity.D<-st$Selectivity.D ##DERMO SELECTIVITY
  selectivity.J<-st$selectivity.J #Seed selectivity
  F.error<-st$F.error
  D.error<-st$D.error
  J.error<-st$J.error
  catchab<-st$catchab
  Nzero<-st$Nzero    #Reduce starting population to 1/4 total!
  iseed2<-get.seed()
  h_bound<-st$h_bound
  F_bound<-st$F_bound
  #spawn a population
  ##### set up the recovery time specific variables
  Nstar_vec<-st$Nstar_vec  
  Ndelta_vec<-st$Ndelta_vec
  sigk<-st$sigk #stdev in growth parameter k
  sigL<-st$sigL #stdev in growth parameter Linfinity
  sigt0<-st$sigt0
  interval <- st$interval


  iseed2<-get.seed()
  #spawn a population
  #st.pop.matrix<-read.csv(paste(dir.data,"st_pop_sample.csv",sep=""))
}

  
#######################################################################################################
#set up done - run the sims...
#######################################################################################################
  # Start the clock!
  #ptm <- proc.time()
  rez<-c()
  nmonths<-600+burnin
  burnin2<-400 #first 200 years of seed fishery, dermo begins at 400
  update<-100
  update2<-update
  ## create a matrix to save surface area data into
  SA.data <- matrix(data = c(1,1,1,1),nrow=1)
  Bed <- "SR" #### Change to SR, MM, or HM, so that the correct shell planting parm is used
  ### Start with a smaller number than dan just to test it out
  ## n<-(50001)   
  n<-3
  
  ploteach <- T
  printsimdat <- F
  
  for (i in 2:n){
    print(i)
    #i<-2
    if(i==update){print(i);update<-update+update2}
    
    extink<-0
    extink100<-0
    ## These parameters will be the same for all regions
    #F.target<-runif(1,0.001,0.2) 
    F.target<-(runif(1,0.0,0.30))/12 #test a set value
    D.target<-(runif(1,0.0,0.25))/12 
    J.target<-(runif(1,0.0,0.25))/12 #seed
    
    if(Bed == "MM"){
    #shellplant MM
    planted.hectares<-runif(1, 474.68, 3433.49)
    shell.plant<-c(planted.hectares, rep(0, times=239))
    shell.plant.mat<-matrix(NA,nrow=800, ncol=nages)
    for(i in c(1:800)){
      shell.plant.mat[i,]<-rep(0,240)
    }
    for(i in seq(1,800,by=12)){
      shell.plant.mat[i,]<-shell.plant
    }
    shell.plant<-shell.plant.mat
    }
    
    if(Bed == "HM"){
      #shellplant HM
      planted.hectares<-runif(1, 890.49, 64441.22)
      shell.plant<-c(planted.hectares, rep(0, times=239))
      shell.plant.mat<-matrix(NA,nrow=800, ncol=nages)
      for(i in c(1:800)){
        shell.plant.mat[i,]<-rep(0,240)
      }
      for(i in seq(1,800,by=12)){
        shell.plant.mat[i,]<-shell.plant
      }
      shell.plant<-shell.plant.mat
    }
    
    if(Bed == "SR"){
      #shellplant SR
      planted.hectares<-runif(1, 171.44, 1240.10)
      shell.plant<-c(planted.hectares, rep(0, times=239))
      shell.plant.mat<-matrix(NA,nrow=800, ncol=nages)
      for(i in c(1:800)){
        shell.plant.mat[i,]<-rep(0,240)
      }
      for(i in seq(1,800,by=12)){
        shell.plant.mat[i,]<-shell.plant
      }
      shell.plant<-shell.plant.mat
    }
    
        
    F.dat<-rep(F.target,nmonths)
    D.dat<-c(rep(0,burnin2),rep(D.target,(nmonths-burnin2))) #account for burn in period
    J.dat<-c(rep(0,burnin),rep(J.target,(nmonths-burnin))) #seed fishery
    
    #F.error<-runif(1,0.0,0.2*F.target)  #implementation error ~cv of 0.5
    F.error<-runif(1,0.0,0*F.target) 
    D.error<-runif(1,0.0,0*D.target) 
    J.error<-runif(1,0.0,0*J.target)
    
    ac<-runif(1,0,0) #degree of autocorrelation in management error 
    man.error<-runif(1,0,0) #degree of management error (assessment error) 
    ass.int<- 3
    #choose a biomass threshold reference point
    pntr_Nstar<-round(runif(1,1,length(Nstar_vec)),0)
    Nstar<-Nstar_vec[pntr_Nstar] 
    #choose a biomass target reference point
    pntr_Ndelta<-round(runif(1,1,length(Ndelta_vec)),0)
    Ndelta<-Ndelta_vec[pntr_Ndelta]  
    Nthresh<-Nzero*Nstar #calculate B theshold as a percentage of virgin B (for the 4 regions)
    Ntarg<-Nthresh+(Ndelta*Nzero)  #The B target is then a function of the threshold 
    #these reference points are sim specific but will hold for all 4 populations  
    ## These things should be different for each region
    ##############################
    ## Just need one population for this run            ##
    ##############################
    # M<-runif(1,Mlo,Mhi)
    steepness<-runif(1,hlo,h_bound)  #reduced range of steepnesses based on He et al 2006    
    # len<-change.growthS(nages=nages) #SVAtoSNE Pars are in function as default
    #selectivity.F<-change.selx.F(len) #need to recalculate selx at age with a change in growth
    #Selectivity.D<-change.selx.D(len)
    begin.sa<-sa.func(nmonths=1000,nages=nages,M=M,a.wt=a.wt,b.wt=b.wt,Mat=Mat,N0=Nzero,D.B=D.B,
                      D.C=D.C,C.L=C.L,C.C=C.C,C.B=C.B,steep=steepness,Reca=rec.a.sa,Recb=rec.b.sa,printsimdat = printsimdat)
    SA.1<-begin.sa$SA1.1
    st.pop2<-c()  
    st.pop2<-stock.recruit(SSN1=Nzero,M=M,Wt=Wt,Mat=Mat,steep=steepness
                           ,N0=Nzero,type=type,Rscale=Rscale,SA.1=SA.1,rec.a.sa=rec.a.sa,rec.b.sa=rec.b.sa)  
    med.recr<-1
    med.recr.in<-1
    st.pop2<-round(st.pop2,digits=0) ### SO IT'S TREATED AS NUMBERS. IT GOES TO DECIMAL PLACES BECAUSE OF MORT ETC
    recs<-miniSpawn(st.pop2,nages,Mat,Wt,M,F.dat,D.dat,J.dat,selectivity.F,Selectivity.D,selectivity.J)
    recs<-round(recs,digits=0)
    st.pop<-as.vector((recs)/(sum(recs))*(Nzero)) ## REMOVED WT IN SUM(RECS*WEIGHT)
    st.pop<-round(st.pop,digits=0)
    #this is now equilibrium starting population in numbers - this requires a burn in period without F!
    ## Populations will be managed separately
    #iseed2<-get.seed(),seed1=iseed2
    type=ifelse(runif(1)<0.5,1,1) #pick a SR curve type
    
    
    yt0<-Spawn(nmonths=nmonths,nages=nages,Mat=Mat,Wt=Wt,Recsd=rvar
               ,steepness=steepness,Nzero=round(Nzero,digits=0),M=M,D.B=D.B,D.C=D.C,C.C=C.C,C.L=C.L,C.B=C.B,a.wt=a.wt,b.wt=b.wt,wt.1=wt.1,saatage=saatage
               ,shell.plant=shell.plant,rec.a.sa=rec.a.sa,rec.b.sa=rec.b.sa,F.dat=F.dat,D.dat=D.dat,J.dat=J.dat,catchab=catchab
               ,st.pop=st.pop,selectivity.I=selectivity.I,selectivity.F=selectivity.F,Selectivity.D=Selectivity.D,selectivity.J=selectivity.J
               ,F.error=F.error,D.error=D.error,J.error=J.error,I.error=I.error,type=type,Plot1=F,Unfished=F
               ,Nthresh=round(Nthresh,digits=0),Ntarg=round(Ntarg,digits=0),ac=ac,man.error=man.error,Rscale=Rscale
               ,ass.int=ass.int,F.target=F.target,D.target=D.target,J.target=J.target)  
    
    
     st.pop.print<-data.frame("Recruits"=melt(yt0$Simulated.population$Recruits),
                              "SSN"=melt(yt0$Simulated.population$SSN))
           
     # 
     if(i==1) {
       app<-F
       c.name<-T
     } else {
       app<-T
       c.name<-F
     }	
    # #might as well always write out the results, just in case something goes wrong!
     output<-T  #set to true for writing each simulation result to a file 
     if(output){
       write.table(st.pop.print,file=paste(dir.graphics,outputname,"_startpop.csv",sep="")
                   ,sep=","
                   ,row.names=F
                   ,col.names=c.name
                   ,append=app
                   ,quote=F
       )
       
       #   
       #   
     } 
     age.try <- yt0$Age.comp
     age1 <- dplyr::filter(age.try,Month == max(age.try$Month))
     age1 <- dplyr::select(age1, -Month)
     age1 <- melt(age1)  
     age1$age <- c(1:240)
     age1.2<-dplyr::filter(age1,age > 33)
     prop <- sum(age1.2$value)/sum(age1$value)
     
     surfaceareadat<-data.frame("Simulation"=i,
                               "Month"=yt0$surfacearea$month,
                             "BSA"=yt0$surfacearea$BSA,
                             "Plant"=yt0$surfacearea$Plant,
                             "InitPlant"=shell.plant[1,1],
                             "CSA"=yt0$surfacearea$CSA,
                             "LSA"=yt0$surfacearea$LSA,
                             "Mort"=yt0$surfacearea$Mort,
                             "Seed"=yt0$Simulated.population$J,
                             "Dermo"=yt0$Simulated.population$D,
                             "F.target"=yt0$surfacearea$F.target,
                             "F.seed"=yt0$Simulated.population$J,
                             "NAA"=yt0$surfacearea$NAA,
                             "Recruits"=yt0$surfacearea$Recruits,
                             "Catch"=yt0$Simulated.population$Catch,
                             "Catch_N"=yt0$Simulated.population$Catch.N,
                             "Ratio"=yt0$surfacearea$LSA/(yt0$surfacearea$LSA+yt0$surfacearea$CSA+yt0$surfacearea$BSA),
                             "PropOver63"=prop)
    
    # Some Quick Plots of Numbers and Surface Area
     if(ploteach == T){
    par(mfrow=c(2,4))
    title<-paste("Dermo",(round(surfaceareadat$Dermo[202],digits=3)*12), ", F",(round(surfaceareadat$F.target,digits=3)*12),", Seed",(round(surfaceareadat$Seed[202],digits=3)*12),", Plant",round(surfaceareadat$InitPlant,digits = 3))
    title<-title[1]
    plot(yt0$surfacearea$NAA,main="NAA")
    plot(yt0$surfacearea$LSA,main="LSA")
    plot(yt0$surfacearea$Plant,main="Plant")
    plot(yt0$surfacearea$BSA,main="BSA")
    plot(yt0$surfacearea$CSA,main="CSA")
    #plot(surfaceareadat$Ratio,main="Ratio") #turned off to make space for plant
    plot(surfaceareadat$Catch_N ~ yt0$surfacearea$LSA,main="Catch v LSA")
    plot(surfaceareadat$Catch_N ~ yt0$surfacearea$BSA,main="Catch v BSA")
    plot(surfaceareadat$Catch_N ~ yt0$surfacearea$CSA,main="Catch v CSA")
    mtext(title[1],side = 3, line = -1.5, outer = T)
     }
    
    ## reset plotting pars
    par(mfrow=c(1,1))
    
     if(i==1) {
      app<-F
      c.name<-T
    } else {
      app<-T
      c.name<-F
    }	
    #c.name<-c("Month","BSA","CSA","LSA","Mort","Dermo","F.target","NAA","Recruits")
    #might as well always write out the results, just in case something goes wrong!
    output<-T  #set to true for writing each simulation result to a file 
    if(output){
      write.table(surfaceareadat,file=paste(dir.graphics,outputname,".csv",sep="")
                  ,sep=","
                  ,row.names=F
                  ,col.names=c.name
                  ,append=app
                  ,quote=F
      )
      
      #   
      #   
    } 
  }

  
  
#### Generate Heat Maps
  ## First- get meter squared based on which reef you're using 
  MM_MeterSquared <- 1.6475090e7
  SR_MeterSquared <- 0.5140455e7
  HM_MeterSquared <- 2.6702761e7
  
  Bed_Area <- SR_MeterSquared ## select which bed you ran, to give per square meter values
  
  data <- read.csv(paste(dir.graphics,outputname,".csv",sep=""), header = F)
  colnames(data)<-c("Simulation","Month","BSA","Plant","Init_Plant","CSA","LSA","Mort","Seed","Dermo","F.target","F.seed","NAA","Recruits","Catch","Catch_N","Ratio","PropOver63")
  firstdat.data<-dplyr::filter(data,Month==1) ## so you can see how many simulations you ran
  nrow(firstdat.data)
  
  data$Dermo <- data$Dermo * 12 ## Convert back to annual % loss
  data$F.target <- data$F.target *12 ## Convert back to annual % loss
  data$Sum_F<-data$Dermo+data$F.target # Sum Dermo and Fishing mortality
  
  ## Get KG cultch per square meters 
  data <- dplyr::filter(data,Month == 800) # how many simulations made it to the terminal Month
  data$CSA_m2<-data$CSA*10000 #convert cultch from hectares to square meters
  data$CSA.vol <- data$CSA_m2/0.111 # convert cultch ESA from m2 to Liters using coversion 1L = 0.111m2
  data$CSA_kg<-data$CSA.vol * (0.696/1) #convert to kg using 0.696kg per 1 liter
  data$CSA_kg_m2 <- data$CSA_kg / Bed_Area #size of in meter squared
  
  ## Number of live animals per square meter (note, includes ALL animals, including spat- for a rough calculation of animals >63mm, multiply NAA by 0.27)
  data$livenpersqm<-((data$NAA)/Bed_Area)

  a<-levelplot(Catch_N~F.target + Dermo, data,panel=panel.levelplot.points,cex=1, pch = 16,
               col.regions=heat.colors(25),main=list("Catch (#s)"),
               xlab=list("Fishing Mortality"),
               ylab=list("Dermo Mortality"),
               par.settings=list(axis.text=list()))+
    layer_(panel.2dsmoother(...))
  
  ## plot of kg / m2
  
  b<-levelplot(CSA_kg_m2~F.target + Dermo, data,panel=panel.levelplot.points,cex=1, pch = 16,
               col.regions=heat.colors(25),main=list("KG of Cultch / Square Meter"),
               xlab=list("Fishing Mortality"),
               ylab=list("Dermo Mortality"),
               par.settings=list(axis.text=list()))+
    layer_(panel.2dsmoother(...))
  
  ## plot of N oysters / m2
  
  c<-levelplot(livenpersqm~F.target + Dermo, data,panel=panel.levelplot.points,cex=1,pch = 16,
               col.regions=heat.colors(25),main=list("# Live Oysters > 3in / Square Meter"),
               xlab=list("Fishing Mortality"),
               ylab=list("Dermo Mortality"),
               par.settings=list(axis.text=list()))+
    layer_(panel.2dsmoother(...))
  


  d <- levelplot(Ratio~F.target + Dermo, data,panel=panel.levelplot.points,cex=1,pch = 16,
               col.regions=heat.colors(25),main=list("Live Ratio at Medium Mortality"),
               xlab=list("Fishing Mortality",cex=1.2),
               ylab=list("Dermo Mortality",cex=1.2))+
    layer_(panel.2dsmoother(...))
  
  e <- levelplot(Catch_N~(F.seed *12) + (Dermo + F.target), data,panel=panel.levelplot.points,cex=1,pch = 16,
                 col.regions=heat.colors(25),main=list("Catch at Medium Mortality"),
                 xlab=list("Seed Fishing Mortality",cex=1.2),
                 ylab=list("Fishing and Dermo Mortality",cex=1.2))+
    layer_(panel.2dsmoother(...))
  
  e  
  plot_grid(a,b,c,d,nrow=2)
  
  