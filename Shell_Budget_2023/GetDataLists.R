## Set up the reef parameters

## Medium Mortality 
setupmmrun <- function(){  
  ## First step is to draw population parameters from a random distribution 
  ## or set them up on your own if you prefer...
  nmonths<-600    # Months of simulated population 
  nages<-240        # Number of Ages- 240 Months specified here, eg an oyster can live to 20 years
  
  #### GENERATE OR SPECIFY LENGTH AT AGE 
  len<- vector(length=nages) # create vector for forloop to insert lengths-at-age to
  for(i in c(1:nages)){
    age<-i/12 # Von bert length at age parameterized for years, so divide i timestep by 12 to represent months 
    len[i]<-125*(1-exp((-0.24)*((age)+.2))) # parameters for Shell Rock from Kraeter et al. 
  }
  
  
  #### GENERATE OR SPECIFY MATURITY AT AGE/LENGTH
  len.1<-as.matrix(len) #create matrix of lengths 
  v<-as.matrix(x=c(1:nages),nrow=nages,ncol=1) # create matrix with dimensions for ages and months- so you can specify
  ## maturity at age in each Month. Useful if you think maturity may vary over time.
  for (i in c(1:nages)){
    x<-len.1[i,] 
    if(x<35){ ## so here determine that if length is < 35mm, oyster is not mature 
      y<-0
    }
    else {
      y<-1/12 # then say that all oysters larger than 35mm are mature. 
      # Divide maturity by 12 to account for monthly recruitment rather than annual
    }
    v[i,]<-y
  }
  Mat<-as.numeric(v)
  
  
  #### SPECIFY OR GENERATE WEIGHT AT LENGTH
  Wt<-(exp(-6.45e-5)*(len)^2.0) # alometric equation for weight-at-length
  # Not used in remaining code because biomass is not a unit of measure here. Still helpful to have, though. 
  
  #### SET STOCK-RECRUITMENT PARAMETERS
  med.recr<-NULL #median recruitment for constant or purely random recruitment
  rvar<-runif(1,0,0) #standard deviation of the recruitment error parameter 
  Rscale<-2   #Dummy parameter for these purposes
  type<-2 #Dummy parameter previously used to determine which Stock-Recruit curve to use. 
  ## This code only implements the Shell Area-Recruit relationship, so this parameter has no actual control 
  steepness<-0.7 #Dummy Parameter for these purposes
  
  
  #### SPECIFY OR GENERATE NATURAL MORTALITY 
  Mlo=0.10/12 #Minimum boundary for natural (non-disease) mortality
  Mhi=0.12/12 #Maximum boundary for natural (non-disease) mortality 
  v<-as.matrix(x=c(1:nages),nrow=nages,ncol=1) # As with maturity, mortality is handled as a matrix 
  # such that mortality could be set to vary annually
  for (i in c(1:nages)){
    x<-len.1[i,]
    if(x<35){
      y<-0.1439115  ## specify elevated mortality for oysters < 35mm (increaesd predation, settlement limitations, etc)
    }
    else {
      y<-(runif(1,Mlo,Mhi)) # individuals > 35mm and < 156 months old, mortality is randomly selected from the
      # min and max mortality previously specified
    }
    if(i>156){
      y<-0.04166 #Above 156 months of age, mortality elevates again 
    }
    #v[i,]<-y
    if(i==240){ # terminal age is 240 months, assign 100% mortality to this age group
      y<-1
    }
    v[i,]<-y
  }
  M.N<-as.numeric(v) # Save as natural mortality 
  M<-M.N
  
  
  ## Disarticulation rates
  D.B<-0.154032 # disarticulation rate of boxes at half life of 4.5 months
  D.C<-0.2/12 # decay rate of cultch at half life for 4 years then converted to month
  C.C<-0.40 # correction rate for cultch SHELL ROCK
  C.L<-1.0#*0.82 # correction factor for live animals SHELL ROCK
  C.B<-0.7#*0.82 # correction factor for boxes SHELL ROCK
  a.wt<-2.71#^(1/12)
  b.wt<-0.71#^(1/12)
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
  #end shellplant
  wt.1<-a.wt*(len.1^b.wt)
  saatage<-(len.1*wt.1*2*0.8)/(1e10) # conversion of len/wid at age to surface area at age in hectares
  #saatage[1:25]<-0
  rec.a.sa<-1.28e-2 # a for SA bev holt
  rec.b.sa<-1.40e-2 # b for SA bev holt
  shell.correct<-211 # correction factor for shell rock from Hemeon et al. 2021
  F.options<-c(0,0)
  F.dat<-rep(0,nmonths) #vector of F in each Month
  F.options.2<-c(0,0)
  D.dat<-rep(F.options.2,nmonths)
  F.options.3<-c(0,0) #seed fishery
  J.dat<-rep(F.options.3,nmonths)
  ### SPECIFCYING F OF 0.05 OR 0.2
  catch<-NULL #vector of catch (weight) in each Month
  #st.pop<-runif(nages,10000,Nzero) #population in Month one
  selectivity.I<-rep(1,nages) #selx for index 
  #I.error<-runif(1,0.0,0.5) #standard deviation of the index error parameter
  I.error<-0.0
  #can have either domed selx
  #selectivity.F<-1-(1/(1+exp(-7.63+(0.105*len)))) #selx for fishery from assessment (2013) S values.     
  #or asymptotic
  v<-as.matrix(x=c(1:nages),nrow=nages,ncol=1)
  for(i in c(1:nages)){
    x<-len.1[i,]
    if(x<63){
      y<-0
    }
    else {
      y<-1.0
    }
    v[i,]<-y
  }
  S.F<-as.numeric(v)
  selectivity.F<-c(S.F)
  #NOW DERMO SELECTIVITY
  v<-as.matrix(x=c(1:nages),nrow=nages,ncol=1)
  for (i in c(1:nages)){
    x<-len.1[i,]
    if(x<40){
      y<-0
    }
    else {
      y<-1
    }
    v[i,]<-y
  }
  S.D<-as.numeric(v)
  Selectivity.D<-c(S.D)
  
  #seed fishery selectivity
  v<-as.matrix(x=c(1:nages),nrow=nages,ncol=1)
  for (i in c(1:nages)){
    x<-len.1[i,]
    if(x<40){
      y<-0
    }
    else {
      y<-1
    }
    v[i,]<-y
  }
  S.J<-as.numeric(v)
  selectivity.J<-c(S.J)
  #selx for fishery from assessment (2013) GBK values.     ### THESE ARE STILL FROM HENNEN
  #F.error<-runif(1,0.0,0.25) #standard deviation of the fishery error parameter
  F.error<-0.0
  D.error<-0.0
  J.error<-0.0
  catchab<-0
  Nzero<-2344180224 ### N individuals at shell rock
  h_bound<-0.7 #upper bound on steepness parameter for Bev-Holt
  hlo<-0.7 #min steepness
  F_bound<-0 #upper bound on F
  F_bound.D<-0
  sigk<-0 #stdev in growth parameter k : sigk=0.0075,sigL=0.0075,sigt0=1.62 ### Changed to rate for 12 months
  sigL<-0 #stdev in growth parameter Linfinity
  sigt0<-0
  burnin<-200
  burnin2<-400
  ##### set up the recovery time specific variables
  Nstar_vec<-seq(0.05,0.5,0.05)  
  Ndelta_vec<-seq(0.05,0.5,0.1)
  
  interval<-160 #interval for big recruit classes not used for Bev Holt recruitment 
  #catchab<-runif(1,0.1,0.8) #index catchability
  return(list("nages"=nages,"nmonths"=nmonths,"len"=len,"Mat"=Mat,"Wt"=Wt,"rvar"=rvar
              ,"med.recr"=med.recr,"steepness"=steepness,"Nzero"=Nzero,"Rscale"=Rscale
              ,"type"=type,"M"=M,"D.B"=D.B,"D.C"=D.C,"C.C"=C.C,"C.L"=C.L,"C.B"=C.B,"shell.plant"=shell.plant,"a.wt"=a.wt,"b.wt"=b.wt,"wt.1"=wt.1,"saatage"=saatage,"rec.a.sa"=rec.a.sa,"rec.b.sa"=rec.b.sa,"shell.correct"=shell.correct
              ,"Mlo"=Mlo,"Mhi"=Mhi,"F.dat"=F.dat,"D.dat"=D.dat,"J.dat"=J.dat,"catch"=catch,"selectivity.I"=selectivity.I
              ,"I.error"=I.error,"selectivity.F"=selectivity.F,"Selectivity.D"=Selectivity.D,"selectivity.J"=selectivity.J,"F.error"=F.error,"D.error"=D.error,"J.error"=J.error,"catchab"=catchab
              ,"h_bound"=h_bound,"F_bound"=F_bound,"Nstar_vec"=Nstar_vec,"Ndelta_vec"=Ndelta_vec
              ,"burnin"=burnin,"burnin2"=burnin2,"sigk"=sigk,"sigL"=sigL,"sigt0"=sigt0,"interval"=interval,"hlo"=hlo))
}
MidMort <- setupmmrun()
save(MidMort, file = "C:/Users/Admin/Documents/Presentation/OysterSeminarPresentationBase/Oyster_MSE/MidMortDEBay.RData")

### High Mortality 
setuphmrun<-function(){  
  ## First step is to draw population parameters from a random distribution 
  ## or set them up on your own if you prefer...
  nmonths<-600    #round(runif(1,20,50),0) #months of simulated data
  nages<-240               #round(runif(1,200,300),0) #ages
  len<- vector(length=nages)
  for(i in c(1:nages)){
    age<-i/12
    len[i]<-140*(1-exp((-0.23)*((age)+.2))) # values from Kraueter for "new beds" high mortality 
  }
  #len<- c(7.82751,9.76421,11.66888,13.54209,15.38432,17.19611,18.97795,20.73034,22.45377,24.14871,25.81564,27.455,29.06729,30.65291,32.21233,33.74598,35.25428,36.73764,38.19649,39.63123,41.04225,42.42995,43.79471,45.13691,46.45693,47.75514,49.03188,50.28753,51.52241,52.73689,53.9313,55.10596,56.2612,57.39736,58.51473,59.61363,60.69437,61.75725,62.80257,63.8306,64.84164,65.83597,66.81387,67.7756,68.72143,69.65163,70.56646,71.46618,72.35101,73.22122,74.07704,74.91873,75.74649,76.56059,77.36122,78.14862,78.923,79.68459,80.43359,81.1702,81.89465,82.60712,83.30781,83.99693,84.67464,85.34116,85.99666,86.64133,87.27534,87.89888,88.51211,89.1152,89.70832,90.29165,90.86532,91.42952,91.98439,92.53009,93.06676,93.59458,94.11366,94.62417,95.12624,95.62,96.10562,96.5832,97.05289,97.51482,97.9691,98.41588,98.85529,99.28741,99.71241,100.13038,100.54144,100.9457,101.34328,101.7343,102.11884,102.49703,102.86897,103.23477,103.59451,103.94832,104.29626,104.63847,104.97501,105.306,105.63152,105.95165,106.26649,106.57613,106.88065,107.18013,107.47467,107.76434,108.04922,108.32939,108.60493,108.87592,109.14243,109.40453,109.6623,109.91582,110.16512,110.41032,110.65148,110.88864,111.12187,111.35126,111.57686,111.79872,112.01691,112.23151,112.44255,112.65011,112.85423,113.05499,113.25243,113.44659,113.63755,113.82536,114.01006,114.1917,114.37035,114.54604,114.71883,114.88876,115.05589,115.22025,115.3819,115.54087,115.69721,115.85098,116.0022,116.15092,116.29718,116.44102,116.58249,116.72163,116.85845,116.99302,117.12536,117.25552,117.38352,117.50941,117.63322,117.75498,117.87473,117.99251,118.10832,118.22224,118.33427,118.44444,118.55279,118.65936,118.76415,118.86723,118.96859,119.06828,119.16632,119.26275,119.35758,119.45084,119.54256,119.63276,119.72147,119.80872,119.89452,119.9789,120.0619,120.14352,120.22379,120.30273,120.38037,120.45673,120.53182,120.60567,120.67831,120.74974,120.81998,120.88908,120.95702,121.02385,121.08957,121.15421,121.21777,121.28028,121.34177,121.40224,121.4617,121.52018,121.5777,121.63426,121.6899,121.74461,121.79841,121.85133,121.90337,121.95455,122.00489,122.0544,122.10308,122.15096,122.19805,122.24436,122.28991,122.3347,122.37876,122.42208,122.46469,122.50661,122.54781,122.58834,122.6282,122.66741,122.70596,122.74387,122.78117,122.81784)
  
  # This is for SVAtoSNE
  #len[1:10]<-len[11]/10*seq(1,10)+0.1  #making growth linear below 10 yrs.  
  #len.1<-as.matrix((125*(1-exp((-0.24)*((seq(1,nages)/12)+.2))))) # remember to change the functions in here
  len.1<-as.matrix(len)
  v<-as.matrix(x=c(1:nages),nrow=nages,ncol=1)
  for (i in c(1:nages)){
    x<-len.1[i,]
    if(x<35){
      y<-0 ### divide by the number of months when individuals are less than 35mm
    }
    else {
      y<-1/12
    }
    v[i,]<-y
  }
  Mat<-as.numeric(v)
  # Mat<-c(Mat,rep(0.1,2640)) ## Making it to where they are only mature during one month of the year. 
  ###THIS DOESN'T END UP GETTING USED. DECIDED BY LENGTH LATER IN CODE. 
  Wt<-(exp(-6.45e-5)*(len)^2.0) #vector of weight at age (NEFSC 2011) - convert g to 1000 mt
  ### exp(-16.65938) derived from log(5.82e-8) reported by Pine et al. 
  med.recr<-NULL #median recruitment for constant or purely random recruitment
  rvar<-runif(1,0,0) #standard deviation of the recruitment error parameter (assessment estimates ~0.14 for S)
  Rscale<-2   #100000000/1.5 #not used for Bev-Holt
  type<-2 #1=BH 2=SR-BH 3=oq type - see documentation for this project or comments in Spawn code   
  steepness<-0.7 #this should be a random variable to test.... the starting value is from the assessment
  #round(runif(1,1,2),0) #1=Bev-Holt, 2=Ricker
  #if(type==1) steepness<-runif(1,0.2,0.99)  #stock recruit parameter for Bev-Holdt type recruitment
  #if(type==2) steepness<-runif(1,0.2,2.0)  #stock recruit parameter for Ricker type recruitment
  #Nzero<-runif(1,100000,10000000)  #stock recruit parameter - mainly affects the asymptote of the B-H        #runif(1,0.1,0.5) #natural mortality
  Mlo=0.10/12 #### ALL RATES, THINGS IN EXPONENTIALS THAT ARE SOMETHING TIMES TIME, IS DIVIDED BY 12 TO GIVE A MONTHLY CALCULATION
  Mhi=0.12/12  ### HAVE NOT BEEN UPDATED ACCORDING TO OYSTER DATA
  v<-as.matrix(x=c(1:nages),nrow=nages,ncol=1)
  for (i in c(1:nages)){
    x<-len.1[i,]
    if(x<35){
      y<-0.14391156### divide by the number of months when individuals are less than 35mm
      ## this makes it to where 70% of individuals die over this period.
    }
    else {
      y<-(runif(1,Mlo,Mhi))
    }
    if(i>156){
      y<-0.5/12 ## think about dividing this by 12?
    }
    #v[i,]<-y
    if(i==240){
      y<-1
    }
    v[i,]<-y
  }
  M.N<-as.numeric(v)
  ### NOW MAKE M THE SUM OF NATURAL AND DERMO MORTALITY
  #  M.N<-c(rep(0.1,240),M.N,rep(0.1,2400)) ## making it to where universal mortality occurs immediately after spawning
  M<-M.N
  D.B<-0.154032 # disarticulation rate of boxes at half life of 4.5 months
  D.C<-0.2/12 # decay rate of cultch at half life for 4 years then converted to month
  C.C<-0.21 # correction rate for cultch SHELL ROCK
  C.L<-1.0 # correction factor for live animals SHELL ROCK
  C.B<-0.73 # correction factor for boxes  High mortality= combined between cultch and box to account for interior/exterior of shell 
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
  #end shellplant
  a.wt<-2.71
  b.wt<-0.71
  wt.1<-a.wt*(len.1^b.wt)
  saatage<-(len.1*wt.1*2*0.8)/(1e10) # conversion of len/wid at age to surface area at age in hectares
  #saatage[1:25]<-0
  rec.a.sa<-4.23e-03 # a for SA bev holt
  rec.b.sa<-2.89e-03 # b for SA bev holt
  shell.correct<-275 # correction factor for shell rock from Hemeon et al. 2021
  F.options<-c(0,0)
  F.dat<-rep(0,nmonths) #vector of F in each month
  D.options<-c(0,0)
  D.dat<-rep(D.options,nmonths)
  J.options<-c(0,0)
  J.dat<-rep(J.options,nmonths)
  ### SPECIFCYING F OF 0.05 OR 0.2
  catch<-NULL #vector of catch (weight) in each month
  #st.pop<-runif(nages,10000,Nzero) #population in month one
  selectivity.I<-rep(1,nages) #selx for index 
  #I.error<-runif(1,0.0,0.5) #standard deviation of the index error parameter
  I.error<-0.0
  #can have either domed selx
  #selectivity.F<-1-(1/(1+exp(-7.63+(0.105*len)))) #selx for fishery from assessment (2013) S values.     
  #or asymptotic
  v<-as.matrix(x=c(1:nages),nrow=nages,ncol=1)
  for(i in c(1:nages)){
    x<-len.1[i,]
    if(x<63){
      y<-0
    }
    else {
      y<-1.0
    }
    v[i,]<-y
  }
  S.F<-as.numeric(v)
  selectivity.F<-c(S.F)
  #NOW DERMO SELECTIVITY
  v<-as.matrix(x=c(1:nages),nrow=nages,ncol=1)
  for (i in c(1:nages)){
    x<-len.1[i,]
    if(x<40){
      y<-0
    }
    else {
      y<-1
    }
    v[i,]<-y
  }
  S.D<-as.numeric(v)
  Selectivity.D<-c(S.D)
  
  #Seed fishery selectivity 
  v<-as.matrix(x=c(1:nages),nrow=nages,ncol=1)
  for (i in c(1:nages)){
    x<-len.1[i,]
    if(x<40){
      y<-0
    }
    else {
      y<-1
    }
    v[i,]<-y
  }
  S.J<-as.numeric(v)
  selectivity.J<-c(S.J)
  #selx for fishery from assessment (2013) GBK values.     ### THESE ARE STILL FROM HENNEN
  #F.error<-runif(1,0.0,0.25) #standard deviation of the fishery error parameter
  F.error<-0.0
  D.error<-0.0
  J.error<-0.0
  catchab<-0
  Nzero<-6044180224 ### N individuals at MM- increased by order of 10 for MM
  h_bound<-0.47 #upper bound on steepness parameter for Bev-Holt
  hlo<-0.47 #min steepness
  F_bound<-0 #upper bound on F
  D_bound<-0
  sigk<-0 #stdev in growth parameter k : sigk=0.0075,sigL=0.0075,sigt0=1.62 ### Changed to rate for 12 months
  sigL<-0 #stdev in growth parameter Linfinity
  sigt0<-0
  burnin<-200
  burnin2<-400
  ##### set up the recovery time specific variables
  Nstar_vec<-seq(0.05,0.5,0.05)  
  Ndelta_vec<-seq(0.05,0.5,0.1)
  
  interval<-160 #interval for big recruit classes not used for Bev Holt recruitment 
  #catchab<-runif(1,0.1,0.8) #index catchability
  return(list("nages"=nages,"nmonths"=nmonths,"len"=len,"Mat"=Mat,"Wt"=Wt,"rvar"=rvar
              ,"med.recr"=med.recr,"steepness"=steepness,"Nzero"=Nzero,"Rscale"=Rscale
              ,"type"=type,"M"=M,"D.B"=D.B,"D.C"=D.C,"C.C"=C.C,"C.L"=C.L,"C.B"=C.B,"shell.plant"=shell.plant,"a.wt"=a.wt,"b.wt"=b.wt,"wt.1"=wt.1,"saatage"=saatage,"rec.a.sa"=rec.a.sa,"rec.b.sa"=rec.b.sa,"shell.correct"=shell.correct
              ,"Mlo"=Mlo,"Mhi"=Mhi,"F.dat"=F.dat,"D.dat"=D.dat,"J.dat"=J.dat,"catch"=catch,"selectivity.I"=selectivity.I
              ,"I.error"=I.error,"selectivity.F"=selectivity.F,"Selectivity.D"=Selectivity.D,"selectivity.J"=selectivity.J,"F.error"=F.error,"D.error"=D.error,"J.error"=J.error,"catchab"=catchab
              ,"h_bound"=h_bound,"F_bound"=F_bound,"Nstar_vec"=Nstar_vec,"Ndelta_vec"=Ndelta_vec
              ,"burnin"=burnin,"burnin2"=burnin2,"sigk"=sigk,"sigL"=sigL,"sigt0"=sigt0,"interval"=interval,"hlo"=hlo)  
  )    
}

HighMort <- setuphmrun()
save(HighMort, file = "C:/Users/Admin/Documents/Presentation/OysterSeminarPresentationBase/Oyster_MSE/HighMortDEBay.RData")

## Shell Rock 
setupsrrun<-function(){  
  ## First step is to draw population parameters from a random distribution 
  ## or set them up on your own if you prefer...
  nmonths<-600    #round(runif(1,20,50),0) #months of simulated data
  nages<-240               #round(runif(1,200,300),0) #ages
  len<- vector(length=nages)
  for(i in c(1:nages)){
    age<-i/12
    len[i]<-125*(1-exp((-0.24)*((age)+.2)))
  }
  #len<- c(7.82751,9.76421,11.66888,13.54209,15.38432,17.19611,18.97795,20.73034,22.45377,24.14871,25.81564,27.455,29.06729,30.65291,32.21233,33.74598,35.25428,36.73764,38.19649,39.63123,41.04225,42.42995,43.79471,45.13691,46.45693,47.75514,49.03188,50.28753,51.52241,52.73689,53.9313,55.10596,56.2612,57.39736,58.51473,59.61363,60.69437,61.75725,62.80257,63.8306,64.84164,65.83597,66.81387,67.7756,68.72143,69.65163,70.56646,71.46618,72.35101,73.22122,74.07704,74.91873,75.74649,76.56059,77.36122,78.14862,78.923,79.68459,80.43359,81.1702,81.89465,82.60712,83.30781,83.99693,84.67464,85.34116,85.99666,86.64133,87.27534,87.89888,88.51211,89.1152,89.70832,90.29165,90.86532,91.42952,91.98439,92.53009,93.06676,93.59458,94.11366,94.62417,95.12624,95.62,96.10562,96.5832,97.05289,97.51482,97.9691,98.41588,98.85529,99.28741,99.71241,100.13038,100.54144,100.9457,101.34328,101.7343,102.11884,102.49703,102.86897,103.23477,103.59451,103.94832,104.29626,104.63847,104.97501,105.306,105.63152,105.95165,106.26649,106.57613,106.88065,107.18013,107.47467,107.76434,108.04922,108.32939,108.60493,108.87592,109.14243,109.40453,109.6623,109.91582,110.16512,110.41032,110.65148,110.88864,111.12187,111.35126,111.57686,111.79872,112.01691,112.23151,112.44255,112.65011,112.85423,113.05499,113.25243,113.44659,113.63755,113.82536,114.01006,114.1917,114.37035,114.54604,114.71883,114.88876,115.05589,115.22025,115.3819,115.54087,115.69721,115.85098,116.0022,116.15092,116.29718,116.44102,116.58249,116.72163,116.85845,116.99302,117.12536,117.25552,117.38352,117.50941,117.63322,117.75498,117.87473,117.99251,118.10832,118.22224,118.33427,118.44444,118.55279,118.65936,118.76415,118.86723,118.96859,119.06828,119.16632,119.26275,119.35758,119.45084,119.54256,119.63276,119.72147,119.80872,119.89452,119.9789,120.0619,120.14352,120.22379,120.30273,120.38037,120.45673,120.53182,120.60567,120.67831,120.74974,120.81998,120.88908,120.95702,121.02385,121.08957,121.15421,121.21777,121.28028,121.34177,121.40224,121.4617,121.52018,121.5777,121.63426,121.6899,121.74461,121.79841,121.85133,121.90337,121.95455,122.00489,122.0544,122.10308,122.15096,122.19805,122.24436,122.28991,122.3347,122.37876,122.42208,122.46469,122.50661,122.54781,122.58834,122.6282,122.66741,122.70596,122.74387,122.78117,122.81784)
  
  # This is for SVAtoSNE
  #len[1:10]<-len[11]/10*seq(1,10)+0.1  #making growth linear below 10 yrs.  
  #len.1<-as.matrix((125*(1-exp((-0.24)*((seq(1,nages)/12)+.2))))) # remember to change the functions in here
  len.1<-as.matrix(len)
  v<-as.matrix(x=c(1:nages),nrow=nages,ncol=1)
  for (i in c(1:nages)){
    x<-len.1[i,]
    if(x<35){
      y<-0 ### divide by the number of months when individuals are less than 35mm
    }
    else {
      y<-1/12
    }
    v[i,]<-y
  }
  Mat<-as.numeric(v)
  # Mat<-c(Mat,rep(0.1,2640)) ## Making it to where they are only mature during one month of the year. 
  ###THIS DOESN'T END UP GETTING USED. DECIDED BY LENGTH LATER IN CODE. 
  Wt<-(exp(-6.45e-5)*(len)^2.0) #vector of weight at age (NEFSC 2011) - convert g to 1000 mt
  ### exp(-16.65938) derived from log(5.82e-8) reported by Pine et al. 
  med.recr<-NULL #median recruitment for constant or purely random recruitment
  rvar<-runif(1,0,0) #standard deviation of the recruitment error parameter (assessment estimates ~0.14 for S)
  Rscale<-2   #100000000/1.5 #not used for Bev-Holt
  type<-2 #1=BH 2=SR-BH 3=oq type - see documentation for this project or comments in Spawn code   
  steepness<-0.7 #this should be a random variable to test.... the starting value is from the assessment
  #round(runif(1,1,2),0) #1=Bev-Holt, 2=Ricker
  #if(type==1) steepness<-runif(1,0.2,0.99)  #stock recruit parameter for Bev-Holdt type recruitment
  #if(type==2) steepness<-runif(1,0.2,2.0)  #stock recruit parameter for Ricker type recruitment
  #Nzero<-runif(1,100000,10000000)  #stock recruit parameter - mainly affects the asymptote of the B-H        #runif(1,0.1,0.5) #natural mortality
  Mlo=0.10/12 #### ALL RATES, THINGS IN EXPONENTIALS THAT ARE SOMETHING TIMES TIME, IS DIVIDED BY 12 TO GIVE A MONTHLY CALCULATION
  Mhi=0.12/12  ### HAVE NOT BEEN UPDATED ACCORDING TO OYSTER DATA
  v<-as.matrix(x=c(1:nages),nrow=nages,ncol=1)
  for (i in c(1:nages)){
    x<-len.1[i,]
    if(x<35){
      y<-0.1439115 ### divide by the number of months when individuals are less than 35mm
      ## this makes it to where 70% of individuals die over this period.
    }
    else {
      y<-(runif(1,Mlo,Mhi))
    }
    if(i>156){
      y<-0.5/12 ## think about dividing this by 12?
    }
    #v[i,]<-y
    if(i==240){
      y<-1
    }
    v[i,]<-y
  }
  M.N<-as.numeric(v)
  ### NOW MAKE M THE SUM OF NATURAL AND DERMO MORTALITY
  #  M.N<-c(rep(0.1,240),M.N,rep(0.1,2400)) ## making it to where universal mortality occurs immediately after spawning
  M<-M.N
  D.B<-0.154032 # disarticulation rate of boxes at half life of 4.5 months
  D.C<-0.2/12 # decay rate of cultch at half life for 4 years then converted to month
  C.C<-0.40 # correction rate for cultch SHELL ROCK
  C.L<-1.0#*0.82 # correction factor for live animals SHELL ROCK
  C.B<-0.7#*0.82 # correction factor for boxes SHELL ROCK
  a.wt<-2.71#^(1/12)
  b.wt<-0.71#^(1/12)  
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
  #end shellplant
  wt.1<-a.wt*(len.1^b.wt)
  saatage<-(len.1*wt.1*2*0.8)/(1e10) # conversion of len/wid at age to surface area at age in hectares
  #saatage[1:25]<-0
  rec.a.sa<-1.28e-2 # a for SA bev holt
  rec.b.sa<-1.40e-2 # b for SA bev holt
  shell.correct<-211 # correction factor for shell rock from Hemeon et al. 2021
  F.options<-c(0,0)
  F.dat<-rep(0,nmonths) #vector of F in each month
  F.options.2<-c(0,0)
  D.dat<-rep(F.options.2,nmonths)
  F.options.3<-c(0,0)
  J.dat<-rep(F.options.3,nmonths)
  ### SPECIFCYING F OF 0.05 OR 0.2
  catch<-NULL #vector of catch (weight) in each month
  #st.pop<-runif(nages,10000,Nzero) #population in month one
  selectivity.I<-rep(1,nages) #selx for index 
  #I.error<-runif(1,0.0,0.5) #standard deviation of the index error parameter
  I.error<-0.0
  #can have either domed selx
  #selectivity.F<-1-(1/(1+exp(-7.63+(0.105*len)))) #selx for fishery from assessment (2013) S values.     
  #or asymptotic
  v<-as.matrix(x=c(1:nages),nrow=nages,ncol=1)
  for(i in c(1:nages)){
    x<-len.1[i,]
    if(x<63){
      y<-0
    }
    else {
      y<-1.0
    }
    v[i,]<-y
  }
  S.F<-as.numeric(v)
  selectivity.F<-c(S.F)
  #NOW DERMO SELECTIVITY
  v<-as.matrix(x=c(1:nages),nrow=nages,ncol=1)
  for (i in c(1:nages)){
    x<-len.1[i,]
    if(x<40){
      y<-0
    }
    else {
      y<-1
    }
    v[i,]<-y
  }
  S.D<-as.numeric(v)
  Selectivity.D<-c(S.D)
  
  #Seed Fishery SELECTIVITY
  v<-as.matrix(x=c(1:nages),nrow=nages,ncol=1)
  for (i in c(1:nages)){
    x<-len.1[i,]
    if(x<40){
      y<-0
    }
    else {
      y<-1
    }
    v[i,]<-y
  }
  S.J<-as.numeric(v)
  selectivity.J<-c(S.J)
  #selx for fishery from assessment (2013) GBK values.     ### THESE ARE STILL FROM HENNEN
  #F.error<-runif(1,0.0,0.25) #standard deviation of the fishery error parameter
  F.error<-0.0
  D.error<-0.0
  J.error<-0.0
  catchab<-0
  Nzero<-2344180224 ### N individuals at shell rock
  h_bound<-0.7 #upper bound on steepness parameter for Bev-Holt
  hlo<-0.7 #min steepness
  F_bound<-0 #upper bound on F
  F_bound.2<-0
  sigk<-0 #stdev in growth parameter k : sigk=0.0075,sigL=0.0075,sigt0=1.62 ### Changed to rate for 12 months
  sigL<-0 #stdev in growth parameter Linfinity
  sigt0<-0
  burnin<-200
  burnin2<-400
  ##### set up the recovery time specific variables
  Nstar_vec<-seq(0.05,0.5,0.05)  
  Ndelta_vec<-seq(0.05,0.5,0.1)
  
  interval<-160 #interval for big recruit classes not used for Bev Holt recruitment 
  #catchab<-runif(1,0.1,0.8) #index catchability
  return(list("nages"=nages,"nmonths"=nmonths,"len"=len,"Mat"=Mat,"Wt"=Wt,"rvar"=rvar
              ,"med.recr"=med.recr,"steepness"=steepness,"Nzero"=Nzero,"Rscale"=Rscale
              ,"type"=type,"M"=M,"D.B"=D.B,"D.C"=D.C,"C.C"=C.C,"C.L"=C.L,"C.B"=C.B,"shell.plant"=shell.plant,"a.wt"=a.wt,"b.wt"=b.wt,"wt.1"=wt.1,"saatage"=saatage,"rec.a.sa"=rec.a.sa,"rec.b.sa"=rec.b.sa,"shell.correct"=shell.correct
              ,"Mlo"=Mlo,"Mhi"=Mhi,"F.dat"=F.dat,"D.dat"=D.dat,"J.dat"=J.dat,"catch"=catch,"selectivity.I"=selectivity.I
              ,"I.error"=I.error,"selectivity.F"=selectivity.F,"Selectivity.D"=Selectivity.D,"selectivity.J"=selectivity.J,"F.error"=F.error,"D.error"=D.error,"J.error"=J.error,"catchab"=catchab
              ,"h_bound"=h_bound,"F_bound"=F_bound,"Nstar_vec"=Nstar_vec,"Ndelta_vec"=Ndelta_vec
              ,"burnin"=burnin,"burnin2"=burnin2,"sigk"=sigk,"sigL"=sigL,"sigt0"=sigt0,"interval"=interval,"hlo"=hlo)  
  )    
}

ShellRock <- setupsrrun()
save(ShellRock, file = "C:/Users/Admin/Documents/Presentation/OysterSeminarPresentationBase/Oyster_MSE/ShellRockDEBay.RData")

