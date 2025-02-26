###################################################################################
# The modification here for OQ will be to include a possibility of a pulse recruitment
#that will be more likely to occur with > spawning stock, but otherwise random in time
###################################################################################
Spawn<-function(
  nages=1  #number of ages recquired
  ,nmonths=1 #number of Months required
  ,Mat=NULL #vector of maturity at age, or a matrix with nmonths X nages
  ,Wt=NULL #vector of weight at age, or a matrix with nmonths X nages
  ,med.recr.in=NULL #median recruitment for constant or purely random recruitment,
  #can be a vector of length nmonths
  ,Recsd=0 #standard deviation of the recruitment error parameter
  ,steepness=NULL  #stock recruit parameter
  ,Nzero=NULL  #stock recruit parameter - mainly affects the asymptote of the B-H
  ,Rscale=NULL #scale parameter for quahog type recruitment
  ,type=NULL     #1 for Bev-Holt, 2 for Ricker, NULL for random
  ,M=NULL #natural mortality can be a constant or vector with nages cells, or a 
  #matrix with nmonths X nages
  ,D.B=NULL # Disarticulation of Boxes
  ,D.C=NULL # Decay rate of cultch
  ,C.C=NULL # Correction factor for cultch
  ,C.B=NULL # Correction factor for boxes
  ,shell.plant=NULL #Shell Plant Addition 
  ,C.L=NULL # Correction factor for live animals
  ,a.wt=NULL # weight conversion from length
  ,b.wt=NULL # weight conversion parm b from length
  ,wt.1=NULL # weight conversion 
  ,saatage=NULL #surface area at age for live and boxes
  ,rec.a.sa=NULL # a for SA bev holt
  ,rec.b.sa=NULL # b for SA bev holt
  ,shell.correct=NULL # correction factor to bring BH curve to zero
  ,F.dat=NULL #vector of F in each Month
  ,D.dat=NULL # vector of Dermo in each Month
  ,J.dat=NULL # vector of seed fishery in each month
  ,catch=NULL #vector of catch (weight) in each Month
  ,st.pop=NULL #population in Month one
  ,selectivity.I=NULL #selx for index can be a vector of selectivity by age, 
  #with different selectivities by age for each Month
  ,I.error=0.0 #standard deviation of the index error parameter
  ,selectivity.F=NULL #selx for fishery can be a vector of selectivity by age, 
  #with different selectivities by age for each Month   
  ,Selectivity.D=NULL # selx for dermo
  ,selectivity.J=NULL #selx for seed fishery
  ,F.error=0.0 #standard deviation of the catch error parameter
  ,D.error=0.0 # erorr associated with dermo prevalence
  ,J.error=0.0 #error for seed fishery
  ,catchab=NULL #index catchability
  ,ac=0.0 #autocorrelation for management error
  ,man.error=0.0 #management error - will alter the way management actions occur
  ,ass.int=1 #assessment interval
  ,Plot1=FALSE #plot the simulated population?
  ,Unfished=FALSE #calculate unfished population parameters (SSN and total N)?
  ,Nthresh=NULL #Threshold for triggering recovery mode
  ,Ntarg=NULL #Target for for triggering recovery mode
  ,F.target=NULL #F target for managed fishery
  ,D.target=NULL # Mortality to Dermo
  ,J.target=NULL #seed mortality
  ,seed1=NULL #random number seed - can be held for replication of a specific recruitment sequence
  ,inputs=NULL #This accepts a list made by a previous run of spawn.  If you want to rerun
  #a sequence, use this (it can include some altered arguments, you just have to manipulate
  #them directly i.e. testrun$inputs$M<-0.25).
) {
  ################################################################################
  #       Function to generate simulated populations with prespecified           #
  # parameters                                                                   #
  #______________________________________________________________________________#
  #	Code modified from JJD 4/3/13
  #_______________________________________________________________________________
  #
  #Store inputs or fill in values from previous runs 
  if(length(inputs)>0) {
    nages<-inputs$nages;nmonths<-inputs$nmonths;Mat<-inputs$Mat;len<-len;
    Wt<-inputs$Wt;med.recr.in<-inputs$med.recr.in;Recsd<-inputs$Recsd;
    steepness<-inputs$steepness;Nzero<-inputs$Nzero;type<-inputs$type;
    M<-inputs$M;D.B<-inputs$D.B;D.C<-inputs$D.C;C.C<-inputs$C.C;C.L<-inputs$C.L;C.B<-inputs$C.B;shell.plant<-inputs$shell.plant;a.wt<-inputs$a.wt;b.wt<-inputs$b.wt;wt.1<-inputs$wt.1;saatage<-inputs$saatage;
    rec.a.sa<-inputs$rec.a.sa;rec.b.sa<-inputs$rec.b.sa;F.dat<-inputs$F.dat;D.dat<-inputs$D.dat;catch<-inputs$catch;st.pop<-inputs$st.pop;
    selectivity.I<-inputs$selectivity.I;I.error<-inputs$I.error;
    selectivity.F<-inputs$selectivity.F;Selectivity.D<-inputs$Selectivity.D;selectivity.J<-inputs$selectivity.J
    F.error<-inputs$F.error;D.error<-inputs$D.error;J.error<-inputs$J.error
    catchab<-inputs$catchab;Plot1<-inputs$Plot1;Unfished<-inputs$Unfished;
    seed1<-inputs$seed1;ac<-inputs$ac;man.error<-inputs$man.error;
    ass.int<-inputs$ass.int;Ntarg<-inputs$Ntarg;Nthresh<-inputs$Nthresh; F.target<-inputs$F.target; D.target<-inputs$D.target; J.target<-inputs$J.target;
  } else {inputs<-list("nages"=nages,"nmonths"=nmonths,"Mat"=Mat,"Wt"=Wt,"len"=len
                       ,"med.recr.in"=med.recr.in,"Recsd"=Recsd,"steepness"=steepness,"Nzero"=Nzero,"Rscale"=Rscale
                       ,"type"=type,"M"=M,"D.B"=D.B,"D.C"=D.C,"C.C"=C.C,"C.L"=C.L,"C.B"=C.B,"shell.plant"=shell.plant,"a.wt"=a.wt,"b.wt"=b.wt,"wt.1"=wt.1,"saatage"=saatage,"rec.a.sa"=rec.a.sa,"rec.b.sa"=rec.b.sa
                       ,"F.dat"=F.dat,"D.dat"=D.dat,"J.dat"=J.dat,"catch"=catch,"st.pop"=st.pop,"selectivity.I"=selectivity.I
                       ,"I.error"=I.error,"selectivity.F"=selectivity.F,"Selectivity.D"=Selectivity.D,"selectivity.J"=selectivity.J,"F.error"=F.error,"D.error"=D.error,"J.error"=J.error,"catchab"=catchab,"ac"=ac
                       ,"man.error"=man.error,"ass.int"=ass.int,"Ntarg"=Ntarg,"Nthresh"=Nthresh,"F.target"=F.target, "D.target"=D.target,"J.target"=J.target
                       ,"Plot1"=Plot1,"Unfished"=Unfished,"seed1"=seed1)}
  
  M<-fill(M,nmonths,nages)
  Wt<-fill(Wt,nmonths,nages)
  Mat<-fill(Mat,nmonths,nages)
  ## for some reason not converted properly, so change to numeric
  saatage<-as.numeric(saatage)
  saatage<-fill(saatage,nmonths,nages)
  selectivity.I<-fill(selectivity.I,nmonths,nages)
  selectivity.F<-fill(selectivity.F,nmonths,nages)
  Selectivity.D<-fill(Selectivity.D,nmonths,nages)
  selectivity.J<-fill(selectivity.J,nmonths,nages)
  med.recr.in<-NULL
  med.recr<-NULL
  
  if(length(med.recr.in)>0 & length(med.recr.in)<nmonths) med.recr.in<-rep(med.recr.in,nmonths)
  #error trap section
  if(length(F.dat)!=nmonths & length(catch)!=nmonths) {
    print("ERROR Length of F or catch vector must be = to number of months!!");return(NULL)}
  if(length(D.dat)!=nmonths & length(catch)!=nmonths) {
    print("ERROR Length of D or catch vector must be = to number of months!!");return(NULL)}
  if(length(Mat)<nages) {
    print("ERROR Length of maturity vector must be = to number of ages!!");return(NULL) }  
  if(length(Wt)<nages) {
    print("ERROR Length of weight vector must be = to number of ages!!"); return(NULL)}
  if(length(F.dat)<1 & length(catch)<1) {print("No catch specified!"); return(NULL)}  
  if(length(D.dat)<1 & length(catch)<1) {print("No catch specified!"); return(NULL)} 
  if(length(J.dat)<1 & length(catch)<1) {print("No catch specified!"); return(NULL)}  
  if(length(st.pop)!=nages) {
    print("ERROR Length of starting population vector must be = to number of ages!!");return(NULL)}
  if(length(selectivity.I)<nages) {
    print("ERROR Length of index selectivity vector must be = to number of ages!!");return(NULL)}
  if(length(selectivity.F)<nages) {
    print("ERROR Length of index selectivity vector must be = to number of ages!!");return(NULL)}
  if(length(Selectivity.D)<nages) {
    print("ERROR Length of index selectivity vector must be = to number of ages!!");return(NULL)}
  if(length(selectivity.J)<nages) {
    print("ERROR Length of index selectivity vector must be = to number of ages!!");return(NULL)}
  if((length(med.recr.in)>0 & length(steepness>0)) | (length(med.recr.in)>0 & length(Nzero>0))) {
    print("ERROR! Cannot specify a median recruitment and stock recruit parameters!!");return(NULL)}
  #add section for invalid values!!!!
  if(length(type)==0) {print("ERROR !! Must specify a valid stock recruit type [type: 0=random, 1=Bev-Holdt, 2=Ricker]")}
  M<-ifelse(M<0,0,M)
  if(Recsd<0) {print("ERROR !! Negative recruitment error sd [Recsd] not allowed")}
  
  #_______________________________________________________________________________  
  #this needs to go after the matrix declaration section  
  
  SSN<-c() #SSN in each month
  SSN_m<-c() #apparent SSN in each month (includes assessment error)
  FAA<-matrix(NA,nrow=nmonths,ncol=nages)  #F at age in each month
  FAA.D<-matrix(NA,nrow=nmonths,ncol=nages) #Dermo mort at age in each month
  FAA.J<-matrix(NA,nrow=nmonths,ncol=nages) #seed fishery mort at age in each month
  ZAA<-matrix(NA,nrow=nmonths,ncol=nages) #total mortality at age in each month
  NAA<-matrix(NA,nrow=nmonths,ncol=nages)  #numbers at age in each month
  SSNAA<-matrix(NA,nrow=nmonths,ncol=nages) #SSN at age in each month
  CAA<-matrix(NA,nrow=nmonths,ncol=nages) #catch at age in each month
  CSAA.1<-matrix(NA,nrow=nmonths,ncol=nages)#cultch at age in each month
  CSAA<-matrix(NA,nrow=nmonths,ncol=1)# summed cultch in each month
  BAA<-matrix(NA,nrow=nmonths,ncol=nages)#boxes at age in each month
  BSA<-matrix(NA,nrow=nmonths,ncol=nages) # box surface area at age in each month
  BSA.plant<-matrix(NA, nrow=nmonths, ncol=nages) #planted shell surface area
  SAAA<-matrix(NA,nrow=nmonths,ncol=nages)# surface area at age in each month
  SAMAA<-matrix(NA,nrow=nmonths,ncol=nages) #mortality that attributes to surface area at age
  BSA_tot<-matrix(NA,nrow=nmonths,ncol=1) # summed box surface area in each month
  plant_tot<-matrix(NA, nrow=nmonths, ncol=1) # planted shell total surface area
  LSA_tot<-matrix(NA,nrow=nmonths,ncol=1)# summed live surface area in each month
  NAA_tot<-matrix(NA,nrow=nmonths,ncol=1)# summed numbers at age in each month
  CSAA_tot<-matrix(NA,nrow=nmonths,ncol=1)#summed cultch in each month
  TotN<-c() #sum of true abundance at age in each month
  TotCatch<-c() #vector of total  catch each month
  TotYieldB<-c() #vector of yield each month (CAA * Wt)
  IAA<-matrix(NA,nrow=nmonths,ncol=nages) #index obs at age
  TotI<-c() #Index observation each month (sum of IAA)
  FullF<-c() #Full each month
  FUllD<-c()
  FullJ<-c()
  recruits<-c()
  MSY1<-c()
   #TotYield<-catch
  NAA[1,]<-st.pop #population in month one
  recruits[1]<-st.pop[1]
  SSN[1]=sum(NAA[1,]*Mat[1,])      #SSN month one ## REMOVED WEIGHT IN * WT[1,]
  ## get some starting values from the sa.function
  #begin.sa<-sa.func(nmonths=1000,nages=nages,M=M,a.wt=a.wt,b.wt=b.wt,Mat=Mat,N0=Nzero,D.B=D.B,
  #                 D.C=D.C,C.L=C.L,C.C=C.C,steep=steepness,Reca=rec.a.sa,Recb=rec.b.sa)
  
  
  BAA[1,]=begin.sa$BAA.1 # boxes at age
  BSA[1,]=begin.sa$BSA.1 #box surface area at age
  BSA_tot[1,1] = sum(BSA[1,1:nages])
  BSA.plant<-shell.plant
  plant_tot[1,1]=sum(BSA.plant[1,1:nages]*C.B)
  CSAA.1[1,]=begin.sa$CSA # make a matrix so it tracks through time but won't actually operate as such 
  #print(CSAA[1,])
  NAA[1,]=begin.sa$NAA.1 # numers at age- this is used instead of the previous specification so they've all been burned in the same
  ## NEED TO HAVE A STARTING SA
  SA.1<-begin.sa$SA1.1
  ###NEED TO INSERT SOMETHING ABOUT CULTCH. MIGHT NEED TO CONVERT TO SAA THOUGH FIRST
  
  #SSNAA[1,1]<-NAA[1,1]
  needF<-TRUE
  if(length(F.dat)>0){
    FAA<-matrix(rep(F.dat,nages),ncol=nages,nrow=nmonths,byrow=FALSE)*selectivity.F 
    FullF<-F.dat
    needF<-FALSE  
  }
  
  if(needF){
    FullF[1]<-Solve.F(selectivity=selectivity.F,M=M[1,],N.age=NAA[1,],catch1=catch[1]) ## REMOVED WEIGHT
  }
  FAA[1,]<-FullF[1]*selectivity.F[1,] #Fishing selectivity
  
  
  #### ADD DERMO MORTALITY 
  needF<-TRUE
  if(length(D.dat)>0){
    FAA.D<-matrix(rep(D.dat,nages),ncol=nages,nrow=nmonths,byrow=FALSE)*Selectivity.D 
    FullD<-D.dat
    needF<-FALSE  
  }
  
  if(needF){
    FullF[1]<-Solve.F(selectivity=Selectivity.D,M=M[1,],N.age=NAA[1,],catch1=catch[1]) ## REMOVED WEIGHT
  }
  FAA.D[1,]<-FullD[1]*Selectivity.D[1,] #Fishing selectivity
  
  #### ADD SEED MORTALITY 
  needF<-TRUE
  if(length(J.dat)>0){
    FAA.J<-matrix(rep(J.dat,nages),ncol=nages,nrow=nmonths,byrow=FALSE)*Selectivity.D 
    FullJ<-J.dat
    needF<-FALSE  
  }
  
  if(needF){
    FullF[1]<-Solve.F(selectivity=selectivity.J,M=M[1,],N.age=NAA[1,],catch1=catch[1]) ## REMOVED WEIGHT
  }
  FAA.J[1,]<-FullJ[1]*selectivity.J[1,] #Fishing selectivity
  
  ZAA[1,]<-FAA[1,]+M[1,]+FAA.D[1,]+FAA.J[1,]   #total mort 
  SAMAA[1,]<-M[1,]+FAA.D[1,] # mortality that contributes to shell mass ## IS THIS REALLY JUST BAA BUT WITH F INCLUDED
  
  
  
  
  IAA[1,]<-NAA[1,]*catchab*selectivity.I[1,] #index at age in month one
  TotI[1]<-sum(IAA[1,]) #index for the month
  TotN[1]<-sum(NAA[1,])
  
  CAA[1,]<-((FAA[1,])/(FAA.D[1,]+FAA[1,]+M[1,]))*NAA[1,]*(1-exp(-1*(FAA[1,]+FAA.D[1,]+M[1,]))) #Catch at age in each month
  CAA[1,]<-ifelse(CAA[1,]<1 || is.finite(CAA[1,])==FALSE,0,CAA[1,])  
  
  TotCatch[1]<-sum(((FAA[1,])/(FAA[1,]+FAA.D[1,]+M[1,]))*NAA[1,]*(1-exp(-1*(FAA[1,]+FAA.D[1,]+M[1,]))))  #Total Catch each month
  if(TotCatch[1]<1 || is.finite(TotCatch[1])==FALSE) TotCatch[1]<-0 
  TotYieldB[1]<-sum(((FAA[1,])/(FAA.D[1,]+FAA[1,]+M[1,]))*NAA[1,]*(1-exp(-1*(FAA.D[1,]+FAA[1,]+M[1,]))))  #Total Y1eld each month ### REMOVED WEIGHT SO WE GET N
  if(TotYieldB[1]<1 || is.finite(TotYieldB[1])==FALSE) TotYieldB[1]<-0   
  
  if(type==3) med.recr.in[1]=NAA[1,1]
  MSY1<-MSY(steep=steepness,N0=Nzero,Mat=Mat[1,]
            ,M=M[1,],selectivity.F=selectivity.F[1,],type=type,med.recr=med.recr.in[1]) ## REMOVED WEIGHT
  MSY2<-MSY.2(steep=steepness,N0=Nzero,Mat=Mat[1,]
              ,M=M[1,],Selectivity.D=Selectivity.D[1,],type=type,med.recr=med.recr.in[1])### SHOULD WE ADD MORTALITY FROM OTHER FISHERY IN THIS CALCULATION?
  
  if(length(seed1)>0){
    seed.t<-seed1 #Sets random number seed
  } else seed.t=NULL
  
  #recovery testing variables 
  old.ac<-0 #no auto correlated errors in first month of simulation
  recovery_time<-0
  recovery<-FALSE
  recovery_start<-0
  Nofish_time<-0
  Nofish<-FALSE
  ######################################################################################################
  for (i in 2:nmonths) #month loop;
  {
    #   if(!is.null(seed.t)) {
    #    seed.t<-seed.t+69 #alter the seed for each month!
    #   set.seed(seed.t)     
    #}
    
    #if(length(steepness)>0 & length(Nzero)>0) {
    # NAA[i,1]=NAA[i-1,1]*add.error(Recsd)  #bias correction is handled in the function e.g. exp(rnorm(1,0,Recsd)-b.correct)
    # NAA[i,(2:(nages-1))]=NAA[i-1,(1:(nages-2))]*exp(-ZAA[(i-1),(1:(nages-2))])  #Fill abundance at age matrix on the diagonals
    # NAA[i,nages]=(NAA[(i-1),(nages-1)]*exp(-ZAA[(i-1),(nages-1)]))+(NAA[(i-1),nages]*exp(-ZAA[(i-1),nages])) #plus group
    # SSNAA[i,]=NAA[i,]*Mat[i,] #SSN at age in each month ## REMOVED WEIGHT
    # 
    # ## SAMAA specified earlier as M + F2
    # BAA[i,(1:nages)]=(NAA[(i-1),(1:nages)]*(1-exp(-SAMAA[(i-1),(1:nages)]))*2*C.L)+(BAA[(i-1),(1:nages)]*exp(-D.B))# multiplied by 2 to get interior and exterior shells 
    # 
    # ## Part 1- so individuals that would have died and turned to boxes.... 
    # ## Part 2- Boxes from last month that would contribute to new boxes and disarticulation rate
    # BSA[i,(1:nages)]=BAA[i,1:nages]*saatage[i,1:nages]*C.B
    # ## Part 3- convert to surface area
    # ## Now do it for cultch
    # ## Part 1- take boxes that disarticulated last month
    # ## Part 2- take cultch that hasn't decayed yet
    # #CSAA[i,]=sum((BSA[(i-1),(1:nages)]*(1-exp(-D.B))))+(CSAA[(i-1)]*exp(-D.C))
    # CSAA.1[i,1:nages]<-(BAA[(i-1),(1:nages)]*(1-exp(-D.B)))*saatage[i-1,1:nages]*C.B*C.C
    # 
    # #CSAA[i,]<-sum(BSA[i-1,1:nages]*(1-exp(-D.B))) ## Cultch is equal to what shell has disarticulated
    # CSAA[i,]<-sum(CSAA.1[i,1:nages])
    # CSAA[i,]<-CSAA[i,]+((sum(CSAA[i-1,]))*exp(-D.C))-(sum(BSA[i,(1:nages)])*0.18)-(sum(NAA[i,1:nages]*saatage[i,1:nages])*C.L*0.18)
    # 
    # 
    # ## NOW GET TOTAL SURFACE AREA AVAILABLE WITH CORRECTION FACTORS 
    # SA.1<-sum(NAA[i,1:nages]*saatage[i,1:nages])*C.L # don't forget to correct for live surface area
    # 
    # SA.1<-SA.1 + sum(CSAA[i,]) + sum(BSA[i,1:nages])
    #if(i==2){print(SA.1)}
    #if(i==100){print(SA.1)}
    med.rec<-stock.recruit(SSN[i-1],steep=steepness,N0=Nzero,Mat=Mat[i,],M=M[i,],type=type,Rscale=Rscale,
                           SA.1=SA.1,rec.a.sa=rec.a.sa,rec.b.sa=rec.b.sa,shell.correct=shell.correct) ## REMOVED WEIGHT
    
    #} else {
    # med.rec<-med.recr.in[i]
    
    #  } #otherwise jsut use the constant recruitment par from user
    NAA[i,1]=med.rec*add.error(Recsd)  #bias correction is handled in the function e.g. exp(rnorm(1,0,Recsd)-b.correct)
    NAA[i,(2:(nages-1))]=NAA[i-1,(1:(nages-2))]*exp(-ZAA[(i-1),(1:(nages-2))])  #Fill abundance at age matrix on the diagonals
    NAA[i,nages]=(NAA[(i-1),(nages-1)]*exp(-ZAA[(i-1),(nages-1)]))+(NAA[(i-1),nages]*exp(-ZAA[(i-1),nages])) #plus group
    SSNAA[i,]=NAA[i,]*Mat[i,] #SSN at age in each month ## REMOVED WEIGHT
    
    ## SAMAA specified earlier as M + F2
    BAA[i,(1:nages)]=(NAA[(i-1),(1:nages)]*(1-exp(-SAMAA[(i-1),(1:nages)]))*2*C.L)+(BAA[(i-1),(1:nages)]*exp(-D.B))#multiplied by two to account for interior and exterior of boxes
    ## Part 1- so individuals that would have died and turned to boxes.... 
    ## Part 2- Boxes from last month that would contribute to new boxes and disarticulation rate
    BSA[i,(1:nages)]=BAA[i,1:nages]*saatage[i,1:nages]*C.B
    
    ## now calculate the destruction of planted shell
    BSA.plant[i,(1:nages)]=(BSA.plant[(i-1),(1:nages)]*exp(-D.B))+shell.plant[i, (1:nages)] #DB being disart. rate
    ## Part 3- convert to surface area
    ## Now do it for cultch
    ## Part 1- take boxes that disarticulated last month
    ## Part 2- take cultch that hasn't decayed yet
    #CSAA[i,]=sum((BSA[(i-1),(1:nages)]*(1-exp(-D.B))))+(CSAA[(i-1)]*exp(-D.C))
    
    #CSAA.1[i,1:nages]<-(BAA[i-1,1:nages]*saatage[i-1,1:nages]*(1-exp(-D.B))*C.B*C.C)
    #CSAA <- (sum(CSAA.1[i,1:nages]))+((CSAA[i-1,])*exp(-D.C))-(sum(BSA.1[j-1,1:nages])*0.18)-(sum(NAA.1[j,1:nages]*saatage.1[j,1:nages])*C.L*0.18)
    
    
    CSAA.1[i,1:nages]<-(BAA[i-1,1:nages]*saatage[i-1,1:nages]*(1-exp(-D.B))*C.B*C.C)+((CSAA.1[i-1,1:nages])*exp(-D.C))+(BSA.plant[(i-1),(1:nages)]*(1-exp(-D.B)))
    #CSAA[i,]<-sum(BSA[i-1,1:nages]*(1-exp(-D.B))) ## Cultch is equal to what shell has disarticulated
    CSAA[i,]<-sum(CSAA.1[i,1:nages])-(sum(NAA[i,1:nages]*saatage[i,1:nages])*C.L*0.18)-(sum(BSA[i,(1:nages)])*0.18)
    #CSAA[i,]<-CSAA[i,]-(sum(NAA[i,1:nages]*saatage[i,1:nages])*C.L*0.18)-(sum(BSA[i,(1:nages)])*0.18)#   ## NOW GET TOTAL SURFACE AREA AVAILABLE WITH CORRECTION FACTORS 
    SA.1<-sum(NAA[i,1:nages]*saatage[i,1:nages])*C.L # don't forget to correct for live surface area
    
    SA.1<-SA.1 + sum(CSAA[i,]) + sum(BSA[i,1:nages]) + sum(BSA.plant[i,1:nages]*C.B)
    
    #my.list<-list(CSA,BAA.1[1999,],BSA.1[1999,],NAA.1[2000,],SA1.1)
    #names(my.list)<-c("CSA","BAA.1","BSA.1","NAA.1","SA1.1")
    #return(my.list)
    
    
    if(needF){ #trigger the iterative F solution if catch is provided instead of F
      FullF[i]<-Solve.F(selectivity=selectivity.F[i,],selectivity.2=Selectivity.D[i,],M=M[i,],N.age=NAA[i,],catch1=catch[i]) ### REMOVED WEIGHT
    }
    
    err<-add.error(F.error)    
    
    if(needF){ #trigger the iterative F solution if catch is provided instead of F
      FullD[i]<-Solve.D(selectivity=Selectivity.D[i,],selectivity.2=selectivity.F[i,],M=M[i,],N.age=NAA[i,],catch1=catch[i]) ### REMOVED WEIGHT
    }
    
    err.D<-add.error(D.error)    
    
    if(needF){ #trigger the iterative F solution if catch is provided instead of F
      FullJ[i]<-Solve.J(selectivity=selectivity.J[i,],selectivity.2=selectivity.J[i,],M=M[i,],N.age=NAA[i,],catch1=catch[i]) ### REMOVED WEIGHT
    }
    
    err.J<-add.error(J.error)  
    
    FAA[i,]<-FullF[i]*selectivity.F[i,]*err
    FAA.D[i,]<-FullD[i]*Selectivity.D[i,]*err.D
    FAA.J[i,]<-FullJ[i]*selectivity.J[i,]*err.J
    ZAA[i,]<-FAA[i,]+M[i,]+FAA.D[i,]+FAA.J[i,]
    SAMAA[i,]<-M[i,]+FAA.D[i,]
    
    CAA[i,]<-((FAA[i,])/(FAA.D[i,]+FAA[i,]+M[i,]))*NAA[i,]*(1-exp(-1*(FAA.D[i,]+FAA[i,]+M[i,]))) #Catch at age in each month
    ifelse(CAA[i,]<1 || is.finite(CAA[i,])==FALSE,0,CAA[i,])  
    
    
    TotCatch[i]<-sum(((FAA[i,])/(FAA[i,]+FAA.D[i,]+M[i,]))*NAA[i,]*(1-exp(-1*(FAA[i,]+FAA.D[i,]+M[i,])))) #*err  #Total Catch each month
    if(TotCatch[i]<1 || is.finite(TotCatch[i])==FALSE) TotCatch[i]<-0 
    TotYieldB[i]<-sum(((FAA[i,])/(FAA[i,]+FAA.D[i,]+M[i,]))*NAA[i,]*(1-exp(-1*(FAA.D[i,]+FAA[i,]+M[i,])))) #*err  #Total Yield each month ## REMOVED WEIGHT
    if(TotYieldB[i]<=0 || is.finite(TotYieldB[i])==FALSE) TotYieldB[i]<-0 
    #Extinction flag if catch this month is more than the population last month... or if the number of animals is < 1
    
    #### DIAGNOSTICS ####
    #print(c(TotCatch[i],sum(NAA[i,]),FullF[i]))   
    
    if(sum(NAA[i-1,])<=(TotCatch[i]) | sum(NAA[i,])<1) {
      print(paste("EXTINCTION in month",i))
      
      Catch.at.age<-data.frame("Month"=c(1:i),CAA[(1:i),])
      names(Catch.at.age)[2:(nages+1)]<-paste("Age_",c(1:nages),sep="")
      
      Age.comp<-data.frame("Month"=c(1:i),NAA[(1:i),]/sum(NAA[(1:i),]))
      names(Age.comp)[2:(nages+1)]<-paste("Age_",c(1:nages),sep="")
      
      rez<-data.frame("Month"=seq(1,i),"Total_N"=TotN[1:i],"SSN"=SSN[1:i]
                      ,"Recruits"=recruits[1:i],"Index"=TotI[1:i],"Catch"=TotYieldB[1:i],"Catch.N"=TotCatch[1:i]
                      ,"F"=FullF[1:i],"D"=FullD[1:i],"J"=FullJ[1:i],"MSY"=rep(MSY1[1,1],i),"MSY.D"=rep(MSY2[1,1],i),"F_MSY"=rep(MSY1[1,2],i),"F_MSY.D"=rep(MSY2[1,1],i))   #,"OQ_MSY"=MSY2[1:i,1],"OQ_F_MSY"=MSY2[1:i,2]
      #Total planted shell surface area
      for(j in c(1:nmonths)){
        plant_tot[j]<-sum(BSA.plant[j,1:nages]*C.B)
      }
      
      # Total box surface area
      for(j in c(1:nmonths)){
        BSA_tot[j]<-sum(BSA[j,1:nages])
      }
      # Total Live Surface area
      for(j in c(1:nmonths)){
        LSA_tot[j]<-sum(NAA[j,1:nages]*saatage[j,1:nages])*C.L 
      }
      # Total number of animals 
      for(j in c(1:nmonths)){
        NAA_tot[j]<-sum(NAA[j,1:nages])
      }
      # Total number of animals 
      for(j in c(1:nmonths)){
        CSAA_tot[j]<-(CSAA[j,])
        
      }
      
      
      
      #BSA_tot<-rbind(BSA_tot,BSA_new)
      surfacearea<- data.frame("month"=c(1:i),"BSA" = BSA_tot[1:i], "Plant" = plant_tot[1:i], "CSA" = CSAA_tot[1:i],"LSA"=LSA_tot[1:i],
                               "Mort"=Mlo,"Dermo"=D.target,"Seed"=J.target,"F.target"=F.target,"NAA"=NAA_tot[1:i],"Recruits"=NAA[1:i,1])
      #surfacearea$CSA<-surfacearea$CSA*C.C #correct for cultch to avoid doing it before and messing up n dimensions
      
      if(Unfished) {
        uf1<-Unfished.pop(nages=nages,nmonths=nmonths,Mat=Mat,med.recr=med.recr
                          ,Recsd=Recsd,steepness=steepness,Nzero=Nzero,type=type,M=M,D.B=D.B,D.C=D.C,C.L=C.L,C.C=C.C,C.B=C.B,st.pop=st.pop
                          ,selectivity.I=selectivity.I,catchab=catchab)
        #return 
        return(list("Simulated.population"=rez,"Unfished.Status"=uf1,"Catch.at.age"=Catch.at.age
                    ,"Age.comp"=Age.comp,"inputs"=inputs,"surfacearea"=surfacearea))    
      } else return(list("Simulated.population"=rez,"Catch.at.age"=Catch.at.age,"Age.comp"=Age.comp
                         ,"inputs"=inputs,"Rec_start"=recovery_start,"Rec_time"=recovery_time,"Nofish_time"=Nofish_time,"surfacearea"=surfacearea))      
      
    }
    
    recruits[i]<-NAA[i,1]
    IAA[i,]=NAA[i,]*catchab*selectivity.I[i,]   #index at age each month  #TEST THIS AGAINST A LAG OF ABUNDANCE
    TotI[i]=sum(IAA[i,])*add.error(I.error)    #total index for the month  
    
    SSN[i]=sum(na.omit(SSNAA[i,])) #Calculate true total SSN in each month (sum across ages)
    if (SSN[i]<=0) {
      SSN[i]=0 #Once SSN is less 1 SSN is set to zero and the stock collapses.
    } else {
      #      MSY1<-rbind(MSY1,MSY(SSN[i],steep=steepness,N0=Nzero,Mat=Mat[i,],selectivity.F=selectivity.F[i,]
      #	  ,Wt=Wt[i,],M=M[i,],type=type,med.recr=median(recruits)))
    }  
    
    ####### Modifications for testing management - include assessment error (implementation error already built in)
    # Only assess the population if there is an assessment this month!
    
  } #end month loop
  
  
  if(Plot1==TRUE) {plot.SPAWN(rowSums(NAA)[1:i],recruits[1:i],SSN[1:i],SSN_m[1:i],Ntarg
                              ,Nthresh,TotYieldB[1:i],st.pop,nages,NAA[i,],FullF,FullD,FullJ,steepness)}
  
  
  Catch.at.age<-data.frame("Month"=c(1:i),CAA[(1:i),])
  names(Catch.at.age)[2:(nages+1)]<-paste("Age_",c(1:nages),sep="")
  
  Age.comp<-data.frame("Month"=c(1:i),NAA[(1:i),]/rowSums(NAA[(1:i),]))
  names(Age.comp)[2:(nages+1)]<-paste("Age_",c(1:nages),sep="")
  
  #return results  
  rez<-data.frame("month"=seq(1,length(TotI)),"Total_N"=TotN[1:i],"SSN"=SSN[1:i]
                  ,"Recruits"=recruits[1:i],"Index"=TotI[1:i],"Catch"=TotYieldB[1:i],"Catch.N"=TotCatch[1:i]
                  ,"F"=FullF[1:i],"D"=FullD[1:i],"J"=FullJ[1:i],"MSY"=rep(MSY1[1,1],i),"MSY.D"=rep(MSY2[1,1],i),"F_MSY"=rep(MSY1[1,2],i),"F_MSY.D"=rep(MSY2[1,2],i),
                  "BSA"=(sum(BSA[i,])*C.B),"Plant"=(sum(BSA.plant[i,])*C.B),"CSAA" = (CSAA[i,]*C.C),"LSA"=(sum(NAA[i,]*saatage[i,])*C.L))   #,"OQ_MSY"=MSY2[1:i,1],"OQ_F_MSY"=MSY2[1:i,2]
  
  # Total Plant surface area
  for(j in c(1:nmonths)){
    plant_tot[j]<-sum(BSA.plant[j,1:nages])*C.B
  }
  
  # Total box surface area
  for(j in c(1:nmonths)){
    BSA_tot[j]<-sum(BSA[j,1:nages])
  }
  # Total Live Surface area
  for(j in c(1:nmonths)){
    LSA_tot[j]<-sum(NAA[j,1:nages]*saatage[j,1:nages])*C.L 
  }
  # Total number of animals 
  for(j in c(1:nmonths)){
    NAA_tot[j]<-sum(NAA[j,1:nages])
  }
  # Total number of animals 
  for(j in c(1:nmonths)){
    CSAA_tot[j]<-sum(CSAA[j,])
  }
  
  
  #BSA_tot<-rbind(BSA_tot,BSA_new)
  surfacearea<- data.frame("months"=seq(1,length(TotI)), "BSA" = BSA_tot[1:i], "Plant"=plant_tot[1:i], "CSA" = (CSAA_tot[1:i]),"LSA"=LSA_tot[1:i],
                           "Mort"=Mlo,"Dermo"=F.target,"Seed"=F.target,"F.target"=F.target,"NAA"=NAA_tot[1:i],"Recruits"=NAA[1:i,1])
  # surfacearea$CSA<-surfacearea$CSA*C.C # correct for cultch here to not mess with numebr of dimensions in df
  
  if(Unfished) {
    uf1<-Unfished.pop(nages=nages,nmonths=nmonths,Mat=Mat,med.recr=med.recr
                      ,Recsd=Recsd,steepness=steepness,Nzero=Nzero,type=type,M=M,D.B=D.B,D.C=D.C,C.L=C.L,C.C=C.C,C.B=C.B,st.pop=st.pop
                      ,selectivity.I=selectivity.I,catchab=catchab)
    #return 
    return(list("Simulated.population"=rez,"Unfished.Status"=uf1,"Catch.at.age"=Catch.at.age
                ,"Age.comp"=Age.comp,"inputs"=inputs,"surfacearea"=surfacearea))    
  } else return(list("Simulated.population"=rez,"Catch.at.age"=Catch.at.age,"Age.comp"=Age.comp
                     ,"inputs"=inputs,"Rec_start"=recovery_start,"Rec_time"=recovery_time,"Nofish_time"=Nofish_time,"surfacearea"=surfacearea))
}
#################################################################################

#################################################################################
Solve.F<-function(F1=0.1,selectivity=NULL,M=NULL,Wt=NULL,N.age=NULL,catch1=NULL)
{
  for (n in 1:20) 
  {
    delta=F1*0.0001
    fofFa=sum(((F1*selectivity)/(F1*selectivity+M))*N.age(1-exp(-1*(F1*selectivity+M))))
    fprimeFhi=sum((((F1+delta)*selectivity)/((F1+delta)*selectivity+M))*N.age*(1-exp(-1*((F1+delta)*selectivity+M))))
    fprimeFlo=sum((((F1-delta)*selectivity)/((F1-delta)*selectivity+M))*N.age*(1-exp(-1*((F1-delta)*selectivity+M))))	        
    fofF=fofFa-catch1      
    if(!is.finite(fofF)) {print(paste("ERROR fofFa = ",fofFa," catch1 = ",catch1,sep=""))}
    fprimeF=(fprimeFhi-fprimeFlo)/(2.0*delta)
    F1=F1-(fofF/fprimeF)
    if(fofF<abs(0.0001)) {break}
    if (F1<0) F1=0.5*(F1+(fofF/fprimeF))
    if (F1>3) F1=3
  } #end newton raphson for loop
  return(F1)
}

Solve.D<-function(F2=0.1,selectivity=NULL,M=NULL,Wt=NULL,N.age=NULL,catch1=NULL)
{
  for (n in 1:20) 
  {
    delta=F2*0.0001
    fofFa=sum(((F2*selectivity)/(F2*selectivity+M))*N.age(1-exp(-1*(F2*selectivity+M))))
    fprimeFhi=sum((((F2+delta)*selectivity)/((F2+delta)*selectivity+M))*N.age*(1-exp(-1*((F2+delta)*selectivity+M))))
    fprimeFlo=sum((((F2-delta)*selectivity)/((F2-delta)*selectivity+M))*N.age*(1-exp(-1*((F2-delta)*selectivity+M))))	        
    fofF=fofFa-catch1      
    if(!is.finite(fofF)) {print(paste("ERROR fofFa = ",fofFa," catch1 = ",catch1,sep=""))}
    fprimeF=(fprimeFhi-fprimeFlo)/(2.0*delta)
    F2=F2-(fofF/fprimeF)
    if(fofF<abs(0.0001)) {break}
    if (F2<0) F2=0.5*(F2+(fofF/fprimeF))
    if (F2>3) F2=3
  } #end newton raphson for loop
  return(F2)
}
###################################################################################


###################################################################################
max.ypr<-function(F,D,J,M,Mat,selectivity.F,N0,type,med.recr,steep,SA.1,rec.a.sa,rec.b.sa){ 
  #function to be minimized for MSY calculation
  nvec<-c()
  nvec[1]<-1 #set the population to one so you don't have to rescale
  Faa<-F*selectivity.F
  Z<-M+Faa #total mortality 
  for(i in 2:(length(Mat)-1)){nvec[i]<-nvec[i-1]*exp(-Z[i])} #ages up to oldest   
  nvec[length(Mat)]<-nvec[(length(Mat)-1)]*exp(-Z[length(Mat)])*(1/(1-exp(-Z[length(Mat)]))) #plus group!
  SSNR<-sum(Mat*nvec)  
  
  nvec0<-c()
  nvec0[1]<-1 #set the population to one so you don't have to rescale
  for(i in 2:(length(Mat)-1)){nvec0[i]<-nvec0[i-1]*exp(-M[i])}   
  nvec0[length(Mat)]<-nvec0[(length(Mat)-1)]*exp(-M[length(Mat)])*(1/(1-exp(-M[length(Mat)]))) #plus group!
  SSNR0<-sum(Mat*nvec0)    
  
  
  if(type==1) { #Beverton-Holt stock recruit
    #### CHANGE THIS TO BE SURFACE AREA
    
    R0<-((1/SSNR)*N0)    
    Reca<-(N0/R0)*((1-steep)/(4*steep)) #Bev-Holt alpha parm (See Mangel et al 2010 Fish and Fisheries)
    Recb<-(5*steep-1)/(4*steep*R0)  #Bev-Holt Beta    
    SSNhat<-(SSNR-Reca)/Recb
    Rhat<-SSNhat/(Reca+Recb*SSNhat)  
  } else if(type==2){ #Bev Holt SA
    beta<-rec.b.sa
    alpha<-rec.a.sa
    SSNhat<-(SA.1-alpha)/beta
    alpha<-rec.a.sa
    beta<-rec.b.sa
    Rhat<-(((alpha*(SA.1-shell.correct))/(1+(beta*(SA.1-shell.correct))))*1e10)/12# recall we have to do this here because maturity isn't bein used to limit recruitment
  } else (Rhat<-med.recr) #Random recruitment - just use the median for the month being evaluated
  #### PUT IN SOMETHING HERE FOR TYPE===3 TO REPRESENT THE NEW CURVE
  yield<-sum((Faa/Z)*nvec*(1-exp(-1*Z)))
  
  #Diagnostics
  #  print(paste("F=",F," R*=",Rhat," Yield=",yield," MSY=",yield*Rhat)) #," S*=",SSNhat
  return(yield*Rhat)
}
###################################################################################

###################################################################################
#SSN1<-SSN[1];steep=steepness;N0=Nzero;Mat=Mat[1,];Wt=Wt[1,];M=M[1,];type=type;med.recr=med.recr.in[1];selectivity.F=selectivity.F[1,];

MSY<-function(M=NULL,Mat=NULL,steep=NULL,selectivity.F=NULL
              ,N0=NULL,type=NULL,med.recr=NULL)
{      
  #We need to calculate a few things in here - given steepness and  Nzero - the starting SSN
  #we need to derive Rzero which is a function of several arguments...
  #use the resulting alpha and beta to get the formula for the equilibrium (i.e. placing you 
  #somewhere on the equilibrium yield to F curve). This goes into an optimization over F 
  #that maximizes the yield per recruit times the recruitment implied by the equilibrium spawning stock.
  # uses max.ypr<-function(F,M,Wt,Mat,nvec,selectivity.F,type,med.recr,Rhat)      
  MSY<-optimize(f=max.ypr,interval=c(0,5),M=M,Mat=Mat,selectivity.F=selectivity.F,N0=Nzero,type=type,med.recr=med.recr.in[1],steep=steepness,SA.1=SA.1,rec.a.sa=rec.a.sa,rec.b.sa=rec.b.sa,maximum=TRUE) # REMOVED WT 
  return(data.frame("MSY"=MSY$objective,"F.MSY"=MSY$maximum))
}

### MSY ATTRIBUTED TO DERMO

MSY.2<-function(M=NULL,Mat=NULL,steep=NULL,Selectivity.D=NULL
                ,N0=NULL,type=NULL,med.recr=NULL)
{      
  #We need to calculate a few things in here - given steepness and  Nzero - the starting SSN
  #we need to derive Rzero which is a function of several arguments...
  #use the resulting alpha and beta to get the formula for the equilibrium (i.e. placing you 
  #somewhere on the equilibrium yield to F curve). This goes into an optimization over F 
  #that maximizes the yield per recruit times the recruitment implied by the equilibrium spawning stock.
  # uses max.ypr<-function(F,M,Wt,Mat,nvec,selectivity.F,type,med.recr,Rhat)      
  MSY<-optimize(f=max.ypr,interval=c(0,5),M=M,Mat=Mat,selectivity.F=Selectivity.D,N0=Nzero,type=type,med.recr=med.recr.in[1],steep=steepness,SA.1=SA.1,rec.a.sa=rec.a.sa,rec.b.sa=rec.b.sa,maximum=TRUE) # REMOVED WT 
  return(data.frame("MSY.2"=MSY$objective,"F.MSY.2"=MSY$maximum))
}
###################################################################################

###################################################################################
add.error<-function(evar){
  ## the result of this function should be multiplied by your variable to add error ##
  #add log normal error to a value evar is the variance in your error term
  b.correct=(evar^2)/2 #correct the bias in logging the normal pars
  return(exp(rnorm(1,0,evar)-b.correct))
}
################################################################################ 

###################################################################################
get.error.ac<-function(ac,prev.err,evar){
  ## add an autocorrelated error term to a value method from Deroba & Bence 2012.
  ac.err<-ac*prev.err+sqrt(1-ac^2)*rnorm(1,0,evar)
  ## the result of this function can now be used as the error term in add.error.ac
  return(ac.err)
}
################################################################################

###################################################################################
add.error.ac<-function(ac.err,evar){
  ## the result of this function should be multiplied by your variable to add error ##
  #add log normal error to a value evar is the variance in your error term
  b.correct=(evar^2)/2 #correct the bias in logging the normal pars
  return(exp(ac.err-b.correct))
}
################################################################################ 

###################################################################################
fill<-function(x,n,m){
  # a function to make a vector out of a constant and a matrix out of a vector
  if(length(x)==1) x<-rep(x,m) #make a vector out of a constant    
  if(is.matrix(x)==FALSE) return(matrix(rep(x,n),ncol=m,nrow=n,byrow=TRUE)) #make a matrix out of a vector
  #A little ticky when you are running the unfished algorith, and need to get 1000 months out
  #of a time varying vector - this will hold the last month's values constant for the remaining months 
  if(dim(x)[1]<n) { 
    return(rbind(x,matrix(rep(x[dim(x)[1],],(n-dim(x)[1])),ncol=m,nrow=(n-dim(x)[1]),byrow=TRUE))) 
  }
  
  if(dim(x)[1]==n & dim(x)[2]==m) return(x)  
}
###################################################################################

###################################################################################
Unfished.pop<-function(
  nages=NULL  #number of ages recquired
  ,nmonths=NULL
  ,Mat=NULL #vector of maturity at age
  ,Wt=NULL #vector of weight at age
  ,med.recr=NULL #median recruitment for constant or purely random recruitment
  ,Recsd=NULL #standard deviation of the recruitment error parameter
  ,steepness=NULL      #stock recruit parameterFD
  ,Nzero=NULL        #stock recruit parameter
  ,type=NULL     #1 for Bev-Holt, 2 for Ricker, NULL for random
  ,M=NULL #natural mortality
  ,D.B=NULL # Disarticulation of boxes
  ,D.C=NULL # Decay of cultch 
  ,C.L=NULL # Correction factor for live
  ,C.C=NULL # Correction factor for cultch
  ,C.B=NULL # Correction factor for boxes
  ,st.pop=NULL #population in month one
  ,selectivity.I=NULL #selx for index   
  ,catchab=NULL #index catchability
) {
  #_______________________________________________________________________________
  #       Function to generate stable population size with prespecified          #
  # 		parameters                                                        #
  #______________________________________________________________________________#
  #	Code modified from JJD 4/3/13
  #_______________________________________________________________________________
  
  M<-fill(M,nmonths,nages)
  Mat<-fill(Mat,nmonths,nages)
  selectivity.I<-fill(selectivity.I,nmonths,nages)
  selectivity.F<-fill(selectivity.F,nmonths,nages)
  Selectivity.D<-fill(Selectivity.D,nmonths,nages)
  selectivity.J<-fill(selectivity.J,nmonths,nages)
  
  SSN<-c() #SSN in each month
  NAA<-matrix(st.pop,nrow=1,ncol=nages)  #numbers at age in each month
  NAA.old<-matrix(st.pop,nrow=1,ncol=nages)  #numbers at age in each month
  SSNAA<-matrix(NA,nrow=1,ncol=nages) #SSN at age in each month
  TotN<-c() #sum of true abundance at age in each month
  
  SSN.old=sum(NAA[1,]*Mat[1,])      #SSN month one  
  TotN.old<-sum(NAA[1,])
  
  for (j in 2:1000) #month loop;
  {
    if(j<=nmonths) i<-j
    if(j>nmonths) i<-nmonths
    if(length(steepness)>0 & length(Nzero)>0) {
      med.rec<-stock.recruit(SSN.old,steep=steepness,N0=Nzero,Mat=Mat[i,],M=M[i,],type=type,Rscale=Rscale,
                             SA.1=SA.1,rec.a.sa=rec.a.sa,rec.b.sa=rec.b.sa,shell.correct=shell.correct)
    } else {
      med.rec<-med.recr[i]
    } #otherwise jsut use the constant recruitment par from user	
    NAA[1,1]=med.rec*add.error(Recsd)  #bias correction is handled in the function e.g. exp(rnorm(1,0,Recsd)-b.correct)
    NAA[1,(2:(nages-1))]=NAA.old[1,(1:(nages-2))]*exp(-M[(i-1),(1:(nages-2))])  #Fill abundance at age matrix on the diagonals
    NAA[1,nages]=(NAA.old[1,(nages-1)]*exp(-M[(i-1),(nages-1)]))+(NAA.old[1,nages]*exp(-M[(i-1),nages])) #plus group      
    SSNAA[1,]=NAA[1,]*Mat[i,]#SSN at age in each month
    
    ## SAMAA specified earlier as M + F2
    BAA[i,(1:nages)]=(NAA[(i-1),(1:nages)]*(1-exp(-M[(i-1),(1:nages)]))*2*C.L)+(BAA[(i-1),(1:nages)]*exp(-D.B)) # instead of SAMAA use M because no fishing in this scenario PS multiplied by 2 to get interior and exterior of shell 
    ## Part 1- so individuals that would have died and turned to boxes.... 
    ## Part 2- Boxes from last month that would contribute to new boxes and disarticulation rate
    BSA[i,(1:nages)]=BAA[i,1:nages]*saatage[i,1:nages]*C.B
    ## Part 3- convert to surface area
    ## Now do it for cultch
    ## Part 1- take boxes that disarticulated last month
    ## Part 2- take cultch that hasn't decayed yet
    CSAA.1[i,1:nages]<-(BAA[i-1,1:nages]*saatage[i-1,1:nages]*(1-exp(-D.B))*C.B*C.C)+((CSAA.1[i-1,1:nages])*exp(-D.C))
    #CSAA[i,]<-sum(BSA[i-1,1:nages]*(1-exp(-D.B))) ## Cultch is equal to what shell has disarticulated
    CSAA[i,]<-sum(CSAA.1[i,1:nages])-(sum(NAA[i,1:nages]*saatage[i,1:nages])*C.L*0.18)-(sum(BSA[i,(1:nages)])*0.18)
    #CSAA[i,]<-CSAA[i,]-(sum(NAA[i,1:nages]*saatage[i,1:nages])*C.L*0.18)-(sum(BSA[i,(1:nages)])*0.18)#   ## NOW GET TOTAL SURFACE AREA AVAILABLE WITH CORRECTION FACTORS 
    ## NOW GET TOTAL SURFACE AREA AVAILABLE WITH CORRECTION FACTORS 
    SA.1<-sum(NAA[i,1:nages]*saatage[i,1:nages]) # don't forget to correct for live surface area
    SA.1<-SA.1+(sum(BSA[i,1:nages])) # add box surfacea area, no correction
    SA.1<-SA.1+(CSAA[i,]) # add cultch surface area with correction for cultch 
    
    
    
    SSN=sum(SSNAA[1,]) #Calculate true total SSN in each month (sum across ages)
    if (SSN<1) SSN=0 #Once SSN is less 1 SSN is set to zero and the stock collapses.
    TotN=sum(NAA[1,])  #Total true abundance in each month (sum across ages)
    if((TotN-TotN.old)<1) {#equilibrium reached
      return(list("Unfished.N"=TotN,"Unfished.SSN"=SSN,"Unfished.Age.Comp"=(NAA[1,]/TotN)))
      #return(data.frame("Unfished.N"=TotN,"Unfished.SSN"=SSN,"Unfished.Age.Comp"=(NAA[1,]/TotN)))
    }
    TotN.old<-TotN
    SSN.old<-SSN
    NAA.old<-NAA
  } #end month loop
  #if you got here, no stability within 1000 months
  return(data.frame("Unfished.N"=NA,"Unfished.SSN"=NA,"Unfished.Age.Comp"=NA))
}
#################################################################################
###################################################################################
Start.virgin.pop<-function(
  nages=NULL  #number of ages recquired
  ,nmonths=NULL
  ,Mat=NULL #vector of maturity at age
  ,Wt=NULL #vector of weight at age
  ,med.recr=NULL #median recruitment for constant or purely random recruitment
  ,Recsd=NULL #standard deviation of the recruitment error parameter
  ,steepness=NULL      #stock recruit parameter
  ,Nzero=NULL        #stock recruit parameter
  ,type=NULL     #1 for Bev-Holt, 2 for Ricker, NULL for random
  ,M=NULL #natural mortality
  ,D.B=NULL # Box disarticulation rate
  ,D.C=NULL # Cultch decay rate 
  ,C.L=NULL # Correction factor for live
  ,C.C=NULL # Correction factor for cultch
  ,C.B=NULL # Correction factor for boxes
  ,st.pop=NULL #population in month one
  ,selectivity.I=NULL #selx for index   
  ,catchab=NULL #index catchability
) {
  #_______________________________________________________________________________
  #       Function to generate stable population size with prespecified          #
  # 		parameters                                                        #
  #______________________________________________________________________________#
  #	Code modified from JJD 4/3/13
  #_______________________________________________________________________________
  
  M<-fill(M,nmonths,nages)
  Mat<-fill(Mat,nmonths,nages)
  selectivity.I<-fill(selectivity.I,nmonths,nages)
  selectivity.F<-fill(selectivity.F,nmonths,nages)
  Selectivity.D<-fill(Selectivity.D,nmonths,nages)
  selectivity.J<-fill(selectivity.J,nmonths,nages)
  SSN<-c() #SSN in each month
  NAA<-matrix(st.pop,nrow=1,ncol=nages)  #numbers at age in each month
  NAA.old<-matrix(st.pop,nrow=1,ncol=nages)  #numbers at age in each month
  SSNAA<-matrix(NA,nrow=1,ncol=nages) #SSN at age in each month
  TotN<-c() #sum of true abundance at age in each month
  
  SSN.old=sum(NAA[1,]*Mat[1,])      #SSN month one  
  TotN.old<-sum(NAA[1,])
  
  for (j in 2:1000) #month loop;
  {
    if(j<=nmonths) i<-j
    if(j>nmonths) i<-nmonths
    if(length(steepness)>0 & length(Nzero)>0) {
      med.rec<-stock.recruit(SSN.old,steep=steepness,N0=Nzero,Mat=Mat[i,],M=M[i,],type=type,Rscale=Rscale,
                             SA.1=SA.1,rec.a.sa=rec.a.sa,rec.b.sa=rec.b.sa,shell.correct=shell.correct)
      print(med.rec)
    } else {
      med.rec<-med.recr
      print(med.rec)
    } #otherwise jsut use the constant recruitment par from user	
    
    NAA[1,1]=med.rec*add.error(Recsd)  #bias correction is handled in the function e.g. exp(rnorm(1,0,Recsd)-b.correct)
    NAA[1,(2:(nages-1))]=NAA.old[1,(1:(nages-2))]*exp(-M[(i-1),(1:(nages-2))])  #Fill abundance at age matrix on the diagonals
    NAA[1,nages]=(NAA.old[1,(nages-1)]*exp(-M[(i-1),(nages-1)]))+(NAA.old[1,nages]*exp(-M[(i-1),nages])) #plus group      
    SSNAA[1,]=NAA[1,]*Mat[i,] #SSN at age in each month
    SSNAA[1,]=NAA[1,]*Mat[i,]#SSN at age in each month
    
    ## SAMAA specified earlier as M + F2
    BAA[i,(1:nages)]=(NAA[(i-1),(1:nages)]*(1-exp(-M[(i-1),(1:nages)]))*2*C.L)+(BAA[(i-1),(1:nages)]*exp(-D.B)) # instead of SAMAA use M because no fishing in this scenario PS multiplied by 2 to get interior and exterior of shell 
    ## Part 1- so individuals that would have died and turned to boxes.... 
    ## Part 2- Boxes from last month that would contribute to new boxes and disarticulation rate
    BSA[i,(1:nages)]=BAA[i,1:nages]*saatage[i,1:nages]*C.B
    ## Part 3- convert to surface area
    ## Now do it for cultch
    ## Part 1- take boxes that disarticulated last month
    ## Part 2- take cultch that hasn't decayed yet
    CSAA.1[i,1:nages]<-(BAA[i-1,1:nages]*saatage[i-1,1:nages]*(1-exp(-D.B))*C.B*C.C)+((CSAA.1[i-1,1:nages])*exp(-D.C))
    #CSAA[i,]<-sum(BSA[i-1,1:nages]*(1-exp(-D.B))) ## Cultch is equal to what shell has disarticulated
    CSAA[i,]<-sum(CSAA.1[i,1:nages])-(sum(NAA[i,1:nages]*saatage[i,1:nages])*C.L*0.18)-(sum(BSA[i,(1:nages)])*0.18)
    #CSAA[i,]<-CSAA[i,]-(sum(NAA[i,1:nages]*saatage[i,1:nages])*C.L*0.18)-(sum(BSA[i,(1:nages)])*0.18)#   ## NOW GET TOTAL SURFACE AREA AVAILABLE WITH CORRECTION FACTORS 
    ## NOW GET TOTAL SURFACE AREA AVAILABLE WITH CORRECTION FACTORS 
    SA.1<-sum(NAA[i,1:nages]*saatage[i,1:nages])*C.L # don't forget to correct for live surface area
    SA.1<-SA.1+(sum(BSA[i,1:nages])) # add box surfacea area, no correction
    SA.1<-SA.1+(CSAA[i,]) # add cultch surface area with correction for cultch 
    
    
    SSN=sum(SSNAA[1,]) #Calculate true total SSN in each month (sum across ages)
    if (SSN<1) SSN=0 #Once SSN is less 1 SSN is set to zero and the stock collapses.
    TotN=sum(NAA[1,])  #Total true abundance in each month (sum across ages)
    if((TotN-TotN.old)<1) {#equilibrium reached
      return(list("Unfished.N"=TotN,"Unfished.SSN"=SSN,"Unfished.Age.Comp"=(NAA[1,]/TotN)))
      #return(data.frame("Unfished.N"=TotN,"Unfished.SSN"=SSN,"Unfished.Age.Comp"=(NAA[1,]/TotN)))
    }
    TotN.old<-TotN
    SSN.old<-SSN
    NAA.old<-NAA
  } #end month loop
  #if you got here, no stability within 1000 months
  return(list(data.frame("Unfished.N"=NA,"Unfished.SSN"=NA,"Unfished.Age.Comp"=NA),"Virgin.pop"=NAA))
}
#################################################################################

################################################################################
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


################################################################################
test.arg<-function(arg1){
  if(!(length(arg1)>0)){print(str(arg1))}
  return((length(arg1)>0)) #return false if null argument
}
################################################################################
################################################################################
reversibly.standardize<-function(x){
  m1<-mean(na.omit(x))
  s1<-sd(na.omit(x))
  vec<-(na.omit(x)-m1)/s1
  return(data.frame("st.x"=vec,"mean"=m1,"sd"=s1))
}
################################################################################
################################################################################
performance<-function(dat1,advice){
  #track some useful NSTOW performance metrics
  status<-dat1$SSN[length(dat1$SSN)]/dat1$Unfished.SSN
  if(length(status)==0) status=0.0
  rez<-data.frame("Stability"=abs(sum(diff(dat1$Index)))/sum(dat1$Index)  
                  ,"Status"=status
                  ,"Advice"=advice$advice,"cor"=advice$correlation,"Extinct"=0)
  return(rez)
}
################################################################################
fill.maturity<-function(n){
  #fill a maturity vector
  set.seed(get.seed()) 
  mat<-rbeta(n,2,5)
  return(cumsum(sort(mat))/sum(mat))
}
################################################################################
fill.weight<-function(n){
  #fill a maturity vector
  set.seed(get.seed()) 
  wt<-rbeta(n,2,5)
  return(cumsum(sort(wt)))
}
################################################################################
fill.selx<-function(n){
  #fill a maturity vector
  set.seed(get.seed())  
  slx<-rbeta(n,2,5)
  return(sort(slx))
}
################################################################################
################################################################################
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
################################################################################
################################################################################
################################################################################
# logit and inverse logit
logitx<-function(x){
  if(x>=1 | x<0) {
    return("MUST BE BETWEEN 0 & 1 FOOL!")
  } else {return(log((x/(1-x))))}
}
elogitx<-function(x){
  return((exp(x)/(1+exp(x))))
}
################################################################################

################################################################################
plot.SPAWN<-function(TotN,recruits,SSN,SSN_m,Ntarg,Nthresh,TotYieldB,st.pop,nages,NAA
                     ,FullF,steepness){
  {
    #print(TotN)
    #plot.new()
    resize.win(Width=15, Height=10)
    l.out<-matrix(c(1,1,2,3),2,2, byrow=T)
    nf <- layout(l.out)
    #layout.show(nf)
    par(mar=c(4,5,4,4)+0.1,cex.axis=1.5,cex.lab=2,cex.main=2)
    
    plot(TotN~seq(1,length(TotN))
         ,type="l"
         ,lwd=3
         ,ylim=c(0,1.2*max(na.omit(TotN)))
         ,ylab="Number"
         ,xlab="Month"
         ,main="Number of animals and recruits"
    )
    par(new=TRUE)
    ytop<-1.2*max(na.omit(recruits))    
    plot(recruits~seq(1,length(recruits))
         ,lwd=2
         ,type="l"
         ,col="red"
         ,axes=FALSE
         ,ylab=""
         ,xlab=""
         ,ylim=c(0,ytop)
    )
    axis(4,pretty(c(0,ytop)))
    mtext("Recruits",side=4,line=3,cex=1.5)
    legend(
      x="topright"
      ,legend=c("Total N","Recruits",paste("Steepness = ",round(steepness,3),sep=""))
      ,lwd=c(3,2,0)
      ,lty=c(1,1,1)
      ,col=c("black","red","white")
      ,bty="n"
    )
    
    #plot of SSN and yield    
    plot(SSN~seq(1,length(SSN))
         ,type="l"
         ,lwd=3
         ,ylim=c(0,max(na.omit(SSN_m)))
         ,ylab="SSN"
         ,xlab="Month"
         ,main="Spawning stock numbers and yield"
    )
    lines(SSN_m~seq(1,length(SSN_m))
          ,lty=1
          ,lwd=2
          ,col="gray80"
    )
    abline(h=Ntarg,col="green",lty=2,lwd=2)
    abline(h=Nthresh,col="red",lty=2,lwd=2)    
    par(new=T)
    ytop<-1.2*max(na.omit(TotYieldB))
    plot(TotYieldB~seq(1,length(TotYieldB))
         ,lwd=2
         ,col="blue"
         ,type="l"
         ,axes=FALSE
         ,ylab=""
         ,xlab=""
         ,ylim=c(0,ytop)
    )
    axis(side=4,at=pretty(c(0,ytop)))
    mtext("Yield",side=4,line=3,cex=1.5)
    legend(
      x="topright"
      ,legend=c("True SSN","Apparent SSN",paste("Yield at F = ",round(mean(FullF),3),sep="")
                ,"Biomass target","Biomass threshold")
      ,lwd=c(3,2,2,2,2)
      ,lty=c(1,1,1,2,2)
      ,col=c("black","gray80","blue","green","red")
      ,bty="n"      
    )
    
    plot((st.pop/sum(st.pop))~seq(1:nages)
         ,lwd=1
         ,type="l"
         ,col="#00000099"
         ,ylab="Frequency"
         ,xlab="Age"
         ,main="Age Distribution"
    )
    par(new=T)
    plot((NAA/sum(NAA))~seq(1:nages)
         ,lwd=1
         ,type="l"
         ,col="#CC000090"
         ,axes=FALSE
         ,ylab=""
         ,xlab=""
         ,main=""
    )  
    legend(
      x="topright"
      ,legend=c("Starting","Ending")
      ,lwd=c(3,2)
      ,lty=c(1,1)
      ,col=c("#00000099","#CC000090")
      ,bty="n"      
    )
    
  } #plot condition  
  
  
}

##########################################################################################################































































