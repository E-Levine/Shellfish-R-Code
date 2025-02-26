#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
change.growthS<-function(nages,sigk=0.04,sigL=30,sigt0=0){
  #function to change the growth curve a bit given the number of ages and parameters for the 
  #stdev of the change in L infinity and k
  #len<-95.46*(1-exp(-0.051*(seq(1,nages)+.041))) #derived from data length in mm average over time
  len<-(125+rnorm(1,0,sigL))*(1-exp((-0.24+rnorm(1,0,sigk))*(seq(1,nages)+(0.2+rnorm(1,0,sigt0))))) 
  return(len)
}
#Surfclam specific functions
change.growthN<-function(nages,sigk=0.04,sigL=30,sigt0=0){
  #function to change the growth curve a bit given the number of ages and parameters for the 
  #stdev of the change in L infinity and k
  #len<-95.46*(1-exp(-0.675*(seq(1,nages)+.49))) #derived from data length in mm average over time
  len<-(125+rnorm(1,0,sigL))*(1-exp((-0.24+rnorm(1,0,sigk))*(seq(1,nages)+(0.2+rnorm(1,0,sigt0))))) 
  return(len)
}



### Selectivity 
change.selx.F<-function(len){
  #function to calculate the selectivity at age given the len at age vector  
  selectivity.F<-1-(1/(1+exp(-7.63+(0.105*len)))) #selx for fishery from Thorarinsdottir and Jacobson, 2005.     
  selectivity.F[1]<-0  #Need to do this or large recruitment classes will swamp last Months total population 
  return(selectivity.F)      
}

change.selx.D<-function(len){
  #function to calculate the selectivity at age given the len at age vector  
  Selectivity.D<-1-(1/(1+exp(-7.63+(0.105*len)))) #selx for fishery from Thorarinsdottir and Jacobson, 2005.     
  Selectivity.D[1]<-0  #Need to do this or large recruitment classes will swamp last Months total population 
  return(Selectivity.D)      
}

change.selx.J<-function(len){
  #function to calculate the selectivity at age given the len at age vector  
  Selectivity.J<-1-(1/(1+exp(-7.63+(0.105*len)))) #selx for fishery from Thorarinsdottir and Jacobson, 2005.     
  Selectivity.J[1]<-0  #Need to do this or large recruitment classes will swamp last Months total population 
  return(Selectivity.J)      
}




################################################################################### 
stock.recruit<-function(SSN1,M=NULL,Wt=NULL,Mat=NULL,steep=NULL
                        ,N0=NULL,type=NULL,Rnoise=0.2,Rscale=NULL,SA.1=NULL,rec.a.sa=NULL,rec.b.sa=NULL,shell.correct=NULL)
{      
  #We need to calculate a few things in here - given steepness and  Nzero the starting SSN
  #we need to derive Rzero which is a function of several arguments...
  if(F){  #check inputs
    print(c("SSN1=",SSN1))
    print(c("M="  ,M))
    print(c("Wt=",Wt))
    print(c("Mat=",Mat))
    print(c("steep=",steep))
    print(c("N0=",N0))
    print(c("type=",type))
    print(c("Rnoise=",Rnoise))
    print(c("Rscale=",Rscale) )
    print(c("SA.1"=SA.1))
    print(c("rec.a.sa"=rec.a.sa))
    print(c("shell.correct"=shell.correct))
    print(c("rec.b.sa"=rec.b.sa))
  }
  interval = 160
  
  dim.var<-dim(Mat)[2]  
  if(is.null(dim.var)) dim.var<-length(Mat) #Mat can be an array or vector...
  #print(c(dim.var,dim(Mat)[2]))
  if(length(M)==1) M<-fill(M,n=dim.var,m=1)
  
  #preallocation saves some time
  nvec<-rep(1,dim.var)     
  for(i in 2:(dim.var-1)){nvec[i]<-nvec[i-1]*exp(-M[i-1])}   
  
  nvec[dim.var]<-nvec[(dim.var-1)]*exp(-M[dim.var-1])*(1/(1-exp(-M[dim.var]))) #plus group!
  SSNR<-sum(Mat*nvec) ### REMOVED WEIGHT FROM THIS AS WE'RE JUST DOING IT BASED ON N WHICH IS RECORDED IN N VEC  
  
  #print(c("SSNR=",SSNR))    
  #print(c("nvec=",nvec))
  alpha<-rec.a.sa
  beta<-rec.b.sa
  #return(((((SA.1-shell.correct)/(alpha+beta*SA.1)))*1e9)/12)# divided by twelve because recruitment is done at each time step but this is an annual value
  return((((alpha*(SA.1-211))/(1+(beta*(SA.1-211))))*1e9)/12)# multipled by 1e9 to get out of hectares and back to n
  
}
###################################################################################
#testing the proper scale of recruitment
#resize.win()
#evar<-2.
#b.correct=(evar^2)/2
#temp<-exp(rnorm(1000,0,evar)-b.correct)
#plot(temp)
#print(c("med=",round(median(temp),3),"mean=",round(mean(temp),3),"% < 0.02 =",length(temp[temp<0.02])/1000))
###################################################################################
add.errorR<-function(evar){
  ## the result of this function should be multiplied by your variable to add error ##
  #add log normal error to a value evar is the variance in your error term
  b.correct=(evar^2)/2 #correct the bias in logging the normal pars
  return(exp(rnorm(1,0,evar)-b.correct))
}
################################################################################  
miniSpawn<-function(recs=NULL #number to start the cohort with.  A function to generate the equilibrium age comp
                    ,nages=1  #number of ages recquired
                    ,Mat=NULL #vector of maturity at age, or a matrix with nmonths X nages
                    ,Wt=NULL #vector of weight at age, or a matrix with nmonths X nages
                    ,M=NULL #natural mortality can be a constant or vector with nages cells, or a 
                    #matrix with nmonths X nages
                    ,F.dat=NULL #vector of F in each Month
                    ,D.dat=NULL # vector of Dermo F in each Month
                    ,J.dat=NULL # vector of Seed F
                    ,selectivity.F=NULL #selx for fishery can be a vector of selectivity by age, 
                    ,Selectivity.D=NULL #selx for Dermo 
                    ,selectivity.J=NULL #selx for Seed fishery
){    
  #concise version of spawn to take a recruit class and drive it out according to the mortality schedule - a cohort maker
  M<-fill(M,nmonths,nages)
  Wt<-fill(Wt,nmonths,nages)  
  Mat<-fill(Mat,nmonths,nages)
  NAA<-matrix(NA,nrow=1,ncol=nages)  #numbers at age in each Month
  FAA<-matrix(NA,nrow=1,ncol=nages)  #F at age in each Month
  FAA.D<-matrix(NA,nrow=1,ncol=nages)
  FAA.J<-matrix(NA,nrow=1,ncol=nages)
  ZAA<-matrix(NA,nrow=1,ncol=nages) #total mortality at age in each Month
  NAA[1,1]<-recs #population in Month one
  FAA[1,]<-F.dat[1]*selectivity.F #Fishing selectivity
  FAA.D[1,]<-D.dat[1]*Selectivity.D
  FAA.J[1,]<-J.dat[1]*selectivity.J
  ZAA[1,]<-FAA[1,]+FAA.D[1,]+FAA.J[1,]+M[1,]  
  for(i in 2:(nages)) {NAA[1,(i)]=NAA[1,(i-1)]*exp(-ZAA[1,(i-1)])}
  #run out the plus group a bit
  for(i in 1:50) {NAA[1,nages]=(NAA[1,(nages-1)]*exp(-ZAA[1,(nages-1)]))+(NAA[1,nages]*exp(-ZAA[1,nages]))} #plus group
  
  return(NAA)
}
################################################################################  
get.SR.pars<-function(M=NULL,Wt=NULL,Mat=NULL,steep=NULL
                      ,N0=NULL,type=NULL,SA.1.1=NULL,rec.a.sa=NULL,rec.b.sa=NULL){
  
  dim.var<-dim(Mat)[2]  
  if(is.null(dim.var)) dim.var<-length(Mat) #Mat can be an array or vector...
  #print(c(dim.var,dim(Mat)[2]))
  if(length(M)==1) M<-fill(M,n=dim.var,m=1)
  
  nvec<-c()
  nvec[1]<-1 #set the population to one so you don't have to rescale
  for(i in 2:(dim.var-1)){nvec[i]<-nvec[i-1]*exp(-M[i-1])}
  nvec[dim.var]<-nvec[(dim.var-1)]*exp(-M[(dim.var-1)])*(1/(1-exp(-M[dim.var]))) #plus group!
  SSNR<-sum(Mat*nvec)  # REMOVED WEIGHTS O WE ONLY DO IT BASED ON N
  
  if(type==1) { #Beverton-Holt stock recruit   
    R0<-((1/SSNR)*N0)     
    Reca<-(N0/R0)*((1-steep)/(4*steep)) #Bev-Holt alpha parm (See Mangel et al 2010 Fish and Fisheries)
    Recb<-(5*steep-1)/(4*steep*R0)  #Bev-Holt Beta    
    return(data.frame("R0"=R0,"alpha"=Reca,"beta"=Recb))
  } else if(type==2){ #SA Bev Holt
    alpha<-rec.a.sa
    beta<-rec.b.sa
    return(data.frame("SA1.1"=NULL,"alpha"=alpha,"beta"=beta)    )
  }    
}
################################################################################  



################################################################################### 
sa.func<-function(nmonths,nages,M,a.wt,b.wt,Mat,N0,D.B,D.C,C.L,C.C,C.B,steep,Reca,Recb,printsimdat=F){
  M<-fill(M,nmonths,nages)
  nvec0<-matrix(NA,nrow=1,ncol=240)
  nvec0[1,1]<-1 #set the population to one so you don't have to rescale
  for(i in 2:(length(Mat)-1)){nvec0[1,i]<-nvec0[1,i-1]*exp(-M[1,i-1])}   
  nvec0[length(Mat)]<-nvec0[(length(Mat)-1)]*exp(-M[length(Mat)])*(1/(1-exp(-M[length(Mat)]))) #plus group!
  
  #interval<-set.up.SC.run()$interval #frequency of large Month classes
  
  nvec0<-as.vector(nvec0)
  nvec2<-nvec0/(sum(nvec0)) ## proportion
  nvec2<-nvec2*N0 #extrapolate to get number of animals...
  fill.vec<-as.numeric(saatage*nvec2) 
  #N0<-fill.vec/saatage
  ## Get surface area of live animals
  SAV.m<-fill(fill.vec,nmonths,nages)
  SAV<-sum(saatage*nvec2) ## pretend they're all live to initiate it, so it's just whatever it is. 
  
  SAV<- SAV*C.L # add correction factor to live animals
  R0<-N0 ## set R0 equal to number of spawning individuals at beginning
  
  
  #Reca<-(SAV/R0)*((1-steep)/(4*steep)) #Bev-Holt alpha parm (See Mangel et al 2010 Fish and Fisheries)
  #Recb<-(5*steep-1)/(4*steep*R0)  #Bev-Holt Beta  
  #reca<-0.0040
  #recb<-0.0013
  
  #rec<-((SA1/(Reca+Recb*SA1)))
  rec<-nvec2[1]
  
  NAA<-fill(nvec2,nmonths,nages)
  M<-fill(M,nmonths,nages)
  
  ## POP Month 2
  NAA[2,(2:(nages-1))]=NAA[1,(1:(nages-2))]*exp(-M[1,(1:(nages-2))])  #Fill abundance at age matrix on the diagonals
  ## plus group
  NAA[2,length(Mat)]<-NAA[1,(length(Mat)-1)]*exp(-M[1,length(Mat)])*(1/(1-exp(-M[1,length(Mat)]))) #plus group!
  NAA[2,1]<-rec ## KEEP THESE IN NUMBERS, NOW GENERATE THE SURFACE AREA
  ### SOMETHING IS GOING ON HERE
  
  saatage<-as.numeric(saatage)
  saatage<-fill(saatage,nmonths,nages)
  
  BAA<-matrix(NA,nrow=nmonths,ncol=nages)#boxes at age in each Month
  BAA[1,(1:nages)]<-NAA[1,1:nages]*(1-exp(-M[1,1:nages]))*2*C.L # multiplied by 2 to get interior and exterior shells 
  
  BSA<-matrix(NA,nrow=nmonths,ncol=nages)
  BSA[1,1:nages]<-BAA[1,1:nages]*saatage[1,1:nages]*C.B
  
  SA.1<-matrix(NA,nrow=nmonths,ncol=nages)
  SA.1[2,1:nages]<-NAA[2,1:nages]*saatage[2,1:nages]
  
  SA1<-sum(sum(SA.1[2,1:nages]))*C.L
  SA1<-SA1+(sum(sum(BSA[1,1:nages])))
  
  
  #Reca<-(SAV/R0)*((1-steep)/(4*steep)) #Bev-Holt alpha parm (See Mangel et al 2010 Fish and Fisheries)
  #Recb<-(5*steep-1)/(4*steep*R0)  #Bev-Holt Beta  
  #reca<-0.0040
  #recb<-0.0013
  
  #rec<-((SA1/(Reca+Recb*SA1)))
  rec<-nvec2[1]
  
  
  ## POP Month 3
  # Live Animals
  NAA[3,(2:(nages-1))]=NAA[2,(1:(nages-2))]*exp(-M[2,(1:(nages-2))])  #Fill abundance at age matrix on the diagonals
  ## plus group
  NAA[3,length(Mat)]<-NAA[2,(length(Mat)-1)]*exp(-M[2,length(Mat)])*(1/(1-exp(-M[2,length(Mat)]))) #plus group!
  NAA[3,1]<-rec ## KEEP THESE IN NUMBERS, NOW GENERATE THE SURFACE AREA
  
  
  ## Boxes
  BAA[2,(1:nages)]<-NAA[2,1:nages]*(1-exp(-M[2,1:nages]))*2*C.L # multiplied by 2 to get interior and exterior shells 
  
  ## Now add old boxes
  BAA[2,1:nages]<-BAA[2,1:nages]+(BAA[1,1:nages]*exp(-D.B))
  
  ## Box Surface Area
  BSA[2,1:nages]<-BAA[2,1:nages]*saatage[2,1:nages]*C.B
  
  # Cultch
  CSA.1<-matrix(NA,nrow=nmonths,ncol=nages)
  CSA.1[1,1:nages]<-BAA[1,1:nages]*saatage[1,1:nages] * (1-exp(-D.B))*C.B*C.C
  CSA<-(sum(CSA.1[1,1:nages]))-(sum(BSA[2,1:nages])*0.18)-(sum(NAA[3,1:nages]*saatage[3,1:nages])*C.L*0.18) ## Cultch is equal to what shell has disarticulated
  
  # Get total surface area
  SA.1[3,1:nages]<-NAA[3,1:nages]*saatage[3,1:nages]*C.L
  SA1<-sum(SA.1[3,1:nages])
  ## Add Cultch Correction Factor
  SA1<-SA1+(sum(BSA[2,1:nages]))+(CSA)
  
  
  #Reca<-(SAV/R0)*((1-steep)/(4*steep)) #Bev-Holt alpha parm (See Mangel et al 2010 Fish and Fisheries)
  #Recb<-(5*steep-1)/(4*steep*R0)  #Bev-Holt Beta  
  #reca<-0.0040
  #recb<-0.0013
  
  #rec<-((SA1/(Reca+Recb*SA1)))
  rec<-nvec2[1]
  
  ## POP Month 4
  NAA[4,(2:(nages-1))]=NAA[3,(1:(nages-2))]*exp(-M[3,(1:(nages-2))])  #Fill abundance at age matrix on the diagonals
  ## plus group
  NAA[4,length(Mat)]<-NAA[3,(length(Mat)-1)]*exp(-M[3,length(Mat)])*(1/(1-exp(-M[3,length(Mat)]))) #plus group!
  
  NAA[4,1]<-rec ## KEEP THESE IN NUMBERS, NOW GENERATE THE SURFACE AREA
  
  
  ## Boxes
  BAA[3,1:nages]<-NAA[3,1:nages]*(1-exp(-M[3,1:nages]))*2*C.L # multiplied by 2 to get interior and exterior shells 
  
  
  ## Now add old boxes
  BAA[3,1:nages]<-BAA[3,1:nages]+(BAA[2,1:nages]*exp(-D.B))
  
  ## Box Surface Area
  BSA[3,1:nages]<-BAA[3,1:nages]*saatage[3,1:nages]*C.B
  
  # Cultch
  CSA.1[2,1:nages]<-(BAA[2,1:nages]*saatage[2,1:nages]*(1-exp(-D.B))*C.B*C.C)+ (CSA.1[1,1:nages])*exp(-D.C) ## Latter part was already corrected for cultch correction factor from last Month#add newly disarticulated boxes to cultch
  CSA <- (sum(CSA.1[2,1:nages]))-(sum(BSA[3,1:nages])*0.18)-(sum(NAA[4,1:nages]*saatage[4,1:nages])*C.L*0.18)
  
  SA.1[4,1:nages]<-NAA[4,1:nages]*saatage[4,1:nages]*C.L
  SA1<-sum(SA.1[4,1:nages])
  SA1<-SA1+(sum(BSA[3,1:nages]))+(CSA)
  
  
  
  #Reca<-(SAV/R0)*((1-steep)/(4*steep)) #Bev-Holt alpha parm (See Mangel et al 2010 Fish and Fisheries)
  #Recb<-(5*steep-1)/(4*steep*R0)  #Bev-Holt Beta  
  #reca<-0.0040
  #recb<-0.0013
  
  #rec<-((SA1/(Reca+Recb*SA1)))
  rec<-nvec2[1]
  
  ## POP Month 5-1000
  for(j in 5:1000){
    
    NAA[j,(2:(nages-1))]=NAA[j-1,(1:(nages-2))]*exp(-M[j-1,(1:(nages-2))])  #Fill abundance at age matrix on the diagonals
    ## plus group
    NAA[j,length(Mat)]<-NAA[j-1,(length(Mat)-1)]*exp(-M[j-1,length(Mat)])*(1/(1-exp(-M[j-1,length(Mat)]))) #plus group!
    NAA[j,1]<-rec ## KEEP THESE IN NUMBERS, NOW GENERATE THE SURFACE AREA
    
    
    ## Boxes
    BAA[j-1,1:nages]<-NAA[j-1,1:nages]*(1-exp(-M[j-1,1:nages]))*2*C.L # multiplied by 2 to get interior and exterior shells 
    
    # Now add old boxes
    BAA[j-1,1:nages]<-BAA[j-1,1:nages]+(BAA[j-2,1:nages]*exp(-D.B))
    
    # Box surface area
    BSA[j-1,1:nages]<-BAA[j-1,1:nages]*saatage[j-1,1:nages]*C.B
    
    
    # Cultch
    CSA.1[j-2,1:nages]<-((BAA[j-2,1:nages]*saatage[j-2,1:nages]*(1-exp(-D.B)))*C.B*C.C)+ (CSA.1[j-3,1:nages])*exp(-D.C)#add newly disarticulated boxes to cultch
    CSA <- (sum(CSA.1[j-2,1:nages]))-(sum(BSA[j-1,1:nages])*0.18)-(sum(NAA[j,1:nages]*saatage[j,1:nages])*C.L*0.18)
    
    
    # Get total surface area
    SA.1[j,1:nages]<-NAA[j,1:nages]*saatage[j,1:nages]*C.L
    SA1<-sum(SA.1[j,1:nages])
    if(printsimdat){
    if(j==6){
      print(j)
      print(CSA)
      print(sum(BSA[j-1,1:nages]))
      print(SA1)
      print(rec)
      print(sum(NAA[j-1,1:nages]))
    }
    if(j==999){
      print(j)
      print(CSA)
      print(sum(BSA[j-1,1:nages])*C.B)
      print(SA1)
      print(rec)
      print(sum(NAA[j-1,1:nages]))
    }
    }
    ## Add Cultch Correction Factor
    SA1<-SA1+(sum(BSA[j-1,1:nages]))+(CSA)
    
    
    #Reca<-(SAV/R0)*((1-steep)/(4*steep)) #Bev-Holt alpha parm (See Mangel et al 2010 Fish and Fisheries)
    #Recb<-(5*steep-1)/(4*steep*R0)  #Bev-Holt Beta  
    #reca<-0.0040
    #recb<-0.0013
    
    #rec<-((SA1/(Reca+Recb*SA1)))
    rec<-nvec2[1]
    #print(rec)
    
    #NAA[(j+1),1]<-rec # go ahead and add rec to next one... maybe... 
  }
  
  
  New<-matrix(NA,nrow=1000,ncol=240)
  NAA.1<-rbind(NAA,NAA)
  BAA.1<-rbind(BAA,BAA)
  BSA.1<-rbind(BSA,BSA)
  CSA.2<-rbind(CSA.1,CSA.1)
  M.1<-rbind(M,M)
  SA.1.1<-rbind(SA.1,SA.1)
  saatage.1<-rbind(saatage,saatage)
  SAV.m.1<-rbind(SAV.m,SAV.m)
  rec<-NAA[1,1]
  
  ## POP Month 5-1000
  for(j in c(1001:2000)){
    
    
    NAA.1[j,(2:(nages-1))]=NAA.1[j-1,(1:(nages-2))]*exp(-M.1[j-1,(1:(nages-2))])  #Fill abundance at age matrix on the diagonals
    ## plus group
    NAA.1[j,length(Mat)]<-NAA.1[j-1,(length(Mat)-1)]*exp(-M.1[j-1,length(Mat)])*(1/(1-exp(-M.1[j-1,length(Mat)]))) #plus group!
    NAA.1[j,1]<-rec ## KEEP THESE IN NUMBERS, NOW GENERATE THE SURFACE AREA
    
    
    ## Boxes
    BAA.1[j-1,1:nages]<-NAA.1[j-1,1:nages]*(1-exp(-M.1[j-1,1:nages]))*2*C.L # multiplied by 2 to get interior and exterior shells 
    ## SHOULD THESE BE -2 TO GET THOSE INDIVDIUALS FROM LAST Month THAT ARE GROWING UP OR DYING 
    # Now add old boxes
    BAA.1[j-1,1:nages]<-BAA.1[j-1,1:nages]+(BAA.1[j-2,1:nages]*exp(-D.B))
    
    # Box surface area
    BSA.1[j-1,1:nages]<-BAA.1[j-1,1:nages]*saatage.1[j-1,1:nages]*C.B
    
    # Cultch
    CSA.2[j-2,1:nages]<-(BAA.1[j-2,1:nages]*saatage.1[j-2,1:nages]*(1-exp(-D.B))*C.B*C.C)+ (CSA.2[j-3,1:nages])*exp(-D.C)#add newly disarticulated boxes to cultch
    CSA <- (sum(CSA.2[j-2,1:nages]))-(sum(BSA.1[j-1,1:nages])*0.18)-(sum(NAA.1[j,1:nages]*saatage.1[j,1:nages])*C.L*0.18)
    
    
    # Get total surface area
    SA.1.1[j,1:nages]<-NAA.1[j,1:nages]*saatage.1[j,1:nages]*C.L
    SA1.1<-sum(SA.1.1[j,1:nages])
    ## Add Cultch Correction Factor
    SA1.1<-SA1.1+(sum(BSA.1[j-1,1:nages]))+(CSA)
    
    # convert back to mm2?
    #SA1.2<-SA1.1*10e10
    if(printsimdat){
    if(j==1002){
      print(j)
      print(CSA)
      print((sum(CSA.2[j-2,1:nages])))
      print(sum(BSA.1[j-1,1:nages]))
      print(SA1.1)
      #print(sum(BAA.1[j-2,1:nages]*saatage.1[j-2,1:nages]*(1-exp(-D.B))*C.B*C.C))
      #print(sum(CSA.2[j-3,1:nages])*exp(-D.C))
    }
    if(j==1900){
      print(j)
      print(CSA)
      print((sum(CSA.2[j-2,1:nages])))
      print(sum(BSA.1[j-1,1:nages]))
      print(SA1.1)
      #print(sum(BAA.1[j-2,1:nages]*saatage.1[j-2,1:nages]*(1-exp(-D.B))*C.B*C.C))
      #print(sum(CSA.2[j-3,1:nages])*exp(-D.C))
    }
    if(j==1999){
      print(j)
      print(CSA)
      print((sum(CSA.2[j-2,1:nages])))
      print(sum(BSA.1[j-1,1:nages]))
      print(SA1.1)
      #print(sum(BAA.1[j-2,1:nages]*saatage.1[j-2,1:nages]*(1-exp(-D.B))*C.B*C.C))
      #print(sum(CSA.2[j-3,1:nages])*exp(-D.C))
    }
    }
    
    #Reca<-(SAV/R0)*((1-steep)/(4*steep)) #Bev-Holt alpha parm (See Mangel et al 2010 Fish and Fisheries)
    #Recb<-(5*steep-1)/(4*steep*R0)  #Bev-Holt Beta  
    
    rec<-(((rec.a.sa*(SA1.1-shell.correct))/(1+(rec.b.sa*(SA1.1-shell.correct))))*1e9)/12
    #print(rec)
    #rec<-nvec2[1]
    
    #NAA[(j+1),1]<-rec # go ahead and add rec to next one... maybe... 
  }
  
  my.list<-list(((CSA.2[1998,])),BAA.1[1999,],BSA.1[1999,],NAA.1[2000,],SA1.1)
  names(my.list)<-c("CSA","BAA.1","BSA.1","NAA.1","SA1.1")
  return(my.list)
}

#try<-sa.func(nmonths=1000,nages=nages,M=M,a.wt=a.wt,b.wt=b.wt,Mat=Mat,N0=Nzero,D.B=D.B,
#             D.C=D.C,C.L=C.L,C.C=C.C,C.B=C.B,steep=steepness,Reca=rec.a.sa,Recb=rec.b.sa)
###################################################################################












