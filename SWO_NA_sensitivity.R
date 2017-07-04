##><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><><><
## Stock Assessment Input File for BSPSP
## Developed by Henning Winker & Felipe Carvalho (Cape Town/Hawaii)  
##><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>


# required packages
#Just in case
library(gplots)
library(coda)
#library(jagsUI) ### THIS IS NEW
library(rjags)
library(R2jags)
library("fitdistrplus")
#Just in case
#detach("package:jagsUI", unload=TRUE)
#detach("package:R2jags", unload=TRUE)
#dev.off()

################################################
# Define basic object so that they exist:
meanCPUE = Projection = pucl = save.all = FALSE  
KOBE.plot =Biplot= CPUE.plot= FALSE
################################################


# Set Working directory file where to store the results
File = "C:/Work/Research/LargePelagics/ASSESSMENTS/JABBA_development"
# Set working directory for JABBA R source code
JABBA.file = "C:/Work/Research/GitHub/JABBA"
# Date subdirectory
# Set Assessment
assessment = "SWO_SA"


# Choose Sceanrio name for creating a seperate folCder

drops = 7
repl = 2
add = 2
# Drop CPUE observations backward
Scenarios = c("All",paste0("repl",c(1:2)),paste0("add",c(1:2)), paste0("drop",1:(drops))) 

H_Hmsy.sen = B_Bmsy.sen = SR = NULL



for(s in 1:length(Scenarios)){

Scenario = Scenarios[s] 

# Choose model type: 
# 1: Schaefer
# 2: Schaefer with Depensation (CMSY: Froese et al. 2016)  
# 3: Fox
# 4: Fox with Depensation (CMSY: Froese et al. 2016)
# 5: Pella-Tomlinsson  
  Model = 3 #ifelse(s>2,5,1)  # model
   
  Mod.names = c("Schaefer","Schaefer.RecImp","Fox","Fox.RecImp","Pella")[Model]
  
  
  if(Model==5){
  BmsyK = 0.4
  # find inflection point
  ishape = NULL
  # Find shape for  SBmsytoK 
  ishape = seq(0.1,10,0.001)
  
  check.shape =((ishape)^(-1/(ishape-1))-BmsyK)^2
  
  #  Set shape (> 0, with 1.001 ~ Fox and 2 = Schaefer)
  shape =  ishape[check.shape==min(check.shape)] 
  } else {shape=FALSE}
  
  #  Set shape (> 0, with 1.001 ~ Fox and 2 = Schaefer)
  setwd(paste(File))
  dir.create(paste(assessment),showWarnings =FALSE)
  # Load assessment data
  cpue = read.csv(paste0(assessment,"/cpue",assessment,".csv"))#
  se =  read.csv(paste0(assessment,"/se",assessment,".csv"))# use 0.001 if not available 
  catch = read.csv(paste0(assessment,"/catch",assessment,".csv"))
  
  names(cpue)
  names(catch)
  
  #--------------------------------------------------
  # option to exclude CPUE time series or catch year
  #--------------------------------------------------
  for(i in 2:ncol(se)){
  se[,i] = ifelse(se[,i]<0.15,se[,i]+0.1,se[,i])
  }
  
  # Exclude TAI
  if(s!=5){
  cpue = cpue[,-c(10:12)]
  se = se[,-c(10:12)]
  }
  if(s!=4){
  cpue = cpue[,-c(2)]
  se = se[,-c(2)]
  }
  
  if(s>5){ 
    cpue = cpue[,-c(s-2)]
    se = se[,-c(s-2)]
  }
  
  names(cpue)
  ncol(catch)
  ncol(cpue)
  ncol(se)
  
  #------------------------------------------------
  # mean and CV and sd for unfished biomass K (SB0)
  #------------------------------------------------
  mu.K = 200000; CV.K = 1; sd.K=sqrt(log(CV.K^2+1)) 
  K.pr = c(mu.K,sd.K)
  
  #-----------------------------------------------------------
  # mean and CV and sd for Initial depletion level P1= SB/SB0
  #-----------------------------------------------------------
  # Set the initial depletion prior B1/K 
  # To be converted into a lognormal prior (with upper bound at 1.1)
  #BK1.prior = "lnorm"
  # or to be converted into a Beta prior
  #BK1.prior = "beta"
  
  BK1.prior= "lnorm"
  # specify as mean and CV 
  mu.psi = 1 
  CV.psi = 0.05 
  
  #--------------------------------------------------------------
  # Determine estimation for catchability q and observation error 
  #--------------------------------------------------------------
  # Assign q to CPUE
  sets.q = 1:(ncol(cpue)-1) 
  if(s==2)sets.q = c(1,2,2,3,4,5,6)
  if(s==3)sets.q = c(1,2,3,4,4,5,6)
  
  
  
  # Series
  #sets.var = rep(1,ncol(cpue)-1)# estimate the same additional variance error
  sets.var = 1:(ncol(cpue)-1) # estimate individual additional variace
  if(s==2) sets.var= c(1,2,2,3,4,5,6)
  if(s==3) sets.var= c(1,2,2,3,4,5,6)
  
  # As option for data-weighing
  # minimum additional observation error for each variance set (optional choose 1 value for both)
  min.obsE = c(0.1) # Important if SE.I is not availble
  
  # Use SEs for abudance indices (TRUE/FALSE)
  SE.I = TRUE
  
  #><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>
  # Prior specification for Models 1-4, i.e. Schaefer, Fox
  #><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>
  
  #----------------------------------------------------
  # Determine r prior
  #----------------------------------------------------
  # The option are: 
  # a) Specifying a lognormal prior 
  # b) Specifying a resiliance category after Froese et al. (2016; CMSY)
  # Resilience: "Very low", "Low", "Medium", High" (requires r.range = TRUE)
  
  # use range or mean/stdev
  r.range = FALSE # Set to false for mean/stdev specifications
  
  #r.prior = "Low"
  r.prior = c(0.42,0.37) # as range with upper and lower bound of lnorm prior
  
  #><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>>
  # Process Error
  #><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>>
  #Estimate set sigma.proc == True
  sigma.proc = TRUE
  if(sigma.proc == TRUE) igamma = c(4,0.01) # standard process with mu.s = 0.05, CV = 0.23
  #sigma.proc = 0.1 IF Fixed: typicallly 0.05-0.15 (see Ono et al. 2012)
  
  #><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>>
  # Optional: Do TAC Projections
  #><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>>
  Projection = FALSE # Switch on by Projection = TRUE 
  
  # Check final year catch
  catch[nrow(catch),]
  
  # Set range for alternative TAC projections
  TACs = c(0,2500,5000,9000,10000,12000,15000) #example
  
  # Set number of projections years
  pyrs = 10
  
  #><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>
  # Execute model and produce output
  #><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>
  
  # set TRUE if PUCL rules should be applied
  pucl = FALSE
  
  # Option to produce standard KOBE plot
  KOBE.plot = TRUE
  Biplot=FALSE
  
  save.all=FALSE 
  
  # MCMC settings
  ni <- 5000 # Number of iterations
  nt <- 1 # Steps saved
  nb <- 1000 # Burn-in
  nc <- 2 # number of chains
  nsaved = (ni-nb)/nt*nc
  
  # Run model (BSPSPexe file must be in the same working directory)
  source(paste0(JABBA.file,"/JABBAv3.r"))
 
  # Save Residuals
  SR = rbind(SR,data.frame(S=s,R=as.numeric(Resids)))
  B_Bmsy.sen = cbind(B_Bmsy.sen,apply(posteriors$BtoBmsy,2,median))
  H_Hmsy.sen = cbind(H_Hmsy.sen,apply(posteriors$HtoHmsy,2,median))
  if(s==1){sRMSE = RMSE} else {sRMSE = c(sRMSE,RMSE)}
  
  }# THE END

dir.create(paste0(assessment,"/Sensetivity"),showWarnings = F)

Legend = c("S3(Base-Case)","comb.SPA","comb.JPN","+ BRA1","+ TAI","- BRA2","- SPA1","- SPA2","- JPN1","- JPN2","- URU","- ZA")

# Prediction Validation
Par = list(mfrow=c(1,2),mai=c(0.2,0.7,0,.15),omi = c(0.3,0.,0.2,0) + 0.1,mgp=c(2,0.5,0), tck = -0.02,cex=0.8)
png(file = paste0(assessment,"/Sensetivity/Sen_",assessment,".png"), width = 9, height = 4, 
    res = 200, units = "in")
par(Par)


plot(years,years,type="n",ylim=c(0,max(H_Hmsy.sen)),ylab="H/Hmsy",xaxt="n")
for(i in 1:(length(Scenarios))){lines(years,H_Hmsy.sen[,i],col=rainbow((length(Scenarios)))[i],lwd=2)}
lines(years,H_Hmsy.sen[,1],col=rainbow((length(Scenarios)))[1],lwd=2)
legend("topleft",paste0(Legend," (",sRMSE,"%)"),col=rainbow((length(Scenarios))),lwd=2,bty="n",cex=0.8)
abline(h=1,lty=2)
axis(1,labels=T,cex=0.9)

plot(years,years,type="n",ylim=c(0,max(B_Bmsy.sen)),ylab="B/Bmsy")
for(i in 1:(length(Scenarios))){lines(years,B_Bmsy.sen[,i],col=rainbow(length(Scenarios))[i],lwd=2)}
lines(years,B_Bmsy.sen[,1],col=rainbow(length(Scenarios))[1],lwd=2)
abline(h=1,lty=2)
mtext(paste("Year"), side=1, outer=T, at=0.5,line=1,cex=0.9)
dev.off()

  