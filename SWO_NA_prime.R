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
KOBE.plot =Biplot= CPUE.plot= TRUE
################################################

# Set Working directory file where to store the results
File = "C:/Assessments/ICCAT"
# Set working directory for JABBA R source code
JABBA.file = "C:/Assessments/JABBA"
# Set Assessment
assessment = "SWO_NA"


# Choose Sceanrio name for creating a seperate folCder

Scenarios = c(paste0("s",1:100))
# Base-Case S1: Fitted to SPA, JPN and ZA LL CPUE assuming Fox 
# and B1950/K = 0.95, CV = 0.05  
# S2: as S1 but introducing time blocks in 1999 for SPA and 2006 for JPN
# S3: as Sw but with lnorm prior B1950/K = 1, CV = 0.05

for(s in 1:2){
#s=1
Scenario = Scenarios[s] 

# Choose model type: 
# 1: Schaefer
# 2: Schaefer with Depensation (CMSY: Froese et al. 2016)  
# 3: Fox
# 4: Fox with Depensation (CMSY: Froese et al. 2016)
# 5: Pella-Tomlinsson  
 Model = ifelse(s==1,1,3) # 1 #ifelse(s>2,5,1)  # model
   
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
  # Load assessment data
  cpue = read.csv(paste0(assessment,"/cpue",assessment,"_Combined.csv"))#
  se =  read.csv(paste0(assessment,"/se",assessment,"_Combined.csv"))# use 0.001 if not available 
  catch = read.csv(paste0(assessment,"/catch",assessment,".csv"))
  
  names(cpue)
  names(catch)
  
  #--------------------------------------------------
  # option to exclude CPUE time series or catch year
  #--------------------------------------------------
 # manipulation of the data input files
  #Drop Canada early CPUE and SE data, which is the 2nd column 
  #if(s==2){
  #  cpue=cpue[,-2]
  #  se=se[,-2]
  #}
  
    cpue=cpue[,c(1,13)]
    se=se[,c(1,13)]
  
  
  #check it was dropped
  names(cpue)
  #------------------------------------------------------
  # Option use mean CPUE from state-space cpue averaging
  #-----------------------------------------------------
  # Produce CPUE plot average plot
  # If there is a single time series, make this FALSE
  CPUE.plot = FALSE#TRUE
  
  #if(s>4){meanCPUE = TRUE} else{meanCPUE = FALSE}
  
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
  mu.psi = 0.85 
  CV.psi = 0.1 
  
  #--------------------------------------------------------------
  # Determine estimation for catchability q and observation error 
  #--------------------------------------------------------------
  # Assign q to CPUE
  sets.q = 1:(ncol(cpue)-1) 
  #if(s==1)sets.q = c(1,2,2,3,3,4,5)
  
  # Series
  #sets.var = rep(1,ncol(cpue)-1)# estimate the same additional variance error
  sets.var = 1:(ncol(cpue)-1) # estimate individual additional variace
  #if(s==1) sets.var= c(1,2,2,2,3,4,5)
  
  # As option for data-weighing
  # minimum additional observation error for each variance set (optional choose 1 value for both)
  min.obsE = c(0.2) # Important if SE.I is not availble
  
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
  r.prior = c(0.424,0.4) # as range with upper and lower bound of lnorm prior
  
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
  Projection = TRUE # Switch on by Projection = TRUE 
  
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
  Biplot=TRUE
  
  save.all=FALSE
  save.kobe.posterior = FALSE 
  
  # MCMC settings
  ni <- 50000 # Number of iterations
  nt <- 10 # Steps saved
  nb <- 10000 # Burn-in
  nc <- 2 # number of chains
  nsaved = (ni-nb)/nt*nc
  
  # Run model (BSPSPexe file must be in the same working directory)
  source(paste0(JABBA.file,"/JABBAv3.r"))
 
  }# THE END
  



 