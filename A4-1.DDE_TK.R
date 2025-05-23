library(tidyverse)
library(deSolve) 
library(flextable)
library(ggplot2)
library(dplyr)

# Cleaning dataset for measured serum DDE 
# (not published in Github)
# source("ENDORUUP_clean.R")

# Parameter Setting -------------------------------------------------------
# time parameter (day)
starttime <- 100  # to maintain steady-state
midtime1 <- 3*30 + starttime
midtime2 <- 6*30 + starttime
stoptime <- 12*30 + starttime
timepoint = c(starttime, midtime1, midtime2, stoptime)
dtout <- 1 ### resolution of output time

Timesbefore <- seq(0, starttime, dtout)   
Times1 <- seq(starttime, midtime1, dtout) 
Times2 <- seq(midtime1, midtime2, dtout)
Times3 <- seq(midtime2, stoptime, dtout)
Totaltime = c(Timesbefore[-1],Times1[-1],Times2[-1],Times3)

# Biological parameters ---------------------------------------------------
# 0. Anthropometric values 
SEX = 1                       # Sex (0 = male; 1 = female)
AGE = 41                    	# age (year) 
hct = SEX*0.38 + (1-SEX)*0.45 # female average, hematocrit
BW_n = 61.3                   # normal weight (kg) from control group
BMI_n = 21.9                  # normal BMI (kg/m2) from control group
BH = sqrt(BW_n/BMI_n)*100     # body height (cm) 
BWchange = c(109.8, 98.6, 93.1, 87.8)   #body weight (kg) at each time point
# Body weight to be continuous function
BWloess <- loess(BWchange~timepoint, span=1)
BWtotal <- predict(BWloess, new=data.frame(timepoint=Totaltime))
BWtotal <- ifelse(is.na(BWtotal), BWchange[1], BWtotal)
BWfunc = approxfun(Totaltime, BWtotal , rule = 2)

# PBPK model and parameters
# 1. Body surface [cm2]
BS_n = BW_n^0.5150 * BH^(0.4220)*234.9      # at normal body weight

# 2. Organ Volume [L]
# 2-1. Blood [L]
VB = (13.1*BH + 18.05*BW_n-480)*(1-hct*0.91)/0.5723/1000      #plasma volume in male
# 2-2. Rest of the body when normal body weight [L]
VL = 0.0501*(BW_n^0.78)                                        #Liver
VPPT_skin = (AGE>20)*((1.834 * (BS_n/10^4)) + 7.850E-2*(BS_n/10^4)^1.049) + 
  (AGE<=20)*(-9.356E-5 - 2.151E-5*AGE - 5.058E-1*(BS_n/10^4) + 1.134E-6*AGE^2 + 0.117*AGE*(BS_n/10^4) - 1.673E-5*(BS_n/10^4)^2 ) 		      #Skin
VPPT_sm = (AGE>18)*(6.780*(BS_n/10^4)^1.629 - 1.492E-3*AGE + 3.580) +
  (AGE<=18)*(1.629E-1*BW_n + 2.603E-2*BH + 4.551E-1*AGE - 3.332)      #Skeletal muscles
VPPT_heart=1.017E-7 *(BH^0.6862 * BW_n^0.3561 * 242.7)^1.420 
VPPT_ms = VPPT_sm+VPPT_skin          		#muscles and skin
VPPT=VPPT_skin +VPPT_ms+VPPT_heart
VRPT = (AGE>18)*((2.331E-3*AGE + 0.1253*BW_n^0.8477 + BH^0.3821 - 4.725) - VL) +
  (AGE<=18)*((2.515E-2*AGE + 7.619*(BW_n^2/BH)^0.8477 - 6.098) - VL) #richly perfused tissue

# 2-3. Adipose tissue [L]
VATSATR = 0.39*SEX + 0.84*(1-SEX) # VAT/SAT volume ratio (Kaess et al., 2013)
# volume of subcutaneous AT (SAT) is fixed throughout the weight change
VAT_n = 0.91*BW_n - VB - VL - VPPT - VRPT
VAT = 0.91*BWtotal - VB - VL - VPPT - VRPT
VSAT = VAT_n*(1/(1+VATSATR))
VVATtotal = VAT - VSAT
VVATfunc = approxfun(Totaltime, VVATtotal , rule = 2) 
#Slope of visceral adipose tissue changes
DVAT <- c(diff(VVATtotal) / dtout, NA)
DVAT[length(DVAT)] = DVAT[length(DVAT)-1]
DVVATfunc = approxfun(Totaltime, DVAT , rule = 2)

#Fraction of tissue volumes that is blood!
FVV_AT = 0.02 # Vascular volume fraction of adipose tissue (%)
FVV_PPT = 0.01 # Vascular volume fraction of poorly perfused tissue (%)
# fat compartment blood volume
VSATV = VSAT*FVV_AT
VSATT = VSAT-VSATV
VPPTV = VPPT*FVV_PPT
VPPTT = VPPT-VPPTV
# volume of vascular VAT is fixed throughout the weight change
VVATV = (VAT_n - VSAT)*FVV_AT

# Chemical specific parameters -------------------------------------
#Fraction
FLAT=74.1   #lipid contents of adipose tissue
FWAT=21.2 #water contents of adipose tissue
FLB=0.6       #lipid contents of blood
FWB=79.0 #water contents of blood
FLL=4.6        #lipid contents of liver
FWL=74.5    #water contents of liver
FLRPT=3.68
FWRPT=78.1

FLPPT_ms= 4.2   #lipid contents of muscles and skin
FWPPT_ms=74.1 #water contents of muscles and skin
FLPPT_heart=1.0 #lipid contents of heart
FWPPT_heart=73.0 #water contents of heart
Kow=10^6.51 #Kow of DDE

FLPPT=(FLPPT_ms*VPPT_ms+FLPPT_heart*VPPT_heart)/(VPPT_ms+VPPT_heart) #lipid contents of PPT
FWPPT=(FWPPT_ms*VPPT_ms+FWPPT_heart*VPPT_heart)/(VPPT_ms+VPPT_heart)         	#water contents of PPT

#partition coefficient to blood
PAT=(Kow*FLAT+FWAT)/(Kow*FLB+FWB)  #partition coefficient of adipose tissue:blood
PL=(Kow*FLL+FWL)/(Kow*FLB+FWB)    #partition coefficient of liver:blood
PPPT=(Kow*FLPPT+FWPPT)/(Kow*FLB+FWB)  #partition coefficient of poorly perfused tissue:blood
PRPT=(Kow*FLRPT+FWRPT)/(Kow*FLB+FWB)  #partition coefficient of richly perfused tissue :blood

#Permeability constants (L/d/kg tissue) (Permeation area cross products)
PAATC = 0.297*24 # Adipose tissue tissue permeability constant 
PASAT = PAATC * VSAT
PAPPTC = 0.944*24 # Poorly perfused tissue permeability constant
PAPPT = PAPPTC * VPPT # Fat:blood permeability surface area coefficients

CLint=0.0156840*24   #Liver instrinsic clearance (L/d/kg-liver)
Eh=(CLint*VL)/(CLint*VL+(1.00*60*24*VL))        # Elimination rate (unitless)

QL = 24*1.00*60*VL 

# 3. PBPK model ----------------------
pbpkmodel <- function(Time, State, Para) {
  with(as.list(c(State, Para)), {
    VVAT = VVATfunc(Time)
    VVATT = VVAT - VVATV
    DVVAT = DVVATfunc(Time)
    ###################################################################################################################
    #Blood flow (L/d)
    QC = 24*15.048*BWfunc(Time)^0.7609    #cardiac output
    QVAT = 24*0.03*60*VVAT
    QSAT = 24*2.6*0.9/100*60*VSAT         #Virtanen et al., 2002
    QL = 24*1.00*60*VL 
    QPPT = 24*(1.80* (VPPT_sm) + 57.60 * VPPT_heart +  9.0*(VPPT_skin))
    QRPT = QC - (QPPT+QVAT+QSAT+QL)
    ###################################################################################################################
    ## Concentration of the chemical in vein compartment
    CVL = AL/(VL * PL) # con' of chem in liver  / PC of plasma: liver
    CVVATV = AVATV/VVATV  # Blood (Vein) of Visceral Adipose Tissue
    CVVATT = AVATT/VVATT  # Tissue of Visceral Adipose Tissue
    CVSATV = ASATV/VSATV  # Blood (Vein) of Subcutaneous Adipose Tissue
    CVSATT = ASATT/VSATT  # Tissue of Subcutaneous Adipose Tissue
    CVRPT = ARPT/(VRPT * PRPT) # con' of chem in RPT    / PC of plasma: richly perfused tissue
    CVPPTT = APPTT/VPPTT
    CVPPTV = APPTV/VPPTV
    ###################################################################################################################      
    ## Blood compartment
    #con' of chemical in the vein
    CV = (((QL+QVAT) * CVL + QSAT * CVSATV + QRPT * CVRPT + QPPT * CVPPTV ) / QC)
    CA = AA/VB # con' of chem in artery = amount of chem in artery / volume of blood
    RA = QC * (CV - CA)  # rate of change in amount of chem in tissue of blood compartment
    dAA = RA
    ###################################################################################################################    
    ## Liver
    RL = QL*CA + QVAT*CVVATV - (QL+QVAT)*CVL - Eh*(QL*CA + QVAT*CVVATV) # rate of change in amount of the chem in liver
    dAL = RL + Intake*BWfunc(Time)
    CL=AL/VL
    ###################################################################################################################
    ###################################################################################################################
    ## Adipose tissue compartment
    # Visceral Adipose tissue
    PAVAT = PAATC * VVAT # Fat:blood permeability surface area coefficients (L/h)
    RVATT = PAVAT * CVVATV - (PAVAT*CVVATT)/PAT + (CVVATT*DVVAT)
    dAVATT = RVATT 
    RVATV = QVAT * (CA - CVVATV) - PAVAT * CVVATV + (PAVAT*CVVATT)/PAT - (CVVATT*DVVAT)
    dAVATV = RVATV
    AVATtotal = AVATT + AVATV
    CVAT = AVATtotal/VVAT
    
    # Subcutaneous Adipose tissue
    RSATT = PASAT * CVSATV - (PASAT*CVSATT)/PAT
    dASATT = RSATT 
    RSATV = QSAT *(CA - CVSATV) - PASAT*CVSATV + (PASAT*CVSATT)/PAT
    dASATV = RSATV
    ASATtotal = ASATT + ASATV
    CSAT = ASATtotal/VSAT
    
    CAT = (AVATtotal+ASATtotal)/(VVAT+VSAT)
    ###################################################################################################################
    ## RPT of body compartment
    RRPT = QRPT * (CA - CVRPT)
    dARPT = RRPT
    CRPT = ARPT/VRPT
    ###################################################################################################################
    ## PPT of body compartment
    ### later consider that PPT should be two compartments and iclude RPPT in paras
    ## PPT of body compartment
    RPPTT = PAPPT*CVPPTV - (PAPPT*CVPPTT)/PPPT
    dAPPTT = RPPTT
    RPPTV = QPPT * (CA - CVPPTV) - PAPPT*CVPPTV + (PAPPT*CVPPTT)/PPPT
    dAPPTV = RPPTV
    APPTtotal = APPTT + APPTV
    CPPT =APPTtotal/VPPT
    
    list(c(dAA, dAL, dAVATV, dAVATT, dASATV, dASATT, dARPT, dAPPTV, dAPPTT),
         c(AT=CAT,Plasma=CA,Liver=CL,RPT=CRPT,PPT=CPPT,SAT=CSAT,VAT=CVAT))
  })
}

# Back calculated intake rate
# Obtained from backcalcualted DDE intake rates
# load("DDE_intake_ENDORUP.RData")

# Or Fixed intake rate from literature (Gunderson 1995)
Intake = (0.013*SEX + 0.0181*(1-SEX)) # ug/kg/d 

bloodtotal <- as.data.frame(matrix(nrow = length(Totaltime), 
                                   ncol = nrow(DDE_all)))

for (i in 1:nrow(DDE_all))  {
  plas_DDE_initial = DDE_all[i, 2]
  Intake = intake_rate[i]
  
  # initial unit: [ug]
  statebefore=c(AA = plas_DDE_initial*VB, 
                AL = (plas_DDE_initial*(1-Eh)+Intake*BWfunc(0)/QL)*PL*VL,
                AVATV = plas_DDE_initial*VVATV,
                AVATT = plas_DDE_initial*PAT*(VVATfunc(0) - VVATV),
                ASATV = plas_DDE_initial*VSATV,
                ASATT = plas_DDE_initial*PAT*VSATT,
                ARPT = plas_DDE_initial*PRPT*VRPT,
                APPTV = plas_DDE_initial*VPPTV,
                APPTT = plas_DDE_initial*PPPT*VPPTT) 
  
  # PBPK output-----------
  outall <- ode(y = statebefore,
                times = Totaltime, 
                func = pbpkmodel, 
                parms = 0) 
  predstotal <- as.data.frame(outall)
  
  bloodtotal[, i] <- predstotal$Plasma
  print(paste("Participant", i, "done"))
}

colnames(bloodtotal) <- DDE_all$ID
