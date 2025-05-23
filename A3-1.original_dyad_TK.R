# Appendix 3-1
# Original TK model from Verner et al. (2016)
# Translated by Jung Y. (University of California, Irvine)

library(deSolve) 
library(dplyr)
library(tidyverse)

Age_conception = 30                     # Maternal Age at conception (years)
Total_bf = 4.1/12                       # Total duration of breastfeeding (years)
Total_pp = 12                           # Total postpartum duration for model predictions (years)
Age_delivery = Age_conception+0.769     # Maternal Age at delivery (years)
Weight_prepregnancy_M = 70              # Pre-pregnancy maternal body weight (kg)
Sex = 0                                 # Sex of child (Female=0; Male=1)

Age_span = seq(from=0,to=100,by=0.1)
Age_gestational = c(0.0, 0.5, 1.0, 2.0, 5.0, 10.0, 13.0, 16.0, 20.0, 100.0)+0.769

# Bodyweight trajectory
table_bw_females = list("years"=c(0.0, 0.217, 0.353, 0.463, 0.551, 0.639, Age_gestational), #gestational age
                        "y"= c(0.0, 0.0, 0.26, 0.69, 1.25, 2.02, 3.4, 7.2, 9.5, 12.1, 18.0, 33.0, 46.0, 54.0, 58.0, 58.0))
table_bw_females = data.frame(approx(table_bw_females$years,table_bw_females$y, xout=round(Age_span,2), rule=2))
BW_F_func <- approxfun(table_bw_females$x, table_bw_females$y, rule = 2)

table_bw_males = list("years"=c(0.0, 0.217, 0.353, 0.463, 0.551, 0.639, Age_gestational),
                      "y"=c(0.0, 0.0, 0.26, 0.69, 1.25, 2.02, 3.5, 7.8, 10.3, 12.7, 18.0, 32.0, 46.0, 61.0, 71.0, 71.0))
table_bw_males = data.frame(approx(table_bw_males$years,table_bw_males$y, xout=round(Age_span,1), rule=2))
BW_M_func <- approxfun(table_bw_males$x, table_bw_males$y, rule = 2) 

# Adult fraction of time spent indoors [unitless] (US EPA, 2011)
q_dust_m = 1/10^3           # Dust ingestion rate [g/day]; indoor settled dust only
f_ind_m = 1159/60/24        # Mean daily time spent indoors at the age of 18-64 y.o.
F_uptake = 0.8              # The gastrointestinal uptake factor (Gebbink et al. 2015; Y. Wu et al., 2020 )

# Children Dust ingestion rate by age [g/day]
# indoor settled dust only
table_q_dust_child <- list("years"=c(0.0, 0.5, 1.0, 6.0, 20.0),
                           "y"= c(0, 30, 60, 60, 30)/10^3)
table_q_dust_child <- data.frame(approx(table_q_dust_child$years,table_q_dust_child$y, 
                                        xout=round(seq(from=0,to=20,by=0.1),1), rule=2))
q_dust_child_func <- approxfun(table_q_dust_child$x, table_q_dust_child$y, rule = 2)

# Children fraction of time spent indoors [unitless] (US EPA, 2011)
table_f_ind_child <- list("years"=round(c(0.0, 1/12, 3/12, 6/12, 1, 2, 3, 6, 11, 16),1),
                          "y"= c(1440, 1432, 1414, 1301, 1353, 1316, 1278, 1244, 1260, 1248)/60/24)
table_f_ind_child <- data.frame(approx(table_f_ind_child$years,table_f_ind_child$y, 
                                       xout=round(seq(from=0,to=20,by=0.1),1), rule=2))
f_ind_child_func <- approxfun(table_f_ind_child$x, table_f_ind_child$y, rule = 2)

# Simulation time
# Age of mother (years; by=1 hour)
Age_M =(Age_conception+seq(from=0, to=0.769+Total_pp, by=1/365))    

preds_orig <- data.frame(Age_M-Age_conception)
for (k in 1:3) {
  Compound = k                            # Compound (1:PFOS, 2:PFOA, 3:PFHxS)
  
  #******************************************************************************************************
  # Compound
  PFOS = 1*(Compound == 1)
  PFOA = 1*(Compound == 2)
  PFHxS = 1*(Compound == 3)
  
  #Dust concentration (ng/g)
  PFOS_dust = 59
  PFOA_dust = 37.1
  PFHxS_dust = 5.55
  # Initial PFAS level in indoor dust (ug/g)
  dust_conc = (PFOS_dust*PFOS + PFOA_dust*PFOA + PFHxS_dust*PFHxS)/1000
  
  # Daily intake through diet (ng/kg/d)
  PFOS_diet = 0.58
  PFOA_diet = 0.18
  PFHxS_diet = 0.08
  # Daily intake through diet (ug/kg/d)
  diet_edi = (PFOS_diet*PFOS + PFOA_diet*PFOA + PFHxS_diet*PFHxS)/1000
  
  # Volume of distribution (l/kg)
  PFOS_VD_BW = 0.230    # PFOS
  PFOA_VD_BW = 0.170    # PFOA
  PFHxS_VD_BW = 0.213   # PFHxS
  VD_BW = PFOS*PFOS_VD_BW+PFOA*PFOA_VD_BW+PFHxS*PFHxS_VD_BW
  
  # Plasma:milk level ratio
  PFOS_PMilk = 0.0138   # PFOS
  PFOA_PMilk = 0.0577   # PFOA
  PFHxS_PMilk = 0.0140  # PFHxS
  PMilk = PFOS*PFOS_PMilk+PFOA*PFOA_PMilk+PFHxS*PFHxS_PMilk
  
  # Cord:maternal plasma level ratio
  PFOS_P_CM = 0.454     # PFOS
  PFOA_P_CM = 0.783     # PFOA
  PFHxS_P_CM = 0.556    # PFHxS
  P_CM = PFOS*PFOS_P_CM+PFOA*PFOA_P_CM+PFHxS*PFHxS_P_CM
  
  # Half-lives (years)
  PFOS_HL = 5.4         # PFOS
  PFOA_HL = 3.8         # PFOA
  PFHxS_HL = 8.5        # PFHxS
  # Half-lives (day)
  HALF_LIFE = (PFOS*PFOS_HL+PFOA*PFOA_HL+PFHxS*PFHxS_HL)*365 
  
  # placental transfer
  ktrans1 = P_CM # Mother->fetus placental transfer  (l/d)
  ktrans2 = 1.0  # Fetus->mother placental transfer  (l/d)
  #*****************************************************************
  
  pbpkmodel <- function(Time, State, Parmeters) {
    with(as.list(c(State, Paras)), {
      #Switches for Time Events    
      #Conception
      if (Age_conception-Time/365>0) {            
        IO_CONCEPTION = 0
      } else {
        IO_CONCEPTION = 1
      }
      #Delivery
      if (Age_delivery-Time/365>0) {              
        IO_DELIVERY = 0
      } else {
        IO_DELIVERY = 1
      }
      #Breasfeeding
      if (Total_bf-(Time/365-Age_delivery)>0 & (Time/365-Age_delivery)>0 ) {  
        IO_BF = 1
      } else {
        IO_BF = 0
      }
      
      # Age of child (year)
      Age_C = IO_DELIVERY*(Time/365-Age_delivery)
      # Gestational age of child (year)
      Age_C_gest =Time/365-Age_conception
      
      # Breastfeeding end
      if (Total_bf-Age_C>0){                   
        IO_END_TOTAL_BF = 0
      } else {
        IO_END_TOTAL_BF = 1
      }
      # Child age > 6 m.o
      if (0.500-Age_C>0){                      .
        IO_CHILD_6m = 0 
      } else {
        IO_CHILD_6m = 1
      }
      # Child age > 12 m.o
      if (1.000-Age_C>0){
        IO_CHILD_12m = 0
      } else {
        IO_CHILD_12m = 1
      }
      # Child age > 30 m.o
      if (2.500-Age_C>0){ 
        IO_CHILD_30m = 0
      } else {
        IO_CHILD_30m = 1
      }
      #**********************************************
      # weight (kg)
      BW_C = (IO_CONCEPTION)*((1-Sex)*BW_F_func(Age_C_gest)+Sex*BW_M_func(Age_C_gest))
      BW_M = (Weight_prepregnancy_M/BW_F_func(Age_conception))*BW_F_func(Time/365)
      #***********************************************
      # Volume of distribution (l)
      VD_M = VD_BW*BW_M   
      VD_C = VD_BW*BW_C+1E-03 
      
      # Concentrations (ug/L)
      C_MOTHER = A_MOTHER/VD_M
      C_CHILD = A_CHILD/VD_C
      
      # Placental diffusion
      # Mother->fetus placental transfer rate (ug/day)
      DIFF_M_C = (IO_CONCEPTION-IO_DELIVERY)*ktrans1*C_MOTHER 
      # Fetus->mother placental transfer rate (ug/day)
      DIFF_C_M = (IO_CONCEPTION-IO_DELIVERY)*ktrans2*C_CHILD 
      #******************************************************************************************************
      # Breast milk transfer (g/d)
      VOLUME_MILK_1ST_YEAR = BW_C*(-0.312*(Age_C*365)+157.7)
      VOLUME_MILK_2ND_YEAR = (-0.763*Age_C*365 + 720.3)
      # Breast milk consumption rate (l/d)
      VOLUME_MILK = IO_BF*(IO_DELIVERY-IO_END_TOTAL_BF)* 
        ((1-IO_CHILD_12m)*VOLUME_MILK_1ST_YEAR 
         +(IO_CHILD_12m-IO_CHILD_30m)*VOLUME_MILK_2ND_YEAR)/1000 
      
      # Concentration in breast milk (ug/l)
      C_MILK = PMilk*C_MOTHER           
      # Mother-child lactational transfer (ug/d)
      TRANSFER_BF = VOLUME_MILK*C_MILK  
      #******************************************************************************************************
      # Elimination
      # Elimination from the maternal compartment (ug/day)
      ELIMINATION_M = C_MOTHER*VD_M*log(2)/HALF_LIFE 
      # Elimination from the child compartment (ug/day)
      ELIMINATION_C = C_CHILD*VD_C*log(2)/HALF_LIFE
      #******************************************************************************************************
      # Dosing
      INTAKE_M =
        # Maternal PFAS dietary intake (ug/day)
        diet_edi*BW_M +                              
        # Maternal PFAS dust intake (ug/day)
        dust_conc*q_dust_m*f_ind_m*F_uptake                        
      INTAKE_C = 
        IO_CHILD_6m*(
            # Child PFAS dietary intake (ug/day)
            diet_edi*BW_C + 
            # Child PFAS dust intake (ug/day)
            dust_conc*q_dust_child_func(Age_C)*f_ind_child_func(Age_C)*F_uptake)         
      
      #******************************************************************************************************
      # Mass-balance differential equations
      # Maternal compartment
      # Rate of amount in maternal compartment (ug/day)
      RA_MOTHER = INTAKE_M-ELIMINATION_M-TRANSFER_BF-DIFF_M_C+DIFF_C_M 
      # Child compartment
      # Rate of amount in child compartment (ug/day)
      RA_CHILD = INTAKE_C+TRANSFER_BF-ELIMINATION_C+DIFF_M_C-DIFF_C_M 
      
      ## Mass balance
      list(c(RA_MOTHER,RA_CHILD), 
           c(serum_M=C_MOTHER, serum_C=C_CHILD,Age_C=Age_C) )
    })
  }
  
  # For initial conditions, assume mom is at steady state at time of conception
  # (ug/day/kg)/(/day)/(l/kg) = ug/L
  C_SS <- (diet_edi+(dust_conc*q_dust_m*f_ind_m*F_uptake)/Weight_prepregnancy_M) / 
    (log(2)/HALF_LIFE) / VD_BW
  State <- c(A_MOTHER=C_SS*VD_BW*Weight_prepregnancy_M, A_CHILD=0) #[ug]
  Paras <- c(ktrans1, ktrans2, HALF_LIFE)
  
  # PBPK output
  # input unit: initial state (ng), time (hour)
  # output unit: serum [ug/L; ng/mL], age [year]
  out <- ode(y = State,
             times = Age_M*365, 
             func = pbpkmodel, 
             parms = Paras) 
  
  # extract results
  out1=as.data.frame(out)[,5]
  preds_orig = cbind(preds_orig, out1)
}
colnames(preds_orig) <- c("time", "PFOS","PFOA","PFHxS","PFNA")

preds_orig <- preds_orig %>% 
  filter(time > 0.769) %>% 
  mutate(time = (time-0.769))
