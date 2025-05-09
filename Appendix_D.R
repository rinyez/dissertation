#---------------------------------------
# Load necessary library
library(deSolve)
library(lubridate)
library(dplyr)
library(tidyverse)
# ---------------------------
# Require

# 1. Greene et al., (2023)
# Parameters
# Define intake rates and body weights at specific ages
age_days <- 
  c(0, 30, 90, 180, 365, c(2, 3, 6, 11, 16, 21, 30, 40, 50, 60)*365)
water_intake <- 
  c(240, 290, 186, 151, 119, 67, 45, 41, 31, 31, 47, 44, 43, 42, 42)  # mL/kg/day
breastmilk_intake <- 
  c(220, 190, 150, 130, 130,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0)  # mL/kg/day
body_weight <- 
  c(3.6, 3.8, 7, 8.9, 10.5, 
    13.4, 18.6, 30.7, 56.8, 71.4, 72.5, 74.5, 78.5, 80.7, 80.7)  # kg

# Linear interpolation of intake rates and body weight
# Water intake rate: input [year]
water_intake_rate_func <- 
  approxfun(age_days/365, water_intake, rule=2)
# Breastmilk intake rate: input [day]
breastmilk_intake_rate <- 
  approx(age_days, breastmilk_intake, xout = seq(0, 80*365, by = 1))$y
breastmilk_intake_rate[1] <- 0*breastmilk_intake_rate[1]
breastmilk_intake_rate[2] <- 0.13*breastmilk_intake_rate[2]
breastmilk_intake_rate[3] <- 0.26*breastmilk_intake_rate[3]
breastmilk_intake_rate[4] <- 0.76*breastmilk_intake_rate[4]
# Body weight of a child: input [year]
BW_child_func <- approxfun(age_days/365, body_weight, rule=2)

Vd <- 0.36                       # Volume of distribution in L/kg
half_life <- 902                 # Elimination half-life in days
placental_transfer = 0.83        # Placental transfer factor
breastmilk_transfer = 0.068      # Breastmilk transfer factor

# A function to calculate serum concentration
EC_func <- 
  function(previous_mass, intake_today, Vd, BW, half_life) {
    k <- log(2) / half_life  # Elimination rate constant
    serum_conc <- ((previous_mass + intake_today) / (Vd * BW)) * exp(-k)
  return(serum_conc)
}

# A main function to simulate serum concentrations over time
EC_mother_func <- 
  function(intake_data, Vd, half_life,
           BW_mom, initial_maternal_serum, 
           breastmilk_intake_rate, BW_child) {
    
    n <- length(intake_data)
    serum_conc <- numeric(n)
    
    # Initial mass in the system
    serum_conc[1] <- initial_maternal_serum
    previous_mass <- 
      (initial_maternal_serum * Vd * BW_mom[1]) - initial_mass_infant # ug
    
    for (i in 2:n) {
      BW_today <- BW_mom[i]
      if (i <= 365) { 
        intake_today <- 
          intake_data[i] - (serum_conc[i-1]) * (breastmilk_transfer) * (breastmilk_intake_rate[i-1]) / 1000 * BW_child[i-1]
      } else { 
        intake_today <- intake_data[i]
      }
      serum_conc[i] <- 
        EC_func(previous_mass, intake_today, Vd, BW_today, half_life)
      # Update previous mass for the next day's calculation
      previous_mass <- serum_conc[i] * Vd * BW_today
    }
    
    return(serum_conc)
  }

# Function to estimate breastfed child's serum concentrations
EC_BF_child_func <- 
  function(intake_data, maternal_serum_conc, 
           breastmilk_intake_rate, 
           body_weights, Vd, half_life) {
    n <- length(maternal_serum_conc)
    infant_serum_conc <- numeric(n)
    
    # Initial mass in the infant at birth
    previous_mass_infant <- initial_mass_infant
    
    for (i in 1:n) {
      # Daily intake
      if (i <= 365) { 
        # From breastmilk [ug]
        intake_today <- maternal_serum_conc[i] * breastmilk_transfer * breastmilk_intake_rate[i] / 1000 *body_weights[i]
      } else { 
        # From drinking water after 1 year [ug]
        intake_today <- intake_data[i]
      }
      infant_serum_conc[i] <- EC_func(previous_mass_infant, intake_today, 
                                                                      Vd, body_weights[i], half_life)
      # Update previous mass for the next day's calculation
      previous_mass_infant <- infant_serum_conc[i] * Vd * body_weights[i]
    }
    
    return(infant_serum_conc)
  }

# Function to convert date format (m/d/Y) into years in decimals 
# (e.g., year 2010.5)
# This can be skipped if the date is provided in demicimal-year format
decimal_year_func <- function(date) {
  year <- year(date)
  
  days_in_year <- ifelse(leap_year(date), 366, 365)
  decimal <- (yday(date) - 1) / days_in_year
  decimal_year <- year + decimal
  return(decimal_year)
}

EC_Mother_Greene = data.frame()
EC_Child_PK_Greene = data.frame()
for (q in 1:dim(EC_Mother)[1]) {
  ptnumb=EC_Mother[q,1]
  
  # A demographic info of pt
  demogpt = filter(demog, ID==ptnumb)
  demogpt["DOB"]=as.Date(paste0(demogpt[1,4],"/6/15"))         #Assumed DOB as July 1
  
  # Initial serum concentration of mother from the existing estimation
  EC_initial = EC_Mother %>% filter(ID==ptnumb)
  
  # 3. Check how many parities the pt has
  childinfopt = as.data.frame(filter(child, Mother_ID==ptnumb))    
  # Number of pregnancy in one mother
  maxpreg = nrow(childinfopt)                               
  
  for (r in 1: maxpreg) {
    if (r == maxpreg) {
      # Set simulation time from the birth of a first child till visitdate
      simulday = seq(as.Date(childinfopt[r,12]),as.Date(demogpt$visitdate,format ="%m/%d/%Y"), by='days')
    } else {
      simulday = seq(as.Date(childinfopt[r,12]),as.Date(childinfopt[r+1,12])-1, by='days')
    } 
    # Number of years after 1/1/1950
    simulday_FR_dec = decimal_year_func(simulday)-1950
    # Mom age in decimal year
    simulday_age_mom_dec = decimal_year_func(simulday)-decimal_year_func(demogpt$DOB)
    # Child age in decimal year
    simulday_age_child_dec = decimal_year_func(simulday) - decimal_year_func(as.Date(childinfopt[r,12]))
    
    # ---------------------------
    # Body weight [kg]
    BWpt = filter(BW, ID==ptnumb)
    BW_mom_func <- approxfun((1:58), BWpt[2:59],rule = 2)
    BW_mom <- BW_mom_func(simulday_FR_dec)
    BW_child <- BW_child_func(simulday_age_child_dec)
    # Water concentration [ug/L]
    watconcpt <- filter(watconc, ID==ptnumb)
    water_conc_func <- approxfun((1:58),watconcpt[2:59], rule = 2)
    water_conc <- water_conc_func(simulday_FR_dec)
    # Water intake rate [mL/kg-bw/day]
    water_intake_rate_mom <- water_intake_rate_func(simulday_age_mom_dec)
    water_intake_rate_child <- water_intake_rate_func(simulday_age_child_dec)
    # daily intake [ug/day]
    water_edi_mom <- water_conc*water_intake_rate_mom*BW_mom/1000
    water_edi_child <- water_conc*water_intake_rate_child*BW_child/1000
    
    # Initial concentration of maternal serum [ug/L]
    if (r == 1) {
      initial_maternal_serum <- water_edi_mom[1] / (Vd * BW_mom[1] * log(2) / half_life)
    } else {
      initial_maternal_serum <- maternal_serum_conc[length(maternal_serum_conc)]
    }
    
    # Initial mass in the infant at birth [ug]
    initial_mass_infant <- initial_maternal_serum * placental_transfer * Vd * BW_child[1]
    
    # Simulate maternal serum concentrations
    maternal_serum_conc <- 
      EC_mother_func(water_edi_mom, Vd, half_life, BW_mom, 
                     initial_maternal_serum, breastmilk_intake_rate,BW_child)
    
    # Estimate breastfed child serum concentrations
    infant_bf_serum_conc <- 
      EC_BF_child_func(water_edi_child, maternal_serum_conc, 
                       breastmilk_intake_rate, BW_child, Vd, half_life)
    infant_serum_visitdate = infant_bf_serum_conc[length(infant_bf_serum_conc)]
    EC_Child_PK_Greene = rbind(EC_Child_PK_Greene, c(as.numeric(childinfopt[r,c(1,8)]),infant_serum_visitdate))
  }  
  
  # Select serum concentrations at serum collection (visit) year
  maternal_serum_visitdate = maternal_serum_conc[length(maternal_serum_conc)]
  EC_Mother_Greene = rbind(EC_Mother_Greene, c(as.numeric(demogpt[,c(1,3)]), maternal_serum_visitdate))
  
  print( paste("This is participant", q) )
}

colnames(EC_Mother_Greene) <- c("ID", "measure", "estimate")
colnames(EC_Child_PK_Greene) <- c("ID", "measure", "estimate")