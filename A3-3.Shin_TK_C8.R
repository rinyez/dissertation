# Appendix C.
# A TK model from Shin et al. (2011)
# Recoded in R by Jung Y. (University of California, Irvine)

# Maternal serum concentrations were estimated from published R codes
# Please direct to Zhu et al. (2022)
# https://doi.org/10.1016/j.envres.2022.112892

# This code is for esimating child serum concentations
# as described in Shin et al. (2011)

# Need a dataframe with 
# 1) Child ID linked with 2) their mother ID
# 3) Year of birth of a child and
# 4) Sex of a child
EC_Child_fix = 
  child %>% filter(Child_ID %in% EC_Child_PK_all$ID) %>% 
  select(Child_ID,Mother_ID,Child_DOB_year,Child_gender) %>% 
  mutate(Child_DOB_year = round(Child_DOB_year,0)) %>% 
  mutate(Child_gender = Child_gender-1) %>%   # Sex of child (Female=0; Male=1)
  mutate(ID = Mother_ID) %>%
  left_join(select(demog, c("ID","visitdate")), by="ID")

# One-compartment TK parameters
VD_BW_fix = 0.198                   # Volume of distribution (l/kg)
HALF_LIFE_fix = 3.5                 # Half-lives (year)
KElim_fix = log(2)/(HALF_LIFE_fix)  #Elimination rate constant (/y)

W1 <- array(0, dim = c(4178, 58, 58), 
            dimnames = list(EC_Child_fix$Child_ID, paste0("Year", 1951:2008), paste0("Year", 1951:2008)))
ChildAge = array(0, dim = c(4178, 58), 
                 dimnames = list(EC_Child_fix$Child_ID, paste0("Year", 1951:2008)))
ChildWeight = array(0, dim = c(4178, 58), 
                    dimnames = list(EC_Child_fix$Child_ID, paste0("Year", 1951:2008)))

for (i in 1:nrow(W1)) {
  sex=as.numeric(EC_Child_fix[i,4])
  
  # r: year for exposure 
  for (r in 1:58){ 
    if(r>=(EC_Child_fix[i,3]-1950)) {
      ChildAge[i,r] <- (r+1950-as.numeric(EC_Child_fix[i,3]))
      
      # j: year for biomarker prediction
      for (j in (max(EC_Child_fix[i,3], 1951) - 1950):r) { 
        
        # Child body weight using body weight trajectory function
        ChildWeight[i,j] <- 
          ((1-sex)*BW_F_func(ChildAge[i,j]+0.769)+sex*BW_M_func(ChildAge[i,j]+0.769))
        # Weight term for each year and each child
        W1[i, j, r] <- 
          ((1 - exp(-KElim_fix)) /(KElim_fix * VD_BW_fix * ChildWeight[i,j]))* exp(-KElim_fix * (r-j))
      }      
    }
  }
}

# Water EDI of each child
wat_edi_chi = EC_Child_fix[,c(1,5)] %>% 
  left_join(wat_edi, by="ID") %>% 
  ungroup() %>% 
  select(-1,-2)

# Air EDI of each child
air_edi_chi = EC_Child_fix[,c(1,5)] %>%
  left_join(air_edi, by="ID") %>% 
  ungroup() %>% 
  select(-1,-2) 

# Total dose [ug/year]
# Water ingestion and air inhalation for > 1 y.o.
ChildDose = (wat_edi_chi/q_wat_m_wo*q_wat_child_func(ChildAge) + 
               air_edi_chi/q_air_m*q_air_child_func(ChildAge))*
  (ChildAge>1) 

EC_Child_fix_1 = EC_Child_fix$Child_ID
for (i in 1:58) {
  EC_Child_fix_1 <- cbind(EC_Child_fix_1, apply(W1[, 1:58, i] * ChildDose, 1, sum))
}

# Add serum concentration at 1 y.o.
# Applying empirical ratio of 1 y.o. serum to maternal serum from the C8 health study
# 1.27 (Shin et al. 2011)
EC_newborn = EC_Child_fix[,c(1,5)] %>% 
  left_join(EC_Mother, by="ID") %>% 
  ungroup() %>% 
  select(-1,-2)
EC_Child_fix_1 = EC_Child_fix_1[,-1]+
  ((ChildAge!=0)*(apply((EC_newborn)*(ChildAge==1)*1.27,1,sum))*exp(-KElim_fix * (ChildAge-1)))

# Add the background serum concentration
# Interpolated NHANES serum levels
for (i in 1:nrow(W1)) {
  for (r in 1:58){
    if(r>=(EC_Child_fix[i,3]-1950+2)){
      EC_Child_fix_1[i, r] <- (EC_Child_fix_1[i, r] + 0.11 *r) * (r < 49) +  # Year < 1999
        (EC_Child_fix_1[i, r] + 5.2 - 0.33*(r-49))* (r >= 49)  # Year >= 1999
    }
  }
}

# Select serum concentrations at serum collection (visit) year
visit = as.numeric(str_sub(EC_Child_fix$visitdate,start = -4))-1950
EC_Child_fix$C8.fixed <- sapply(1:length(visit), function(i) EC_Child_fix_1[i, visit[i]])
