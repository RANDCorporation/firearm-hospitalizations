########################################################################################
### AssembleInjuryData.r ###

## last update: 02/05/21
## authors:     Andrew R.Morral
##
## description: Reads in and organizes data



#################### Functions ###############################

### Create missing value index
miss.index = function(u){
  missIndex = ifelse( is.na(u) , 1, 0)
  missIndex = sapply( 1:length(missIndex), function(n) missIndex[n]*sum(missIndex[1:n]))
  return(as.matrix(missIndex))
}

### Impute missing completeness values linearly between known values, or to max of nearest values
imputer = function(u){
  count = sum(!is.na(u))
  ind = !is.na(u)
  if (count==1) u = rep(u[which(ind),length(u)])
  if (count>1) {
    instance = which(ind)
    numInst = length(instance)
    for (i in 1:length(u)){
      if (i<instance[1] & is.na(u[i])) u[i] = u[instance[1]]   
      if (i>instance[numInst] & is.na(u[i])) u[i] = u[instance[numInst]]
      if (i>instance[1] & i<instance[numInst] & is.na(u[i])) {
        bracketL = which(instance<i)[length(which(instance<i))]
        bracketH = which(instance>i)[1]
        u[i] = u[instance[bracketL]] + (i-instance[bracketL])*(u[instance[bracketH]]-u[instance[bracketL]])/(instance[bracketH]-instance[bracketL]+1)
      }
    }
  }
  return(u)
}


################### Assemble data #################################################################
##### Read in and merge discharge data, covariates, and completeness data
x = read.csv("RAND_State_Hospitalizations_Firearm_Injury_02052021.csv", stringsAsFactors = FALSE)
x = x[order(x$State, x$Year),]
x = x[ x$State != "District of Columbia", ]     

state_index = as.numeric(factor(x$State))
num_state = max(state_index)
year_index = as.numeric(factor(x$Year))
num_year = max(year_index)
state_year_index = sapply(unique(x$State),function(u) {which((x$Year %in% 2000:2016) & x$State==u)})

### Impute completeness, replace completeness with 1 if discharges missing ##
x$Completeness = x$Completeness/100
for(i in 1:num_state){
  state_completeness = x$Completeness[state_year_index[,i]]
  #2015 was a strange year for completeness because of introduction of ICD10. Instead of imputing
  #missing values based on this year's value, we strip it out, impute missing, then return the
  #2015 value (including missing values) to the completeness variable.
  temp_2015 = state_completeness[16]    
  state_completeness[16] = NA
  state_completeness = imputer(state_completeness)
  state_completeness[16] = temp_2015
  x$Completeness[state_year_index[,i]]= state_completeness
}
x$Completeness = ifelse(is.na(x$pcGD_I16),1,x$Completeness)
CompMiss = miss.index(x$Completeness)
### Replace remaining missing values with 999, so they will be estimated in the model
Completeness = as.matrix(ifelse(is.na(x$Completeness),999,x$Completeness))

### Assemble FA discharge information
Discharges = matrix(x$pcGD_I16 *10000,ncol=1)
DissMiss = miss.index(Discharges) #Index of missing discharge values
Discharges[is.na(Discharges)] = 999 #replace missing with 999
num_data = length(Discharges)

str(x)

### Other inputs
popinv = 1E6/sapply(unique(x$State),function(u) {x$Population[x$State==u]})
x$yearL = scale(x$Year)
x$yearSq = x$yearL^2

### Create normalized covariate matrix and collect matrix info
p.vars = c("smth_Inc2000","smth_pov","smth_pcTFnoS", 
           "smth_HFR","smth_pct_urban","smth_uninsured","pct_female","smth_pct_HS",
           "smth_pctBachelor","smth_pcCops","smth_pctWNH", "NIS_pct_Dxgun",
           "smth_pcVioGunCrime","yearL","yearSq")   
num_pred = length(p.vars)
Pred = scale(x[,p.vars])
PredMiss = apply(Pred,2,miss.index)
Pred = apply(Pred,2,function(u) ifelse(is.na(u),999,u))  ### Replace missing values again with 999


### Name all included interaction terms in which predictors have missing values.
### Interactions are constructed within the Stan model so that missing covariate values can be handled
### correctly
int1 = c(which(p.vars == "smth_pctWNH"),which(p.vars == "smth_pcVioGunCrime"))
int2 = c(which(p.vars == "smth_pct_urban"),which(p.vars == "smth_pcVioGunCrime"))
int3 = c(which(p.vars == "smth_HFR"),which(p.vars == "smth_pct_HS"))
int4 = c(which(p.vars == "smth_HFR"),which(p.vars == "smth_pctBachelor"))
int5 = c(which(p.vars == "smth_pcTFnoS"),which(p.vars == "smth_pct_urban"))
int_mat = rbind(int1,int2,int3,int4,int5)
num_int = dim(int_mat)[1]