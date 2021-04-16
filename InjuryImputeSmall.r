
###Copyright (C) 2021 The RAND Corporation###
###This code is licensed under the GPLv3 or later###

### InjuryImputeSource.r ###

## last update: 06/31/20
## authors:     Andrew R.Morral
##
## description: Imputes state firearm injury hospital discharges for missing years 
##              and tests model fit against held out observations

library(rstan)
options(mc.cores=parallel::detectCores())
rstan_options(auto_write = TRUE)

source("AssembleInjuryData.r")

### Assemble all Stan inputs
stan.list = list(Discharges=Discharges,DissMiss=DissMiss,num_data=num_data,
                 state_index=state_index,num_state=num_state,year_index=year_index,
                 num_year=num_year,num_pred=num_pred,Pred=Pred,PredMiss=PredMiss,
                 state_year_index=state_year_index,int_mat=int_mat, num_int=num_int,popinv=popinv, 
                 completeness=Completeness,CompMiss=CompMiss)

inits.list1 = list(beta0=1,betaPred=rep(0,num_pred), a=2, b=1, covPrior=1,
                   Diss_impute=rep(1,max(DissMiss)), Pred_Impute=rep(0,sum(apply(PredMiss,2,max))), 
                   gamma=diag(rep(1,num_pred)), se=.1, etasq=.2, rhosq=2, etasq_c=.2, rhosq_c=1)

#############  Run Stan ###############################################################################
sm<-stan_model("InjuryModel5.9.stan")
fit.1 <- sampling(sm, data = stan.list, iter = 10000, chains = 4, seed=12346,
                  init=list(inits.list1, inits.list1,inits.list1, inits.list1),
                  control=list( max_treedepth = 11)) #adapt_delta=.80,
save(fit.1, file="InjuryNatural.Rdat")
post = extract(fit.1)
save(post,file="postI.Rdat")
