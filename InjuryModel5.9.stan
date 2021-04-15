//    InjuryModel5.9.stan
//
//    Last update: September 7, 2020
//    Authors: Andrew R Morral
//
//    Description: Full Bayes imputation of missing state-year firearm injury hospital  
//			   discharges, with imputation of missing predictors of those discharges 
//
//
//    Note on missing value indices: Each missing value index below has the same length 
//    as the data. The index is all 0's, except for those positions where the indexed 
//    variable has a missing value. In those positions, sequential positive numbers are 
//    inserted. So a missing variable index length 10 with 4 missing values might look 
//    like: (0,1,0,0,2,0,3,4,0,0).


functions {
  //replace missing values (999) in data matrix with parameters  
  matrix merge_missingM( int[,] miss_indexes , matrix x_obs , vector x_miss ) {
  	int N = dims(x_obs)[1];
  	int P = dims(x_obs)[2];
  	int counter = 1;
  	matrix [N,P] merged;
  	merged = x_obs;
  	for ( i in 1:P ){
  	  if (sum(miss_indexes[,i])>0){
  	    for ( j in 1:N ){
  	      if( miss_indexes[j,i] > 0 ) {
  	        merged[j,i] = x_miss[ counter ];
  	        counter = counter + 1;
  	      }
  	    }	  
  	  }
  	}
  	return merged; 
  }
}


data { 
  int<lower=0> num_data; 						  //Number of observations in dataset
  int<lower=0> num_pred;						  //Number of predictor variables
  int<lower=0> num_year;						  //Number of years of data
  int<lower=0> num_state;						  //Number of states
  int<lower=0> num_int;							  //Number of interaction terms
  
  matrix[num_data,1] Discharges; 			//Outcome variable
  int DissMiss[num_data,1];						//Index of missing outcome variables
  matrix[num_data,num_pred] Pred; 		//Matrix with predictors
  int PredMiss[num_data,num_pred]; 		//Predictor missing values indices array
  int year_index[num_data];						//Year index
  int state_index[num_data];					//State index
  int state_year_index[num_year,num_state];		//State-year index
  int<lower=0> int_mat[num_int,2];			   	  //Matrix of interacted variables w/missing
  matrix<lower=0>[num_year,num_state] popinv; //Inverse of population
  matrix<lower=0>[num_data,1] completeness;		//Completeness of injury coding
  int CompMiss[num_data,1];						//Completeness missing value index
} 

transformed data{
  vector[num_pred] mu = rep_vector( 0, num_pred);   //Prior for mean predictor value
  vector[num_year] mu_c = rep_vector( .9, num_year);//Prior for mean completeness value
  int<lower=0> miss_count = 0;
  real<lower=0> yrs[num_year];
  
  for (i in 1:num_year)
    yrs[i] = i;									      //Year number

  for (i in 1:num_pred)
    miss_count += max(PredMiss[,i]);	//Count of missing predictors
}

parameters { 
  real<lower=0> beta0;								//Discharge model intercept
  vector[num_pred] betaPred;					//Predictor coefficients
  vector[num_int] betaInt;						//Interactions coefficients
  real<lower=0> hyper;							  //hyperparam for betas priors
  vector<lower=0,upper=1>[num_state] d;			//Proportion uncoded injuries that are FA
  real<lower=0> a;								    //hyperparam for beta dist for d
  real<lower=0> b;								    //hyperparam for beta dist for d
  vector<lower=0>[max(DissMiss[,1])] Diss_impute;	//Outcome missing value params
  vector[miss_count] Pred_impute;			//Predictor missing value params
  vector<lower=0,upper=1>[max(CompMiss[,1])] Comp_impute; //Completeness missing val params
  corr_matrix[num_pred] gamma;        //Prior corr'n between predictors
  real<lower=0> etasq; 							  //max correlation parameter: discharges
  real<lower=0> rhosq; 							  //decay rate: discharge corr's
  real<lower=0> etasq_c; 						  //max correlation parameter: completeness
  real<lower=0> rhosq_c; 						  //decay rate: completeness corr's
  real<lower=0> covPrior;						  //hyperprior for pred covariance shape 
  real<lower=0> se;								    //param controlling pop effects on cov mat
}

transformed parameters{							
  matrix[num_data,1] Diss_hat;				      //Prior for missing discharge values
  matrix<lower=0>[num_data,1] Diss_merge;		//Discharge w/missing val parameters merged
  matrix<lower=0,upper=1>[num_data,1] Comp_merge; //Completeness w/missing val params merged
  matrix[num_data,num_pred] Pred_merge;			//Pred matrix with missing val params merged
  matrix[num_data,num_int] interactions;		//Interaction terms
  
  
  												            //Merge missing value parameters:
  Diss_merge = merge_missingM( DissMiss, Discharges, Diss_impute );
  Comp_merge = merge_missingM( CompMiss, completeness, Comp_impute );
  Pred_merge = merge_missingM( PredMiss, Pred, Pred_impute );

  for ( i in 1:num_int){						  //Create interaction terms
    interactions[,i] = Pred_merge[,int_mat[i,1]] .* Pred_merge[,int_mat[i,2]];
    interactions[,i] = interactions[,i] - mean(interactions[,i]);
  }
  
  for (i in 1:num_data){              //Model of expected observed discharge rate
    Diss_hat[i] = ( beta0 + Pred_merge[i,]*betaPred + interactions[i,]*betaInt ) *     
  			 ( 1- d[state_index[i]] * (1 - Comp_merge[i]) ) ; 	
  }
}
 
model { 
  matrix[num_year,num_year] L_K;      // Cholesky decomposition of discharge covariance matrix
  matrix[num_year,num_year] L_K_c;    // Cholesky decomposition of completeness covariance matrix
  matrix[num_year,num_year] SIGMA = cov_exp_quad( yrs[], etasq, rhosq);	//Discharge covariance matrix
  matrix[num_year,num_year] SIGMA_c = cov_exp_quad( yrs[], etasq_c, rhosq_c);	//Completeness covariance matrix
  vector[num_year] Diag;              // Vector used to retain diagonal of SIGMA
  
  hyper ~ normal(0,.5);								// hyperparam for beta priors
  beta0 ~ normal(0,hyper);						// model intercept
  betaPred ~ normal(0,hyper);					// coefficients main effects		
  betaInt ~ normal(0,hyper);					// coefficients interaction terms
  etasq ~ std_normal(); 							// prior on max correlation between years
  rhosq ~ inv_gamma(5,5); 						// prior on decay rate
  etasq_c ~ std_normal(); 						// prior on max correlation between years
  rhosq_c ~ inv_gamma(5,5); 					// prior on decay rate
  gamma ~ lkj_corr(covPrior);					// prior for predictor covariance matrix
  covPrior ~ normal(1,2);             // prior for lkj_corr shape parameter
  se ~ normal(0,.5);                  // prior for inv_pop effect on discharge variance
  d ~ beta(a,b);                      // prior for effect of completeness on observed d/c rate
  a ~ gamma(2,1);                     // hyperprior for d
  b ~ gamma(2,1);                     // hyperprior for d
  
  for ( m in 1:num_data )             // likelihood/prior for predictor observations/paramters
    Pred_merge[m,] ~ multi_normal( mu, gamma );  	//centered at 0 b/c normalized
  
  L_K_c = cholesky_decompose(SIGMA_c);
  Diag = diagonal(SIGMA);
  for ( k in 1:num_state){ 
  	for ( m in 1:num_year )					
		SIGMA[m,m] = Diag[m] + se * popinv[m,k];
  	L_K = cholesky_decompose(SIGMA);
  	Diss_merge[state_year_index[,k],1] ~   // likelihood/prior for discharge obs/params
  		multi_normal_cholesky( Diss_hat[state_year_index[,k],1], L_K );
  	Comp_merge[state_year_index[,k],1] ~   // likelihood/prior for completeness obs/params
  		multi_normal_cholesky( mu_c, L_K_c );
  }
}

