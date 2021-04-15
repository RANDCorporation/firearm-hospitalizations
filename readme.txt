
Project Description

This Github repository includes statistical programs and data used to construct estimates of state firearm injury hospitalizations described in “Inpatient Hospitalizations for Firearm Injury: Estimating State-Level Rates from 2000 to 	2016” by Rosanna Smart, Samuel Peterson, Terry L. Schell, Rose Kerber, and Andrew R. Morral. RAND Corporation: Santa Monica, TL-A243-3, February, 2021. https://www.rand.org/pubs/tools/TLA243-3.html


Methods used

Stan

Technologies
R with package rstan

Getting started

1. Download "RAND_State_Hospitalizations_Firearm_Injury_02052021.csv" and the accompanying data dictionary
2. Download FullPost.csv for the full set of posterior draws from the stan model
3. Download scripts AssembleInjuryData.r, InjuryImputeSmall.r, and InjuryModel5.9.stan
4. Running InjuryImputeSmall.r will run the other 2 scripts

Other notes:

InjuryModel5.9.stan
-STAN model that produces FullPost.csv

FullPost.csv
-All posterior draws from STAN model

RAND_State_Hospitalizations_Firearm_Injury_02052021.csv
-includes set of covariates and 50 draws from STAN full posterior

AssembleInjuryData.r
-reads in "RAND_State_Hospitalizations_Firearm_Injury_02052021.csv"
-imputes completeness

InjuryImputeSmall.r
-runs AssembleInjuryData.r and InjuryModel5.9.stan

Project members:

Rosanna Smart (rsmart@rand.org), Samuel Peterson (speterso@rand.org), Terry L. Schell (tschell@rand.org), Rose Kerber (rkerber@rand.org), and Andrew R. Morral (morral@rand.org)

Suggested citation for this repository:

Smart, Rosanna, Peterson, Samuel, Schell, Terry L., Kerber, Rose, and Andrew Morral, “Inpatient Hospitalizations for Firearm Injury: Estimating State-Level Rates from 2000 to 2016” GitHub, RAND Corporation Repository, last updated 24 March 2021. As of March 24, 2021
