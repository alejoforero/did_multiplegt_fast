////////////////////////////////////////////////////////////////////////////////
///// Program #4: performs outcome change residualisation, and requests ////////
///// computation of all point estimates asked by user /////////////////////////
////////////////////////////////////////////////////////////////////////////////

capture program drop did_multiplegt_estim2
program did_multiplegt_estim2, eclass
	version 12.0
	syntax varlist(min=4 numeric) [if] [in] [, THRESHOLD_stable_treatment(real 0) trends_nonparam(varlist numeric) trends_lin(varlist numeric) controls(varlist numeric) counter(varlist numeric) placebo(integer 0) dynamic(integer 0) breps(integer 50) cluster(varlist numeric) bootstrap_rep(integer 0) switchers(string) LONGdiff_placebo FIRSTdiff_placebo count_switchers_tot count_switchers_contr robust_dynamic discount(real 1) drop_larger_lower always_trends_nonparam always_trends_lin]


tempvar tag_obs group_incl

if `bootstrap_rep'==0{
gen tag_switchers_contr_XX=0
}
	
if "`trends_lin'" !=""|"`controls'" !=""{

sum d_cat_group_XX, meanonly
local D_min=r(min)
local D_max=r(max)

forvalue d=`=`D_min''/`=`D_max'' {

global cond_increase "abs(diff_stag_d_XX)>`threshold_stable_treatment'&increase_d_XX==1&lag_d_cat_group_XX==`d'&diff_stag_d_XX!=."
global cond_stable "abs(diff_stag_d_XX)<=`threshold_stable_treatment'&ever_change_d_XX==0&lag_d_cat_group_XX==`d'"
global cond_decrease "abs(diff_stag_d_XX)>`threshold_stable_treatment'&increase_d_XX==0&lag_d_cat_group_XX==`d'&diff_stag_d_XX!=."

sum diff_y_XX if $cond_increase, meanonly
scalar n_increase=r(N)
sum diff_y_XX if $cond_stable, meanonly
scalar n_stable=r(N)
sum diff_y_XX if $cond_decrease, meanonly
scalar n_decrease=r(N)

// Assessing if the residualization needs to be done for that treatment value
if ("`switchers'"==""&n_stable>0&(n_increase>0|n_decrease>0))|("`switchers'"=="in"&n_stable>0&n_increase>0)|("`switchers'"=="out"&n_stable>0&n_decrease>0) {

sum diff_y_XX if $cond_stable, meanonly

// Assessing if too many controls
if r(N)>=1{
scalar too_many_controls_temp=(r(N)<wordcount("`controls'"))
scalar too_many_controls=max(too_many_controls,too_many_controls_temp)
}

///////////////// Regression of diff_y on controls and FEs of trends_lin, and storing coefficients

cap drop FE*

if "`trends_lin'"!=""{

capture $noisily reghdfe diff_y_XX `controls' [aweight=`counter'] if $cond_stable, absorb(FE1=`trends_lin' FE2=time_XX) resid keepsingletons

}

if "`trends_lin'"==""&"`trends_nonparam'"==""&"`controls'"!=""{

capture $noisily reghdfe diff_y_XX `controls' [aweight=`counter'] if $cond_stable, absorb(FE1=time_XX) resid keepsingletons

}

if "`trends_nonparam'"!=""&"`controls'"!=""{

capture $noisily reghdfe diff_y_XX `controls' [aweight=`counter'] if $cond_stable, absorb(FE1=trends_var_XX) resid keepsingletons

}

capture matrix didmgt_B = e(b)

// Patching if not enough observations in this regression
if _rc!=301{

if "`trends_lin'"!=""{

gen `tag_obs'=e(sample)
bys `trends_lin': egen `group_incl'=max(`tag_obs') 

fcollapse (mean) FE_1=FE1, by(`trends_lin') merge

sum FE_1 [aweight=`counter'] if lag_d_cat_group_XX==`d'&`group_incl'==1

** AJOUT MELITINE **************************************************************
if "`always_trends_lin'"==""{
	replace FE_1=r(mean) if lag_d_cat_group_XX==`d'&`group_incl'==0 
	// what does the r(mean) do exactly here ?
}

if "`always_trends_lin'"!=""{
	drop if lag_d_cat_group_XX==`d'&`group_incl'==0
}

********************************************************************************

if `bootstrap_rep'==0{
replace tag_switchers_contr_XX=1 if `group_incl'==1&($cond_increase | $cond_decrease) 
}
}

// Creating variables with controls coefficients
local j = 0
foreach var of local controls {
local j = `j' + 1
gen coeff`j' = didmgt_B[1,`j']
}

///////////////////// Residualizing outcome changes

// Current outcome FD, for instantaneous effect estimation
if "`trends_lin'"!=""{
replace diff_y_XX=diff_y_XX-FE_1 if lag_d_cat_group_XX==`d'
}

local j=0
foreach var of local controls{
local j=`j'+1
replace diff_y_XX=diff_y_XX-coeff`j'*ZZ_cont`j' if lag_d_cat_group_XX==`d'
}

// Lagged outcome FD, for FD placebo estimation

** MODIFICATION MELITINE ** "`firstdiff_placebo'"!="" remplace "`longdiff_placebo'"==""
if "`placebo'"!="0"&"`firstdiff_placebo'"!=""{

forvalue i=1/`=`placebo''{

if "`trends_lin'"!=""{
replace diff_y_lag`i'_XX=diff_y_lag`i'_XX-FE_1 if lag_d_cat_group_XX==`d'
}

local j=0
foreach var of local controls{
local j=`j'+1
replace diff_y_lag`i'_XX=diff_y_lag`i'_XX-coeff`j'*ZZ_cont_lag`i'_`j' if lag_d_cat_group_XX==`d'
}

}

}

// Lagged outcome long diff, for long diff placebo estimation

** MODIFICATION MELITINE ** "`firstdiff_placebo'"=="" remplace "`longdiff_placebo'"!=""
if "`placebo'"!="0"&"`firstdiff_placebo'"==""{

forvalue i=1/`=`placebo''{

if "`trends_lin'"!=""{
replace ldiff_y_`i'_lag_XX=ldiff_y_`i'_lag_XX-FE_1*`i' if lag_d_cat_group_XX==`d'
}

local j=0
foreach var of local controls{
local j=`j'+1
replace ldiff_y_`i'_lag_XX=ldiff_y_`i'_lag_XX-coeff`j'*ZZ_cont_ldiff_`i'_lag_`j' if lag_d_cat_group_XX==`d'
}

}

}

// Lead outcome long diff, for dynamic effect estimation
if "`dynamic'"!="0"{

forvalue i=1/`=`dynamic''{

if "`trends_lin'"!=""{
replace ldiff_y_for`i'_XX=ldiff_y_for`i'_XX-FE_1*(`i'+1) if lag_d_cat_group_XX==`d'
}

local j=0
foreach var of local controls{
local j=`j'+1
replace ldiff_y_for`i'_XX=ldiff_y_for`i'_XX-coeff`j'*ZZ_cont_ldiff_for`i'_`j' if lag_d_cat_group_XX==`d'
}

}

}

cap drop `group_incl' `tag_obs'
cap drop FE*
cap drop coeff*

/// End of patch if not enough observations in the regression with controls
}

/// End of condition assessing if residualization needed for that treatment value
}

// End of the loop over values of D_cat
}

// End of the if condition assessing if trends_lin or controls requested
}


// Counting time periods

sum time_XX, meanonly
local max_time=r(max)

// Estimating the instantaneous effect

*Running did_multiplegt_core
did_multiplegt_core2, threshold_stable_treatment(`threshold_stable_treatment') trends_nonparam(`trends_nonparam') trends_lin(`trends_lin') controls(`controls') placebo(`placebo') d_cat_group(d_cat_group_XX) lag_d_cat_group(lag_d_cat_group_XX) diff_d(diff_d_XX) diff_y(diff_y_XX) counter(`counter') time(time_XX) group_int(group_XX) max_time(`max_time') counter_placebo(0) counter_dynamic(0) bootstrap_rep(`bootstrap_rep') switchers(`switchers') `robust_dynamic' discount(`discount') trends_lin(`trends_lin') `drop_larger_lower' `always_trends_nonparam' `always_trends_lin' `firstdiff_placebo' `longdiff_placebo'

*Collecting point estimate and number of observations

scalar effect_0_2=effect_XX
scalar N_effect_0_2=N_effect
scalar N_switchers_effect_0_2=N_switchers
if "`count_switchers_tot'"!=""{
scalar N_switchers_effect_0_tot_2=N_switchers_tot2
}
if "`count_switchers_contr'"!=""&("`trends_lin'"!=""|"`trends_nonparam'"!=""){
scalar N_switchers_effect_0_contr_2=N_switchers_contr2
} 
scalar denom_DID_ell_0=denom_DID_ell_XX
scalar denom_delta_0=denom_XX

// If first difference placebos requested, estimate them and number of observations used in that estimation

** MODIFICATION MELITINE ** "`firstdiff_placebo'"!="" remplace "`longdiff_placebo'"==""
if "`placebo'"!="0"&"`firstdiff_placebo'"!=""{

tempvar cond_placebo  
gen `cond_placebo'=1

*Looping over the number of placebos requested

forvalue i=1/`=`placebo''{

*Replacing FD of outcome by lagged FD of outcome, FD of controls by lagged FD of controls, and excluding from placebo observations whose lagged FD of treatment non 0. 

// Note: the line below is superfluous if the robust_dynamic option is specified, because then 
// only (g,t)s with diff_stag_d_XX>`threshold_stable_treatment', meaning those changing treatment for the first time at t satisfy "cond_increase_t" and "cond_decrease_t" anyways.
// But that line plays a role if the robust_dynamic option is not specified so it is important to keep it. 
replace `cond_placebo'=0 if abs(diff_d_lag`i'_XX)>`threshold_stable_treatment'

preserve

replace diff_y_XX=diff_y_lag`i'_XX

if "`controls'" !=""{
local j=0
foreach var of local controls{
local j=`j'+1
replace `var'=ZZ_cont_lag`i'_`j'
}
}

*If no observation satisfy `cond_placebo'==1, set N_placebo_`i'_2 to 0

sum diff_y_XX if `cond_placebo'==1

if r(N)==0{

scalar N_placebo_`i'_2=.
scalar placebo_`i'_2=.
scalar N_switchers_placebo_`i'_2=.
if "`count_switchers_tot'"!=""{
scalar N_switchers_placebo_`i'_tot_2=.
}
if "`count_switchers_contr'"!=""&("`trends_lin'"!=""|"`trends_nonparam'"!=""){
scalar N_switchers_placebo_`i'_contr_2=.
} 
}

*Otherwise, run did_multiplegt_core

else{

did_multiplegt_core2 if `cond_placebo'==1, threshold_stable_treatment(`threshold_stable_treatment') trends_nonparam(`trends_nonparam') trends_lin(`trends_lin') controls(`controls') placebo(`placebo') d_cat_group(d_cat_group_XX) lag_d_cat_group(lag_d_cat_group_XX) diff_d(diff_d_XX) diff_y(diff_y_XX) counter(`counter') time(time_XX) group_int(group_XX) max_time(`max_time') counter_placebo(`i') counter_dynamic(0) bootstrap_rep(`bootstrap_rep') switchers(`switchers') `robust_dynamic' discount(`discount') trends_lin(`trends_lin') `drop_larger_lower' `always_trends_nonparam' `always_trends_lin' `firstdiff_placebo' `longdiff_placebo'

*Collecting point estimate and number of observations

scalar placebo_`i'_2=effect_XX
scalar N_placebo_`i'_2=N_effect
scalar N_switchers_placebo_`i'_2=N_switchers
if "`count_switchers_tot'"!=""{
scalar N_switchers_placebo_`i'_tot_2=N_switchers_tot2
}
if "`count_switchers_contr'"!=""&("`trends_lin'"!=""|"`trends_nonparam'"!=""){
scalar N_switchers_placebo_`i'_contr_2=N_switchers_contr2
} 

}

restore

*End of the loop on the number of placebos
}

*End of the condition assessing if the computation of placebos was requested by the user
}


// If long-difference placebos requested, estimate them and number of observations used in that estimation

** MODIFICATION MELITINE ** "`firstdiff_placebo'"=="" remplace "`longdiff_placebo'"!=""
if "`placebo'"!="0"&"`firstdiff_placebo'"==""{

tempvar cond_placebo
gen `cond_placebo'=1

*Looping over the number of placebos requested

forvalue i=1/`=`placebo''{

*Replacing FD of outcome by long diff of outcome, and creating variable to exclude from placebo observations whose lead FD of treatment non 0. 

if `i'>1{
replace `cond_placebo'=0 if abs(diff_stag_d_for`=`i'-1'_XX)>`threshold_stable_treatment'
}

preserve

replace diff_y_XX=ldiff_y_`i'_lag_XX
if `i'>1{
replace `counter'=counter_F`=`i'-1'_XX
}
if "`controls'" !=""{
local j=0
foreach var of local controls{
local j=`j'+1
replace `var'=ZZ_cont_ldiff_`i'_lag_`j'
}
}

*If no observation satisfy `cond_placebo'==1, set N_placebo_`i'_2 to 0

sum diff_y_XX if `cond_placebo'==1

if r(N)==0{
scalar N_placebo_`i'_2=.
scalar placebo_`i'_2=.
scalar N_switchers_placebo_`i'_2=.
if "`count_switchers_tot'"!=""{
scalar N_switchers_placebo_`i'_tot_2=.
}
if "`count_switchers_contr'"!=""&("`trends_lin'"!=""|"`trends_nonparam'"!=""){
scalar N_switchers_placebo_`i'_contr_2=.
} 
}

*Otherwise, run did_multiplegt_core

else{

did_multiplegt_core2 if `cond_placebo'==1, threshold_stable_treatment(`threshold_stable_treatment') trends_nonparam(`trends_nonparam') trends_lin(`trends_lin') controls(`controls') placebo(`placebo') d_cat_group(d_cat_group_XX) lag_d_cat_group(lag_d_cat_group_XX) diff_d(diff_d_XX) diff_y(diff_y_XX) counter(`counter') time(time_XX) group_int(group_XX) max_time(`max_time') counter_placebo(`i') counter_dynamic(0) bootstrap_rep(`bootstrap_rep') switchers(`switchers') `robust_dynamic' discount(`discount') trends_lin(`trends_lin') `drop_larger_lower' `always_trends_nonparam' `always_trends_lin' `firstdiff_placebo' `longdiff_placebo'

*Collecting point estimate and number of observations

scalar placebo_`i'_2=-effect_XX
scalar N_placebo_`i'_2=N_effect
scalar N_switchers_placebo_`i'_2=N_switchers
if "`count_switchers_tot'"!=""{
scalar N_switchers_placebo_`i'_tot_2=N_switchers_tot2
}
if "`count_switchers_contr'"!=""&("`trends_lin'"!=""|"`trends_nonparam'"!=""){
scalar N_switchers_placebo_`i'_contr_2=N_switchers_contr2
} 
}

restore

*End of the loop on the number of placebos
}

*End of the condition assessing if the computation of long-diff placebos was requested by the user
}


// If dynamic effects requested, estimate them and number of observations used in that estimation

if "`dynamic'"!="0"{

tempvar cond_dynamic
gen `cond_dynamic'=1

*Looping over the number of placebos requested

forvalue i=1/`=`dynamic''{

*Replacing FD of outcome by long diff of outcome, and creating variable to exclude from placebo observations whose lead FD of treatment non 0. 

replace `cond_dynamic'=0 if abs(diff_stag_d_for`i'_XX)>`threshold_stable_treatment'

preserve

replace diff_y_XX=ldiff_y_for`i'_XX
replace `counter'=counter_F`i'_XX

replace diff_d_XX=ldiff_d_for`i'_XX

if "`controls'" !=""{
local j=0
foreach var of local controls{
local j=`j'+1
replace `var'=ZZ_cont_ldiff_for`i'_`j'
}
}

*If no observation satisfy `cond_dynamic'==1, set N_effect_`i'_2 to 0

sum diff_y_XX if `cond_dynamic'==1

if r(N)==0{

scalar N_effect_`i'_2=.
scalar N_switchers_effect_`i'_2=.
scalar effect_`i'_2=.
if "`count_switchers_tot'"!=""{
scalar N_switchers_effect_`i'_tot_2=.
}
if "`count_switchers_contr'"!=""&("`trends_lin'"!=""|"`trends_nonparam'"!=""){
scalar N_switchers_effect_`i'_contr_2=.
} 
}

*Otherwise, run did_multiplegt_core

else{


did_multiplegt_core2 if `cond_dynamic'==1, threshold_stable_treatment(`threshold_stable_treatment') trends_nonparam(`trends_nonparam') trends_lin(`trends_lin') controls(`controls') placebo(`placebo') d_cat_group(d_cat_group_XX) lag_d_cat_group(lag_d_cat_group_XX) diff_d(diff_d_XX) diff_y(diff_y_XX) counter(`counter') time(time_XX) group_int(group_XX) max_time(`max_time') counter_placebo(0) counter_dynamic(`i') bootstrap_rep(`bootstrap_rep') switchers(`switchers') `robust_dynamic' discount(`discount') trends_lin(`trends_lin') `drop_larger_lower' `always_trends_nonparam' `always_trends_lin' `firstdiff_placebo' `longdiff_placebo'

*Collecting point estimate and number of observations

scalar effect_`i'_2=effect_XX
scalar N_effect_`i'_2=N_effect
scalar N_switchers_effect_`i'_2=N_switchers
if "`count_switchers_tot'"!=""{
scalar N_switchers_effect_`i'_tot_2=N_switchers_tot2
}
if "`count_switchers_contr'"!=""&("`trends_lin'"!=""|"`trends_nonparam'"!=""){
scalar N_switchers_effect_`i'_contr_2=N_switchers_contr2
} 
scalar denom_DID_ell_`i'=denom_DID_ell_XX
scalar denom_delta_`i'=denom_XX
}

restore

drop diff_stag_d_for`i'_XX 

*End of the loop on the number of dynamic effects
}

*End of the condition assessing if the computation of dynamic effects was requested by the user
}

end

////////////////////////////////////////////////////////////////////////////////
///// Program #5: performs computation of all individual point estimates ///////
////////////////////////////////////////////////////////////////////////////////

capture program drop did_multiplegt_core2
program did_multiplegt_core2
	version 12.0
	syntax [if] [in] [, THRESHOLD_stable_treatment(real 0) trends_nonparam(varlist numeric) trends_lin(varlist numeric) controls(varlist numeric) placebo(integer 0) d_cat_group(varlist numeric) lag_d_cat_group(varlist numeric) diff_d(varlist numeric) diff_y(varlist numeric) counter(varlist numeric) time(varlist numeric) group_int(varlist numeric) max_time(integer 0) counter_placebo(integer 0) counter_dynamic(integer 0) bootstrap_rep(integer 0) switchers(string) robust_dynamic discount(real 1) trends_lin(varlist numeric) drop_larger_lower always_trends_nonparam always_trends_lin FIRSTdiff_placebo LONGdiff_placebo]

tempvar diff_y_res1 diff_y_res2 group_incl treatment_dummy tag_obs tag_switchers counter_tot counter_switchers tag_switchers_tot counter_switchers_tot tag_switchers_contr counter_switchers_contr

preserve

// Selecting the sample

	if "`if'" !=""{
	keep `if'
	}
	
// Drop if diff_y_XX missing, to avoid that those observations are used in estimation

drop if diff_y_XX==.
	
// Creating residualized first diff outcome if trends_nonparam specified in estimation

if "`trends_nonparam'" !=""{

sum d_cat_group_XX, meanonly
local D_min=r(min)
local D_max=r(max)

forvalue d=`=`D_min''/`=`D_max'' {

global cond_increase "abs(diff_stag_d_XX)>`threshold_stable_treatment'&increase_d_XX==1&lag_d_cat_group_XX==`d'&diff_stag_d_XX!=."

global cond_stable "abs(diff_stag_d_XX)<=`threshold_stable_treatment'&ever_change_d_XX==0&lag_d_cat_group_XX==`d'"
global cond_decrease "abs(diff_stag_d_XX)>`threshold_stable_treatment'&increase_d_XX==0&lag_d_cat_group_XX==`d'&diff_stag_d_XX!=."

sum diff_y_XX if $cond_stable, meanonly

/////////////// Regression of diff_y on time FEs interacted with trends_nonparam variable, and computation of residuals

cap drop FE*
capture $noisily reghdfe diff_y_XX  [aweight=`counter'] if $cond_stable, absorb(FE1=trends_var_XX) resid keepsingletons

// Patching if not enough observations in this regression
capture matrix didmgt_B = e(b)
if _rc!=301{ // r(301) in an error for "last estimate not found"

gen `tag_obs'=e(sample)
bys trends_var_XX: egen `group_incl'=max(`tag_obs')

fcollapse (mean) FE=FE1 , by(trends_var_XX) merge

gen `diff_y_res1'  = diff_y_XX - FE
gen constant = didmgt_B[1,1]
replace `diff_y_res1' = `diff_y_res1' - constant

cap drop FE constant

** AJOUT MELITINE **************************************************************

if "`always_trends_nonparam'"==""{
$noisily reg diff_y_XX i.time_XX [aweight=`counter'] if $cond_stable, $no_header_no_table
predict `diff_y_res2', r
replace diff_y_XX=`diff_y_res2' if lag_d_cat_group_XX==`d'&`group_incl'==0
}

replace diff_y_XX=`diff_y_res1' if lag_d_cat_group_XX==`d'&`group_incl'==1


if "`always_trends_nonparam'"!=""{
drop if lag_d_cat_group_XX==`d'&`group_incl'==0
}

********************************************************************************

if `bootstrap_rep'==0{
replace tag_switchers_contr_XX=1 if `group_incl'==1&($cond_increase | $cond_decrease) 
}

** RETRAIT MELITINE - plus de `diff_y_res2' dans le drop, exécuté après
drop `diff_y_res1' `group_incl' `tag_obs'


** AJOUT MELITINE **************************************************************
if "`always_trends_nonparam'"==""{
drop `diff_y_res2'
}
********************************************************************************

}
// End of the loop over values of D_cat
}

// End of the if condition assessing if trends_nonparam included in estimation
}
 
// Treatment effect

// Initializing estimate, weight, and variable to count observations used in estimation
scalar effect_XX=0
scalar N_effect =0
scalar N_switchers=0
scalar N_switchers_tot2=0
scalar N_switchers_contr2=0
scalar denom_XX=0
scalar denom_DID_ell_XX=0
gen `tag_obs'=0
gen `tag_switchers'=0
gen `tag_switchers_tot'=0

$noisily di "Computing DIDM"

// Looping over time periods
forvalue t=`=`counter_placebo'+2'/`=`max_time'-`counter_dynamic''{

// Determining the min and max value of group of treatment at t-1

sum lag_d_cat_group_XX if time_XX==`t', meanonly

local D_min=r(min)
local D_max=r(max)



// Ensuring that there are observations with non missing lagged treatment

if `D_min'!=.&`D_max'!=.{

// Looping over possible values of lag_D at time t
forvalue d=`=`D_min''/`=`D_max'' {

// Defining conditions for groups where treatment increased/remained stable/decreased between t-1 and t

// Note: If robust_dynamic option specified, 
// cond_increase_t=those:  
// 1) whose treatment changes for first time at t (t=t-\ell in paper): abs(diff_stag_d_XX)>`threshold_stable_treatment'
// 2) whose treatment cost is higher than if they had kept status quo treatment: increase_d_XX==1
// Same thing for cond_decrease_t
// cond_stable_t= those: 
// 1) whose treatment does not change at t, 
// 2) whose treatment has not changed before t (was implied by 1 in staggered designs, not the case anymore, hence the addition of ever_change_d_XX==0)
// those two conditions are sufficient to select the right control groups when we estimate the instantaneous treatment effect, or first difference placebos,
// but if we are running this to estimate a dynamic effect at time t+\ell we need a third condition: 
// 3) the ``if `cond_dynamic'==1'' condition when we call the command in program #4 above ensures we do not have obs whose treatment changed for the first time somehwere between t+1 and t+\ell. 
// If robust_dynamic option not specified, ever_change_d_XX==0, diff_stag_d_XX=diff_d_XX, and increase_d_XX=(diff_d_XX>0) so the new conditions below are equal to the old ones commented above.
 
global cond_increase_t "abs(diff_stag_d_XX)>`threshold_stable_treatment'&increase_d_XX==1&lag_d_cat_group_XX==`d'&time_XX==`t'&diff_stag_d_XX!=."
global cond_stable_t "abs(diff_stag_d_XX)<=`threshold_stable_treatment'&ever_change_d_XX==0&lag_d_cat_group_XX==`d'&time_XX==`t'"
global cond_decrease_t "abs(diff_stag_d_XX)>`threshold_stable_treatment'&increase_d_XX==0&lag_d_cat_group_XX==`d'&time_XX==`t'&diff_stag_d_XX!=."

// Counting number of units in each supergroup
sum d_cat_group_XX if $cond_increase_t, meanonly
scalar n_increase=r(N)
sum d_cat_group_XX if $cond_stable_t, meanonly
scalar n_stable=r(N)
sum d_cat_group_XX if $cond_decrease_t, meanonly
scalar n_decrease=r(N)

// If there are units whose treatment increased and units whose treatment remained stable, estimate corresponding DID, 
// increment point estimate and weight, and tag observations used in estimation

if "`switchers'" !="out"{
if `bootstrap_rep'==0{
replace `tag_switchers_tot'=1 if $cond_increase_t
}
if n_increase*n_stable>0 {
gen `treatment_dummy' =($cond_increase_t)

if `bootstrap_rep'==0{
replace `tag_obs'=1 if (($cond_increase_t)|($cond_stable_t))
replace `tag_switchers'=1 if $cond_increase_t
}

$noisily reg diff_y_XX `treatment_dummy' [aweight=`counter'] if ($cond_increase_t)|($cond_stable_t), $no_header_no_table
sum `counter' if $cond_increase_t, meanonly
scalar effect_XX=effect_XX+_b[`treatment_dummy']*r(N)*r(mean)*(`discount'^`t')
$noisily reg diff_d_XX `treatment_dummy' [aweight=`counter'] if ($cond_increase_t)|($cond_stable_t), $no_header_no_table
sum `counter' if $cond_increase_t, meanonly
scalar denom_XX=denom_XX+_b[`treatment_dummy']*r(N)*r(mean)*(`discount'^`t')
scalar denom_DID_ell_XX=denom_DID_ell_XX+r(N)*r(mean)*(`discount'^`t')
drop `treatment_dummy' 
}
}

// If there are units whose treatment decreased and units whose treatment remained stable, estimate corresponding DID, 
// increment point estimate and weight, and tag observations used in estimation

if "`switchers'"!="in"{
if `bootstrap_rep'==0{
replace `tag_switchers_tot'=1 if $cond_decrease_t
}
if n_decrease*n_stable>0 {
gen `treatment_dummy' =($cond_decrease_t)

if `bootstrap_rep'==0{
replace `tag_obs'=1 if (($cond_decrease_t)|($cond_stable_t))
replace `tag_switchers'=1 if $cond_decrease_t
}

$noisily reg diff_y_XX `treatment_dummy' [aweight=`counter'] if ($cond_decrease_t)|($cond_stable_t), $no_header_no_table
sum `counter' if $cond_decrease_t, meanonly
scalar effect_XX=effect_XX-_b[`treatment_dummy']*r(N)*r(mean)*(`discount'^`t')
$noisily reg diff_d_XX `treatment_dummy' [aweight=`counter'] if ($cond_decrease_t)|($cond_stable_t), $no_header_no_table
sum `counter' if $cond_decrease_t, meanonly
scalar denom_XX=denom_XX-_b[`treatment_dummy']*r(N)*r(mean)*(`discount'^`t')
scalar denom_DID_ell_XX=denom_DID_ell_XX+r(N)*r(mean)*(`discount'^`t')
drop `treatment_dummy' 
}
}

// End of loop on recat treatment values at t-1 
}

// End of condition ensuring that there are observations with non missing lagged treatment
}

// End of loop on time
}

if "`robust_dynamic'"==""{
scalar effect_XX=effect_XX/denom_XX
}

if "`robust_dynamic'"!=""{
scalar effect_XX=effect_XX/denom_DID_ell_XX
}

if `bootstrap_rep'==0{
egen `counter_tot'=total(`counter') if `tag_obs'==1
sum `counter_tot', meanonly
scalar N_effect=r(mean)
egen `counter_switchers'=total(`counter') if `tag_switchers'==1
sum `counter_switchers', meanonly
scalar N_switchers=r(mean)
egen `counter_switchers_tot'=total(`counter') if `tag_switchers_tot'==1
sum `counter_switchers_tot', meanonly
scalar N_switchers_tot2=r(mean)
egen `counter_switchers_contr'=total(`counter') if tag_switchers_contr_XX==1&`tag_switchers'==1
sum `counter_switchers_contr', meanonly
scalar N_switchers_contr2=r(mean)
}

restore

end
