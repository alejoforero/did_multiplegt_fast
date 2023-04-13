// Generate a complete panel of 300 units observed in 15 periods
clear all
timer clear
set seed 10
global T = 15
global I = 3e2

set obs `=$I*$T'
gen i = int((_n-1)/$T )+1 					// unit id
gen t = mod((_n-1),$T )+1					// calendar period
tsset i t

// Randomly generate treatment rollout years uniformly across Ei=10..16 (note that periods t>=16 would not be useful since all units are treated by then)
gen Ei = ceil(runiform()*7)+$T -6 if t==1	// year when unit is first treated
bys i (t): replace Ei = Ei[1]
gen K = t-Ei 								// "relative time", i.e. the number periods since treated (could be missing if never-treated)
gen D = K>=0 & Ei!=. 						// treatment indicator

// Generate the outcome with parallel trends and heterogeneous treatment effects
gen tau = cond(D==1, (t-12.5), 0) 			// heterogeneous treatment effects (in this case vary over calendar periods)
gen eps = rnormal()							// error term
gen Y = i + 3*t + tau*D + eps 				// the outcome (FEs play no role since all methods control for them)
//save five_estimators_data, replace



// Estimation with did_multiplegt of de Chaisemartin and D'Haultfoeuille (2020)


loc dynamic 5
loc placebo 5
loc breps 10
loc seed 1234


did_multiplegt Y i t D, 		robust_dynamic dynamic(`dynamic') placebo(`placebo') breps(`breps') seed(`seed')
graph export "original.png", replace
matrix variancemgt_1=e(didmgt_variances)
matrix variance_1=e(variances)

did_multiplegt_fast Y i t D, 	robust_dynamic dynamic(`dynamic') placebo(`placebo') breps(`breps') seed(`seed') maxprocessors(16)
graph export "fast.png", replace
matrix variancemgt_2=e(didmgt_variances)
matrix variance_2=e(variances)

matrix variance_difference=variance_1-variance_2
matrix variancemgt_difference=variancemgt_1-variancemgt_2

matrix list variance_difference
matrix list variancemgt_difference

matrix list variance_1
matrix list variance_2



