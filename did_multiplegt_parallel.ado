******************************************************************************
*Define the program for looop
cap program drop did_multiplegt_parallel
program did_multiplegt_parallel
	version 12.0


syntax  [if] [in] , varlist(string) [THRESHOLD_stable_treatment(real 0) trends_nonparam(string) trends_lin(string) controls(string) counter(string) placebo(integer 0) dynamic(integer 0) breps(integer 50) cluster(string) bootstrap_rep(integer 0) switchers(string) LONGdiff_placebo FIRSTdiff_placebo count_switchers_tot count_switchers_contr robust_dynamic discount(real 1) drop_larger_lower always_trends_nonparam always_trends_lin seed(integer 0) paralleloffset(integer 100) ] maxprocessors(integer) filename(string)
	

loc i=$task_id

sleep `=mod(`i',`maxprocessors')*`paralleloffset''

scalar too_many_controls=0
 
*END NEWW

loc workers=`maxprocessors'
loc all=`breps'
di as red "Worker" `i'
loc even=floor(`all'/`workers')
di as result "Even work is " `even'

loc residual=mod(`all',`workers')
di "Residual is "  `residual'

loc assigned=`even'+(`residual'>=`i')
di "Assigned is "  `assigned'

loc initial=(`i'-1)*`even'+1
di "Initial wihout adjustment is "`initial'
di "Adjustment is "cond(`residual'>=`i'-1,`i'-1,`residual')
loc initial=`initial'+cond(`residual'>=`i'-1,`i'-1,`residual')
di "Initial after adjustment "`initial'

loc final=`initial'+`assigned'-1
di "Final is "`final'

us "`filename'", clear
**Solve an issue with a previous temp variable
rename __* ___*
loc counter="_`counter'"
di as red "`seed'"

forvalues y=`initial'/`final'{
	di as red "seed is `seed', iteration seed is `=`seed'+`y''"
*NEWW
	preserve
		///////////////
		set rng mt64
		set seed `=`seed'+`y''
		*g seed=`=`seed'+`y''
		*g cluster="`cluster'"
		*g filename="`filename'"
		*g random="`c(seed)'"
		*g current="`c(rng_current)'"

		bsample, cluster(`cluster')
		*sa sample_parallel_`y', replace

		//Indicate that program will run bootstrap replications

		local bootstrap_rep=1


		did_multiplegt_estim2 `varlist', threshold_stable_treatment(`threshold_stable_treatment') trends_nonparam(`trends_nonparam') trends_lin(`trends_lin') controls(`controls') counter(`counter') placebo(`placebo') dynamic(`dynamic') breps(`breps') cluster(`cluster') bootstrap_rep(`bootstrap_rep') switchers(`switchers') `longdiff_placebo' `firstdiff_placebo' `count_switchers_tot' `count_switchers_contr' `robust_dynamic' discount(`discount') `drop_larger_lower' `always_trends_nonparam' `always_trends_lin'


		// Put results into a matrix 

		matrix didmgt_bootstrap_`y'=effect_0_2
		forvalue j=1/`dynamic'{
		matrix didmgt_bootstrap_`y'=didmgt_bootstrap_`y',effect_`j'_2
		}
		forvalue j=1/`placebo'{
		matrix didmgt_bootstrap_`y'=didmgt_bootstrap_`y',placebo_`j'_2
		}

		sleep `=mod(`i',`maxprocessors')*`paralleloffset''

		clear

		matsave didmgt_bootstrap_`y', p("`c(tmpdir)'") replace
	restore
}
//end of forvalues	
	
end



