*! efa v1.0.2
** Author: Michael W. Gruszczynski, University of Nebraska, Lincoln, NE, USA
**			(mikegruz@huskers.unl.edu)
program define efa, rclass
	version 10.0
	
syntax varlist(max=1 numeric) ,           /*
    */ [INDicator(varlist)                /*
    */ CLEAR                        	/* Required to run
    */ TINTerval(varlist)                 /* Declare time var
    */ WINdows(numlist)		          /* Time window declaration
    */ RAWcount(numlist)	          /* Raw frame count threshold
    */ SALience				  /* Create file of salience scores
    */ FILE(string)			  /* Denote save file names
    */ LOADings				  /*
    */ TAILS			          /* Include tails in analysis
    */ SAVE				  /* Save adjusted factor loadings
    */ LABELS                               /* Save variable labels
    */ DROPzeros                            /* Drop variables with no contrib
    */ COUNT(varlist)                         /* Count variable
    */ TOLerance(numlist max =1 >0 <1) *] /* Factor loading tolerance */
		
		
    if "`tinterval'" != "" {
        local tinterval = "`tinterval'"
    }
		
    else {
        di as err "must declare time var using {cmd:tinterval}"
        exit 198
    }
		
    if "`clear'" == "" {
        di as err "{cmd:evolfact} destroys data; cannot be " /*
        */ "specified without {cmd:clear}"
        exit 198
    }
		
    if "`windows'" != "" {
        local windows = "`windows'"
    }
		
    else {
        local windows = 5
    }
		
    if "`tolerance'" != "" {
        local tolerance = "`tolerance'"
    }
		
    else {
        local tolerance = .90
    }
		
    if "`rawcount'" == "" {
        local rawcount = 5
    }
		
    else {
        local rawcount = "`rawcount'"
    }
		    
    if "`verbose'" != "" {
        local verbose = "noisily"
    }
    
    else {
        local verbose = "quietly"
    }
    
    if "`labels'" != "" {
        quietly levelsof `varlist', local(levels)
        quietly labelbook `varlist'
        local labels = r(names)
        foreach val of local levels {
                local `val' : label `varlist' `val'
                }
    }
    
    else {
        local value = 0
    }
    
    if "`count'" != "" {
        local stub = "`count'"
    }
    
    else {
        local stub = "_"
    }
					
    /* Create indicator vars for each frame type and reshape
     * by year
     */
		
    tempfile contracted salience allcounts countloads
		
    capture {
        
        if "`count'" == "" {
            contract `tinterval' `varlist', freq(`stub')
            drop if `varlist' == .
        }
        
        else {
            /* contract `tinterval' `varlist' `stub', freq(_) */
            drop if `stub' == .
            /* drop _freq */
        }
        
        save "`contracted'"
        keep `varlist'
        contract `varlist'
        drop _freq
        save "`salience'"
        
        use "`contracted'", clear
        reshape wide `stub', i(`tinterval') j(`varlist')
        recode `stub'* (.=0)
        save "`allcounts'"
    }
		
        quietly sum `tinterval'
        local span = `windows' - 1
        local timemin = r(min) + `span'
        local timemax = r(max)
     
    * Dynamic Factor Analyses of Upper Tail Windows
     
    if "`tails'" != "" {
        local tspan = `span' - 1
        local lbound = `timemax' - `tspan'
        forvalues utail = `timemax'(-1)`lbound' {
            local tempspan = `timemax' - `utail'
            if `tempspan' > 1 {
                use "`allcounts'", clear
                quietly keep if `tinterval' >= `utail' & `tinterval' <= `timemax'
                
                foreach var of varlist `stub'* {
                    quietly sum `var'
                    if r(sum) <= `rawcount' {
                        drop `var'
                    }
                }
                
            capture {
                factor `stub'*, pcf
                rotate, varimax
                matrix loadings = e(r_L)
                collapse (sum) `stub'*
                gen y = 1
                reshape long `stub', i(y)
                drop y
                rename _j `varlist'
                rename `stub' count_`utail'_`timemax'
														
                svmat loadings, names(col)
                renpfix Factor salience
					
                foreach var of varlist salience* {
                    replace `var' = abs(`var')
                    replace `var' = 0 if `var' < `tolerance'	
                    replace `var' = `var' * count_`utail'_`timemax'
                    rename `var' `var'_`utail'_`timemax'
                }
					
                svmat loadings, names(col)
                renpfix Factor factor
					
                foreach var of varlist factor* {
                    replace `var' = abs(`var')
                    /* replace `var' = 0 if `var' < `tolerance' */	
                    rename `var' `var'_`utail'_`timemax'
                }
					
                merge 1:1 `varlist' using "`salience'", nogen
                save "`salience'", replace
                
                }
            }
        }
    }
    
    * Dynamic Factor Analyses of t-year windows                    
                                                                
    forvalues timestart = `timemax'(-1)`timemin' {
        use "`allcounts'", clear
        local timeend = `timestart' - `span'
        quietly keep if `tinterval' <= `timestart' & `tinterval' >= `timeend'

        foreach var of varlist `stub'* {
            quietly sum `var'
            if r(sum) <= `rawcount' {
                drop `var'
            }
        }
				
        capture {
            factor `stub'*, pcf
            rotate, varimax
            matrix loadings = e(r_L)
            collapse (sum) `stub'*
            gen y = 1
            reshape long `stub', i(y)
            drop y
            rename _j `varlist'
            rename `stub' count_`timeend'_`timestart'
														
            svmat loadings, names(col)
            renpfix Factor salience
					
            foreach var of varlist salience* {
                replace `var' = abs(`var')
                replace `var' = 0 if `var' < `tolerance'	
                replace `var' = `var' * count_`timeend'_`timestart'
                rename `var' `var'_`timeend'_`timestart'
            }
					
            svmat loadings, names(col)
            renpfix Factor factor
					
            foreach var of varlist factor* {
                replace `var' = abs(`var')
                /* replace `var' = 0 if `var' < `tolerance' */	
                rename `var' `var'_`timeend'_`timestart'
            }
					
            merge 1:1 `varlist' using "`salience'", nogen
            save "`salience'", replace
							
        }			
						
				
    }

    * Dynamic Factor Analyses of Lower Tail Windows

    if "`tails'" != "" {
        
        use "`allcounts'", clear
        
        quietly sum `tinterval'
        local tspan = `windows' - 2
        local timemin = r(min)
        local ubound = `timemin' + `tspan'
        forvalues ltail = `ubound'(-1)`timemin' {
            local tempspan = `ltail' - `timemin'
            if `tempspan' > 1 {
                use "`allcounts'", clear
                quietly keep if `tinterval' >= `timemin' & `tinterval' <= `ltail'
                
                foreach var of varlist `stub'* {
                    quietly sum `var'
                    if r(sum) <= `rawcount' {
                        drop `var'
                    }
                }
                
            capture {
                factor `stub'*, pcf
                rotate, varimax
                matrix loadings = e(r_L)
                collapse (sum) `stub'*
                gen y = 1
                reshape long `stub', i(y)
                drop y
                rename _j `varlist'
                rename `stub' count_`timemin'_`ltail'
														
                svmat loadings, names(col)
                renpfix Factor salience
					
                foreach var of varlist salience* {
                    replace `var' = abs(`var')
                    replace `var' = 0 if `var' < `tolerance'	
                    replace `var' = `var' * count_`timemin'_`ltail'
                    rename `var' `var'_`timemin'_`ltail'
                }
					
                svmat loadings, names(col)
                renpfix Factor factor
					
                foreach var of varlist factor* {
                    replace `var' = abs(`var')
                    /* replace `var' = 0 if `var' < `tolerance' */	
                    rename `var' `var'_`timemin'_`ltail'
                }
					
                merge 1:1 `varlist' using "`salience'", nogen
                                
                save "`salience'", replace
                
                }
            }
        }
    }       
    
    if "`labels'" != "" {         
        foreach l of local levels {
            label define `varlist' `l' "``l''", modify
        }
                    
        label values `varlist' `varlist'
    
        quietly save "`salience'", replace
    }
    
 
    quietly {
        recode factor* (.=0)
        egen _tot = rowtotal(factor*)
        drop if _tot == 0
        drop _tot
        mkmat `varlist', matrix(indicator)
        mkmat factor*, matrix(loadings)
        mkmat salience*, matrix(salience)
        return matrix fo_L loadings
        return matrix sal salience
    }
    
    return local cmd efa
    return local time `tinterval'
    return local var `varlist'
    return local value `value'
    return local labels labels
                                                			
		/*if "`save'" != "" {
			use "`counts'", clear
			save counts_`file'.dta
			use "`salience'", clear
			save salience_`file'.dta
			use "`loadings'", clear
			save loadings_`file'.dta
		}*/
		
end
