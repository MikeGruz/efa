*! efaho v0.0.1
** Author: Michael W Gruszczynski, University of Nebraska, Lincoln, NE, USA
**                      (mikegruz@huskers.unl.edu)
program define efaho, rclass
        version 10.0
    
syntax [,                /*
    */ KEEP                                   /* Keep 2nd-order loadings
    */ CLEAR                                  /* Required to run
    */ LABELS                                 /* Keep value labels
    */ PERsist                                /* Save ef persistence
    */ TOLerance(numlist max=1 >0 <1)     *] /* 2nd-order loading tolerance */
    
    if "`clear'" == "" {
        di as err "{cmd:efaho} destroys data; cannot be specified " /*
        */ "without {cmd:clear}"
        exit 198
    }
                
    if "`r(cmd)'" != "efa" {
        di as err "{cmd:efaho} cannot run without first running {cmd:efa}"
        exit 198
    }
    
    if "`tolerance'" != "" {
        local tolerance = "`tolerance'"
        }
        
    else {
        local tolerance = .50
    }
    
    tempfile marker salience factor
    local var = r(var)
        
    matrix load = r(fo_L)
    matrix salience = r(sal)
    keep `var' salience*
    save "`salience'", replace
   
    clear
    
        svmat load, names(col)
        local time = r(time)
        
        factor factor*, pcf
        rotate, varimax
        matrix load2 = e(r_L)
        xpose, clear varname
        drop v*
        svmat load2, names(col)
        
        foreach x of varlist Factor* {
            replace `x' = 0 if `x' < `tolerance'
            /* count if `x' != 0
            if r(N) < 2 {
                drop `x'
                } */
            }
        
        local nfactor = 0               
        foreach x of varlist Factor* {
            replace `x' = 1 if `x' != 0
            local nfactor = `nfactor' + 1
            }
        
        save "`factor'", replace    
                    
        foreach x of varlist Factor* { 
            use "`factor'", clear
            local n = _N
            local list = ""
            
            forvalues obs = 1/`n'  {
                qui sum `x' in `obs'
                if r(sum) != 0 {
                    local a = _varname in `obs'
                    tokenize `a', parse("r""_")
                    local a = "salience`3'`4'`5'`6'`7' "
                    local list "`list' `a'"
                    }
            }  
                          
            use "`salience'", clear
            
            if "`labels'" != "" {
                levelsof `var', local(levels)
            }
            
            keep `var' `list'
            
            if missing("`list'") == 0 { 
                egen tot = rowtotal(salience*)
            
                
            drop if tot == 0
            drop tot
            
            if "`labels'" != "" {
                foreach val of local levels {
                    local `val' : label `var' `val'
                } 
            }
            
            if `=_N' > 0 {
            xpose, clear varname
            
            
            if "`labels'" != "" {
                foreach v of varlist v* {
                    local lab = `v' in 1
                    label var `v' "``lab''"
                }
            }
            
            drop in 1
            recode v* (0=.)
            
            split _varname, parse(_)
            
            gen `time' = _n
            local obs = _N
            
            forval ob = 1/`obs'  {
                local lb = _varname2 in `ob'
                local ub = _varname3 in `ob'
                label define `time' `ob' "`lb'-`ub'", modify
            }
            
            label values `time' `time'
            drop _varname* 
            
            
                save `x', replace
            }
            }
        
        }
        
        use "`factor'", clear
        
        mkmat Factor*, matrix(marker)
        drop Factor*
        save "`marker'", replace
        clear       
        
        svmat salience, names(col)
        
        foreach x of varlist salience* {
            sum `x'
            replace `x' = r(sum)
            }
        
        keep in 1/`nfactor'
        mkmat salience*, matrix(_sal2)
        matrix _sal2 = _sal2'
        drop salience*
        
        use "`marker'", clear
        svmat marker, names(col)
        svmat _sal2, names(col)

        forvalues f = 1/`nfactor' {
            replace Factor`f' = Factor`f'*r`f'
            drop r`f'
            }
            
        split _varname, parse(_)
        destring _varname2 _varname3, replace
        
        if "`persist'" != "" {
            svmat load2
            renpfix load2 persist
            foreach x of varlist persist* {
                replace `x' = 0 if `x' < `tolerance'
            }
            
            collapse (sum) Factor* persist*, by(_varname2 _varname3)
            recode Factor* persist* (0=.)
        }
        
        else {
            collapse (sum) Factor*, by(_varname2 _varname3)
            recode Factor* (0=.)
        }
        
        gen `time' = _n
        
        local n = _N
        
        forvalues ob = 1/`n'  {
            local lb = _varname2 in `ob'
            local ub = _varname3 in `ob'
            label define `time' `ob' "`lb'-`ub'", modify
            }
        
        drop _varname2 _varname3
        label values `time' `time'
        renpfix Factor ef

    
    
end
