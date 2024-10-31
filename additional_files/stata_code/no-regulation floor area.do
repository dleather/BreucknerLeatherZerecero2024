cd "C:\Users\davle\Dropbox (Personal)\Bunching\replication_repository\"


import delimited using "additional_stata_code\new_df_1p25.csv", clear

*DIFFERENT CSV FILE NEEDED FOR EACH FBAR GROUP

*VALUE OF d AND fbar MUST BE ADJUSTED FOR EACH GROUP, AND theta VALUE SPECIFIED

gen d = .1
gen fbar = 1.25
gen flrarea_now = .
replace flrarea_now = far * lotarea if far < fbar - d
replace flrarea_now = far * lotarea if far > fbar + d
replace flrarea_now = fbar * lotarea if far >= fbar - d & far <= fbar + d
collapse (sum) flrarea_now
display flrarea_now

import delimited using "additional_stata_code\new_df_1p25.csv", clear

gen d = .1
gen fbar = 1.25
gen theta = 1.03
gen flrarea_new = .
replace flrarea_new = far * lotarea if far < fbar - d
replace flrarea_new = theta * far * lotarea if far > fbar + d
replace flrarea_new = 0.5* (1 + theta) * fbar * lotarea if far >= fbar - d & far <= fbar + d
collapse (sum) flrarea_new
display flrarea_new

