cd "F:/replication_repository/"
import delimited using "additional_stata_code/new_df_0p5.csv", clear

keep far lotarea
gen flsp = far * lotarea
summarize far lotarea flsp

gen d = 0.04
gen fbar = 0.5

*LABELS OBSERVATIONS AS BEING IN minus, plus, OR b (NEAR FBAR) RANGES

gen lbl = "na"
replace lbl = "minus" if float(far) <= float(2 - d) & float(far) > float(2 - 2*d)
replace lbl = "plus" if float(far) >=  float(2 + d) & float(far) < float(2 + 2*d)
replace lbl = "b" if float(far) > float(2 - d) & float(far) < float(2 + d)

*GENERATES DUMMY INDICATORS FOR THE THREE RANGES

gen Hminus = 0
replace Hminus = 1 if lbl == "minus"
gen Hplus = 0
replace Hplus = 1 if lbl == "plus"

*USES COLLAPSE COMMAND TO SUM UP OBERVATIONS IN THREE RANGES TO GENERATE AREAS 

*Hminus, Hplus

collapse (sum) Hminus Hplus (mean) far flsp lotarea

gen K = lotarea * far

*GENERATES DENSITY HEIGHTS USING AREAS

gen fbar = 0.5 
gen d = 0.04
gen th = 1.117
gen hminus = Hminus/d
gen hplus = Hplus/d

gen m = 0.5 * ((th - 1) * fbar)^2 *((1 + th)/th)
display m

gen M = 0.5 * ((th - 1) * fbar)^2 *((1 + th)/th) * hplus * lotarea
display M

