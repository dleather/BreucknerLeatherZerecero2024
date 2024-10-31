import delimited using "C:\Users\jkbrueck\Documents\STATA\Bunching\csv files\new_df_2p0.csv", clear

keep far lotarea
gen flsp = far * lotarea
summarize far lotarea flsp

gen d = .15
gen fbar = 2.0

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

gen fbar = 2 
gen d = 0.15
gen th = 1.095
gen hminus = Hminus/d
gen hplus = Hplus/d

gen q = 0.5 * ((th - 1) * fbar)^2
display q

gen m = 0.5 * ((th - 1) * fbar)^2 *((1 + th)/th)
display m

gen M = 0.5 * ((th - 1) * fbar)^2 *((1 + th)/th) * hplus * lotarea
display M

