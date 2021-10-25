use "/Users/samharper/Dropbox/work/research/projects/TAMANI project/Data/Impact Eval/dta/hh_female_impact_rpt.dta", clear

by hh1, sort: egen txdel_time = min(cond(txdel == 1, delwave, .))
tab txdel_time
tab hh1 txdel_time, nol
csdid sba_birth txdel , ivar(hh1) time(delwave) gvar(txdel_time)
csdid sba_birth txdel , cluster(hh1) time(delwave) gvar(txdel_time) method(dripw)
estat all
csdid_plot, group(1)
gen g1 = (txdel_time==1)
table g1 delwave if delwave<2 & txdel_time<2, c(mean sba_birth) row col
table g1 delwave if delwave<2, c(mean sba_birth) row col
reg sba_birth g1##delwave if delwave<2, vce(cl hh1)
preserve
collapse (sum) sba = sba_birth (mean) group = txdel_time (count) pop = sba_birth, by(hh1 delwave)
gen p_sba = sba / pop
drop if delwave==.
reshape wide sba group pop p_sba, i(hh1) j(delwave)
reshape long
reshape wide sba pop p_sba, i(hh1 group) j(delwave)
list hh1 group pop* p_sba*
restore
reg sba_birth g1##delwave if delwave<2, vce(cl hh1)
margins g1#delwave
margins r.g1#r.delwave
help concindex
logit sba_birth i.txdel i.hh1 i.delwave, vce(cl hh1)
margins r.txdel
csdid sba_birth txdel , cluster(hh1) time(delwave) gvar(txdel_time) method(dripw)
estat event, window(-2 2)
csdid_plot, style(rarea)
by hh1, sort: egen tx_time = min(cond(txdist == 1, wave, .))
tab hh1 tx_time, nol
logit contracept_f i.txdist i.hh1 i.wave, vce(cl hh1)
margins txdist
margins r.txdist
csdid contracept_f txdist , cluster(hh1) time(wave) gvar(tx_time) method(dripw)
estat all
csdid_plot, group(1)
csdid fecund_marr_f txdist , cluster(hh1) time(wave) gvar(tx_time) method(dripw)
estat all
