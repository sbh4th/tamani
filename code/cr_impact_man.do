capture log close
log using "logs/manuscript_an.txt", text replace

//  program:  cr_manuscript_dta.do
//  task:     creates a stripped down clean dataset for analysis and open source
// 	input:    hh_female_impact_rpt.dta
//  output:   hh_female_impact_man.dta
//  project:  TAMANI
//  author:   Erin Hetherington \ 2021-09-24


//  #0
//  program setup

set linesize 80
clear all
macro drop _all
version 16


//   #1
//   Data cleaning

use "../dta/hh_female_impact_rpt.dta", clear

* create categories for age
egen age3 = cut(wi2), at(15 20 30 50)
gen age_1519 = (age3==15) if age3!=.
gen age_2029 = (age3==20) if age3!=.
gen age_3049 = (age3==30) if age3!=.
label var age3 "age 3 categories"

*create categories for education
gen educ3=eduatt
recode educ3 1=0 2=1 3/5=2
label define educ3 0"less than primary" 1"completed primary" 2"more than primary"
label values educ3 educ3
label var educ3 "Educational Attainment (3 cat)"
gen less_primary = (eduatt==1) if eduatt!=.
gen primary = (eduatt==2) if eduatt!=.
gen more_primary = (eduatt>=3 & eduatt<=5) if eduatt!=.

* ever given birth
label var birth "Have you ever given birth?"
label define noyes 0"no" 1"yes"
label values birth noyes

*married
label values marr_f noyes

*religion
gen christian = (hh10==1) if hh10!=. 
gen islam = (hh10==2) if hh10!=. 
gen other_religion = (hh10==6 | hh10==7) if hh10!=.
rename hh10 religion
recode religion 6/7=3
label define religion 1"Christianity" 2"Islam" 3"other/no religion"
label values religion religion 

*ethnicity
gen nyamwezi = (hh12==1) if hh12!=. 
gen sukuma = (hh12==2) if hh12!=. 
gen waha = (hh12==3) if hh12!=. 
gen other_ethnic = (hh12==96) if hh12!=.
rename hh12 ethnic
recode ethnic 96=4
label define ethnic 1"Nyamwezi" 2"Sukuma" 3"Waha" 4"Other"
label values ethnic ethnic

*wealth
gen livestock = (ha6==1) if ha6!=.
gen land = (ha8==1) if ha8!=.


*reorder districts for easier interpretation
gen district = hh1
label var district "District"
recode district 6=2 2=3 8=4 7=6 3=7 4=8
label define district_reordered 1"Kaliua DC" 2"Urambo DC" 3 "Nzega DC" 4 "Igunga DC" 5"Tabora MC" 6 "Uyui DC" 7"Nzega TC" 8"Sikonge DC"
label values district district_reordered

*label wave
label var wave "Wave of data collection"

* assign value to time of switching to treatment
gen txtime=.
label var txtime "Time of intervention"
replace txtime = 1 if district == 1 | district == 2
replace txtime = 2 if district == 3 | district == 4
replace txtime = 3 if district == 5 | district == 6
replace txtime = 4 if district == 7 | district == 8


//   #2
//   saving stripped dataset

keep sba_birth contracept_f anc4 ///
district txtime txdist txpreg  txdel wave delwave pregwave ///
age3 age_1519 age_2029 age_3049 ///
marr_f birth livestock land ///
educ3 less_primary primary more_primary ///
religion christian islam other_religion ///
ethnic nyamwezi sukuma waha other_ethnic  

save "..\dta\hh_female_impact_man.dta", replace

log close