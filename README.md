# Thesis
Galactic Chemical Evolution of R-process Elements in Cosmological Simulations

My R-Process run commands:

screen -S [screen name]

to exit: ctrl+a+d

to see avail: screen -ls 

to enter: screen -r [screen number]

------------------------------------------------------------------------
step0

emacs -nw pickstars.pro

idl

.rne pickstars.pro

------------------------------------------------------------------------
step1

emacs -nw paint_rp.c

gcc -o paint_rp paint_rp.c tipsyxdr.c Tipsy.c tipsmall_defs3.h -lm   

emacs -nw runpaint_eu_serial.sh

./runpaint_eu_serial.sh

------------------------------------------------------------------------
step2

emacs -nw make_ord_sml_MEG.c

gcc -o make_ord_sml_MEG make_ord_sml_MEG.c tipsyxdr.c Tipsy.c tipsmall_defs3.h -lm

emacs -nw runmake_ord_sml.sh

./runmake_ord_sml.sh

------------------------------------------------------------------------
step3

emacs -nw read_gord_eu_cumu_MEG.c

gcc -o read_gord_eu_cumu_MEG read_gord_eu_cumu_MEG.c tipsyxdr.c Tipsy.c tipsmall_defs3.h -lm

./read_gord_eu_cumu_MEG

------------------------------------------------------------------------
step4

emacs -nw make_ord_eucumu_MEG.c

gcc -o make_ord_eucumu_MEG make_ord_eucumu_MEG.c tipsyxdr.c Tipsy.c tipsmall_defs3.h -lm

emacs -nw runmake_ord_eucumu.sh

./runmake_ord_eucumu.sh

------------------------------------------------------------------------
Step5

emacs -nw rp_stellar_combine_MEG.pro

idl

.rne rp_stellar_combine_MEG.pro

------------------------------------------------------------------------
step6

emacs -nw eu_fe_MEG.pro

idl

.com lookback

.com rtipsy

.com percentile_16

.com percentile_84

.com binning_ave

.com binning

.com histlike

.com sym

.com legend

.rne eu_fe_MEG.pro


------------------------------------------------------------------------
In Eris_AHF/

emacs -nw writeAHF.sh

./writeAHF.sh


outputs:

tipsyfile.step.input


emacs -nw runahf.sh

./runahf.sh


outputs:

tipsyfile.step.log

tipsyfile.step.out

tipsyfile.step.redshift.#.AHF_halos

tipsyfile.step.redshift.#.AHF_particles

tipsyfile.step.redshift.#.AHF_profiles

tipsyfile.step.redshift.#.AHF_substructure


emacs -nw convertahf.pro

emacs -nw new_grp_stat.pro

idl

.rne convertahf.pro


outputs:

tipsyfile.step.amiga.grp

tipsyfile.step.amiga.gtp

tipsyfile.step.amiga.stat
