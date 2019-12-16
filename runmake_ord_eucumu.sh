#!/bin/tcsh -m
#foreach iz(`seq 1 1 10`)        #loop (start inc end)
#foreach iz(`seq 13 1 40`)       #no .gord_eu before 13
#foreach iz(`seq 41 1 70`)
#foreach iz(`seq 71 1 100`)
#foreach iz(`seq 101 1 130`)
#foreach iz(`seq 131 1 160`)
#foreach iz(`seq 161 1 190`)
#foreach iz(`seq 191 1 220`)
#foreach iz(`seq 221 1 250`)
#foreach iz(`seq 251 1 280`)
#foreach iz(`seq 281 1 310`)
#foreach iz(`seq 311 1 340`)
#foreach iz(`seq 341 1 370`)
#foreach iz(`seq 371 1 400`)

#foreach iz(`seq 399 1 400`)
foreach iz(399)

set name = L90Mpc8000_hithres                                         #the sim
set step = `printf "%05d\n" $iz`                                      #which sim
set tipsydir=/mn/stornext/u3/shens/scratch/Eris_data/                 #where the sims are
set iorddir=/mn/stornext/u3/shens/scratch/Eris_data/iordfile/         #where .iord files are
set outdir=/mn/stornext/d17/extragalactic/personal/gabrierg/rprocess/ #where outfile goes


set tipfile=$name.$step                          #the main sim snapshot
set gordeufile=$name.$step.gord_eu_idx1          #step2 output
set euonlyfile=$name.$step.gord_euonly_cumu_idx1 #step3 output
set listfile=$name.$step.iord                    #snapshot with extension .iord
#set outfile=$name.$step.stellar_rp               #output file
set outfile=$name.$step.Xi_star


./make_ord_eucumu_MEG /mn/stornext/u3/shens/scratch/Eris_data/$tipfile /mn/stornext/d17/extragalactic/personal/gabrierg/rprocess_halos/$gordeufile /mn/stornext/d17/extragalactic/personal/gabrierg/rprocess_halos/$euonlyfile /mn/stornext/u3/shens/scratch/Eris_data/iordfile/$listfile /mn/stornext/d17/extragalactic/personal/gabrierg/rprocess_halos/$outfile >> /mn/stornext/d17/extragalactic/personal/gabrierg/rprocess_halos/output.$step.stellar_rp.Xi_star &


#./make_ord_eucumu /mn/stornext/u3/shens/scratch/Eris_data/$tipfile /mn/stornext/d17/extragalactic/personal/gabrierg/rprocess/$gordeufile /mn/stornext/d17/extragalactic/personal/gabrierg/rprocess/$euonlyfile /mn/stornext/u3/shens/scratch/Eris_data/iordfile/$listfile /mn/stornext/d17/extragalactic/personal/gabrierg/rprocess/$outfile >> /mn/stornext/d17/extragalactic/personal/gabrierg/rprocess/output.$step.stellar_rp &


#./step4copy /mn/stornext/u3/shens/scratch/Eris_data/$tipfile /mn/stornext/d17/extragalactic/personal/gabrierg/rprocess_halos/$gordeufile /mn/stornext/d17/extragalactic/personal/gabrierg/rprocess_halos/$euonlyfile /mn/stornext/u3/shens/scratch/Eris_data/iordfile/$listfile /mn/stornext/d17/extragalactic/personal/gabrierg/rprocess_halos/$outfile >> /mn/stornext/d17/extragalactic/personal/gabrierg/rprocess_halos/output.$step.Xi_star &

end 
exit 
