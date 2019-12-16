#!/bin/tcsh -m 
#foreach iz(`seq 1 1 10`)    #loop (start inc end)
#foreach iz(`seq 11 1 40`)
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
foreach iz(`seq 371 1 400`)
#foreach iz(370)

set name = L90Mpc8000_hithres                                   #the sim
set step = `printf "%05d\n" $iz`                                #which sim
set tipsydir=/mn/stornext/u3/shens/scratch/Eris_data/           #where the sims are
set outdir=/mn/stornext/d17/extragalactic/personal/gabrierg/rprocess/ #where fiducial outfile goes
set outdir=/mn/stornext/d17/extragalactic/personal/gabrierg/rprocess_idx2p0/ #
set outdir=/mn/stornext/d17/extragalactic/personal/gabrierg/rprocess_halos/ #
#change paths below


set tipfile=$name.$step
set listfile=$name.$step.sid_merger_0p05 #3 merger files check which yeild from paint_rp.c
set outfile=$name.$step.painteu.gidlist_idx1p0
set logfile=$name.$step.painteu.gidlist_idx1p0.log


#set pbsfile = $name.$step.painteu.pbs
#set pbsfile = $name.$step.gord_eu.pbs
#set pbsfile = $name.$step.gord_eu_cumu_200Myr.pbs
#set pbsfile = $name.$step.smooth.pbs
#echo $pbsfile
#qsub $pbsfile
#echo $tipfile


./paint_rp /mn/stornext/u3/shens/scratch/Eris_data/$tipfile /mn/stornext/d17/extragalactic/personal/gabrierg/rprocess/$listfile /mn/stornext/d17/extragalactic/personal/gabrierg/rprocess_halos/$outfile /mn/stornext/d17/extragalactic/personal/gabrierg/rprocess_halos/$logfile >> /mn/stornext/d17/extragalactic/personal/gabrierg/rprocess_halos/output.$step.sid_merger_0p05 &


end 
exit 
    
