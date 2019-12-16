#!/bin/tcsh -m
#foreach iz(`seq 1 1 12`)     #loop (start inc end)
#foreach iz(`seq 13 1 40`)    #no .painteu.gidlist before 13
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
#foreach iz(1 2 3 4 5 `seq 371 1 400`)
#foreach iz(400)
foreach iz(`seq 399 1 400`)


set name = L90Mpc8000_hithres                                                 #the sim
set step = `printf "%05d\n" $iz`                                            #which sim
set tipsydir=/mn/stornext/u3/shens/scratch/Eris_data/              #where the sims are
set iorddir=/mn/stornext/u3/shens/scratch/Eris_data/iordfile/   #where .iord files are
set outdir=/mn/stornext/d17/extragalactic/personal/gabrierg/rprocess_halos/ #where outfile goes


set tipfile=$name.$step                         #the main sim snapshot
set listfile=$name.$step.iord                   #snapshot with extension .iord
set gridfile=$name.$step.painteu.gidlist_idx1p0 #the output from step 1
set grpfile=$name.$step.amiga.grp               #the output of the halo number
set outfile=$name.$step.gord_eu_idx1.3dplot            #this output
#echo $tipfile
#echo $gridfile
#echo $outfile


#set pbsfile = $name.$step.painteu.pbs
#set pbsfile = $name.$step.gord_eu.pbs
#set pbsfile = $name.$step.gord_eu_cumu_200Myr.pbs
#set pbsfile = $name.$step.smooth.pbs
#echo $pbsfile
#qsub $pbsfile
#echo $tipfile


./make_ord_sml_MEG /mn/stornext/u3/shens/scratch/Eris_data/$tipfile /mn/stornext/u3/shens/scratch/Eris_data/iordfile/$listfile /mn/stornext/d17/extragalactic/personal/gabrierg/rprocess_halos/$gridfile /mn/stornext/d17/extragalactic/personal/gabrierg/Eris_AHF/$grpfile /mn/stornext/d17/extragalactic/personal/gabrierg/rprocess_halos/$outfile >> /mn/stornext/d17/extragalactic/personal/gabrierg/rprocess_halos/output.$step.ord_sml.3dplot &

#./make_ord_sml_MEG /mn/stornext/u3/shens/scratch/Eris_data/$tipfile /mn/stornext/u3/shens/scratch/Eris_data/iordfile/$listfile /mn/stornext/d17/extragalactic/personal/gabrierg/rprocess_halos/$gridfile /mn/stornext/d17/extragalactic/personal/gabrierg/Eris_AHF/$grpfile /mn/stornext/d17/extragalactic/personal/gabrierg/rprocess_halos/$outfile

#./make_ord_sml_test /mn/stornext/u3/shens/scratch/Eris_data/$tipfile /mn/stornext/u3/shens/scratch/Eris_data/iordfile/$listfile /mn/stornext/d17/extragalactic/personal/gabrierg/rprocess/$gridfile /mn/stornext/d17/extragalactic/personal/gabrierg/rprocess/$outfile >> /mn/stornext/d17/extragalactic/personal/gabrierg/rprocess_halos/output.$step.eu_out &

end
exit
