#!/bin/tcsh -m                                                                       
#foreach iz(`seq 1 1 10`)        #loop (start inc end)
#foreach iz(`seq 13 1 20`)       #no .gord_eu before 13
#foreach iz(`seq 21 1 30`) 
#foreach iz(`seq 31 1 40`) 
#foreach iz(`seq 41 1 50`) 
#foreach iz(`seq 51 1 60`) 
#foreach iz(`seq 61 1 70`) 
#foreach iz(`seq 71 1 100`) 
#foreach iz(`seq 81 1 90`) 
#foreach iz(`seq 91 1 100`)                                                         
#foreach iz(`seq 101 1 110`)                                               
#foreach iz(`seq 111 1 120`)                                                      
#foreach iz(`seq 121 1 130`)                                                        
#foreach iz(`seq 131 1 160`)                                                         
#foreach iz(`seq 141 1 150`)                                                        
#foreach iz(`seq 151 1 160`)                                               
#foreach iz(`seq 161 1 170`)                                
#foreach iz(`seq 171 1 180`)                                                    
#foreach iz(`seq 181 1 190`)
#foreach iz(`seq 191 1 220`)                                                        
#foreach iz(`seq 201 1 210`)                                                         
#foreach iz(`seq 211 1 220`)                                                        
#foreach iz(`seq 221 1 230`)                                  
#foreach iz(`seq 231 1 240`)                                                         
#foreach iz(`seq 241 1 250`)                                                      
#foreach iz(`seq 251 1 280`)                                                        
#foreach iz(`seq 261 1 270`)                                                         
#foreach iz(`seq 271 1 280`)                                                        
#foreach iz(`seq 281 1 310`)                
#foreach iz(`seq 291 1 300`)                                                       
#foreach iz(`seq 301 1 310`)                                            
#foreach iz(`seq 311 1 320`)                                 
#foreach iz(`seq 321 1 330`)                                                         
#foreach iz(`seq 331 1 340`)                                                       
#foreach iz(`seq 341 1 350`)                                                       
#foreach iz(`seq 351 1 360`)                                                     
#foreach iz(`seq 361 1 370`)                                                        
#foreach iz(`seq 361 1 380`)                                                     
#foreach iz(`seq 381 1 390`)                                                        
foreach iz(`seq 361 1 400`)


set name = L90Mpc8000_hithres                                         #the sim
set step = `printf "%05d\n" $iz`                                      #which sim
set tipsydir=/mn/stornext/u3/shens/scratch/Eris_data/             #where the sims are
set iorddir=/mn/stornext/u3/shens/scratch/Eris_data/iordfile/   #where .iord files are
set outdir=/mn/stornext/d17/extragalactic/personal/gabrierg/rprocess/ #outfile paths


set tipfile=$name.$step                          #the main sim snapshot              
set gordeufile=$name.$step.gord_eu_idx1          #step2 output                       
set euonlyfile=$name.$step.gord_euonly_cumu_idx1 #step3 output                       
set listfile=$name.$step.iord                    #snapshot with extension .iord      
set outfile=$name.$step.stellar_rp               #output file                       


./make_ord_eucumu_MEG /mn/stornext/u3/shens/scratch/Eris_data/$tipfile /mn/stornext/d17/extragalactic/personal/gabrierg/rprocess_halos/$gordeufile /mn/stornext/d17/extragalactic/personal/gabrierg/rprocess_halos/$euonlyfile /mn/stornext/u3/shens/scratch/Eris_data/iordfile/$listfile /mn/stornext/d17/extragalactic/personal/gabrierg/rprocess_halos/$outfile >> /mn/stornext/d17/extragalactic/personal/gabrierg/rprocess_halos/output.$step.stellar_rp &

end
exit
