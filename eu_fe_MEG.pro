;pro eu_fe_MEG, tipfile, grpfile, eufil e, oxyfile, fefile, dyrunit, ismooth, outfile, output, plottops,filename1, filename2

dir = '/mn/stornext/u3/shens/scratch/Eris_data/'
dir1 = '/mn/stornext/d17/extragalactic/personal/gabrierg/rprocess_halos/'
;dir2 = '/mn/stornext/u3/shens/scratch/Eris_amiga/'
dir2 = '/mn/stornext/d17/extragalactic/personal/gabrierg/Eris_AHF/' 
;dir3 = '/mn/stornext/u3/shens/scratch/rprocess/'
;dir4 = '/mn/stornext/d17/extragalactic/personal/gabrierg/Eris_AHF/'
name = 'L90Mpc8000_hithres.'
step = [400]
tipfile = dir + name + string(step[0], format='(i5.5)')
grpfile = dir2 + name + string(step[0], format='(i5.5)')  + '.amiga.grp' ;Eris_amiga
;oxyfile = tipfile + '.OxMassFrac'
;fefile = tipfile + '.FeMassFrac'
posCenter = dir2 + name + string(step[0], format='(i5.5)')  + '.amiga.stt'
oxyfile = dir1 + name + 'oxy_smooth_128'      ; in rprocess_halos dir
fefile = dir1 + name + 'fe_smooth_128'        ; in rprocess_halos dir
eufile = dir1 + name + 'seu_sm_0p01_128'      ; step5 output smoothing
ns_eufile = dir1 + name + 'seu_nosm_0p01_128' ; step5 output no smoothing
dwarfile = dir1 + name + 'seu_dwarf_0p01_128' ; step5 output with DWARF FRACTION
dir_dwarfile = dir1 + name + 'seu_dir_dwarf_0p01_128' ; step5 output with DIRECT DWARF FRACTION 

;dsecunit = 1.2236e+18
;dyrunit = dsecunit/3600./24./365.
dyrunit = 38767.4e6

outfile = dir1 + 'Eris_host_abundance_smooth_idx1_350pc.dat'
filename1 = "1Age_dist_lowZ_smooth.eps"
filename2 = "2Age_FeH_smooth_128.eps"
filename3 = "3Age_EuFe_host_idx1.eps"
filename4 = "4FeH_EuFe_host_idx1.eps"
filename5 = "5FeH_EuFe_dir_dhost_idx1.eps"
filename6 = "6FeH_EuFe_dhost_idx1.eps"
filename7 = "7FeH_EuFe_sat_idx1.eps"
filename8 = "8FeH_EuFe_dsat_idx1.eps"
filename9 = "9FeH_EuFe_dir_dsat_idx1.eps"
filename10 = "10Age_EuH_host_idx1.eps"
filename11 = "11Age_EuH_dhost_idx1.eps"
filename12 = "12Age_EuH_dir_dhost_idx1.eps"
filename13 = "13Age_EuH_sat_idx1.eps"
filename14 = "14Age_EuH_dsat_idx1.eps"
filename15 = "15Age_EuH_dir_dsat_idx1.eps"
filename16 = "16Age_FeH_sat_128.eps"
filename17 = "17Age_dwarfrac_128.eps"
filename18 = "18Age_seu_128.eps"
filename19 = "19Age_seu_dwarf_128.eps"
filename20 = "20Age_dwarfrac_sat_128.eps"
filename21 = "21Age_seu_sat_128.eps"
filename22 = "22Age_seu_dwarf_sat_128.eps"
filename23 = "23FeH_dwarfrac_128.eps"
filename24 = "24FeH_seu_128.eps"
filename25 = "25FeH_seu_dwarf_128.eps"
filename26 = "26FeH_dwarfrac_sat_128.eps"
filename27 = "27FeH_seu_sat_128.eps"
filename28 = "28FeH_seu_dwarf_sat_128.eps"
filename29 = "29Dwarfrac_EuH_idx1.eps"
filename30 = "30Dwarfrac_EuH_dhost_idx1.eps"
;2P
filename31 = "31Dwarfrac_EuH_sat_idx1.eps"
filename32 = "32Dwarfrac_EuH_dsat_idx1.eps"
filename33 = "33Dir_Dwarfrac_EuH_idx1.eps"
filename34 = "34Dir_Dwarfrac_EuH_dhost_idx1.eps"
filename35 = "35Dir_Dwarfrac_EuH_sat_idx1.eps"
filename36 = "36Dir_Dwarfrac_EuH_dsat_idx1.eps"
filename37 = "37Age_EuFe_dhost_idx1.eps"
filename38 = "38Age_EuFe_dir_dhost_idx1.eps"
filename39 = "39Age_EuFe_sat_idx1.eps"
filename40 = "40Age_EuFe_dsat_idx1.eps"
filename41 = "41Age_EuFe_dir_dsat_idx1.eps"
;3P
filename42 = "42gXY2000.eps"
filename43 = "43gXY1000.eps"
filename44 = "44gXY500.eps"
filename45 = "45gXY200.eps"
filename46 = "46gXY100.eps"
filename47 = "47gXZ2000.eps"
filename48 = "48gXZ1000.eps"
filename49 = "49gXZ500.eps"
filename50 = "50gXZ200.eps"
filename51 = "51gXZ100.eps"
filename52 = "52sXY500.eps"
filename53 = "53sXY200.eps"
filename54 = "54sXY100.eps"
filename55 = "55sXZ500.eps"
filename56 = "56sXZ200.eps"
filename57 = "57sXZ100.eps"
;4P
filename58 = "58EuFe_EuH_host_idx1.eps"
filename59 = "59EuFe_EuH_dhost_idx1.eps"
filename60 = "60EuFe_EuH_sat_idx1.eps"
filename61 = "61EuFe_EuH_dsat_idx1.eps"
;retall
;--------------------------------------------------------------------;
;                                BLOCK1                              ;
print, "                         block1"
ismooth = 1
output = 1
;plottops = 0 ;plot
plottops = 1 ;save
ifloor = 1
icheckh = 0
;goto, plot
;goto, plot2
;goto, plot3
;goto, plot4
;goto, plot5
;goto, plot6
if (icheckh eq 1) then begin 
   file = dir1 + name + 'hsmooth_combine'
   close, 1
   nstot = 0L
   openr, 1, file
   readf, 1, nstot 
   print, 'nstot' , nstot
   timestep = intarr(nstot)
   hsm = fltarr(nstot)
   readf, 1, timestep 
   readf, 1, hsm 
   close, 1
endif 
;goto, plot

if (ismooth eq 1) then begin 
   label=['NSNS merger only + mixing 128']
endif else begin
   label=['NSNS merger only + no mixing']
endelse

print, "                         block1.1"
FeH_solar = 7.45-12
OH_solar = 8.66-12
EuH_solar = 0.52 -12 

Ofe_solar = OH_solar - FeH_solar
Eufe_solar = EuH_solar - FeH_solar 
EuO_solar = EuH_solar - OH_solar 

ageofuniverse = lookback(1000000) ;in yrs 
Ah  = 1.00794
Ao  = 15.9994
Afe = 55.845
Aeu = 151.964 
meu_mo_solar = 10^(EuO_solar)*Aeu/Ao 

mfemh_solar = 10^FeH_solar*Afe/Ah
momh_early = 10^(OH_solar+0.4)*Ao/Ah
fe_floor = mfemh_solar * 1e-4 
O_floor = momh_early *1e-4 

rtipsy, tipfile, h, g, d, s
redshift = (1-h.time)/h.time

;--------------------------------------------------------------------;
;                                BLOCK2                              ;
print, "                         block2"
print, "read iron and oxygen"
if (ismooth eq 0) then begin 
   fetot = rbarray(fefile)
   oxytot = rbarray(oxyfile)
   
   gfe = fetot(0:h.nsph-1)
   sfe1 = fetot(h.nsph+h.ndark:h.nbodies-1)
   goxy = oxytot(0:h.nsph-1)
   soxy1 = oxytot(h.nsph+h.ndark:h.nbodies-1)
endif else begin 
   ns = h.nstar
   close, 1
   openr, 1, oxyfile
   ns1 = 0l
   readf, 1, ns1
   print, "ns1, ns", ns1, ns
   soxy1 = fltarr(ns1)
   readf, 1, soxy1
   close, 1
   
   openr, 1, fefile
   ns1 = 0l
   readf, 1, ns1
   print, "ns1, ns", ns1, ns
   sfe1 = fltarr(ns1)
   readf, 1, sfe1
   close, 1
endelse 

if (ifloor eq 1) then begin 
   sfe = sfe1 + fe_floor
   soxy = soxy1 + o_floor
endif else begin 
   sfe = sfe1
   soxy = soxy1 
endelse 

ns = h.nstar

print, "                         block2.1"
print, "Read in eu file"
close, 1
openr, 1, eufile
ns_eu = 0l
readf, 1, ns_eu
print, "ns_eu, ns", ns_eu, ns
seu = fltarr(ns_eu)
seu_0 = fltarr(ns_eu)
seu_per = fltarr(ns_eu)
readf, 1, seu
close, 1

index = where(seu ne 0)
print, "fraction of eu enriched", n_elements(index)*1.0/ns

seu_SN = dblarr(ns)
seu_SN = soxy*meu_mo_solar
seu_tot = seu + seu_SN

print, "                         block2.2"
print, "Read in group file"
close, 1
ntot = 0L
openr, 1, grpfile 
readf, 1, ntot 
grpid = lonarr(ntot)
readf, 1, grpid
close, 1
sgrp = grpid(h.nsph+h.ndark:h.nbodies-1)
help, sgrp
help, hsm

ihost = where(sgrp eq 1) 
if (icheckh eq 1) then begin 
   ihost = where(sgrp eq 1 and hsm le 0.35)
endif 

print, "number of stars in host", n_elements(ihost)
nhost =  n_elements(ihost)

;age = ageofuniverse -s.tform * dyrunit
;ihost = where(sgrp eq 1 and age/1e9 lt 11.)
;print, "number of stars in host", n_elements(ihost)
;nhost =  n_elements(ihost)

print, "                         block2.3"
;HOST PARAMETERS
fh = 0.764
age_host =  ageofuniverse -s(ihost).tform * dyrunit
tform_host =  s(ihost).tform * dyrunit
sFe_host = s(ihost).mass*sfe(ihost)
sOxy_host = s(ihost).mass*soxy(ihost)
seu_host = s(ihost).mass*seu(ihost) 
seu_host_tot = s(ihost).mass*seu_tot(ihost)
smass_host = s(ihost).mass

FeH_host = alog10(sFe_host/(smass_host*fh) *(Ah/Afe))-FeH_solar
OFe_host = alog10(soxy_host/sfe_host *(Afe/Ao))-OFe_solar

index = where(OFe_host gt 0)
nindex = n_elements(index) 
print, "nhost, nindex", nhost, nindex

EuFe_host = alog10(seu_host/sfe_host*(Afe/Aeu))-Eufe_solar
EuFe_host_tot = alog10(seu_host_tot/sfe_host*(Afe/Aeu))-Eufe_solar
EuH_host = alog10(seu_host/(smass_host*fh)*(Ah/Aeu)) - EuH_solar

nFenH = sFe_host/(smass_host*fh) *(Ah/Afe)
nOxynH = soxy_host/(smass_host*fh) *(Ah/Ao)
nEunH = seu_host/(smass_host*fh) *(Ah/Aeu)
nEunH_tot = seu_host_tot/(smass_host*fh) *(Ah/Aeu)

;--------------------------------------------------------------------;
;^^^^^^^^^^^^^^^^^^^^^^^^^^^^MY ADDITIONS^^^^^^^^^^^^^^^^^^^^^^^^^^^^;
print, "                 new addition block I"
print, "Read in non-smoothing eu file"
close, 1
openr, 1, ns_eufile
ns_eu = 0l
readf, 1, ns_eu_nosm
print, "ns_nseu, ns", ns_eu_nosm, ns
seu_nosm = fltarr(ns_eu_nosm)
readf, 1, seu_nosm
close, 1

index_nosm = where(seu_nosm ne 0)
print, "fraction of non-smooth eu enriched", n_elements(index_nosm)*1.0/ns

;^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^;
print, "Read in dwarf file"
close, 1
openr, 1, dwarfile
dwarf = 0l
readf, 1, dwarf
print, "dwarf, ns", dwarf, ns
dwarfeu = fltarr(dwarf) ;array of all dwarf eu values               
readf, 1, dwarfeu
close, 1
dwarfindex = where(dwarfeu ne 0) ;use index for something
print, "fraction of eu enrichment from dwarfs",n_elements(dwarfindex)*1.0/ns ;number of non-zero elements

;retall 
dwarfeu_SN = dblarr(ns)
dwarfeu_SN = soxy*meu_mo_solar
dwarfeu_tot = dwarfeu + dwarfeu_SN

;^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^;     
print, "Read in direct dwarf file"
close, 1
openr, 1, dir_dwarfile
dir_dwarf = 0l
readf, 1, dir_dwarf
print, "dir_dwarf, ns", dir_dwarf, ns
dir_dwarfeu = fltarr(dir_dwarf) ;array of all direct dwarf eu values
readf, 1, dir_dwarfeu
close, 1
dirdwarfindex = where(dir_dwarfeu ne 0) ;use index for something              
print, "fraction of direct eu enrichment from dwarfs",n_elements(dirdwarfindex)*1.0/ns ;number of non-zero elements

;^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^;
;we look at non-host particles
isat = where(sgrp gt 1)

print, "number of stars in satellites", n_elements(isat)
nsat =  n_elements(isat)

;^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^;
print, "                 new addition block II"
;GAS PARAMETERS
ggrp = grpid(h.nsph+h.ndark:h.nbodies-1)
gihost = where(ggrp eq 1)
print, "number of gas in host", n_elements(gihost)
gnhost =  n_elements(gihost)

;gage_host =  ageofuniverse -g(ihost).tform * dyrunit
;gtform_host =  g(ihost).tform * dyrunit
;gFe_host = g(ihost).mass*gfe(ihost)
;gOxy_host = g(ihost).mass*goxy(ihost)
;geu_host = g(ihost).mass*geu(ihost)
;geu_host_tot = g(gihost).mass*geu_tot(gihost)
;gmass_host = g(gihost).mass

;gFeH_host = alog10(gFe_host/(gmass_host*fh) *(Ah/Afe))-FeH_solar
;gOFe_host = alog10(goxy_host/gfe_host *(Afe/Ao))-OFe_solar

;gasindex = where(gOFe_host gt 0)
;gindex = n_elements(gindex)
;print, "ghost, gindex", ghost, gindex

;gEuFe_host = alog10(geu_host/gfe_host*(Afe/Aeu))-Eufe_solar
;gEuFe_host_tot = alog10(geu_host_tot/gfe_host*(Afe/Aeu))-Eufe_solar
;gEuH_host = alog10(geu_host/(gmass_host*fh)*(Ah/Aeu)) - EuH_solar

;gnFenH = gFe_host/(gmass_host*fh) *(Ah/Afe)
;gnOxynH = goxy_host/(gmass_host*fh) *(Ah/Ao)
;gnEunH = geu_host/(gmass_host*fh) *(Ah/Aeu)
;gnEunH_tot = geu_host_tot/(gmass_host*fh) *(Ah/Aeu)

;^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^;
print, "                 new addition block III"
;SATELLITE PARAMETERS
age_sat = ageofuniverse -s(isat).tform * dyrunit
tform_sat = s(isat).tform * dyrunit
sFe_sat = s(isat).mass*sfe(isat)
sOxy_sat = s(isat).mass*soxy(isat)
seu_sat = s(isat).mass*seu(isat)
seu_sat_tot = s(isat).mass*dwarfeu_tot(isat)
smass_sat = s(isat).mass

FeH_sat = alog10(sFe_sat/(smass_sat*fh) *(Ah/Afe))-FeH_solar
OFe_sat = alog10(soxy_sat/sfe_sat *(Afe/Ao))-OFe_solar

mindex = where(OFe_sat gt 0)
nindex_sat = n_elements(mindex)
print, "nsat, nindex_sat", nsat, nindex_sat

EuFe_sat = alog10(seu_sat/sfe_sat*(Afe/Aeu))-Eufe_solar
EuFe_sat_tot = alog10(seu_sat_tot/sfe_sat*(Afe/Aeu))-Eufe_solar
EuH_sat = alog10(seu_sat/(smass_sat*fh)*(Ah/Aeu)) - EuH_solar

nFenH_sat = sFe_sat/(smass_sat*fh) *(Ah/Afe)
nOxynH_sat = soxy_sat/(smass_sat*fh) *(Ah/Ao)
nEunH_sat = seu_sat/(smass_sat*fh) *(Ah/Aeu)
nEunH_sat_tot = seu_sat_tot/(smass_sat*fh) *(Ah/Aeu)

;^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^;
print, "                 new addition block IV"
;make array for smoothed dwarf fraction in host 
seudwarf_host = s(ihost).mass * dwarfeu(ihost)*seu(ihost) ;eu from dwarf present in disk  
;make array for smoothed dwarf fraction in satellites  
seudwarf_sat = s(isat).mass * dwarfeu(isat)*seu(isat) ;eu present in only satellites

;make array for non-smoothed direct dwarf fraction in host        
nseudwarf_host = s(ihost).mass * dir_dwarfeu(ihost)*seu(ihost)
;make array for non-smoothed directdwarf fraction in satellites     
nseudwarf_sat = s(isat).mass * dir_dwarfeu(isat)*seu(isat)

;^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^;
print, "                 new addition block V"
;SOME HOST PARAMETERS WITH DWARF ENRICHMENT
;calc [Eu/H] with dwarf enrichment
EuFe_dhost = alog10(seudwarf_host/sfe_host*(Afe/Aeu)) - Eufe_solar
EuFe_dhost_tot = alog10(seu_host_tot/sfe_host*(Afe/Aeu))-Eufe_solar
EuH_dhost = alog10(seudwarf_host/(smass_host*fh)*(Ah/Aeu)) - EuH_solar

;SOME SATALLITE PARAMETERS WITH DWARF ENRICHMENT 
;calc sat [Eu/H] with dwarf enrichment
EuFe_dsat = alog10(seudwarf_sat/sfe_sat*(Afe/Aeu))-Eufe_solar
EuFe_dsat_tot = alog10(seu_sat_tot/sfe_sat*(Afe/Aeu))-Eufe_solar
EuH_dsat = alog10(seudwarf_sat/(smass_sat*fh)*(Ah/Aeu)) - EuH_solar

;^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^;
print, "                 new addition block VI"
;SOME HOST PARAMETERS WITH DIRECT DWARF ENRICHMENT
;calc [Eu/H] with direct dwarf enrichment
EuFe_dir_dhost = alog10(nseudwarf_host/sfe_host*(Afe/Aeu)) - Eufe_solar
EuFe_dir_dhost_tot = alog10(seu_host_tot/sfe_host*(Afe/Aeu))-Eufe_solar
EuH_dir_dhost = alog10(nseudwarf_host/(smass_host*fh)*(Ah/Aeu)) - EuH_solar

;SOME SATALLITE PARAMETERS WITH DIRECT DWARF ENRICHMENT
;calc sat [Eu/H] with direct dwarf enrichment
EuFe_dir_dsat = alog10(nseudwarf_sat/sfe_sat*(Afe/Aeu))-Eufe_solar
EuFe_dir_dsat_tot = alog10(seu_sat_tot/sfe_sat*(Afe/Aeu))-Eufe_solar
EuH_dir_dsat = alog10(nseudwarf_sat/(smass_sat*fh)*(Ah/Aeu)) - EuH_solar

;^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^;
;CENTRAL POSITION VECTOR
;plot:
print, "Read in position center file"
close, 1
openr, 1, posCenter
row = fltarr(18)
readf, 1, row
print, "Xc, Yc, Zc", row[13], row[14], row[15]
close, 1

;^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^;
print, "                 new addition block VII"
;UNIT CHANGE
kpcunit = 1000
Mpc = 1000000
lunit = 90

Xc_codeunit=(row[13]/lunit)-0.5
Yc_codeunit=(row[14]/lunit)-0.5
Zc_codeunit=(row[15]/lunit)-0.5

;ALL PARTICLES
xgas = (g.pos[0] - Xc_codeunit) * lunit * kpcunit
ygas = (g.pos[1] - Yc_codeunit) * lunit * kpcunit
zgas = (g.pos[2] - Zc_codeunit) * lunit * kpcunit

xstar = (s.pos[0] - Xc_codeunit) * lunit * kpcunit
ystar = (s.pos[1] - Yc_codeunit) * lunit * kpcunit
zstar = (s.pos[2] - Zc_codeunit) * lunit * kpcunit

;HOST
xgas_host = (g(ihost).pos[0] - Xc_codeunit) * lunit * kpcunit
ygas_host = (g(ihost).pos[1] - Yc_codeunit) * lunit * kpcunit
zgas_host = (g(ihost).pos[2] - Zc_codeunit) * lunit * kpcunit

xstar_host = (s(ihost).pos[0] - Xc_codeunit) * lunit * kpcunit
ystar_host = (s(ihost).pos[1] - Yc_codeunit) * lunit * kpcunit
zstar_host = (s(ihost).pos[2] - Zc_codeunit) * lunit * kpcunit

;SATELLITE
xgas_sat = (g(isat).pos[0] - Xc_codeunit) * lunit * kpcunit
ygas_sat = (g(isat).pos[1] - Yc_codeunit) * lunit * kpcunit
zgas_sat = (g(isat).pos[2] - Zc_codeunit) * lunit * kpcunit

xstar_sat = (s(isat).pos[0] - Xc_codeunit) * lunit * kpcunit
ystar_sat = (s(isat).pos[1] - Yc_codeunit) * lunit * kpcunit
zstar_sat = (s(isat).pos[2] - Zc_codeunit) * lunit * kpcunit

;^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^;
msolunit=1.078e17
sat_mass = fltarr(10671) 
sat_eumass = fltarr(10671)
sat_deumass = fltarr(10671)
satgrp = fltarr(10671)
for i=0,10670 do begin
   grpidx = where(sgrp eq i)
   if (grpidx[0] eq -1) then continue
   sat_mass[i] = total(s(grpidx).mass)
   sat_eumass[i] = total(s(grpidx).mass*seu(grpidx))
   sat_deumass[i] = total(s(grpidx).mass*seu(grpidx)*dwarfeu(grpidx))
   print, sat_mass[i]*msolunit, sat_eumass[i]*msolunit, sat_deumass[i]*msolunit
endfor

;^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^;
;--------------------------------------------------------------------;
;                                BLOCK3                              ;
print, "                         block3"
print, "Now Binning"
if (output eq 1) then begin 
   close, 1
   openw, 1, outfile
;printf, 1, "Assuming solar abunfance [Fe/H] = 7.45, [O/H]=8.66, [Eu/H] = 0.52" 
   printf, 1, "Age [yr], n_Fe/n_H, n_O/n_H, n_Eu/n_H (NSNS merger only), n_Eu/n_H (merger+Type II), distance [kpc]"
   printf, 1, nhost 
   for i=0l,nhost-1 do begin 
      printf, 1, age_host[i], nFenH[i], nOxynH[i], nEunH[i], nEunH_tot[i]
   endfor   
   close, 1 
endif 

;--------------------------------------------------------------------;
;                                  AGE                               ;
plot:
logtmin = 8.5
logtmax = 11.5
nbins = 21

;HOST
;age|FeH
binning, tform_host, FeH_host, FeH_host, logtmin, logtmax, nbins, logtbins, median_FeH, percent_up_FeH, percent_down_FeH, number, /log
binning_ave, tform_host, sfe_host, logtmin, logtmax, nbins, logtbins, ave_sfe, rms_sfe, /log
binning_ave, tform_host, smass_host, logtmin, logtmax, nbins, logtbins, ave_smass, rms_smass, /log
ave_FeH = alog10(ave_sfe/(ave_smass*fh) *(Ah/Afe))-FeH_solar

;age|OFe
binning, tform_host, OFe_host, OFe_host, logtmin, logtmax, nbins, logtbins, median_OFe_time, percent_up_OFe_time, percent_down_OFe_time, number, /log
;age|EuFe
binning, tform_host, EuFe_host, seu_host, logtmin, logtmax, nbins, logtbins, median_EuFe_age, percent_up_EuFe_age, percent_down_EuFe_age, number1, /log, /nozero
;age|EuFe dhost
binning, tform_host, EuFe_dhost, seudwarf_host, logtmin, logtmax, nbins, logtbins, median_EuFe_dhost_age, percent_up_EuFe_dhost_age, percent_down_EuFe_dhost_age, number1d, /log, /nozero
;age|EuFe dir_dhost 
binning, tform_host, EuFe_dir_dhost, nseudwarf_host, logtmin, logtmax, nbins, logtbins, median_dir_EuFe_dhost_age, percent_up_EuFe_dir_dhost_age, percent_down_EuFe_dir_dhost_age, number1dir_d, /log, /nozero
;age|EuFe_tot
binning, tform_host, EuFe_host_tot, seu_host_tot, logtmin, logtmax, nbins, logtbins, median_EuFe_age_tot, percent_up_EuFe_age_tot, percent_down_EuFe_age_tot, number1_tot, /log, /nozero
;age|EuH
binning, tform_host, EuH_host, seu_host, logtmin, logtmax, nbins, logtbins, median_EuH_age, percent_up_EuH_age, percent_down_EuH_age, number11, /log, /nozero
;age|EuH dhost
binning, tform_host, EuH_dhost, seudwarf_host, logtmin, logtmax, nbins, logtbins, median_EuH_dhost_age, percent_up_EuH_dhost_age, percent_down_EuH_dhost_age, number11d, /log, /nozero
;age|EuH dir_dhost
binning, tform_host, EuH_dir_dhost, nseudwarf_host, logtmin, logtmax, nbins, logtbins, median_dir_EuH_dhost_age, percent_up_EuH_dir_dhost_age, percent_down_EuH_dir_dhost_age, number11dir_d, /log, /nozero

;age|dwarf fraction
binning, tform_host, dwarfeu(ihost), dwarfeu(ihost), logtmin, logtmax, nbins, logtbins, median_deu_age17, percent_up_deu_age17, percent_down_deu_age17, number17deu, /log, /nozero
;age|dir dwarf fraction
binning, tform_host, dir_dwarfeu(ihost), dir_dwarfeu(ihost), logtmin, logtmax, nbins, logtbins, median_dir_deu_age17, percent_up_dir_deu_age17, percent_down_dir_deu_age17, number17ddeu, /log, /nozero
;goto, plot2
;age|seu 
binning, tform_host, seu_host, seu_host, logtmin, logtmax, nbins, logtbins, median_eu_age, percent_up_eu_age, percent_down_eu_age, number1deu, /log, /nozero
;age|seudwarf         
binning, tform_host, seudwarf_host, seudwarf_host, logtmin, logtmax, nbins, logtbins, median_dseu_age, percent_up_dseu_age, percent_down_dseu_age, number1deu, /log, /nozero

;SAT
;age|FeH
binning, tform_sat, FeH_sat, FeH_sat, logtmin, logtmax, nbins, logtbins, median_FeH_sat, percent_up_FeH_sat, percent_down_FeH_sat, numbers, /log
binning_ave, tform_sat, sfe_sat, logtmin, logtmax, nbins, logtbins, ave_sfe, rms_sfe, /log
binning_ave, tform_sat, smass_sat, logtmin, logtmax, nbins, logtbins, ave_smass, rms_smass, /log
ave_FeH_sat = alog10(ave_sfe/(ave_smass*fh) *(Ah/Afe))-FeH_solar

;age|OFe
binning, tform_sat, OFe_sat, OFe_sat, logtmin, logtmax, nbins, logtbins, median_OFe_time_sat, percent_up_OFe_time_sat, percent_down_OFe_time_sat, numbers, /log
;age|EuFe
binning, tform_sat, EuFe_sat, seu_sat, logtmin, logtmax, nbins, logtbins, median_EuFe_age_sat, percent_up_EuFe_age_sat, percent_down_EuFe_age_sat, number1s, /log, /nozero
;age|EuFe dsat
binning, tform_sat, EuFe_dsat, seudwarf_sat, logtmin, logtmax, nbins, logtbins, median_EuFe_dsat_age, percent_up_EuFe_dsat_age, percent_down_EuFe_dsat_age, number1sd, /log, /nozero
;age|EuFe dir_dsat
binning, tform_sat, EuFe_dir_dsat, nseudwarf_sat, logtmin, logtmax, nbins, logtbins, median_dir_EuFe_dsat_age, percent_up_EuFe_dir_dsat_age, percent_down_EuFe_dir_dsat_age, number1sdir_d, /log, /nozero
;age|EuFe_tot
binning, tform_sat, EuFe_sat_tot, seu_sat_tot, logtmin, logtmax, nbins, logtbins, median_EuFe_age_tot_sat, percent_up_EuFe_age_tot_sat, percent_down_EuFe_age_tot_sat, number1s_tot, /log, /nozero
;age|EuH
binning, tform_sat, EuH_sat, seu_sat, logtmin, logtmax, nbins, logtbins, median_EuH_age_sat, percent_up_EuH_age_sat, percent_down_EuH_age_sat, number11s, /log, /nozero
;age|EuH dsat
binning, tform_sat, EuH_dsat, seudwarf_sat, logtmin, logtmax, nbins, logtbins, median_EuH_dsat_age, percent_up_EuH_dsat_age, percent_down_EuH_dsat_age, number11sd, /log, /nozero
;age|EuH dir_dsat
binning, tform_sat, EuH_dir_dsat, nseudwarf_sat, logtmin, logtmax, nbins, logtbins, median_dir_EuH_dsat_age, percent_up_EuH_dir_dsat_age, percent_down_EuH_dir_dsat_age, number11sdir_d, /log, /nozero
plot2:
;age|dwarf fraction
binning, tform_sat, dwarfeu(isat), dwarfeu(isat), logtmin, logtmax, nbins, logtbins, median_deu_age_sat17, percent_up_deu_age_sat17, percent_down_deu_age_sat17, number17sdeu, /log, /nozero
;age|dir dwarf fraction
binning, tform_sat, dir_dwarfeu(isat), dir_dwarfeu(isat), logtmin, logtmax, nbins, logtbins, median_dir_deu_age_sat17, percent_up_dir_deu_age_sat17, percent_down_dir_deu_age_sat17, number17sddeu, /log, /nozero
;goto, plot3
;age|seu
binning, tform_sat, seu_sat, seu_sat, logtmin, logtmax, nbins, logtbins, median_eu_age_sat, percent_up_eu_age_sat, percent_down_eu_age_sat, number1seu, /log, /nozero
;age|seudwarf
binning, tform_sat, seudwarf_sat, seudwarf_sat, logtmin, logtmax, nbins, logtbins, median_dseu_age_sat, percent_up_dseu_age_sat, percent_down_dseu_age_sat, number1sdseu, /log, /nozero

;--------------------------------------------------------------------;
;                             METALLICITY                            ;
;METALLICITY
zmin = -4.0
zmax = 1.0 
nzbins = 51

;HOST
;FeH|OFe
binning, FeH_host, OFe_host, sOxy_host, zmin, zmax, nzbins, zbins, median_OFe, percent_up_OFe, percent_down_OFe,number_oxy,/nozero 
;FeH|EuFe
binning, FeH_host, EuFe_host, seu_host, zmin, zmax, nzbins, zbins, median_EuFe, percent_up_EuFe, percent_down_EuFe,number,/nozero 
;FeH|EuFe dhost
binning, FeH_host, EuFe_dhost, seudwarf_host, zmin, zmax, nzbins, zbins, median_EuFe_dhost, percent_up_EuFe_dhost, percent_down_EuFe_dhost, numberd, /nozero
;FeH|EuFe dir_dhost
binning, FeH_host, EuFe_dir_dhost, nseudwarf_host, zmin, zmax, nzbins, zbins, median_EuFe_dir_dhost, percent_up_EuFe_dir_dhost, percent_down_EuFe_dir_dhost, numberdir_d, /nozero
;FeH|EuFe_tot
binning, FeH_host, EuFe_host_tot, seu_host_tot, zmin, zmax, nzbins, zbins, median_EuFe1, percent_up_EuFe1, percent_down_EuFe1, number_tot,/nozero
;FeH|dwarf fraction
binning, FeH_host, dwarfeu, dwarfeu, zmin, zmax, nzbins, zbins, median_z_deu, percent_up_z_deu, percent_down_z_deu, number1FeD, /log, /nozero
;FeH|seu
binning, FeH_host, seu_host, seu_host, zmin, zmax, nzbins, zbins, median_z_eu, percent_up_z_eu, percent_down_z_eu, number1Fe, /log, /nozero
;FeH|seudwarf
binning, FeH_host, seudwarf_host, seudwarf_host, zmin, zmax, nzbins, zbins, median_z_dseu, percent_up_z_dseu, percent_down_z_dseu, number1Fe_dseu, /log, /nozero

;SAT
;FeH|OFe
binning, FeH_sat, OFe_sat, sOxy_sat, zmin, zmax, nzbins, zbins, median_OFe_sat, percent_up_OFe_sat, percent_down_OFe_sat, number_oxy_sat, /nozero
;FeH|EuFe
binning, FeH_sat, EuFe_sat, seu_sat, zmin, zmax, nzbins, zbins, median_EuFe_sat, percent_up_EuFe_sat, percent_down_EuFe_sat, numbers, /nozero
;FeH|EuFe dsat
binning, FeH_sat, EuFe_dsat, seudwarf_sat, zmin, zmax, nzbins, zbins, median_EuFe_dsat, percent_up_EuFe_dsat, percent_down_EuFe_dsat, numbersd, /nozero
;FeH|EuFe dir_dsat
binning, FeH_sat, EuFe_dir_dsat, nseudwarf_sat, zmin, zmax, nzbins, zbins, median_EuFe_dir_dsat, percent_up_EuFe_dir_dsat, percent_down_EuFe_dir_dsat, numbersdir_d, /nozero
;FeH|EuFe_tot
binning, FeH_sat, EuFe_sat_tot, seu_sat_tot, zmin, zmax, nzbins, zbins, median_EuFe1_sat, percent_up_EuFe1_sat, percent_down_EuFe1_sat, numbers_tot_sat,/nozero
;FeH|dwarf fraction
binning, FeH_sat, dwarfeu, dwarfeu, zmin, zmax, nzbins, zbins, median_z_deu_sat, percent_up_z_deu_sat, percent_down_z_deu_sat, number1sFeD, /log, /nozero
;FeH|seu
binning, FeH_sat, seu_sat, seu_sat, zmin, zmax, nzbins, zbins, median_z_eu_sat, percent_up_z_eu_sat, percent_down_z_eu_sat, number1sFe, /log, /nozero
;FeH|seudwarf
binning, FeH_sat, seudwarf_sat, seudwarf_sat, zmin, zmax, nzbins, zbins, median_z_dseu_sat, percent_up_z_dseu_sat, percent_down_z_dseu_sat, number1sFe_dseu, /log, /nozero

;--------------------------------------------------------------------;
;                               EUROPIUM                             ;
;HOST
;EuFe|EuH
binning, EuFe_host, EuH_host, seu_host, zmin, zmax, nzbins, zbins, median_EuFe_EuH, percent_up_EuFe_EuH, percent_down_EuFe_EuH, number2, /log, /nozero
;EuFe|EuH dhost
binning, EuFe_host, EuH_dhost, seudwarf_host, zmin, zmax, nzbins, zbins, median_EuFe_EuH_dhost, percent_up_EuFe_EuH_dhost, percent_down_EuFe_EuH_dhost, number2d, /log, /nozero
;Dwarfrac|EuH
binning, dwarfeu, EuH_host, dwarfeu, zmin, zmax, nzbins, zbins, median_dfrac_EuH, percent_up_dfrac_EuH, percent_down_dfrac_EuH, number22, /log, /nozero
;Dwarfrac|EuH dhost
binning, dwarfeu, EuH_dhost, dwarfeu, zmin, zmax, nzbins, zbins, median_dfrac_EuH_dhost, percent_up_dfrac_EuH_dhost, percent_down_dfrac_EuH_dhost, number22d, /log, /nozero

;SAT
;EuFe|EuH
binning, EuFe_sat, EuH_sat, seu_sat, zmin, zmax, nzbins, zbins, median_EuFe_EuH_sat, percent_up_EuFe_EuH_sat, percent_down_EuFe_EuH_sat, number2s, /log, /nozero
;EuFe|EuH dsat
binning, EuFe_sat, EuH_dsat, seudwarf_sat, zmin, zmax, nzbins, zbins, median_EuFe_EuH_dsat, percent_up_EuFe_EuH_dsat, percent_down_EuFe_EuH_dsat, number2sd, /log, /nozero
;Dwarfrac|EuH
binning, dwarfeu, EuH_sat, dwarfeu, zmin, zmax, nzbins, zbins, median_dfrac_EuH_sat, percent_up_dfrac_EuH_sat, percent_down_dfrac_EuH_sat, number22, /log, /nozero
;Dwarfrac|EuH dhost
binning, dwarfeu, EuH_dsat, dwarfeu, zmin, zmax, nzbins, zbins, median_dfrac_EuH_dsat, percent_up_dfrac_EuH_dsat, percent_down_dfrac_EuH_dsat, number22d, /log, /nozero

;--------------------------------------------------------------------;
;                               POSITION                             ;

;ALL


;HOST


;SAT


;--------------------------------------------------------------------;
;--------------------------------------------------------------------;
print, "                       pre-plot"
plot3:
zmin1 = -5.0
zmax1 = 0.0
nzbins = 51
dz = (zmax1 - zmin1)/(nzbins-1)
zbins1 = findgen(nzbins)*dz + zmin1
index = where(FeH_host lt -2.5 and Ofe_host gt 0)
help, index
hist = histogram(EuH_host(index), binsize = dz, min=zmin1, max=zmax1)
histlike, zbins1, hist, zbins1p, histp 

index = where(FeH_host lt -2.0)
index1 = where(FeH_host lt -2.0 and EuFe_host gt 0.5)

logamin = 8.0
logamax = 10.2
nabins = 23
da = (logamax-logamin)/(nabins-1)
abins = findgen(nabins)*da + logamin 
hist_time = histogram(alog10(tform_host(index)), binsize = da, min = logamin, max=logamax)
hist_time1 = histogram(alog10(tform_host(index1)), binsize = da, min = logamin, max=logamax)

histlike, abins, hist_time, abinsp, hist_timep
histlike, abins, hist_time1, abinsp, hist_timep1

;--------------------------------------------------------------------;
;--------------------------------------------------------------------;
;                                PLOT1                               ;
print, "                         plots"
;retall
;goto, plot:
;goto, plot2
;goto, plot3
;goto, plot4
;goto, plot5
;plot:
if (1 eq 0) then begin
loadct, 39 
if(plottops eq 1) then begin
   entry_device = !d.name
   set_plot,'PS'
   device, bits_per_pixel = 8,/portrait, color = 1,filename =filename1
   tk = 2
endif else begin 
   window, 1, retain = 2
   tk = 1
endelse
;plot, zbins1p, histp, xthick=tk, ythick = tk, charthick = tk, xtitle = '[Eu/H]', ytitle='distribution dN/d[Eu/H]', xrange = [-5, 0], charsize = 1.5, thick = tk
;label = ['[Fe/H] < -2.5, no smooth']
plot, abinsp, hist_timep1, xthick=tk, ythick = tk, charthick = tk, xtitle = 'log (Formation time [yr])', ytitle='Age distribution', xrange = [8, 10.2], charsize = 1.1, thick = tk
;oplot, abinsp, hist_timep1, thick = tk, linestyle = 2
label = ['[Fe/H] < -2.0 and [Eu/Fe] > 0.5, no smooth' ]
;label = ['[Fe/H] < -2.0, no smooth']
;lines = [0, 2]
legend, label, box =0, /top, /left

if (plottops eq 1) then begin 
   device, /close_file
   set_plot, entry_device 
endif
endif

;--------------------------------------------------------------------;
;                                PLOT2                               ;
;plot:
data2=[8.1971463,     -0.41588226, $
       1.0827049,     -0.68736558, $
       2.2585745,     -0.56238463, $
       3.8856139,     -0.35125073, $
       4.3981409,    -0.070412221, $
       4.7890802,     -0.33748231, $
       5.3823209,     -0.32644775, $
       5.7757495,     -0.24863002, $
       6.8323426,     -0.16256979, $
       6.4782772,   -0.0095409624, $
       7.7495623,    -0.068144151, $
       8.5163528,     -0.11276329, $
       8.9205413,    -0.054416710, $
       9.8075567,    -0.054552561, $
       10.124896,    -0.010099460, $
       10.777775,     0.020395469, $
       11.240031,    0.0064178902, $
       11.859208,     0.086982397, $
       12.172605,     0.025744557, $
       12.341887,    -0.054940706, $
       12.898253,     0.070135124, $
       13.233793,     0.056176952, $
       13.890792,     0.028262765, $
       14.185663,      0.10609558, $
       11.633342,     -0.24674578]
nn2 =n_elements(data2)
data2x=data2(0:nn2-1:2)
data2y=data2(1:nn2-1:2)

if (1 eq 0) then begin
loadct, 39
if(plottops eq 1) then begin
   entry_device = !d.name
   set_plot,'PS'
   device, bits_per_pixel = 8,/portrait, color = 1,filename =filename2 ;have from original run
   tk = 2
endif else begin 
   window, 2, retain = 2
   tk = 1
endelse
plot, alog10(tform_host(0:nhost-1:100)), FeH_host(0:nhost-1:100), psym = sym(0), xthick=2, ythick=2, charthick=2, xtitle = 'log (Formation time [yr])', ytitle='[Fe/H]', /nodata, xrange = [8.4, 10.3], xstyle = 1, yrange = [-4, 0.5], charsize = 1.1, Color=cgColor('black'), background=cgColor('white')
oplot,alog10(tform_host(0:nhost-1:100)), FeH_host(0:nhost-1:100),psym=sym(0), color = 250

;plot, alog10(tform_host(index(0:nindex-1:100))), FeH_host(index(0:nindex-1:100)), psym = 3, xthick=tk, ythick = tk, charthick = tk, xtitle = 'log (Formation time [yr])', ytitle='[Fe/H]', /nodata, xrange = [8, 10.5], xstyle = 1, yrange = [-4, 0.5], charsize = 1.5 
;oplot,alog10(tform_host(index(0:nindex-1:100))), FeH_host(index(0:nindex-1:100)),psym = 3, color = 250

oplot, logtbins, median_FeH, thick = tk+1, color = 0
oplot, logtbins, ave_FeH, thick = tk+10, linestyle = 1, color = 0
oplot, logtbins, percent_down_FeH, thick = tk+1, linestyle = 2, color = 0
oplot, logtbins, percent_up_FeH, thick = tk+1, linestyle = 2, color = 0
oplot, data2x, data2y, psym =2, color=85


if (plottops eq 1) then begin 
   device, /close_file
   set_plot, entry_device 
endif
endif

;--------------------------------------------------------------------;
;                                PLOT3                               ;
;plot:
if (1 eq 0) then begin
loadct, 39
if(plottops eq 1) then begin
   entry_device = !d.name
   set_plot,'PS'
   device, bits_per_pixel = 8,/portrait, color = 1,filename = filename3
   tk = 2
endif else begin
   window, 3, retain = 2
   tk = 1
endelse
plot, alog10(tform_host(0:nhost-1:50)), EuFe_host(0:nhost-1:50),psym=sym(0), xthick=2, ythick=2, charthick=2, xtitle = 'log (Formation time [yr])', ytitle='[Eu/Fe]', /nodata, xrange = [8.6, 10.3], xstyle = 1, yrange = [-2, 2], charsize = 1.1, Color=cgColor('black'), background=cgColor('white')
oplot,alog10(tform_host(0:nhost-1:50)), EuFe_host(0:nhost-1:50),psym=sym(0), color = 250, symsize = 0.3

;plot,alog10(tform_host(index(0:nindex-1:100))), EuFe_host(index(0:nindex-1:100)), psym = sym(1), xthick=tk, ythick = tk, charthick = tk, xtitle = 'log (Formation time [yr])', ytitle='[Eu/Fe]', /nodata, xrange = [8.6, 10.5], xstyle = 1, yrange = [-3, 2], charsize = 1.5 
;oplot,alog10(tform_host(index(0:nindex-1:100))),EuFe_host(index(0:nindex-1:100)),psym = sym(1), color = 70, symsize = 0.3

oplot, logtbins, median_EuFe_age, thick = tk+1, color = 0
oplot, logtbins, percent_down_EuFe_age, thick = tk+1, linestyle = 2, color = 0
oplot, logtbins, percent_up_EuFe_age, thick = tk+1, linestyle = 2, color = 0

;label=['merger + SN II, 128 Nbr']
legend, label, box =0, charthick = tk, /bottom, /left

if (plottops eq 1) then begin 
   device, /close_file
   set_plot, entry_device 
endif
endif

;--------------------------------------------------------------------;
;                                PLOT4                               ;
data4=[-2.3619047,     0.046327969, $
      -2.3314286,   -0.0081628013, $
      -2.2196825,      0.39818256, $
      -2.2323809,      0.35210218, $
      -1.9073016,      0.21779437, $
      -1.8692063,      0.22195740, $
      -1.8641270,      0.17167451, $
      -1.8996825,      0.28901799, $
      -1.7473016,      0.22187127, $
      -1.7726984,      0.26378863, $
      -1.7346032,      0.30566112, $
      -1.7193651,      0.35173971, $
      -1.6914286,      0.38523951, $
      -1.6660317,      0.44807068, $
      -1.5441270,  -0.00033917679, $
      -1.5136508,      0.12952747, $
      -1.4882540,      0.18397877, $
      -1.4577778,      0.24261641, $
      -1.4526984,      0.22166312, $
      -1.3638095,      0.37662815, $
      -1.3409524,      0.33052264, $
      -1.3257143,      0.30118229, $
      -1.3282539,      0.27604443, $
      -3.2380952,     -0.48517553, $
      -3.2406349,     -0.49774356, $
      -3.2660317,     -0.30498831, $
      -3.2660317,     -0.30498831, $
      -3.2152381,      0.52458420, $
      -3.1771428,     -0.38047006, $
      -3.1161905,       1.6432286, $
      -3.2203174,       1.5176039, $
      -1.0463492,      0.30517485, $
     -0.99047617,      0.25066613, $
     -0.95238093,      0.20035992, $
     -0.97269839,      0.15428492, $
      -1.0488889,      0.22137781, $
      -1.0768254,      0.26748691, $
      -1.1758730,     0.079009524, $
      -1.1479365,     0.028710489, $
      -1.5771428,      0.35163923, $
      -1.4831746,      0.49403084, $
      -1.4425397,      0.47724237, $
      -1.0742857,      0.50631177, $
      -1.0311111,      0.51885109, $
     -0.99301585,      0.40569576, $
     -0.92698411,      0.43916863, $
     -0.90920633,      0.38468684, $
     -0.89904760,      0.30926071, $
     -0.93206347,      0.25062486, $
      -1.0158730,      0.22554443, $
     -0.91174601,      0.25061051, $
     -0.91936506,      0.33441472, $
     -0.78730157,      0.35108117, $
     -0.85841268,      0.38884089, $
     -0.86349204,      0.41817407, $
     -0.83809522,      0.44329577, $
     -0.79492061,      0.43488539, $
     -0.81777776,      0.38462224, $
     -0.84317458,      0.30503130, $
     -0.84317458,      0.28408159, $
     -0.83301585,      0.22122529, $
     -0.84063490,      0.18771114, $
     -0.69079363,      0.14989580, $
     -0.73396823,      0.23791507, $
     -0.72380950,      0.29656708, $
     -0.69587299,      0.39291599, $
     -0.65015871,      0.41383340, $
     -0.64253966,      0.34259901, $
     -0.61968252,      0.36353257, $
     -0.60952379,      0.39704492, $
     -0.60952379,      0.41799463, $
     -0.60698411,      0.45151237, $
     -0.60190474,      0.46407860, $
     -0.58412696,      0.13725061, $
     -0.58158728,     0.057639932, $
     -0.58158728,   -0.0093991301, $
     -0.57904760,    -0.068060104, $
     -0.53587299,    -0.042950960, $
     -0.54349204,     0.049233134, $
     -0.67555553,     0.049326441, $
     -0.78730157,    -0.017633668, $
      -1.9453968,  -5.5664531e-05, $
      -1.9758730,    -0.020983839, $
      -2.0038095,     0.037695079, $
      -3.0095238,     -0.14176183, $
      -3.4438095,     -0.11212540, $
      -3.5073016,     -0.16654978, $
     -0.47492061,    -0.030424201, $
     -0.47492061,    -0.030424201, $
     -0.44698410,     0.019835357, $
     -0.44698410,     0.044975006, $
     -0.47238093,     0.070132598, $
     -0.49015871,     0.095284807, $
     -0.46222220,      0.12878460, $
     -0.42920633,      0.12038139, $
     -0.41396823,     0.074281269, $
     -0.37587299,     0.074254353, $
     -0.37079363,     0.099390412, $
     -0.36571426,      0.12033653, $
     -0.34031744,     0.074229231, $
     -0.33523807,     0.044896053, $
     -0.33523807,   -0.0095731851, $
     -0.35555553,    -0.089167716, $
     -0.37841268,     -0.13105098, $
     -0.41650791,     -0.16035365, $
     -0.38603172,    -0.051436711, $
     -0.36825395,    -0.038879448, $
     -0.36825395,    -0.022119682, $
     -0.30730156,     -0.16881070, $
     -0.17015871,     -0.12281824, $
     -0.23619045,    -0.072492287, $
     -0.27174601,     0.057421018, $
     -0.25142855,     0.074166428, $
     -0.24888887,     0.078354575, $
     -0.22857141,      0.11604969, $
     -0.22857141,      0.12023963, $
     -0.25904760,      0.16216058, $
     -0.27174601,      0.22082873, $
     -0.10412696,      0.17881089, $
    -0.048253945,     0.094972585, $
     -0.10920633,     -0.14381101, $
     -0.11682537,     -0.19408492, $
      -1.1377778,      0.38065839, $
      -1.1504762,      0.39323718, $
      -1.1961905,      0.38488960, $
      -1.2292063,      0.38910287, $
      -1.7803174,      0.56127985, $
      -1.7650793,      0.44814066, $
      -2.0368254,      0.37710366, $
      -2.0622222,      0.41064113, $
      -2.1409524,    -0.029247087]

data4_2=[-0.044692568,      0.29640610, $
    -0.094971903,      0.29649715, $
     -0.16201102,      0.32281498, $
     -0.16759761,      0.33592331, $
     -0.13966465,      0.20489062, $
     -0.11731827,      0.13935910, $
    -0.072625532,     0.060688906, $
    -0.055865754,     0.034462134, $
    -0.067038939,      0.12616984, $
     -0.12849146,      0.16557576, $
     -0.17318420,     0.073969225, $
     -0.24580991,     0.034806118, $
     -0.31284902,      0.12661499, $
     -0.30726243,      0.21829235, $
     -0.22346354,      0.32292627, $
     -0.18435739,      0.34905187, $
     -0.22346354,      0.12645312, $
     -0.27932946,      0.28373281, $
     -0.32960880,      0.46719881, $
     -0.39664791,      0.48041842, $
     -0.37430154,      0.28390480, $
     -0.35754176,      0.23148161, $
     -0.34636858,      0.23146138, $
     -0.44692725,      0.21854528, $
     -0.48603340,      0.36269641, $
     -0.48603340,      0.44128567, $
     -0.49161999,      0.53298325, $
     -0.49720658,      0.53299337, $
     -0.56424569,      0.53311478, $
     -0.54748592,      0.41520054, $
     -0.52513955,      0.20558871, $
     -0.49161999,      0.11384053, $
     -0.52513955,      0.11390124, $
     -0.58100547,     0.074707778, $
     -0.61452503,      0.21884879, $
     -0.64804459,      0.33679339, $
     -0.65363118,      0.42849097, $
     -0.65363118,      0.44158918, $
     -0.68156414,      0.52022903, $
     -0.70949711,      0.57267245, $
     -0.74301666,      0.45484927, $
     -0.72625688,      0.23214935, $
     -0.73743007,      0.15358032, $
     -0.72625688,      0.38932787, $
     -0.77094963,      0.55968553, $
     -0.80446918,      0.66453192, $
     -0.85474852,      0.70391760, $
     -0.87150830,      0.69084974, $
     -0.83240215,      0.42881472, $
     -0.72625688,      0.35003324, $
     -0.60335184,      0.29741782, $
     -0.49161999,      0.34960831, $
     -0.45251384,      0.51981422, $
     -0.54748592,      0.65096832, $
     -0.59776525,      0.74274684, $
      -2.1675978,      0.71939335, $
      -2.1452514,      0.61456720, $
      -2.1787710,      0.56223506, $
      -2.2234637,      0.65400347, $
      -2.2625699,      0.70646713, $
      -2.2849163,      0.70650760, $
      -2.3351956,      0.54942013, $
      -2.3910615,      0.44473562, $
      -2.4078213,      0.45786418, $
      -2.4692738,      0.34009158, $
      -2.5363129,      0.20923089, $
      -2.5754191,      0.11761424, $
      -2.6033521,     0.012879144, $
      -2.4972068,    -0.026607712, $
      -2.4525140,     0.038802401, $
      -2.3407822,    -0.013792783, $
      -2.2402235,     0.051516158, $
      -2.1731844,      0.15618043, $
      -2.3910615,      0.47093204, $
      -2.5418995,      0.53669626, $
      -2.6703912,      0.57622358, $
      -2.6592180,      0.45831946, $
      -2.7430169,      0.34058732, $
      -2.7765364,      0.34064803, $
      -2.8044694,      0.34069861, $
      -2.8659219,      0.38010453, $
      -2.8659219,      0.40630095, $
      -2.7541901,      0.43229503, $
      -2.7150839,      0.40602779, $
      -2.4748604,      0.70685158, $
      -2.3407822,      0.78519803, $
      -2.3016760,      0.83752005, $
      -2.3463688,      0.92928846, $
      -2.3463688,      0.96858309, $
      -2.2960894,       1.0339831, $
      -2.3016760,       1.1256807, $
      -2.4022347,       1.2568449, $
      -2.3519554,       1.3222449, $
      -2.2681565,       1.3220931, $
      -2.1675978,       1.1778307, $
      -2.2011174,       1.0731057, $
      -2.2793297,       1.0863456, $
      -2.2234637,       1.0207533, $
      -2.1731844,      0.92897482, $
      -2.1340783,      0.85031474, $
      -2.0614526,       2.2385935, $
      -2.0782123,       1.9897578, $
      -2.2067040,       1.8983031, $
      -2.2346369,       2.1603179, $
      -2.3072626,       1.8198959, $
      -2.5977655,       1.8597166, $
      -2.6256984,       1.8728654, $
      -2.7486035,     -0.32741127, $
      -2.7486035,     -0.34050948, $
      -2.7262571,     -0.48463026, $
      -2.7430169,     -0.55009096, $
      -2.5698325,     -0.56350280, $
      -2.5474861,     -0.68142716, $
      -1.9273743,       1.8454044, $
      -1.9441341,      0.74518508, $
      -1.9050280,      0.71891784, $
      -1.8659218,      0.60096313, $
      -1.8100559,      0.37819239, $
      -1.7932961,      0.22098352, $
      -1.7877095,      0.28646445, $
      -1.7597765,      0.52218164, $
      -1.7318436,      0.69240779, $
      -1.6927374,      0.86261370, $
      -1.6871508,      0.73162148, $
      -1.7653631,      0.43050429, $
      -1.7597765,      0.37810133, $
      -1.8826816,      0.48310959, $
      -2.1117319,      0.53591723, $
      -2.0614526,      0.37864766, $
      -1.9162011,      0.27359893, $
      -1.4860335,      0.71815905, $
      -1.4245810,      0.67875313, $
      -1.3854748,      0.54770021, $
      -1.3072625,      0.41657647, $
      -1.2513966,      0.33788604, $
      -1.1787709,      0.28536168, $
      -1.0558659,      0.28513910, $
      -1.0279329,      0.40297240, $
      -1.0223463,      0.54704260, $
      -1.1173184,      0.66509848, $
      -1.0391061,      0.73044789, $
     -0.95530719,      0.74339434, $
      -1.8324023,     -0.89229352, $
      -1.8268157,       1.3998831, $
      -1.6424581,       1.1113886, $
      -1.4134078,       1.1109738, $
      -1.3966480,       1.0716488, $
     -0.37988813,      0.34940597, $
     -0.32402221,      0.27071554, $
     -0.25698309,      0.27059413, $
     -0.21229035,      0.32290604, $
     -0.11731827,      0.30963583, $
    -0.039105976,      0.25710135, $
     -0.76536303,      0.28461301, $
     -0.90502785,      0.57302656, $
     -0.67039096,      0.58569984, $
     -0.49720658,      0.54609158, $
     -0.60335184,      0.46769455, $
     -0.63128481,      0.87378965]

nn4 =n_elements(data4)
data4x=data4(0:nn4-1:2)
data4y=data4(1:nn4-1:2)
nn4_2= n_elements(data4_2)
data4_2x=data4_2(0:nn4_2-1:2)
data4_2y=data4_2(1:nn4_2-1:2)

if (1 eq 0) then begin
loadct, 39 
if(plottops eq 1) then begin
   entry_device = !d.name
   set_plot,'PS'
   device, bits_per_pixel = 8,/portrait, color = 1,filename =filename4 ;have from original run
   tk = 2
endif else begin 
   window, 4, retain = 2
   tk = 1
endelse 
plot, FeH_host(0:nhost-1:25), EuFe_host(0:nhost-1:25), psym = sym(0), xthick=2, ythick = 2, charthick = 2, xtitle = '[Fe/H]', ytitle='[Eu/Fe]', /nodata, xrange = [-4, 0.5], charsize = 1.1, yrange = [-3, 2.5], Color=cgColor('black'), background=cgColor('white')
oplot,FeH_host(0:nhost-1:25), EuFe_host(0:nhost-1:25), psym = sym(0), symsize= 0.5, color = 250

;plot, FeH_host(index(0:nindex-1:50)),EuFe_host(index(0:nindex-1:50)), xthick=tk, ythick = tk, charthick = tk, xtitle = '[Fe/H]', ytitle='[Eu/Fe]', /nodata, xrange = [-4, 0.5], charsize = 1.5, yrange = [-4, 4]
;oplot, FeH_host(index(0:nindex-1:50)),EuFe_host(index(0:nindex-1:50)), psym = sym(1),symsize= 0.3, color = 70

oplot, zbins, median_EuFe, thick = tk+2, color = 0
oplot, zbins, percent_up_EuFe, linestyle = 2, thick = tk+1, color = 0
oplot, zbins, percent_down_EuFe, linestyle = 2, thick = tk+1, color = 0
oplot, data4x, data4y, psym =2, color=85
oplot, data4_2x, data4_2y, psym =1, color=135


;label=['']
legend, label, box =0, charthick = tk, /bottom, /left, color = 0

if (plottops eq 1) then begin 
   device, /close_file
   set_plot, entry_device 
endif
endif

;^^^^^^^^^^^^^^^^^^^^^^^^^^^^MY ADDITIONS^^^^^^^^^^^^^^^^^^^^^^^^^^^^;
;--------------------------------------------------------------------;
;                                PLOT5                               ;
;plot:
if (1 eq 0) then begin
loadct, 39
plottops=1 
if(plottops eq 1) then begin
   entry_device = !d.name
   set_plot,'PS'
   device, bits_per_pixel = 8,/portrait, color = 1,filename = filename5 ;same as plot4 but with dir_dwarfrac
   tk = 2
endif else begin
   window, 5, retain = 1
   tk = 1
endelse
plot, FeH_host(0:nhost-1:50), EuFe_dir_dhost(0:nhost-1:50), psym = sym(0), xthick=2, ythick = 2, charthick = 2, xtitle = '[Fe/H]', ytitle='[Eu/Fe] nosmooth', /nodata, xrange = [-4, 0.5], charsize = 1.1, yrange = [-3, 2], Color=cgColor('black'), background=cgColor('white')
oplot,FeH_host(0:nhost-1:50), EuFe_dir_dhost(0:nhost-1:50), psym = sym(0), symsize= 0.5, color = 250
oplot, zbins, median_EuFe, thick = tk+2, color = 85
oplot, zbins, percent_up_EuFe, linestyle = 4, thick = tk+1, color = 85
oplot, zbins, percent_down_EuFe, linestyle = 4, thick = tk+1, color = 85

oplot, zbins, median_EuFe_dir_dhost, thick = tk+2, color = 0
oplot, zbins, percent_up_EuFe_dir_dhost, linestyle = 2, thick = tk+1, color = 0
oplot, zbins, percent_down_EuFe_dir_dhost, linestyle = 2, thick = tk+1, color = 0

;label=['NSNS merger only + mixing 128']
;legend, label, box =0, charthick = tk, /bottom, /left, color = 0

if (plottops eq 1) then begin
   device, /close_file
   set_plot, entry_device
endif
endif
plottops = 1
;--------------------------------------------------------------------;
;                                PLOT6                               ;
;plot:
if (1 eq 0) then begin
loadct, 39
if(plottops eq 1) then begin
   entry_device = !d.name
   set_plot,'PS'
   device, bits_per_pixel = 8,/portrait, color = 1,filename = filename6 ;same as plot4 but with dwarf
   tk = 2
endif else begin
   window, 6, retain = 1
   tk = 1
endelse
plot, FeH_host(0:nhost-1:50), EuFe_dhost(0:nhost-1:50), psym = sym(0), xthick=2, ythick = 2, charthick = 2, xtitle = '[Fe/H]', ytitle='[Eu/Fe]', /nodata, xrange = [-4, 0.5], charsize = 1.1, yrange = [-6, 2], Color=cgColor('black'), background=cgColor('white')
oplot,FeH_host(0:nhost-1:50), EuFe_dhost(0:nhost-1:50), psym = sym(0), symsize= 0.4, color = 250

oplot, zbins, median_EuFe_dhost, thick = tk+2, color = 0
oplot, zbins, percent_down_EuFe_dhost, linestyle = 2, thick = tk+1, color = 0
oplot, zbins, percent_up_EuFe_dhost, linestyle = 2, thick = tk+1, color = 0

;legend, label, box =0, charthick = tk, /bottom, /left, color = 0

if (plottops eq 1) then begin
   device, /close_file
   set_plot, entry_device
endif
endif
;retall
;--------------------------------------------------------------------;
;                                PLOT7                               ;
data7=[-2.4381348,    -0.083270704, $
      -2.4837814,     -0.17278672, $
      -2.5280331,     -0.11912622, $
      -2.4297029,     -0.13695576, $
      -2.4294414,     -0.11011016, $
      -2.8424158,     -0.21777493, $
      0.12880122,     -0.40365610, $
     0.076379088,     -0.26946495, $
     -0.32996993,      0.30295875, $
     -0.45507788,      0.32971841, $
     -0.34866444,      0.22240968, $
     -0.54671974,      0.11489225, $
     -0.59966494,      0.19539221, $
     -0.67095593,      0.23113723, $
     -0.78710892,      0.25790302, $
     -0.78710892,      0.25790302, $
     -0.80536755,      0.22209662, $
     -0.98455399,      0.21302531, $
      -1.0561065,      0.22192474, $
      -1.1635661,      0.22185107, $
      -1.0309852,     0.042972505, $
     -0.93274209,     0.016194432, $
     -0.83388879,     0.052056086, $
     -0.83354009,     0.087850215, $
     -0.69026068,     0.087948432, $
      -1.7636358,      0.21249126, $
      -1.8710953,      0.21241759, $
      -2.2294682,      0.19427499, $
      -2.4883775,      0.27463376, $
      -1.1761824,     -0.15399343, $
     -0.98074242,     -0.31493196, $
     -0.86511249,     -0.39538895, $
      -1.6497494,    -0.046936380, $
      -1.6327112,     -0.13640943, $
      -2.4047559,     -0.33380506, $
      -2.3427686,     -0.40535035, $
      -2.5403008,     -0.45917659, $
      -2.5319561,     -0.52181018, $
      -2.1196307,       1.3576517, $
      0.10507468,    -0.081527351, $
    -0.028552358,    -0.010031170, $
    -0.072978468,     0.025732267, $
    -0.029772825,     -0.13531062, $
     -0.15479360,    -0.099602434, $
     -0.27033636,    -0.010196911, $
     -0.32371743,     0.025560387, $
     -0.34337088,     -0.15342254, $
     -0.26329927,     -0.20705849, $
    -0.039773902,     -0.24269915, $
     -0.27138247,     -0.11757930, $
     -0.52194708,    -0.099854115, $
     -0.58402159,    -0.037257358, $
      -1.9980340,     0.051258072, $
      -1.9893405,     0.024418613, $
      -1.9095305,    -0.056062931, $
      -1.9063049,      0.27503277, $
      -2.0131543,      0.33759883, $
      -1.2598041,      0.45444539, $
      -1.2963213,      0.38283258, $
      -1.3240580,      0.29332884, $
      -1.3524921,      0.13223684, $
      -1.4336970,     0.069541866, $
      -1.5855826,      0.10523164, $
     -0.47656203,    -0.037183696, $
     -0.61097366,    -0.046224307, $
     -0.61149671,    -0.099915501, $
     -0.54942221,     -0.16251226, $
      -2.4897965,     -0.79023546, $
      -2.4543254,     -0.82600503, $
      -1.5216774,      0.23055406, $
      -1.4410827,      0.23060931, $
      -1.4658554,      0.44535567, $
      -1.3403988,      0.45439014, $
      -1.3313566,      0.46334482, $
      -1.9964648,      0.21233165, $
      -2.0961027,     0.095933209, $
      -2.1246239,    -0.074107322, $
      -2.3205870,     0.033140018, $
      -2.4728214,     0.033035663, $
      -2.4881160,      0.30147936, $
       2.1568945,     0.036209303]

nn7 =n_elements(data7)
data7x=data4(0:nn4-1:2)
data7y=data4(1:nn4-1:2)
if (1 eq 0) then begin
loadct, 39
if(plottops eq 1) then begin
   entry_device = !d.name
   set_plot,'PS'
   device, bits_per_pixel = 8,/portrait, color = 1,filename = filename7 ;same as plot4
   tk = 2
endif else begin
   window, 7, retain = 2
   tk = 1
endelse
plot, FeH_sat(0:nsat-1), EuFe_sat(0:nsat-1), psym = sym(0), xthick=2, ythick = 2, charthick = 2, xtitle = '[Fe/H]', ytitle='[Eu/Fe]', /nodata, xrange = [-3, 0.5], charsize = 1.1, yrange = [-3, 2], Color=cgColor('black'), background=cgColor('white')
oplot,FeH_sat(0:nsat-1), EuFe_sat(0:nsat-1), psym = sym(0), symsize= 0.4, color = 250

oplot, zbins, median_EuFe, thick = tk+2, color = 85
;oplot, zbins, percent_down_EuFe, linestyle = 2, thick = tk
;oplot, zbins, percent_up_EuFe, linestyle = 2, thick = tk
oplot, zbins, median_EuFe_sat, thick = tk+2, color = 0
oplot, zbins, percent_down_EuFe_sat, linestyle = 2, thick = tk+1, color = 0
oplot, zbins, percent_up_EuFe_sat, linestyle = 2, thick = tk+1, color = 0
oplot, data7x, data7y, psym =2, color=85

label=['NSNS merger only + mixing 128']
legend, label, box =0, charthick = tk, /bottom, /left, color = 0

if (plottops eq 1) then begin
   device, /close_file
   set_plot, entry_device
endif
endif

;--------------------------------------------------------------------;
;                                PLOT8                               ;
;plot:
if (1 eq 0) then begin
loadct, 39
if(plottops eq 1) then begin
   entry_device = !d.name
   set_plot,'PS'
   device, bits_per_pixel = 8,/portrait, color = 1,filename = filename8 ;same as plot4
   tk = 2
endif else begin
   window, 8, retain = 1
   tk = 1
endelse
plot, FeH_sat(0:nsat-1), EuFe_dsat(0:nsat-1), psym = sym(0), xthick=2, ythick = 2, charthick = 2, xtitle = '[Fe/H]', ytitle='[Eu/Fe]', /nodata, xrange = [-3, 0.5], charsize = 1.1, yrange = [-6, 2], Color=cgColor('black'), background=cgColor('white')
oplot,FeH_sat(0:nsat-1), EuFe_dsat(0:nsat-1), psym = sym(0), symsize= 0.4, color = 250

;oplot, zbins, median_EuFe, thick = tk+2
;oplot, zbins, percent_down_EuFe, linestyle = 2, thick = tk
;oplot, zbins, percent_up_EuFe, linestyle = 2, thick = tk
oplot, zbins, median_EuFe_dhost, thick = tk+2, color = 85
oplot, zbins, median_EuFe_dsat, thick = tk+2, color = 0
oplot, zbins, percent_up_EuFe_dsat, linestyle = 2, thick = tk+1, color = 0
oplot, zbins, percent_down_EuFe_dsat, linestyle = 2, thick = tk+1, color = 0

;legend, label, box =0, charthick = tk, /bottom, /left, color = 0

if (plottops eq 1) then begin
   device, /close_file
   set_plot, entry_device
endif
endif

;--------------------------------------------------------------------;
;                                PLOT9                               ;
;plot:
if (1 eq 0) then begin
loadct, 39
if(plottops eq 1) then begin
   entry_device = !d.name
   set_plot,'PS'
   device, bits_per_pixel = 8,/portrait, color = 1,filename = filename9 ;same as plot4
   tk = 2
endif else begin
   window, 9, retain = 1
   tk = 1
endelse
plot, FeH_sat(0:nsat-1:50), EuFe_dir_dsat(0:nsat-1:50), psym = sym(0), xthick=2, ythick = 2, charthick = 2, xtitle = '[Fe/H]', ytitle='[Eu/Fe] no smooth', /nodata, xrange = [-4, 0.5], charsize = 1.1, yrange = [-4, 4], Color=cgColor('black'), background=cgColor('white')
oplot,FeH_sat(0:nsat-1:50), EuFe_dir_dsat(0:nsat-1:50), psym = sym(0), symsize= 0.4, color = 250

;oplot, zbins, median_EuFe, thick = tk+2
;oplot, zbins, percent_down_EuFe, linestyle = 2, thick = tk
;oplot, zbins, percent_up_EuFe, linestyle = 2, thick = tk

oplot, zbins, median_EuFe_dir_dsat, thick = tk+2, color = 0
oplot, zbins, percent_up_EuFe_dir_dsat, linestyle = 2, thick = tk+1, color = 0
oplot, zbins, percent_down_EuFe_dir_dsat, linestyle = 2, thick = tk+1, color = 0

legend, label, box =0, charthick = tk, /bottom, /left, color = 0

if (plottops eq 1) then begin
   device, /close_file
   set_plot, entry_device
endif
endif

;--------------------------------------------------------------------;
;                                PLOT10                              ;
;plot:
if (1 eq 0) then begin
loadct, 39
if(plottops eq 1) then begin
   entry_device = !d.name
   set_plot,'PS'
   device, bits_per_pixel = 8,/portrait, color = 1,filename = filename10
   tk = 2
endif else begin
   window, 10, retain = 2
   tk = 1
endelse
plot, alog10(tform_host(0:nhost-1:100)), EuH_host(0:nhost-1:100), psym = sym(0), xthick=2, ythick = 2, charthick = 2, xtitle = 'log (ForMation time [yr])', ytitle='[Eu/H]', /nodata, xrange = [8.5, 10.3], xstyle = 1, yrange = [-3, 0.5], charsize = 1.1, Color=cgColor('black'), background=cgColor('white')
oplot,alog10(tform_host(0:nhost-1:100)), EuH_host(0:nhost-1:100),psym=sym(0), symsize= 0.4, color = 250

;oplot, logtbins, median_EuFe_age, thick = tk+1, color = 0
;oplot, logtbins, percent_down_EuFe_age, thick = tk+1, linestyle = 2, color = 0
;oplot, logtbins, percent_up_EuFe_age, thick = tk+1, linestyle = 2, color = 0

oplot, logtbins, median_EuH_age, thick = tk+1, color = 0
oplot, logtbins, percent_up_EuH_age, thick = tk+1, linestyle = 2, color = 0
oplot, logtbins, percent_down_EuH_age, thick = tk+1, linestyle = 2, color = 0

;label=['merger + SN II, 128 Nbr']
;legend, label, box =0, charthick = tk, /bottom, /left

if (plottops eq 1) then begin
   device, /close_file
   set_plot, entry_device
endif
endif


;--------------------------------------------------------------------;
;                                PLOT11                              ;
;plot:
if (1 eq 0) then begin
loadct, 39
if(plottops eq 1) then begin
   entry_device = !d.name
   set_plot,'PS'
   device, bits_per_pixel = 8,/portrait, color = 1,filename = filename11
   tk = 2
endif else begin
   window, 11, retain = 2
   tk = 1
endelse
plot, alog10(tform_host(0:nhost-1:100)), EuH_dhost(0:nhost-1:100), psym = sym(0), xthick=2, ythick = 2, charthick = 2, xtitle = 'log (ForMation time [yr])', ytitle='[Eu/H] with dwarf contamination', /nodata, xrange = [8.5, 10.3], xstyle = 1, yrange = [-5, 0.5], charsize = 1.1, Color=cgColor('black'), background=cgColor('white')
oplot,alog10(tform_host(0:nhost-1:100)), EuH_dhost(0:nhost-1:100),psym=sym(0), symsize= 0.4, color = 250

;oplot, logtbins, median_EuFe_age, thick = tk+1, color = 0
;oplot, logtbins, percent_down_EuFe_age, thick = tk+1, linestyle = 2, color = 0
;oplot, logtbins, percent_up_EuFe_age, thick = tk+1, linestyle = 2, color = 0

oplot, logtbins, median_EuH_dhost_age, thick = tk+1, color = 0
oplot, logtbins, percent_up_EuH_dhost_age, thick = tk+1, linestyle = 2, color = 0
oplot, logtbins, percent_down_EuH_dhost_age, thick = tk+1, linestyle = 2, color = 0

oplot, logtbins, median_EuH_age, thick = tk+1, color = 85
oplot, logtbins, percent_up_EuH_age, thick = tk+1, linestyle = 4, color = 85
oplot, logtbins, percent_down_EuH_age, thick = tk+1, linestyle = 4, color = 85


;label=['merger + SN II, 128 Nbr']
;legend, label, box =0, charthick = tk, /bottom, /left

if (plottops eq 1) then begin
   device, /close_file
   set_plot, entry_device
endif
endif

;--------------------------------------------------------------------;      
;                                PLOT12                              ;
;plot:
if (1 eq 0) then begin
loadct, 39
plottops = 1
if(plottops eq 1) then begin
   entry_device = !d.name
   set_plot,'PS'
   device, bits_per_pixel = 8,/portrait, color = 1,filename = filename12
   tk = 2
endif else begin
   window, 12, retain = 2
   tk = 1
endelse
plot, alog10(tform_host(0:nhost-1:25)), EuH_dir_dhost(0:nhost-1:25),psym=sym(0), xthick=2, ythick = 2, charthick = 2, xtitle = 'log (ForMation time [yr])', ytitle='[Eu/H] with dwarf contamination', /nodata, xrange = [8.3, 10.3], xstyle = 1, yrange = [-3, 0.5], charsize = 1.1, Color=cgColor('black'), background=cgColor('white')
oplot,alog10(tform_host(0:nhost-1:25)), EuH_dir_dhost(0:nhost-1:25),psym=sym(0), symsize= 0.4, color = 250

oplot, logtbins, median_EuFe_age, thick = tk+1, color = 85
oplot, logtbins, percent_up_EuFe_age, thick = tk+1, linestyle = 4, color = 85
oplot, logtbins, percent_down_EuFe_age, thick = tk+1, linestyle = 4, color = 85

oplot, logtbins, median_dir_EuH_dhost_age, thick = tk+1, color = 0
oplot, logtbins, percent_up_EuH_dir_dhost_age, thick = tk+1, linestyle = 2, color = 0
oplot, logtbins, percent_down_EuH_dir_dhost_age, thick = tk+1, linestyle = 2, color = 0

;label=['merger + SN II, 128 Nbr']
;legend, label, box =0, charthick = tk, /bottom, /left

if (plottops eq 1) then begin
   device, /close_file
   set_plot, entry_device
endif
endif
plottops = 1
;--------------------------------------------------------------------;
;                                PLOT13                              ;
;plot:
if (1 eq 0) then begin
loadct, 39
if(plottops eq 1) then begin
   entry_device = !d.name
   set_plot,'PS'
   device, bits_per_pixel = 8,/portrait, color = 1,filename = filename13
   tk = 2
endif else begin
   window, 13, retain = 2
   tk = 1
endelse
plot, alog10(tform_sat(0:nsat-1)), EuH_sat(0:nsat-1),psym=sym(0), xthick=2, ythick = 2, charthick = 2, xtitle = 'log (ForMation time [yr])', ytitle='[Eu/H]', /nodata, xrange = [8.5, 10.3], xstyle = 1, yrange = [-3, 1], charsize = 1.1, Color=cgColor('black'), background=cgColor('white')
oplot,alog10(tform_sat(0:nsat-1)), EuH_sat(0:nsat-1),psym=sym(0), symsize= 0.4, color = 250

oplot, logtbins, median_EuH_age, thick = tk+1, color = 85
oplot, logtbins, median_EuH_age_sat, thick = tk+1, color = 0
oplot, logtbins, percent_up_EuH_age_sat, thick = tk+1, linestyle = 2, color = 0
oplot, logtbins, percent_down_EuH_age_sat, thick = tk+1, linestyle = 2, color = 0

if (plottops eq 1) then begin
   device, /close_file
   set_plot, entry_device
endif
endif

;--------------------------------------------------------------------;
;                                PLOT14                              ;
;plot:
if (1 eq 0) then begin
loadct, 39
if(plottops eq 1) then begin
   entry_device = !d.name
   set_plot,'PS'
   device, bits_per_pixel = 8,/portrait, color = 1,filename =filename14
   tk = 2
endif else begin
   window, 14, retain = 2
   tk = 1
endelse
plot, alog10(tform_sat(0:nsat-1)), EuH_dsat(0:nsat-1),psym=sym(0), xthick=2, ythick = 2, charthick = 2, xtitle = 'log (ForMation time [yr])', ytitle='[Eu/H]', /nodata, xrange = [8.5, 10.3], xstyle = 1, yrange = [-6, 0.5], charsize = 1.1, Color=cgColor('black'), background=cgColor('white')
oplot,alog10(tform_sat(0:nsat-1)), EuH_dsat(0:nsat-1),psym=sym(0), symsize= 0.4, color = 250

oplot, logtbins, median_EuH_dhost_age, thick = tk+1, color = 85
oplot, logtbins, median_EuH_dsat_age, thick = tk+1, color = 0
oplot, logtbins, percent_up_EuH_dsat_age, thick = tk+1, linestyle = 2, color = 0
oplot, logtbins, percent_down_EuH_dsat_age, thick = tk+1, linestyle = 2, color = 0

if (plottops eq 1) then begin
   device, /close_file
   set_plot, entry_device
endif
endif

;--------------------------------------------------------------------;
;                                PLOT15                              ;
;plot:
if (1 eq 0) then begin
loadct, 39
if(plottops eq 1) then begin
   entry_device = !d.name
   set_plot,'PS'
   device, bits_per_pixel = 8,/portrait, color = 1,filename =filename15
   tk = 2
endif else begin
   window, 15, retain = 2
   tk = 1
endelse
plot, alog10(tform_sat(0:nsat-1:100)), EuH_dir_dsat(0:nsat-1:100),psym=sym(0), xthick=2, ythick = 2, charthick = 2, xtitle = 'log (ForMation time [yr])', ytitle='[Eu/H]', /nodata, xrange = [8.5, 10.5], xstyle = 1, yrange = [-4, 0.5], charsize = 1.1, Color=cgColor('black'), background=cgColor('white')
oplot,alog10(tform_sat(0:nsat-1:100)), EuH_dir_dsat(0:nsat-1:100),psym=sym(0), symsize= 0.4, color = 250

oplot, logtbins, median_dir_EuH_dsat_age, thick = tk+1, color = 0
oplot, logtbins, percent_up_EuH_dir_dsat_age, thick = tk+1, linestyle = 2, color = 0
oplot, logtbins, percent_down_EuH_dir_dsat_age, thick = tk+1, linestyle = 2, color = 0

if (plottops eq 1) then begin
   device, /close_file
   set_plot, entry_device
endif
endif

;--------------------------------------------------------------------;
;                                PLOT16                              ;
;plot:
if (1 eq 0) then begin
loadct, 39
if(plottops eq 1) then begin
   entry_device = !d.name
   set_plot,'PS'
   device, bits_per_pixel = 8,/portrait, color = 1,filename =filename16
   tk = 2
endif else begin
   window, 16, retain = 2
   tk = 1
endelse
plot, alog10(tform_sat(0:nsat-1:100)), FeH_sat(0:nsat-1:100),psym=sym(0), xthick=2, ythick=2, charthick=2, xtitle = 'log (Formation time [yr])', ytitle='[Fe/H]', /nodata, xrange = [8.4, 10.3], xstyle = 1, yrange = [-4.5, 0.5], charsize = 1.1, Color=cgColor('black'), background=cgColor('white')
oplot,alog10(tform_sat(0:nsat-1:100)), FeH_sat(0:nsat-1:100),psym=sym(0), color = 250

;oplot, logtbins, median_FeH, thick = tk+1, color = 0
;oplot, logtbins, ave_FeH, thick = tk+10, linestyle = 1, color = 0
;oplot, logtbins, percent_down_FeH, thick = tk+1, linestyle = 2, color = 0
;oplot, logtbins, percent_up_FeH, thick = tk+1, linestyle = 2, color = 0

oplot, logtbins, median_FeH_sat, thick = tk+1, color = 0
oplot, logtbins, ave_FeH_sat, thick = tk+10, linestyle = 1, color = 0
oplot, logtbins, percent_down_FeH_sat, thick = tk+1, linestyle = 2, color = 0
oplot, logtbins, percent_up_FeH_sat, thick = tk+1, linestyle = 2, color = 0

if (plottops eq 1) then begin
   device, /close_file
   set_plot, entry_device
endif
endif

;--------------------------------------------------------------------;
;                                PLOT17                              ;
;plot:
if (1 eq 1) then begin
loadct, 39
if(plottops eq 1) then begin
   entry_device = !d.name
   set_plot,'PS'
   device, bits_per_pixel = 8,/portrait, color = 1,filename =filename17
   tk = 2
endif else begin
   window, 17, retain = 2
   tk = 1
endelse
plot, alog10(tform_host(0:nhost-1:100)), alog10(dwarfeu(ihost(0:nhost-1:100))),psym=sym(0), xthick=2, ythick = 2, charthick = 2, xtitle = 'log (ForMation time [yr])', ytitle='Dwarf fraction', /nodata, xrange = [8.5, 10.3], xstyle = 1, yrange = [-5, 0.1], charsize = 1.1, Color=cgColor('black'), background=cgColor('white')
oplot,alog10(tform_host(0:nhost-1:100)), alog10(dwarfeu(ihost(0:nhost-1:100))),psym=sym(0), symsize= 0.4, color = 250

oplot, logtbins, alog10(median_deu_age17), thick = tk+1, color = 0
oplot, logtbins, alog10(percent_up_deu_age17), thick = tk+1, linestyle = 2, color = 0
oplot, logtbins, alog10(percent_down_deu_age17), thick = tk+1, linestyle = 2, color = 0

if (plottops eq 1) then begin
   device, /close_file
   set_plot, entry_device
endif
endif

;                                PLOT17.2                            ;
;plot:                                                                                
if (1 eq 1) then begin
loadct, 39
plottops = 1
if(plottops eq 1) then begin
   entry_device = !d.name
   set_plot,'PS'
   device, bits_per_pixel = 8,/portrait, color = 1,filename = '17.2Age_dir_dwarfrac_128.eps'
   tk = 2
endif else begin
   window, 0, retain = 2
   tk = 1
endelse
plot, alog10(tform_host(0:nhost-1:10)), alog10(dir_dwarfeu(ihost(0:nhost-1:10))),psym=sym(0), xthick=2, ythick = 2, charthick = 2, xtitle = 'log (ForMation time [yr])', ytitle='Direct dwarf fraction', /nodata, xrange = [8.5, 10.3], xstyle = 1, yrange = [-5, 0.1], charsize = 1.1, Color=cgColor('black'), background=cgColor('white')
oplot,alog10(tform_host(0:nhost-1:10)), alog10(dir_dwarfeu(ihost(0:nhost-1:10))),psym=sym(0), symsize= 0.4, color = 250

oplot, logtbins, alog10(median_dir_deu_age17), thick = tk+1, color = 0
oplot, logtbins, alog10(percent_up_dir_deu_age17), thick = tk+1, linestyle = 2, color = 0
oplot, logtbins, alog10(percent_down_dir_deu_age17), thick = tk+1, linestyle = 2, color = 0

if (plottops eq 1) then begin
   device, /close_file
   set_plot, entry_device
endif
endif
plottops = 1
;--------------------------------------------------------------------;
;                                PLOT18                              ;
;plot:
if (1 eq 0) then begin
loadct, 39
if(plottops eq 1) then begin
   entry_device = !d.name
   set_plot,'PS'
   device, bits_per_pixel = 8,/portrait, color = 1,filename =filename18
   tk = 2
endif else begin
   window, 18, retain = 2
   tk = 1
endelse
plot, alog10(tform_host(0:nhost-1:100)), alog10(seu_host(0:nhost-1:100)),psym=sym(0), xthick=2, ythick = 2, charthick = 2, xtitle = 'log (ForMation time [yr])', ytitle='Eu mass', /nodata, xrange = [8.5, 10.3], xstyle = 1, yrange = [-30, -20], charsize = 1.1, Color=cgColor('black'), background=cgColor('white')
oplot,alog10(tform_host(0:nhost-1:100)), alog10(seu_host(0:nhost-1:100)),psym=sym(0), symsize= 0.4, color = 250

oplot, logtbins, alog10(median_eu_age), thick = tk+1, color = 0
oplot, logtbins, alog10(percent_up_eu_age), thick = tk+1, linestyle = 2, color = 0
oplot, logtbins, alog10(percent_down_eu_age), thick = tk+1, linestyle = 2, color = 0

if (plottops eq 1) then begin
   device, /close_file
   set_plot, entry_device
endif
endif

;--------------------------------------------------------------------;
;                                PLOT19                              ;
;plot:
if (1 eq 0) then begin
loadct, 39
if(plottops eq 1) then begin
   entry_device = !d.name
   set_plot,'PS'
   device, bits_per_pixel = 8,/portrait, color = 1,filename =filename19
   tk = 2
endif else begin
   window, 19, retain = 2
   tk = 1
endelse
plot, alog10(tform_host(0:nhost-1:100)), alog10(seudwarf_host(0:nhost-1:100)),psym=sym(0), xthick=2, ythick = 2, charthick = 2, xtitle = 'log (ForMation time [yr])', ytitle='Eu dwarf mass', /nodata, xrange = [8.5, 10.3], xstyle = 1, yrange = [-35, -20], charsize = 1.1, Color=cgColor('black'), background=cgColor('white')
oplot,alog10(tform_host(0:nhost-1:100)), alog10(seudwarf_host(0:nhost-1:100)),psym=sym(0), symsize= 0.4, color = 250

oplot, logtbins, alog10(median_dseu_age), thick = tk+1, color = 0
oplot, logtbins, alog10(percent_up_dseu_age), thick = tk+1, linestyle = 2, color = 0
oplot, logtbins, alog10(percent_down_dseu_age), thick = tk+1, linestyle = 2, color = 0

if (plottops eq 1) then begin
   device, /close_file
   set_plot, entry_device
endif
endif

;--------------------------------------------------------------------;
;                                PLOT20                              ;
;plot:
if (1 eq 1) then begin
loadct, 39
if(plottops eq 1) then begin
   entry_device = !d.name
   set_plot,'PS'
   device, bits_per_pixel = 8,/portrait, color = 1,filename =filename20
   tk = 2
endif else begin
   window, 20, retain = 2
   tk = 1
endelse
plot, alog10(tform_sat(0:nsat-1)), alog10(dwarfeu(isat)),psym=sym(0), xthick=2, ythick = 2, charthick = 2, xtitle = 'log (ForMation time [yr])', ytitle='Dwarf fraction' , /nodata, xrange = [8.5, 10.3], xstyle = 1, yrange = [-5, 1], charsize = 1.1;, Color=cgColor('black'), background=cgColor('white')
oplot,alog10(tform_sat(0:nsat-1)), alog10(dwarfeu(isat)),psym=sym(0), symsize= 1.4, color = 250

oplot, logtbins, alog10(median_deu_age_sat17), thick = tk+1, color = 0
oplot, logtbins, alog10(percent_up_deu_age_sat17), thick = tk+1, linestyle = 2, color = 0
oplot, logtbins, alog10(percent_down_deu_age_sat17), thick = tk+1, linestyle = 2, color = 0

if (plottops eq 1) then begin
   device, /close_file
   set_plot, entry_device
endif
endif
;retall
;--------------------------------------------------------------------;
;                                PLOT21                              ;
;plot:
if (1 eq 0) then begin
loadct, 39
if(plottops eq 1) then begin
   entry_device = !d.name
   set_plot,'PS'
   device, bits_per_pixel = 8,/portrait, color = 1,filename =filename21
   tk = 2
endif else begin
   window, 21, retain = 2
   tk = 1
endelse
plot, alog10(tform_sat(0:nsat-1:100)), alog10(seu_sat(0:nsat-1:100)),psym=sym(0), xthick=2, ythick = 2, charthick = 2, xtitle = 'log (ForMation time [yr])', ytitle='Eu mass', /nodata, xrange = [8.5, 10.3], xstyle = 1, yrange = [-27, -20], charsize = 1.1, Color=cgColor('black'), background=cgColor('white')
oplot,alog10(tform_sat(0:nsat-1:100)), alog10(seu_sat(0:nsat-1:100)),psym=sym(0), symsize= 0.4, color = 250

oplot, logtbins, alog10(median_eu_age_sat), thick = tk+1, color = 0
oplot, logtbins, alog10(percent_up_eu_age_sat), thick = tk+1, linestyle = 2, color = 0
oplot, logtbins, alog10(percent_down_eu_age_sat), thick = tk+1, linestyle = 2, color = 0

if (plottops eq 1) then begin
   device, /close_file
   set_plot, entry_device
endif
endif

;--------------------------------------------------------------------;
;                                PLOT22                              ;
;plot:
if (1 eq 0) then begin
loadct, 39
if(plottops eq 1) then begin
   entry_device = !d.name
   set_plot,'PS'
   device, bits_per_pixel = 8,/portrait, color = 1,filename =filename22
   tk = 2
endif else begin
   window, 22, retain = 2
   tk = 1
endelse
plot, alog10(tform_sat(0:nsat-1:100)), alog10(seudwarf_sat(0:nsat-1:100)),psym=sym(0), xthick=2, ythick = 2, charthick = 2, xtitle = 'log (ForMation time [yr])', ytitle='Eu dwarf mass', /nodata, xrange = [8.5, 10.3], xstyle = 1, yrange = [-35, -20], charsize = 1.1, Color=cgColor('black'), background=cgColor('white')
oplot,alog10(tform_sat(0:nsat-1:100)), alog10(seudwarf_sat(0:nsat-1:100)),psym=sym(0), symsize= 0.4, color = 250

oplot, logtbins, alog10(median_dseu_age_sat), thick = tk+1, color = 0
oplot, logtbins, alog10(percent_up_dseu_age_sat), thick = tk+1, linestyle = 2, color = 0
oplot, logtbins, alog10(percent_down_dseu_age_sat), thick = tk+1, linestyle = 2, color = 0

if (plottops eq 1) then begin
   device, /close_file
   set_plot, entry_device
endif
endif

;--------------------------------------------------------------------;
;                                PLOT23                              ;
;plot:
if (1 eq 0) then begin
loadct, 39
if(plottops eq 1) then begin
   entry_device = !d.name
   set_plot,'PS'
   device, bits_per_pixel = 8,/portrait, color = 1,filename =filename23
   tk = 2
endif else begin
   window, 23, retain = 2
   tk = 1
endelse
plot, FeH_host(0:nhost-1:50), alog10(dwarfeu(0:nhost-1:50)),psym=sym(0), xthick=2, ythick = 2, charthick = 2, xtitle = '[Fe/H]', ytitle='Dwarf fraction', /nodata, xrange = [-4.2, 0.5], xstyle = 1, yrange = [-15, 1], charsize = 1.1, Color=cgColor('black'), background=cgColor('white')
oplot,FeH_host(0:nhost-1:50), alog10(dwarfeu(0:nhost-1:100)),psym=sym(0), symsize= 0.4, color = 250

oplot, zbins, alog10(median_z_deu), thick = tk+2, color = 0
oplot, zbins, alog10(percent_up_z_deu), linestyle = 2, thick = tk+1, color = 0
oplot, zbins, alog10(percent_down_z_deu), linestyle = 2, thick = tk+1, color = 0

if (plottops eq 1) then begin
   device, /close_file
   set_plot, entry_device
endif
endif

;--------------------------------------------------------------------;
;                                PLOT24                              ;
;plot:
if (1 eq 0) then begin
loadct, 39
if(plottops eq 1) then begin
   entry_device = !d.name
   set_plot,'PS'
   device, bits_per_pixel = 8,/portrait, color = 1,filename =filename24
   tk = 2
endif else begin
   window, 24, retain = 2
   tk = 1
endelse
plot, FeH_host(0:nhost-1:50), alog10(seu_host(0:nhost-1:100)),psym=sym(0), xthick=2, ythick = 2, charthick = 2, xtitle = '[Fe/H]', ytitle='Eu mass', /nodata, xrange = [-4.5, 0.5], xstyle = 1, yrange = [-30, -20], charsize = 1.1, Color=cgColor('black'), background=cgColor('white')
oplot,FeH_host(0:nhost-1:50), alog10(seu_host(0:nhost-1:100)),psym=sym(0), symsize= 0.4, color = 250

oplot, zbins, alog10(median_z_eu), thick = tk+2, color = 0
oplot, zbins, alog10(percent_up_z_eu), linestyle = 2, thick = tk+1, color = 0
oplot, zbins, alog10(percent_down_z_eu), linestyle = 2, thick = tk+1, color = 0

if (plottops eq 1) then begin
   device, /close_file
   set_plot, entry_device
endif
endif

;--------------------------------------------------------------------;
;                                PLOT25                              ;
;plot:
if (1 eq 0) then begin
loadct, 39
if(plottops eq 1) then begin
   entry_device = !d.name
   set_plot,'PS'
   device, bits_per_pixel = 8,/portrait, color = 1,filename =filename25
   tk = 2
endif else begin
   window, 25, retain = 2
   tk = 1
endelse
plot, FeH_host(0:nhost-1:50), alog10(seudwarf_host(0:nhost-1:100)),psym=sym(0), xthick=2, ythick = 2, charthick = 2, xtitle = '[Fe/H]', ytitle='Eu dwarf mass', /nodata, xrange = [-4.5, 0.5], xstyle = 1, yrange = [-40, -20], charsize = 1.1, Color=cgColor('black'), background=cgColor('white')
oplot,FeH_host(0:nhost-1:50), alog10(seudwarf_host(0:nhost-1:100)),psym=sym(0), symsize= 0.4, color = 250

oplot, zbins, alog10(median_z_dseu), thick = tk+2, color = 0
oplot, zbins, alog10(percent_up_z_dseu), linestyle = 2, thick = tk+1, color = 0
oplot, zbins, alog10(percent_down_z_dseu), linestyle = 2, thick = tk+1, color = 0

if (plottops eq 1) then begin
   device, /close_file
   set_plot, entry_device
endif
endif

;--------------------------------------------------------------------;
;                                PLOT26                              ;
;plot:
if (1 eq 0) then begin
loadct, 39
if(plottops eq 1) then begin
   entry_device = !d.name
   set_plot,'PS'
   device, bits_per_pixel = 8,/portrait, color = 1,filename =filename26
   tk = 2
endif else begin
   window, 26, retain = 2
   tk = 1
endelse
plot, FeH_sat(0:nsat-1:50), alog10(dwarfeu(0:nsat-1:50)),psym=sym(0), xthick=2, ythick = 2, charthick = 2, xtitle = '[Fe/H]', ytitle='Dwarf fraction', /nodata, xrange = [-4.5, 0.5], xstyle = 1, yrange = [-8, 1], charsize = 1.1, Color=cgColor('black'), background=cgColor('white')
oplot,FeH_sat(0:nsat-1:50), alog10(dwarfeu(0:nsat-1:50)),psym=sym(0), symsize= 0.4, color = 250

oplot, zbins, alog10(median_z_deu_sat), thick = tk+2, color = 0
oplot, zbins, alog10(percent_up_z_deu_sat), linestyle = 2, thick = tk+1, color = 0
oplot, zbins, alog10(percent_down_z_deu_sat), linestyle = 2, thick = tk+1, color = 0

if (plottops eq 1) then begin
   device, /close_file
   set_plot, entry_device
endif
endif

;--------------------------------------------------------------------;
;                                PLOT27                              ;
;plot:
if (1 eq 0) then begin
loadct, 39
if(plottops eq 1) then begin
   entry_device = !d.name
   set_plot,'PS'
   device, bits_per_pixel = 8,/portrait, color = 1,filename =filename27
   tk = 2
endif else begin
   window, 27, retain = 2
   tk = 1
endelse
plot, FeH_sat(0:nsat-1:50), alog10(seu_sat(0:nsat-1:100)),psym=sym(0), xthick=2, ythick = 2, charthick = 2, xtitle = '[Fe/H]', ytitle='Eu mass', /nodata, xrange = [-4.5, 0.5], xstyle = 1, yrange = [-30, -20], charsize = 1.1, Color=cgColor('black'), background=cgColor('white')
oplot,FeH_sat(0:nsat-1:50), alog10(seu_sat(0:nsat-1:100)),psym=sym(0), symsize= 0.4, color = 250

oplot, zbins, median_z_eu_sat, thick = tk+2, color = 0
oplot, zbins, percent_up_z_eu_sat, linestyle = 2, thick = tk+1, color = 0
oplot, zbins, percent_down_z_eu_sat, linestyle = 2, thick = tk+1, color = 0

if (plottops eq 1) then begin
   device, /close_file
   set_plot, entry_device
endif
endif

;--------------------------------------------------------------------;
;                                PLOT28                              ;
;plot:
if (1 eq 0) then begin
loadct, 39
if(plottops eq 1) then begin
   entry_device = !d.name
   set_plot,'PS'
   device, bits_per_pixel = 8,/portrait, color = 1,filename =filename28
   tk = 2
endif else begin
   window, 28, retain = 2
   tk = 1
endelse
plot, FeH_sat(0:nsat-1:50), alog10(seudwarf_sat(0:nsat-1:50)),psym=sym(0), xthick=2, ythick = 2, charthick = 2, xtitle = '[Fe/H]', ytitle='Eu dwarf mass', /nodata, xrange = [-4.5, 0.5], xstyle = 1, yrange = [-35, -20], charsize = 1.1, Color=cgColor('black'), background=cgColor('white')
oplot,FeH_sat(0:nsat-1:50), alog10(seudwarf_sat(0:nsat-1:50)),psym=sym(0), symsize= 0.4, color = 250

oplot, zbins, median_z_dseu_sat, thick = tk+2, color = 0
oplot, zbins, percent_up_z_dseu_sat, linestyle = 2, thick = tk+1, color = 0
oplot, zbins, percent_down_z_dseu_sat, linestyle = 2, thick = tk+1, color = 0

if (plottops eq 1) then begin
   device, /close_file
   set_plot, entry_device
endif
endif

;--------------------------------------------------------------------;
;                                PLOT29                              ;
;plot:
if (1 eq 0) then begin
loadct, 39
if(plottops eq 1) then begin
   entry_device = !d.name
   set_plot,'PS'
   device, bits_per_pixel = 8,/portrait, color = 1,filename =filename29
   tk = 2
endif else begin
   window, 29, retain = 2
   tk = 1
endelse
plot, alog10(dwarfeu(0:nhost-1:50)), EuH_host(0:nhost-1:100),psym=sym(0), xthick=2, ythick = 2, charthick = 2, xtitle = 'Dwarf fraction', ytitle='[Eu/H]', /nodata, xrange = [-15, 1], xstyle = 1, yrange = [-6, 2], charsize = 1.1, Color=cgColor('black'), background=cgColor('white')
oplot,alog10(dwarfeu(0:nhost-1:50)), EuH_host(0:nhost-1:100),psym=sym(0), symsize= 0.4, color = 250

oplot, zbins, median_dfrac_EuH, thick = tk+2, color = 0
oplot, zbins, percent_up_dfrac_EuH, linestyle = 2, thick = tk+1, color = 0
oplot, zbins, percent_down_dfrac_EuH, linestyle = 2, thick = tk+1, color = 0

if (plottops eq 1) then begin
   device, /close_file
   set_plot, entry_device
endif
endif

;--------------------------------------------------------------------;
;                                PLOT30                              ;
;plot:
if (1 eq 0) then begin
loadct, 39
if(plottops eq 1) then begin
   entry_device = !d.name
   set_plot,'PS'
   device, bits_per_pixel = 8,/portrait, color = 1,filename =filename30
   tk = 2
endif else begin
   window, 30, retain = 2
   tk = 1
endelse
plot, alog10(dwarfeu(0:nhost-1:50)), EuH_dhost(0:nhost-1:100),psym=sym(0), xthick=2, ythick = 2, charthick = 2, xtitle = 'Dwarf fraction', ytitle='[Eu/H] with dwarf contamination', /nodata, xrange = [-15, 1], xstyle = 1, yrange = [-13, 2], charsize = 1.1, Color=cgColor('black'), background=cgColor('white')
oplot,alog10(dwarfeu(0:nhost-1:50)), EuH_dhost(0:nhost-1:100),psym=sym(0), symsize= 0.4, color = 250

oplot, zbins, median_dfrac_EuH_dhost, thick = tk+2, color = 0
oplot, zbins, percent_up_dfrac_EuH_dhost, linestyle = 2, thick = tk+1, color = 0
oplot, zbins, percent_down_dfrac_EuH_dhost, linestyle = 2, thick = tk+1, color = 0

if (plottops eq 1) then begin
   device, /close_file
   set_plot, entry_device
endif
endif
;retall

;--------------------------------------------------------------------;
;--------------------------------------------------------------------;
;                                2PLOT1                              ;
;plot2:
if (1 eq 0) then begin
loadct, 39
if(plottops eq 1) then begin
   entry_device = !d.name
   set_plot,'PS'
   device, bits_per_pixel = 8,/portrait, color = 1,filename =filename31
   tk = 2
endif else begin
   window, 1, retain = 2
   tk = 1
endelse
plot, alog10(dwarfeu(0:nsat-1:50)), EuH_sat(0:nsat-1:100),psym=sym(0), xthick=2, ythick = 2, charthick = 2, xtitle = 'Dwarf fraction', ytitle='[Eu/H]', /nodata, xrange = [-10, 1], xstyle = 1, yrange = [-2, 1], charsize = 1.1, Color=cgColor('black'), background=cgColor('white')
oplot,alog10(dwarfeu(0:nsat-1:50)), EuH_sat(0:nsat-1:100),psym=sym(1), symsize= 0.4, color = 250

oplot, zbins, median_dfrac_EuH_sat, thick = tk+2, color = 0
oplot, zbins, percent_up_dfrac_EuH_sat, linestyle = 2, thick = tk+1, color = 0
oplot, zbins, percent_down_dfrac_EuH_sat, linestyle = 2, thick = tk+1, color = 0

if (plottops eq 1) then begin
   device, /close_file
   set_plot, entry_device
endif
endif

;--------------------------------------------------------------------;
;                                2PLOT2                              ;
if (1 eq 0) then begin
loadct, 39
if(plottops eq 1) then begin
   entry_device = !d.name
   set_plot,'PS'
   device, bits_per_pixel = 8,/portrait, color = 1,filename =filename32
   tk = 2
endif else begin
   window, 2, retain = 2
   tk = 1
endelse
plot, alog10(dwarfeu(0:nsat-1:50)), EuH_dsat(0:nsat-1:50),psym=sym(0), xthick=2, ythick = 2, charthick = 2, xtitle = 'Dwarf fraction', ytitle='[Eu/H]', /nodata, xrange = [-8, 1], xstyle = 1, yrange = [-4, 0.5], charsize = 1.1, Color=cgColor('black'), background=cgColor('white')
oplot,alog10(dwarfeu(0:nsat-1:50)), EuH_dsat(0:nsat-1:50),psym=sym(0), symsize= 0.4, color = 250

oplot, zbins, median_dfrac_EuH_dsat, thick = tk+2, color = 0
oplot, zbins, percent_up_dfrac_EuH_dsat, linestyle = 2, thick = tk+1, color = 0
oplot, zbins, percent_down_dfrac_EuH_dsat, linestyle = 2, thick = tk+1, color = 0

if (plottops eq 1) then begin
   device, /close_file
   set_plot, entry_device
endif
endif

;--------------------------------------------------------------------;
;                                2PLOT3                              ;
;plot:
if (1 eq 0) then begin
loadct, 39
if(plottops eq 1) then begin
   entry_device = !d.name
   set_plot,'PS'
   device, bits_per_pixel = 8,/portrait, color = 1,filename =filename33
   tk = 2
endif else begin
   window, 3, retain = 2
   tk = 1
endelse
plot, alog10(dir_dwarfeu(0:nhost-1:50)), EuH_host(0:nhost-1:100),psym=sym(0), xthick=2, ythick = 2, charthick = 2, xtitle = 'Direct dwarf fraction', ytitle='[Eu/H]', /nodata, xrange = [-10, 1], xstyle = 1, yrange = [-4, 2], charsize = 1.1, Color=cgColor('black'), background=cgColor('white')
oplot,alog10(dir_dwarfeu(0:nhost-1:50)), EuH_host(0:nhost-1:100),psym=sym(0), symsize= 0.4, color = 250

if (plottops eq 1) then begin
   device, /close_file
   set_plot, entry_device
endif
endif

;--------------------------------------------------------------------;
;                                2PLOT4                              ;
;plot: 
if (1 eq 0) then begin
loadct, 39
if(plottops eq 1) then begin
   entry_device = !d.name
   set_plot,'PS'
   device, bits_per_pixel = 8,/portrait, color = 1,filename =filename34
   tk = 2
endif else begin
   window, 4, retain = 2
   tk = 1
endelse
plot, alog10(dir_dwarfeu(0:nhost-1:50)), EuH_dhost(0:nhost-1:50),psym=sym(0), xthick=2, ythick = 2, charthick = 2, xtitle = 'Direct dwarf fraction', ytitle='[Eu/H] with dwarf contamination', /nodata, xrange = [-10, 1], xstyle = 1, yrange = [-14, 2], charsize = 1.1, Color=cgColor('black'), background=cgColor('white')
oplot,alog10(dir_dwarfeu(0:nhost-1:50)), EuH_dhost(0:nhost-1:50),psym=sym(0), symsize= 0.4, color = 250

if (plottops eq 1) then begin
   device, /close_file
   set_plot, entry_device
endif
endif

;--------------------------------------------------------------------;
;                                2PLOT5                              ;
;plot:
if (1 eq 0) then begin
loadct, 39
if(plottops eq 1) then begin
   entry_device = !d.name
   set_plot,'PS'
   device, bits_per_pixel = 8,/portrait, color = 1,filename =filename35
   tk = 2
endif else begin
   window, 5, retain = 2
   tk = 1
endelse
plot, alog10(dir_dwarfeu(0:nsat-1:50)), EuH_sat(0:nsat-1:50),psym=sym(0), xthick=2, ythick = 2, charthick = 2, xtitle = 'Direct dwarf fraction', ytitle='[Eu/H]', /nodata, xrange = [-1, 1], xstyle = 1, yrange = [-4, 2], charsize = 1.1, Color=cgColor('black'), background=cgColor('white')
oplot,alog10(dir_dwarfeu(0:nsat-1:50)), EuH_sat(0:nsat-1:50),psym=sym(0), symsize= 0.4, color = 250

if (plottops eq 1) then begin
   device, /close_file
   set_plot, entry_device
endif
endif

;--------------------------------------------------------------------;
;                                2PLOT6                              ;
;plot:
if (1 eq 0) then begin
loadct, 39
if(plottops eq 1) then begin
   entry_device = !d.name
   set_plot,'PS'
   device, bits_per_pixel = 8,/portrait, color = 1,filename =filename36
   tk = 2
endif else begin
   window, 6, retain = 2
   tk = 1
endelse
plot, alog10(dir_dwarfeu(0:nsat-1:50)), EuH_dsat(0:nsat-1:50),psym=sym(0), xthick=2, ythick = 2, charthick = 2, xtitle = 'Direct dwarf fraction', ytitle='[Eu/H]', /nodata, xrange = [-1, 1], xstyle = 1, yrange = [-14, 2], charsize = 1.1, Color=cgColor('black'), background=cgColor('white')
oplot,alog10(dir_dwarfeu(0:nsat-1:50)), EuH_dsat(0:nsat-1:50),psym=sym(0), symsize= 0.4, color = 250

if (plottops eq 1) then begin
   device, /close_file
   set_plot, entry_device
endif
endif

;--------------------------------------------------------------------;
;                                2PLOT7                              ;
;plot:
if (1 eq 1) then begin
loadct, 39
if(plottops eq 1) then begin
   entry_device = !d.name
   set_plot,'PS'
   device, bits_per_pixel = 8,/portrait, color = 1,filename = filename37
   tk = 2
endif else begin
   window, 7, retain = 2
   tk = 1
endelse
plot, alog10(tform_host(0:nhost-1:50)), EuFe_dhost(0:nhost-1:50),psym=sym(0), xthick=2, ythick=2, charthick=2, xtitle = 'log (Formation time [yr])', ytitle='[Eu/Fe]', /nodata, xrange = [8.6, 10.3], xstyle = 1, yrange = [-6, 1], charsize = 1.1, Color=cgColor('black'), background=cgColor('white')
oplot,alog10(tform_host(0:nhost-1:50)), EuFe_dhost(0:nhost-1:50),psym=sym(0), color = 250, symsize = 0.3

;oplot, logtbins, median_EuFe_age, thick = tk+1, color = 0
;oplot, logtbins, percent_down_EuFe_age, thick = tk+1, linestyle = 2, color = 0
;oplot, logtbins, percent_up_EuFe_age, thick = tk+1, linestyle = 2, color = 0

oplot, logtbins, median_EuFe_dhost_age, thick = tk+1, color = 0
oplot, logtbins, percent_up_EuFe_dhost_age, thick = tk+1, linestyle = 2, color = 0
oplot, logtbins, percent_down_EuFe_dhost_age, thick = tk+1, linestyle = 2, color = 0

;label=['merger + SN II, 128 Nbr']
;legend, label, box =0, charthick = tk, /bottom, /left

if (plottops eq 1) then begin
   device, /close_file
   set_plot, entry_device
endif
endif

;--------------------------------------------------------------------;
;                                2PLOT8                              ;
;plot:
if (1 eq 1) then begin
loadct, 39
if(plottops eq 1) then begin
   entry_device = !d.name
   set_plot,'PS'
   device, bits_per_pixel = 8,/portrait, color = 1,filename = filename38
   tk = 2
endif else begin
   window, 8, retain = 2
   tk = 1
endelse
plot, alog10(tform_host(0:nhost-1:50)), EuFe_dir_dhost(0:nhost-1:50),psym=sym(0), xthick=2, ythick=2, charthick=2, xtitle = 'log (Formation time [yr])', ytitle='[Eu/Fe]', /nodata, xrange = [8.6, 10.3], xstyle = 1, yrange = [-2, 2], charsize = 1.1, Color=cgColor('black'), background=cgColor('white')
oplot,alog10(tform_host(0:nhost-1:50)), EuFe_dir_dhost(0:nhost-1:50),psym=sym(0), color = 250, symsize = 0.3

;oplot, logtbins, median_EuFe_age, thick = tk+1, color = 0
;oplot, logtbins, percent_down_EuFe_age, thick = tk+1, linestyle = 2, color = 0
;oplot, logtbins, percent_up_EuFe_age, thick = tk+1, linestyle = 2, color = 0

oplot, logtbins, median_dir_EuFe_dhost_age, thick = tk+1, color = 0
oplot, logtbins, percent_up_EuFe_dir_dhost_age, thick = tk+1, linestyle = 2, color = 0
oplot, logtbins, percent_down_EuFe_dir_dhost_age, thick = tk+1, linestyle = 2, color = 0

;label=['merger + SN II, 128 Nbr']
;legend, label, box =0, charthick = tk, /bottom, /left

if (plottops eq 1) then begin
   device, /close_file
   set_plot, entry_device
endif
endif

;--------------------------------------------------------------------;
;                                2PLOT9                              ;
;plot:
if (1 eq 1) then begin
loadct, 39
if(plottops eq 1) then begin
   entry_device = !d.name
   set_plot,'PS'
   device, bits_per_pixel = 8,/portrait, color = 1,filename = filename39
   tk = 2
endif else begin
   window, 9, retain = 2
   tk = 1
endelse
plot, alog10(tform_sat(0:nsat-1)), EuFe_sat(0:nsat-1),psym=sym(0), xthick=2, ythick=2, charthick=2, xtitle = 'log (Formation time [yr])', ytitle='[Eu/Fe]', /nodata, xrange = [8.5, 10.3], xstyle = 1, yrange = [-2, 1], charsize = 1.1, Color=cgColor('black'), background=cgColor('white')
oplot,alog10(tform_sat(0:nsat-1)), EuFe_sat(0:nsat-1),psym=sym(0), color = 250, symsize = 0.3

;oplot, logtbins, median_EuFe_age, thick = tk+1, color = 0
;oplot, logtbins, percent_down_EuFe_age, thick = tk+1, linestyle = 2, color = 0
;oplot, logtbins, percent_up_EuFe_age, thick = tk+1, linestyle = 2, color = 0

oplot, logtbins, median_EuFe_age, thick = tk+1, color = 85
oplot, logtbins, median_EuFe_age_sat, thick = tk+1, color = 0
oplot, logtbins, percent_up_EuFe_age_sat, thick = tk+1, linestyle = 2, color = 0
oplot, logtbins, percent_down_EuFe_age_sat, thick = tk+1, linestyle = 2, color = 0

;label=['merger + SN II, 128 Nbr']
;legend, label, box =0, charthick = tk, /bottom, /left

if (plottops eq 1) then begin
   device, /close_file
   set_plot, entry_device
endif
endif

;--------------------------------------------------------------------;
;                               2PLOT10                              ;
;plot:
if (1 eq 1) then begin
loadct, 39
if(plottops eq 1) then begin
   entry_device = !d.name
   set_plot,'PS'
   device, bits_per_pixel = 8,/portrait, color = 1,filename = filename40
   tk = 2
endif else begin
   window, 10, retain = 2
   tk = 1
endelse
plot, alog10(tform_sat(0:nsat-1)), EuFe_dsat(0:nsat-1),psym=sym(0), xthick=2, ythick=2, charthick=2, xtitle = 'log (Formation time [yr])', ytitle='[Eu/Fe]', /nodata, xrange = [8.5, 10.3], xstyle = 1, yrange = [-6, 0.1], charsize = 1.1, Color=cgColor('black'), background=cgColor('white')
oplot,alog10(tform_sat(0:nsat-1)), EuFe_dsat(0:nsat-1),psym=sym(0), color = 250, symsize = 0.3

;oplot, logtbins, median_EuFe_age, thick = tk+1, color = 0
;oplot, logtbins, percent_down_EuFe_age, thick = tk+1, linestyle = 2, color = 0
;oplot, logtbins, percent_up_EuFe_age, thick = tk+1, linestyle = 2, color = 0

oplot, logtbins, median_EuFe_dhost_age, thick = tk+1, color = 85
oplot, logtbins, median_EuFe_dsat_age, thick = tk+1, color = 0
oplot, logtbins, percent_up_EuFe_dsat_age, thick = tk+1, linestyle = 2, color = 0
oplot, logtbins, percent_down_EuFe_dsat_age, thick = tk+1, linestyle = 2, color = 0

;label=['merger + SN II, 128 Nbr']
;legend, label, box =0, charthick = tk, /bottom, /left

if (plottops eq 1) then begin
   device, /close_file
   set_plot, entry_device
endif
endif

;--------------------------------------------------------------------;
;                               2PLOT11                              ;
;plot: 
if (1 eq 0) then begin
loadct, 39
if(plottops eq 1) then begin
   entry_device = !d.name
   set_plot,'PS'
   deviCe, bits_per_pixel = 8,/portrait, color = 1,filename = filename41
   tk = 2
endif else begin
   window, 11, retain = 2
   tk = 1
endelse
plot, alog10(tform_sat(0:nsat-1:50)), EuFe_dir_dsat(0:nsat-1:50),psym=sym(0), xthick=2, ythick=2, charthick=2, xtitle = 'log (Formation time [yr])', ytitle='[Eu/Fe]', /nodata, xrange = [8.5, 10.3], xstyle = 1, yrange = [-3, 2], charsize = 1.1, Color=cgColor('black'), background=cgColor('white')
oplot,alog10(tform_sat(0:nsat-1:50)), EuFe_dir_dsat(0:nsat-1:50),psym=sym(0), color = 250, symsize = 0.3

;oplot, logtbins, median_EuFe_age, thick = tk+1, color = 0
;oplot, logtbins, percent_down_EuFe_age, thick = tk+1, linestyle = 2, color = 0
;oplot, logtbins, percent_up_EuFe_age, thick = tk+1, linestyle = 2, color = 0

oplot, logtbins, median_dir_EuFe_dsat_age, thick = tk+1, color = 0
oplot, logtbins, percent_up_EuFe_dir_dsat_age, thick = tk+1, linestyle = 2, color = 0
oplot, logtbins, percent_down_EuFe_dir_dsat_age, thick = tk+1, linestyle = 2, color = 0

;label=['merger + SN II, 128 Nbr']
legend, label, box =0, charthick = tk, /bottom, /left

if (plottops eq 1) then begin
   device, /close_file
   set_plot, entry_device
endif
endif
;retall

;--------------------------------------------------------------------;
;--------------------------------------------------------------------;
;                               3PLOT1                               ;
;plot3:
if (1 eq 0) then begin
loadct, 39
if(plottops eq 1) then begin
   entry_device = !d.name
   set_plot,'PS'
   device, bits_per_pixel = 8,/portrait, color = 1,filename=filename42
   tk = 2
endif else begin
   window, 1, retain = 2
   tk = 1
endelse
plot, xgas(0:nhost-1:100), ygas(0:nhost-1:100),psym=sym(0), xthick=2, ythick = 2, charthick = 2, xtitle = 'Direct dwarf fraction', ytitle='[Eu/H]', /nodata, xrange = [-2000,2000], xstyle = 1, yrange = [-2000,2000], charsize = 1.1, Color=cgColor('black'), background=cgColor('white')
oplot,xgas(0:nhost-1:100), ygas(0:nhost-1:100),psym=sym(0), symsize= 0.4, color = 250

if (plottops eq 1) then begin
   device, /close_file
   set_plot, entry_device
endif
endif

;--------------------------------------------------------------------;
;                               3PLOT2                               ;
;plot: 
if (1 eq 0) then begin
loadct, 39
if(plottops eq 1) then begin
   entry_device = !d.name
   set_plot,'PS'
   device, bits_per_pixel = 8,/portrait, color = 1,filename =filename43
   tk = 2
endif else begin
   window, 2, retain = 2
   tk = 1
endelse
plot, xgas(0:nhost-1:100), ygas(0:nhost-1:100),psym=sym(0), xthick=2, ythick = 2, charthick = 2, xtitle = 'Direct dwarf fraction', ytitle='[Eu/H]', /nodata, xrange = [-1000,1000], xstyle = 1, yrange = [-1000,1000], charsize = 1.1, Color=cgColor('black'), background=cgColor('white')
oplot,xgas(0:nhost-1:100), ygas(0:nhost-1:100),psym=sym(0), symsize= 0.4, color = 250

if (plottops eq 1) then begin
   device, /close_file
   set_plot, entry_device
endif
endif

;--------------------------------------------------------------------;
;                               3PLOT3                               ;
;plot:
if (1 eq 0) then begin
loadct, 39
if(plottops eq 1) then begin
   entry_device = !d.name
   set_plot,'PS'
   device, bits_per_pixel = 8,/portrait, color = 1,filename =filename44
   tk = 2
endif else begin
   window, 3, retain = 2
   tk = 1
endelse
plot, xgas(0:nhost-1:100), ygas(0:nhost-1:100),psym=sym(0), xthick=2, ythick = 2, charthick = 2, xtitle = 'Direct dwarf fraction', ytitle='[Eu/H]', /nodata, xrange = [-500,500], xstyle = 1, yrange = [-500,500], charsize = 1.1, Color=cgColor('black'), background=cgColor('white')
oplot,xgas(0:nhost-1:100), ygas(0:nhost-1:100),psym=sym(0), symsize= 0.4, color = 250

if (plottops eq 1) then begin
   device, /close_file
   set_plot, entry_device
endif
endif

;--------------------------------------------------------------------;
;                               3PLOT4                               ;
;plot:
if (1 eq 0) then begin
loadct, 39
if(plottops eq 1) then begin
   entry_device = !d.name
   set_plot,'PS'
   device, bits_per_pixel = 8,/portrait, color = 1,filename =filename45
   tk = 2
endif else begin
   window, 4, retain = 2
   tk = 1
endelse
plot, xgas(0:nhost-1:100), ygas(0:nhost-1:100),psym=sym(0), xthick=2, ythick = 2, charthick = 2, xtitle = 'Direct dwarf fraction', ytitle='[Eu/H]', /nodata, xrange = [-200,200], xstyle = 1, yrange = [-200,200], charsize = 1.1, Color=cgColor('black'), background=cgColor('white')
oplot,xgas(0:nhost-1:100), ygas(0:nhost-1:100),psym=sym(0), symsize= 0.4, color = 250

if (plottops eq 1) then begin
   device, /close_file
   set_plot, entry_device
endif
endif

;--------------------------------------------------------------------;
;                               3PLOT5                               ;
;plot:
if (1 eq 0) then begin
loadct, 39
if(plottops eq 1) then begin
   entry_device = !d.name
   set_plot,'PS'
   device, bits_per_pixel = 8,/portrait, color = 1,filename =filename46
   tk = 2
endif else begin
   window, 5, retain = 2
   tk = 1
endelse
plot, xgas(0:nhost-1:100), ygas(0:nhost-1:100),psym=sym(0), xthick=2, ythick = 2, charthick = 2, xtitle = 'Direct dwarf fraction', ytitle='[Eu/H]', /nodata, xrange = [-100,100], xstyle = 1, yrange = [-100,100], charsize = 1.1, Color=cgColor('black'), background=cgColor('white')
oplot,xgas(0:nhost-1:100), ygas(0:nhost-1:100),psym=sym(0), symsize= 0.4, color = 250

if (plottops eq 1) then begin
   device, /close_file
   set_plot, entry_device
endif
endif

;--------------------------------------------------------------------;
;                               3PLOT6                               ;
;plot:
if (1 eq 0) then begin
loadct, 39
if(plottops eq 1) then begin
   entry_device = !d.name
   set_plot,'PS'
   device, bits_per_pixel = 8,/portrait, color = 1,filename =filename47
   tk = 2
endif else begin
   window, 6, retain = 2
   tk = 1
endelse
plot, xgas(0:nhost-1:100), zgas(0:nhost-1:100),psym=sym(0), xthick=2, ythick = 2, charthick = 2, xtitle = 'Direct dwarf fraction', ytitle='[Eu/H]', /nodata, xrange = [-3000,3000], xstyle = 1, yrange = [-2000,2000], charsize = 1.1, Color=cgColor('black'), background=cgColor('white')
oplot,xgas(0:nhost-1:100), zgas(0:nhost-1:100),psym=sym(0), symsize= 0.4, color = 250

if (plottops eq 1) then begin
   device, /close_file
   set_plot, entry_device
endif
endif

;--------------------------------------------------------------------;
;                               3PLOT7                               ;
;plot:
if (1 eq 0) then begin
loadct, 39
if(plottops eq 1) then begin
   entry_device = !d.name
   set_plot,'PS'
   device, bits_per_pixel = 8,/portrait, color = 1,filename =filename48
   tk = 2
endif else begin
   window, 7, retain = 2
   tk = 1
endelse
plot, xgas(0:nhost-1:100), zgas(0:nhost-1:100),psym=sym(0), xthick=2, ythick = 2, charthick = 2, xtitle = 'Direct dwarf fraction', ytitle='[Eu/H]', /nodata, xrange = [-1000,1000], xstyle = 1, yrange = [-1000,1000], charsize = 1.1, Color=cgColor('black'), background=cgColor('white')
oplot,xgas(0:nhost-1:100), zgas(0:nhost-1:100),psym=sym(0), symsize= 0.4, color = 250

if (plottops eq 1) then begin
   device, /close_file
   set_plot, entry_device
endif
endif

;--------------------------------------------------------------------;
;                               3PLOT8                               ;
;plot:
if (1 eq 0) then begin
loadct, 39
if(plottops eq 1) then begin
   entry_device = !d.name
   set_plot,'PS'
   device, bits_per_pixel = 8,/portrait, color = 1,filename =filename49
   tk = 2
endif else begin
   window, 8, retain = 2
   tk = 1
endelse
plot, xgas(0:nhost-1:100), zgas(0:nhost-1:100),psym=sym(0), xthick=2, ythick = 2, charthick = 2, xtitle = 'Direct dwarf fraction', ytitle='[Eu/H]', /nodata, xrange = [-500,500], xstyle = 1, yrange = [-500,500], charsize = 1.1, Color=cgColor('black'), background=cgColor('white')
oplot,xgas(0:nhost-1:100), zgas(0:nhost-1:100),psym=sym(0), symsize= 0.4, color = 250

if (plottops eq 1) then begin
   device, /close_file
   set_plot, entry_device
endif
endif

;--------------------------------------------------------------------;
;                               3PLOT9                               ;
;plot:
if (1 eq 0) then begin
loadct, 39
if(plottops eq 1) then begin
   entry_device = !d.name
   set_plot,'PS'
   device, bits_per_pixel = 8,/portrait, color = 1,filename =filename50
   tk = 2
endif else begin
   window, 9, retain = 2
   tk = 1
endelse
plot, xgas(0:nhost-1:100), zgas(0:nhost-1:100),psym=sym(0), xthick=2, ythick = 2, charthick = 2, xtitle = 'Direct dwarf fraction', ytitle='[Eu/H]', /nodata, xrange = [-200,200], xstyle = 1, yrange = [-200,200], charsize = 1.1, Color=cgColor('black'), background=cgColor('white')
oplot,xgas(0:nhost-1:100), zgas(0:nhost-1:100),psym=sym(0), symsize= 0.4, color = 250

if (plottops eq 1) then begin
   device, /close_file
   set_plot, entry_device
endif
endif

;--------------------------------------------------------------------;
;                              3PLOT10                               ;
;plot:
if (1 eq 0) then begin
loadct, 39
if(plottops eq 1) then begin
   entry_device = !d.name
   set_plot,'PS'
   device, bits_per_pixel = 8,/portrait, color = 1,filename =filename51
   tk = 2
endif else begin
   window, 10, retain = 2
   tk = 1
endelse
plot, xgas(0:nhost-1:100), zgas(0:nhost-1:100),psym=sym(0), xthick=2, ythick = 2, charthick = 2, xtitle = 'Direct dwarf fraction', ytitle='[Eu/H]', /nodata, xrange = [-100,100], xstyle = 1, yrange = [-100,100], charsize = 1.1, Color=cgColor('black'), background=cgColor('white')
oplot,xgas(0:nhost-1:100), zgas(0:nhost-1:100),psym=sym(0), symsize= 0.4, color = 250

if (plottops eq 1) then begin
   device, /close_file
   set_plot, entry_device
endif
endif

;--------------------------------------------------------------------;
;                              3PLOT11                               ;
;plot:
if (1 eq 0) then begin
loadct, 39
if(plottops eq 1) then begin
   entry_device = !d.name
   set_plot,'PS'
   device, bits_per_pixel = 8,/portrait, color = 1,filename =filename52
   tk = 2
endif else begin
   window, 11, retain = 2
   tk = 1
endelse
plot, xstar(0:nhost-1:100), ystar(0:nhost-1:100),psym=sym(0), xthick=2, ythick = 2, charthick = 2, xtitle = 'Direct dwarf fraction', ytitle='[Eu/H]', /nodata, xrange = [-500,500], xstyle = 1, yrange = [-500,500], charsize = 1.1, Color=cgColor('black'), background=cgColor('white')
oplot,xstar(0:nhost-1:100), ystar(0:nhost-1:100),psym=sym(0), symsize= 0.4, color = 250

if (plottops eq 1) then begin
   device, /close_file
   set_plot, entry_device
endif
endif

;--------------------------------------------------------------------;
;                               3PLOT12                              ;
;plot:
if (1 eq 0) then begin
loadct, 39
if(plottops eq 1) then begin
   entry_device = !d.name
   set_plot,'PS'
   device, bits_per_pixel = 8,/portrait, color = 1,filename =filename53
   tk = 2
endif else begin
   window, 12, retain = 2
   tk = 1
endelse
plot, xstar(0:nhost-1:100), ystar(0:nhost-1:100),psym=sym(0), xthick=2, ythick = 2, charthick = 2, xtitle = 'Direct dwarf fraction', ytitle='[Eu/H]', /nodata, xrange = [-200,200], xstyle = 1, yrange = [-200,200], charsize = 1.1, Color=cgColor('black'), background=cgColor('white')
oplot,xstar(0:nhost-1:100), ystar(0:nhost-1:100),psym=sym(0), symsize= 0.4, color = 250

if (plottops eq 1) then begin
   device, /close_file
   set_plot, entry_device
endif
endif

;--------------------------------------------------------------------;
;                               3PLOT13                              ;
;plot:
if (1 eq 0) then begin
loadct, 39
if(plottops eq 1) then begin
   entry_device = !d.name
   set_plot,'PS'
   device, bits_per_pixel = 8,/portrait, color = 1,filename =filename54
   tk = 2
endif else begin
   window, 13, retain = 2
   tk = 1
endelse
plot, xstar(0:nhost-1:100), ystar(0:nhost-1:100),psym=sym(0), xthick=2, ythick = 2, charthick = 2, xtitle = 'Direct dwarf fraction', ytitle='[Eu/H]', /nodata, xrange = [-100,100], xstyle = 1, yrange = [-100,100], charsize = 1.1, Color=cgColor('black'), background=cgColor('white')
oplot,xstar(0:nhost-1:100), ystar(0:nhost-1:100),psym=sym(0), symsize= 0.4, color = 250

if (plottops eq 1) then begin
   device, /close_file
   set_plot, entry_device
endif
endif

;--------------------------------------------------------------------;
;                               3PLOT14                              ;
;plot:
if (1 eq 0) then begin
loadct, 39
if(plottops eq 1) then begin
   entry_device = !d.name
   set_plot,'PS'
   device, bits_per_pixel = 8,/portrait, color = 1,filename =filename55
   tk = 2
endif else begin
   window, 14, retain = 2
   tk = 1
endelse
plot, xstar(0:nhost-1:100), zstar(0:nhost-1:100),psym=sym(0), xthick=2, ythick = 2, charthick = 2, xtitle = 'Direct dwarf fraction', ytitle='[Eu/H]', /nodata, xrange = [-500,500], xstyle = 1, yrange = [-500,500], charsize = 1.1, Color=cgColor('black'), background=cgColor('white')
oplot,xstar(0:nhost-1:100), zstar(0:nhost-1:100),psym=sym(0), symsize= 0.4, color = 250

if (plottops eq 1) then begin
   device, /close_file
   set_plot, entry_device
endif
endif

;--------------------------------------------------------------------;
;                               3PLOT15                              ;
;plot:
if (1 eq 0) then begin
loadct, 39
if(plottops eq 1) then begin
   entry_device = !d.name
   set_plot,'PS'
   device, bits_per_pixel = 8,/portrait, color = 1,filename =filename56
   tk = 2
endif else begin
   window, 15, retain = 2
   tk = 1
endelse
plot, xstar(0:nhost-1:100), zstar(0:nhost-1:100),psym=sym(0), xthick=2, ythick = 2, charthick = 2, xtitle = 'Direct dwarf fraction', ytitle='[Eu/H]', /nodata, xrange = [-200,200], xstyle = 1, yrange = [-200,200], charsize = 1.1, Color=cgColor('black'), background=cgColor('white')
oplot,xstar(0:nhost-1:100), zstar(0:nhost-1:100),psym=sym(0), symsize= 0.4, color = 250

if (plottops eq 1) then begin
   device, /close_file
   set_plot, entry_device
endif
endif

;--------------------------------------------------------------------;
;                              3PLOT16                               ;
;plot:
if (1 eq 0) then begin
loadct, 39
if(plottops eq 1) then begin
   entry_device = !d.name
   set_plot,'PS'
   device, bits_per_pixel = 8,/portrait, color = 1,filename =filename57
   tk = 2
endif else begin
   window, 15, retain = 2
   tk = 1
endelse
plot, xstar(0:nhost-1:100), zstar(0:nhost-1:100),psym=sym(0), xthick=2, ythick = 2, charthick = 2, xtitle = 'Direct dwarf fraction', ytitle='[Eu/H]', /nodata, xrange = [-100,100], xstyle = 1, yrange = [-100,100], charsize = 1.1, Color=cgColor('black'), background=cgColor('white')
oplot,xstar(0:nhost-1:100), zstar(0:nhost-1:100),psym=sym(0), symsize= 0.4, color = 250

if (plottops eq 1) then begin
   device, /close_file
   set_plot, entry_device
endif
endif
;retall

;--------------------------------------------------------------------;
;--------------------------------------------------------------------;
;                               4PLOT1                               ;
;plot4:
if (1 eq 0) then begin
loadct, 39
if(plottops eq 1) then begin
   entry_device = !d.name
   set_plot,'PS'
   device, bits_per_pixel = 8,/portrait, color = 1,filename = filename58
   tk = 2
endif else begin
   window, 1, retain = 1
   tk = 1
endelse
plot, EuFe_host(0:nhost-1:100), EuH_host(0:nhost-1:100), psym = sym(0), xthick=2, ythick = 2, charthick = 2, xtitle = '[Fe/H]', ytitle='[Eu/Fe] >= 0.95', /nodata, xrange = [-4, 0.5], charsize = 1.1, yrange = [-4, 4], Color=cgColor('black'), background=cgColor('white')
oplot,EuFe_host(0:nhost-1:100), EuH_host(0:nhost-1:100), psym = sym(0), symsize= 0.4, color = 250

;oplot, zbins, median_EuFe, thick = tk+2, color = 0
;oplot, zbins, percent_down_EuFe, linestyle = 2, thick = tk+1, color = 0
;oplot, zbins, percent_up_EuFe, linestyle = 2, thick = tk+1, color = 0

oplot, zbins, median_EuFe_EuH, thick = tk+2, color = 0
oplot, zbins, percent_up_EuFe_EuH, linestyle = 2, thick = tk+1, color = 0
oplot, zbins, percent_down_EuFe_EuH, linestyle = 2, thick = tk+1, color = 0

legend, label, box =0, charthick = tk, /bottom, /left, color = 0

if (plottops eq 1) then begin
   device, /close_file
   set_plot, entry_device
endif
endif

;--------------------------------------------------------------------;
;                               4PLOT2                               ;
;plot:
if (1 eq 0) then begin
loadct, 39
if(plottops eq 1) then begin
   entry_device = !d.name
   set_plot,'PS'
   device, bits_per_pixel = 8,/portrait, color = 1,filename = filename59
   tk = 2
endif else begin
   window, 2, retain = 1
   tk = 1
endelse
plot, EuFe_dhost(0:nhost-1:100), EuH_dhost(0:nhost-1:100), psym = sym(0), xthick=2, ythick = 2, charthick = 2, xtitle = '[Fe/H]', ytitle='[Eu/Fe] >= 0.95', /nodata, xrange = [-4, 0.5], charsize = 1.1, yrange = [-4, 4], Color=cgColor('black'), background=cgColor('white')
oplot,EuFe_dhost(0:nhost-1:100), EuH_dhost(0:nhost-1:100), psym = sym(0), symsize= 0.4, color = 250

;oplot, zbins, median_EuFe, thick = tk+2, color = 0
;oplot, zbins, percent_down_EuFe, linestyle = 2, thick = tk+1, color = 0
;oplot, zbins, percent_up_EuFe, linestyle = 2, thick = tk+1, color = 0

oplot, zbins, median_EuFe_EuH_dhost, thick = tk+2, color = 0
oplot, zbins, percent_up_EuFe_EuH_dhost, linestyle = 2, thick = tk+1, color = 0
oplot, zbins, percent_down_EuFe_EuH_dhost, linestyle = 2, thick = tk+1, color = 0

legend, label, box =0, charthick = tk, /bottom, /left, color = 0

if (plottops eq 1) then begin
   device, /close_file
   set_plot, entry_device
endif
endif

;--------------------------------------------------------------------;
;                               4PLOT3                               ;
;plot:
if (1 eq 0) then begin
loadct, 39
if(plottops eq 1) then begin
   entry_device = !d.name
   set_plot,'PS'
   device, bits_per_pixel = 8,/portrait, color = 1,filename = filename60
   tk = 2
endif else begin
   window, 3, retain = 1
   tk = 1
endelse
plot, EuFe_sat(0:nsat-1:100), EuH_sat(0:nsat-1:100), psym = sym(0), xthick=2, ythick = 2, charthick = 2, xtitle = '[Fe/H]', ytitle='[Eu/Fe] >= 0.95', /nodata, xrange = [-4, 0.5], charsize = 1.1, yrange = [-4, 4], Color=cgColor('black'), background=cgColor('white')
oplot,EuFe_sat(0:nsat-1:100), EuH_sat(0:nsat-1:100), psym = sym(0), symsize= 0.4, color = 250

;oplot, zbins, median_EuFe, thick = tk+2, color = 0
;oplot, zbins, percent_down_EuFe, linestyle = 2, thick = tk+1, color = 0
;oplot, zbins, percent_up_EuFe, linestyle = 2, thick = tk+1, color = 0

oplot, zbins, median_EuFe_EuH_sat, thick = tk+2, color = 0
oplot, zbins, percent_up_EuFe_EuH_sat, linestyle = 2, thick = tk+1, color = 0
oplot, zbins, percent_down_EuFe_EuH_sat, linestyle = 2, thick = tk+1, color = 0

legend, label, box =0, charthick = tk, /bottom, /left, color = 0

if (plottops eq 1) then begin
   device, /close_file
   set_plot, entry_device
endif
endif

;--------------------------------------------------------------------;
;                               4PLOT4                               ;
;plot:
if (1 eq 0) then begin
loadct, 39
if(plottops eq 1) then begin
   entry_device = !d.name
   set_plot,'PS'
   device, bits_per_pixel = 8,/portrait, color = 1,filename = filename61
   tk = 2
endif else begin
   window, 4, retain = 1
   tk = 1
endelse
plot, EuFe_dsat(0:nsat-1:100), EuH_dsat(0:nsat-1:100), psym = sym(0), xthick=2, ythick = 2, charthick = 2, xtitle = '[Fe/H]', ytitle='[Eu/Fe] >= 0.95', /nodata, xrange = [-4, 0.5], charsize = 1.1, yrange = [-4, 4], Color=cgColor('black'), background=cgColor('white')
oplot,EuFe_dsat(0:nsat-1:100), EuH_dsat(0:nsat-1:100), psym = sym(0), symsize= 0.4, color = 250

;oplot, zbins, median_EuFe, thick = tk+2, color = 0
;oplot, zbins, percent_down_EuFe, linestyle = 2, thick = tk+1, color = 0
;oplot, zbins, percent_up_EuFe, linestyle = 2, thick = tk+1, color = 0

oplot, zbins, median_EuFe_EuH_dsat, thick = tk+2, color = 0
oplot, zbins, percent_up_EuFe_EuH_dsat, linestyle = 2, thick = tk+1, color = 0
oplot, zbins, percent_down_EuFe_EuH_dsat, linestyle = 2, thick = tk+1, color = 0

legend, label, box =0, charthick = tk, /bottom, /left, color = 0

if (plottops eq 1) then begin
   device, /close_file
   set_plot, entry_device
endif
endif
;retall

;--------------------------------------------------------------------;
;^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^;
;^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^;
plot5:
;PERCENTILE RANGE
heu_idx_per = where(dwarfeu(ihost) gt 0.85)
nn = n_elements(heu_idx_per)

EuFe_host_per = alog10(seu_host(heu_idx_per)/sfe_host(heu_idx_per)*(Afe/Aeu)) - Eufe_solar
EuH_host_per = alog10(seu_host(heu_idx_per)/(smass_host(heu_idx_per)*fh)*(Ah/Aeu)) - EuH_solar

heu_sat_idx_per = where(dwarfeu(isat) gt 0.85)
nns = n_elements(heu_sat_idx_per)

EuFe_sat_per = alog10(seu_sat(heu_sat_idx_per)/sfe_host(heu_sat_idx_per)*(Afe/Aeu)) - Eufe_solar
EuH_sat_per = alog10(seu_sat(heu_sat_idx_per)/(smass_host(heu_sat_idx_per)*fh)*(Ah/Aeu)) - EuH_solar

indexfe = where(FeH_host(heu_idx_per) lt 0.0)
nnfe = n_elements(indexfe)


;MORE BINNING
logtmin85=8.3
logtmax85=10.3
nbins85=21
;age_EuFe
binning, tform_host(heu_idx_per), EuFe_host_per, seu_host(heu_idx_per), logtmin85, logtmax85, nbins85, logtbins, median_EuFe_age_85, percent_up_EuFe_age_85, percent_down_EuFe_age_85, number85, /log, /nozero
;age_EuFe_sat
binning, tform_sat(heu_sat_idx_per), EuFe_sat_per, seu_sat(heu_sat_idx_per), logtmin85, logtmax85, nbins85, logtbins, median_EuFe_age_sat_85, percent_up_EuFe_age_sat_85, percent_down_EuFe_age_sat_85, numbers85, /log, /nozero

zmin85=-4.0
zmax85=1.0
nzbins85=51
;FeH_EuF
binning, FeH_host(heu_idx_per), EuFe_host_per, seu_host(heu_idx_per), zmin85, zmax85, nzbins85, zbins, median_EuFe_85, percent_up_EuFe_85, percent_down_EuFe_85, number185,/nozero
;FeH_EuFe_sat
binning, FeH_sat(heu_sat_idx_per), EuFe_sat_per, seu_sat(heu_sat_idx_per), zmin85, zmax85, nzbins85, zbins, median_EuFe_sat_85, percent_up_EuFe_sat_85, percent_down_EuFe_sat_85, number1s85, /nozero

;CHANGE BINNING FOR PERCENTILE PLOTS
filename62 = "62Age_EuH_host_heu_75_idx1.eps"
filename63 = "63Age_EuFe_host_heu_75_idx1.eps"
filename64 = "64FeH_EuFe_host_heu_75_idx1.eps"
filename65 = "65Age_EuH_sat_heu_75_idx1.eps"
filename66 = "66Age_EuFe_sat_heu_75_idx1.eps"
filename67 = "67FeH_EuFe_sat_heu_75_idx1.eps"
filename68 = "68Age_EuH_host_heu_85_idx1.eps"
filename69 = "69Age_EuFe_host_heu_85_idx1.eps"
filename70 = "70FeH_EuFe_host_heu_85_idx1.eps"
filename71 = "71Age_EuH_sat_heu_85_idx1.eps"
filename72 = "72Age_EuFe_sat_heu_85_idx1.eps"
filename73 = "73FeH_EuFe_sat_heu_85_idx1.eps"
filename74 = "74Age_EuH_host_heu_95_idx1.eps"
filename75 = "75Age_EuFe_host_heu_95_idx1.eps"
filename76 = "76FeH_EuFe_host_heu_95_idx1.eps"
filename77 = "77Age_EuH_sat_heu_95_idx1.eps"
filename78 = "78Age_EuFe_sat_heu_95_idx1.eps"
filename79 = "79FeH_EuFe_sat_heu_95_idx1.eps"

;^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^;
;--------------------------------------------------------------------;
;                               5PLOT5                               ;
if (1 eq 0) then begin
loadct, 39
if(plottops eq 1) then begin
   entry_device = !d.name
   set_plot,'PS'
   device, bits_per_pixel = 8,/portrait, color = 1,filename = filename62
   tk = 2
endif else begin
   window, 5, retain = 1
   tk = 1
endelse

plot, alog10(tform_host(heu_idx_per)), EuH_host_per, psym = sym(0), xthick=2, ythick = 2, charthick = 2, xtitle = 'time', ytitle='[Eu/H] with dwarf >= 0.95', /nodata, xrange = [8.5, 10], charsize = 1.1, yrange = [-4, 4], Color=cgColor('black'), background=cgColor('white')
oplot,alog10(tform_host(heu_idx_per)), EuH_host_per, psym = sym(0), symsize= 0.4, color = 250

;oplot, zbins, median_EuFe, thick = tk+2, color = 0
;oplot, zbins, percent_down_EuFe, linestyle = 2, thick = tk+1, color = 0
;oplot, zbins, percent_up_EuFe, linestyle = 2, thick = tk+1, color = 0

oplot, logtbins, median_EuH_age, thick = tk+1, color = 0
oplot, logtbins, percent_up_EuH_age, thick = tk+1, linestyle = 2, color = 0
oplot, logtbins, percent_down_EuH_age, thick = tk+1, linestyle = 2, color = 0

legend, label, box =0, charthick = tk, /bottom, /left, color = 0

if (plottops eq 1) then begin
   device, /close_file
   set_plot, entry_device
endif
endif

;^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^;
;                               5PLOT6                               ;
if (1 eq 0) then begin
loadct, 39
if(plottops eq 1) then begin
   entry_device = !d.name
   set_plot,'PS'
   device, bits_per_pixel = 8,/portrait, color = 1,filename = filename63
   tk = 2
endif else begin
   window, 6, retain = 1
   tk = 1
endelse

plot, alog10(tform_host(heu_idx_per)), EuFe_host_per, psym = sym(0), xthick=2, ythick = 2, charthick = 2, xtitle = 'time', ytitle='[Eu/Fe] dwarf frac >= 0.95', /nodata, xrange = [8.5, 10], charsize = 1.1, yrange = [-4, 4], Color=cgColor('black'), background=cgColor('white')
oplot,alog10(tform_host(heu_idx_per)), EuFe_host_per, psym = sym(0), symsize= 0.4, color = 250

;oplot, zbins, median_EuFe, thick = tk+2, color = 0                                       
;oplot, zbins, percent_down_EuFe, linestyle = 2, thick = tk+1, color =0                   
;oplot, zbins, percent_up_EuFe, linestyle = 2, thick = tk+1, color = 0                    

oplot, logtbins, median_EuFe_age, thick = tk+1, color = 0
oplot, logtbins, percent_up_EuFe_age, thick = tk+1, linestyle = 2, color = 0
oplot, logtbins, percent_down_EuFe_age, thick = tk+1, linestyle = 2, color = 0

legend, label, box =0, charthick = tk, /bottom, /left, color = 0

if (plottops eq 1) then begin
   device, /close_file
   set_plot, entry_device
endif
endif

;^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^;
;                               5PLOT7                               ;
if (1 eq 0) then begin
loadct, 39
if(plottops eq 1) then begin
   entry_device = !d.name
   set_plot,'PS'
   device, bits_per_pixel = 8,/portrait, color = 1,filename = filename64
   tk = 2
endif else begin
   window, 7, retain = 1
   tk = 1
endelse
                                                                
plot, FeH_host(heu_idx_per), EuFe_host_per, psym = sym(0), xthick=2, ythick = 2, charthick = 2, xtitle = '[Fe/H]', ytitle='[Eu/Fe] dwarf frac >= 0.95', /nodata, xrange = [-4, 1], charsize = 1.1, yrange = [-4, 4], Color=cgColor('black'), background=cgColor('white')
oplot,FeH_host(heu_idx_per), EuFe_host_per, psym = sym(0), symsize= 0.4, color = 250

;oplot, zbins, median_EuFe, thick = tk+2, color = 0                                       
;oplot, zbins, percent_down_EuFe, linestyle = 2, thick = tk+1, color =0                                                         
;oplot, zbins, percent_up_EuFe, linestyle = 2, thick = tk+1, color = 0        
oplot, zbins, median_EuFe, thick = tk+2, color = 0
oplot, zbins, percent_up_EuFe, linestyle = 2, thick = tk+1, color = 0
oplot, zbins, percent_down_EuFe, linestyle = 2, thick = tk+1, color = 0

legend, label, box =0, charthick = tk, /bottom, /left, color = 0

if (plottops eq 1) then begin
   device, /close_file
   set_plot, entry_device
endif
endif

;^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^;
;                               5PLOT8                               ;
if (1 eq 0) then begin
loadct, 39
if(plottops eq 1) then begin
   entry_device = !d.name
   set_plot,'PS'
   device, bits_per_pixel = 8,/portrait, color = 1,filename = filename65
   tk = 2
endif else begin
   window, 8, retain = 1
   tk = 1
endelse
                                                    
plot, alog10(tform_sat(heu_sat_idx_per)), EuH_sat_per, psym = sym(0), xthick=2, ythick = 2, charthick = 2, xtitle = 'time', ytitle='[Eu/H] with dwarf >= 0.95', /nodata, xrange = [8.5, 10], charsize = 1.1, yrange = [-4, 4], Color=cgColor('black'), background=cgColor('white')
oplot,alog10(tform_sat(heu_sat_idx_per)), EuH_sat_per, psym = sym(0), symsize= 0.4, color = 250

;oplot, zbins, median_EuFe, thick = tk+2, color = 0                                       
;oplot, zbins, percent_down_EuFe, linestyle = 2, thick = tk+1, color =0                   
;oplot, zbins, percent_up_EuFe, linestyle = 2, thick = tk+1, color = 0                    

oplot, logtbins, median_EuH_age_sat, thick = tk+1, color = 0
oplot, logtbins, percent_up_EuH_age_sat, thick = tk+1, linestyle = 2, color = 0
oplot, logtbins, percent_down_EuH_age_sat, thick = tk+1, linestyle = 2, color = 0


legend, label, box =0, charthick = tk, /bottom, /left, color = 0

if (plottops eq 1) then begin
   device, /close_file
   set_plot, entry_device
endif
endif


;^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^;
;                               5PLOT9                               ;
if (1 eq 0) then begin
loadct, 39
if(plottops eq 1) then begin
   entry_device = !d.name
   set_plot,'PS'
   device, bits_per_pixel = 8,/portrait, color = 1,filename = filename66
   tk = 2
endif else begin
   window, 9, retain = 1
   tk = 1
endelse

plot, alog10(tform_sat(heu_sat_idx_per)), EuFe_sat_per, psym = sym(0), xthick=2, ythick = 2, charthick = 2, xtitle = 'log (ForMation time [yr])', ytitle='[Eu/Fe] dwarf frac >= 0.85', /nodata, xrange = [8.5, 10.3], charsize = 1.1, yrange = [-2, 4], Color=cgColor('black'), background=cgColor('white')
oplot,alog10(tform_sat(heu_sat_idx_per)), EuFe_sat_per, psym = sym(0), symsize= 0.4, color = 250

;oplot, zbins, median_EuFe, thick = tk+2, color = 0                           
;oplot, zbins, percent_down_EuFe, linestyle = 2, thick = tk+1, color =0          
;oplot, zbins, percent_up_EuFe, linestyle = 2, thick = tk+1, color = 0          

oplot, logtbins, median_EuFe_age, thick = tk+1, color = 155
oplot, logtbins, median_EuFe_age_sat, thick = tk+1, color = 85
oplot, logtbins, percent_up_EuFe_age_sat, thick = tk+1, linestyle = 2, color = 85
oplot, logtbins, percent_down_EuFe_age_sat, thick = tk+1, linestyle = 2, color = 85

oplot, logtbins, median_EuFe_age_sat_85, thick = tk+1, color = 0
oplot, logtbins, percent_up_EuFe_age_sat_85, thick = tk+1, linestyle = 2, color = 0
oplot, logtbins, percent_down_EuFe_age_sat_85, thick = tk+1, linestyle = 2, color = 0

;legend, label, box =0, charthick = tk, /bottom, /left, color = 0

if (plottops eq 1) then begin
   device, /close_file
   set_plot, entry_device
endif
endif

;^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^;
;                               5PLOT10                              ;
if (1 eq 0) then begin
loadct, 39
if(plottops eq 1) then begin
   entry_device = !d.name
   set_plot,'PS'
   device, bits_per_pixel = 8,/portrait, color = 1,filename = filename67
   tk = 2
endif else begin
   window, 10, retain = 1
   tk = 1
endelse

plot, FeH_sat(heu_sat_idx_per), EuFe_sat_per, psym = sym(0), xthick=2, ythick = 2, charthick = 2, xtitle = '[Fe/H]', ytitle='[Eu/Fe] dwarf frac >= 0.85', /nodata, xrange = [-2.9, 0], charsize = 1.1, yrange = [-2.5, 4], Color=cgColor('black'), background=cgColor('white')
oplot,FeH_sat(heu_sat_idx_per), EuFe_sat_per, psym = sym(0), symsize= 0.4, color = 250

;oplot, zbins, median_EuFe, thick = tk+2, color = 0                            
;oplot, zbins, percent_down_EuFe, linestyle = 2, thick = tk+1, color =0          
;oplot, zbins, percent_up_EuFe, linestyle = 2, thick = tk+1, color = 0          

oplot, zbins, median_EuFe_sat_85, thick = tk+2, color = 0
oplot, zbins, percent_up_EuFe_sat_85, linestyle = 2, thick = tk+1, color = 0
oplot, zbins, percent_down_EuFe_sat_85, linestyle = 2, thick = tk+1, color = 0

;legend, label, box =0, charthick = tk, /bottom, /left, color = 0

if (plottops eq 1) then begin
   device, /close_file
   set_plot, entry_device
endif
endif

;--------------------------------------------------------------------;         
;                               5PLOT11                              ;
if (1 eq 0) then begin
loadct, 39
if(plottops eq 1) then begin
   entry_device = !d.name
   set_plot,'PS'
   device, bits_per_pixel = 8,/portrait, color = 1,filename = filename68
   tk = 2
endif else begin
   window, 11, retain = 1
   tk = 1
endelse

plot, alog10(tform_host(heu_idx_per)), EuH_host_per, psym = sym(0), xthick=2, ythick = 2, charthick = 2, xtitle = 'time', ytitle='[Eu/H] with dwarf >= 0.95', /nodata, xrange = [8.5, 10], charsize = 1.1, yrange = [-4, 4], Color=cgColor('black'), background=cgColor('white')
oplot,alog10(tform_host(heu_idx_per)), EuH_host_per, psym = sym(0), symsize= 0.4, color = 250

oplot, logtbins, median_EuH_age, thick = tk+1, color = 0
oplot, logtbins, percent_up_EuH_age, thick = tk+1, linestyle = 2, color = 0
oplot, logtbins, percent_down_EuH_age, thick = tk+1, linestyle = 2, color = 0

legend, label, box =0, charthick = tk, /bottom, /left, color = 0

if (plottops eq 1) then begin
   device, /close_file
   set_plot, entry_device
endif
endif

;^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^;
;                               5PLOT12                              ;
if (1 eq 0) then begin
loadct, 39
if(plottops eq 1) then begin
   entry_device = !d.name
   set_plot,'PS'
   device, bits_per_pixel = 8,/portrait, color = 1,filename = filename69
   tk = 2
endif else begin
   window, 12, retain = 1
   tk = 1
endelse

plot, alog10(tform_host(heu_idx_per)), EuFe_host_per, psym = sym(0), xthick=2, ythick = 2, charthick = 2, xtitle = 'log (ForMation time [yr])', ytitle='[Eu/Fe] dwarf frac >= 0.85', /nodata, xrange = [8.5, 10], charsize = 1.1, yrange = [-2, 2], Color=cgColor('black'), background=cgColor('white')
oplot,alog10(tform_host(heu_idx_per)), EuFe_host_per, psym = sym(0), symsize= 0.4, color = 250

oplot, logtbins, median_EuFe_age, thick = tk+1, color = 85
oplot, logtbins, percent_up_EuFe_age, thick = tk+1, linestyle = 4, color = 85
oplot, logtbins, percent_down_EuFe_age, thick = tk+1, linestyle = 4, color = 85

oplot, logtbins, median_EuFe_age_85, thick = tk+1, color = 0
oplot, logtbins, percent_up_EuFe_age_85, thick = tk+1, linestyle = 2, color = 0
oplot, logtbins, percent_down_EuFe_age_85, thick = tk+1, linestyle = 2, color = 0

;legend, label, box =0, charthick = tk, /bottom, /left, color = 0

if (plottops eq 1) then begin
   device, /close_file
   set_plot, entry_device
endif
endif

;^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^;
;                               5PLOT13                              ;
if (1 eq 0) then begin
loadct, 39
if(plottops eq 1) then begin
   entry_device = !d.name
   set_plot,'PS'
   device, bits_per_pixel = 8,/portrait, color = 1,filename = filename70;'70.0FeH_EuFe_host_heu_85_idx.eps'
   tk = 2
endif else begin
   window, 13, retain = 1
   tk = 1
endelse

;     FeH_host(indexfe), EuFe_host_per(0:nnfe-1)
plot, FeH_host(heu_idx_per), EuFe_host_per, psym = sym(0), xthick=2, ythick = 2, charthick = 2, xtitle = '[Fe/H]', ytitle='[Eu/Fe] dwarf frac >= 0.85', /nodata, xrange = [-4, 1], charsize = 1.1, yrange = [-2.8, 2.8], Color=cgColor('black'), background=cgColor('white')
oplot,FeH_host(heu_idx_per), EuFe_host_per, psym = sym(0), symsize= 0.4, color = 250

oplot, zbins, median_EuFe_85, thick = tk+2, color = 0
oplot, zbins, percent_up_EuFe_85, linestyle = 2, thick = tk+1, color = 0
oplot, zbins, percent_down_EuFe_85, linestyle = 2, thick = tk+1, color = 0

;legend, label, box =0, charthick = tk, /bottom, /left, color = 0

if (plottops eq 1) then begin
   device, /close_file
   set_plot, entry_device
endif
endif

;^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^;        
;                               5PLOT14                              ;            
if (1 eq 0) then begin
loadct, 39
if(plottops eq 1) then begin
   entry_device = !d.name
   set_plot,'PS'
   device, bits_per_pixel = 8,/portrait, color = 1,filename = filename71
   tk = 2
endif else begin
   window, 14, retain = 1
   tk = 1
endelse

plot, alog10(tform_sat(heu_sat_idx_per)), EuH_sat_per, psym = sym(0), xthick=2, ythick = 2, charthick = 2, xtitle = 'time', ytitle='[Eu/H] with dwarf >= 0.95', /nodata, xrange = [8.5, 10], charsize = 1.1, yrange = [-4, 4], Color=cgColor('black'), background=cgColor('white')
oplot,alog10(tform_sat(heu_sat_idx_per)), EuH_sat_per, psym = sym(0), symsize= 0.4, color = 250

oplot, logtbins, median_EuH_age_sat, thick = tk+1, color = 0
oplot, logtbins, percent_up_EuH_age_sat, thick = tk+1, linestyle = 2, color = 0
oplot, logtbins, percent_down_EuH_age_sat, thick = tk+1, linestyle = 2, color = 0

legend, label, box =0, charthick = tk, /bottom, /left, color = 0

if (plottops eq 1) then begin
   device, /close_file
   set_plot, entry_device
endif
endif

;^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^;                                        
;                               5PLOT15                              ;                                         
if (1 eq 0) then begin
loadct, 39
if(plottops eq 1) then begin
   entry_device = !d.name
   set_plot,'PS'
   device, bits_per_pixel = 8,/portrait, color = 1,filename = filename72
   tk = 2
endif else begin
   window, 15, retain = 1
   tk = 1
endelse

plot, alog10(tform_sat(heu_sat_idx_per)), EuFe_sat_per, psym = sym(0), xthick=2, ythick = 2, charthick = 2, xtitle = 'time', ytitle='[Eu/Fe] dwarf frac >= 0.95', /nodata, xrange = [8.5, 10], charsize = 1.1, yrange = [-4, 4], Color=cgColor('black'), background=cgColor('white')
oplot,alog10(tform_sat(heu_sat_idx_per)), EuFe_sat_per, psym = sym(0), symsize= 0.4, color = 250

oplot, logtbins, median_EuFe_age_sat, thick = tk+1, color = 0
oplot, logtbins, percent_up_EuFe_age_sat, thick = tk+1, linestyle = 2, color = 0
oplot, logtbins, percent_down_EuFe_age_sat, thick = tk+1, linestyle = 2, color = 0

legend, label, box =0, charthick = tk, /bottom, /left, color = 0

if (plottops eq 1) then begin
   device, /close_file
   set_plot, entry_device
endif
endif

;^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^;                                                      
;                               5PLOT16                              ;                                                       
if (1 eq 0) then begin
loadct, 39
if(plottops eq 1) then begin
   entry_device = !d.name
   set_plot,'PS'
   device, bits_per_pixel = 8,/portrait, color = 1,filename = filename73
   tk = 2
endif else begin
   window, 16, retain = 1
   tk = 1
endelse

plot, FeH_sat(heu_sat_idx_per), EuFe_sat_per, psym = sym(0), xthick=2, ythick = 2, charthick = 2, xtitle = '[Fe/H]', ytitle='[Eu/Fe] dwarf frac >= 0.95', /nodata, xrange = [-4, 1], charsize = 1.1, yrange = [-4, 4], Color=cgColor('black'), background=cgColor('white')
oplot,FeH_sat(heu_sat_idx_per), EuFe_sat_per, psym = sym(0), symsize= 0.4, color = 250

oplot, zbins, median_EuFe_sat, thick = tk+2, color = 0
oplot, zbins, percent_up_EuFe_sat, linestyle = 2, thick = tk+1, color = 0
oplot, zbins, percent_down_EuFe_sat, linestyle = 2, thick = tk+1, color = 0

legend, label, box =0, charthick = tk, /bottom, /left, color = 0

if (plottops eq 1) then begin
   device, /close_file
   set_plot, entry_device
endif
endif

;--------------------------------------------------------------------;                                                       
;                               5PLOT17                              ;
;plot5
if (1 eq 0) then begin
loadct, 39
if(plottops eq 1) then begin
   entry_device = !d.name
   set_plot,'PS'
   device, bits_per_pixel = 8,/portrait, color = 1,filename = filename74
   tk = 2
endif else begin
   window, 17, retain = 1
   tk = 1
endelse

plot, alog10(tform_host(heu_idx_per)), EuH_host_per, psym = sym(0), xthick=2, ythick = 2, charthick = 2, xtitle = 'time', ytitle='[Eu/H] with dwarf >= 0.95', /nodata, xrange = [8.5, 10], charsize = 1.1, yrange = [-4, 4], Color=cgColor('black'), background=cgColor('white')
oplot,alog10(tform_host(heu_idx_per)), EuH_host_per, psym = sym(0), symsize= 0.4, color = 250

oplot, logtbins, median_EuH_age, thick = tk+1, color = 0
oplot, logtbins, percent_up_EuH_age, thick = tk+1, linestyle = 2, color = 0
oplot, logtbins, percent_down_EuH_age, thick = tk+1, linestyle = 2, color = 0

legend, label, box =0, charthick = tk, /bottom, /left, color = 0

if (plottops eq 1) then begin
   device, /close_file
   set_plot, entry_device
endif
endif

;^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^;                                                            
;                               5PLOT18                              ;                                                             
if (1 eq 0) then begin
loadct, 39
if(plottops eq 1) then begin
   entry_device = !d.name
   set_plot,'PS'
   device, bits_per_pixel = 8,/portrait, color = 1,filename = filename75
   tk = 2
endif else begin
   window, 18, retain = 2
   tk = 1
endelse

plot, alog10(tform_host(heu_idx_per)), EuFe_host_per, psym = sym(0), xthick=2, ythick = 2, charthick = 2, xtitle = 'time', ytitle='[Eu/Fe] dwarf frac >= 0.95', /nodata, xrange = [8.5, 10], charsize = 1.1, yrange = [-4, 4], Color=cgColor('black'), background=cgColor('white')
oplot,alog10(tform_host(heu_idx_per)), EuFe_host_per, psym = sym(0), symsize= 0.4, color = 250

oplot, logtbins, median_EuFe_age, thick = tk+1, color = 0
oplot, logtbins, percent_up_EuFe_age, thick = tk+1, linestyle = 2, color = 0
oplot, logtbins, percent_down_EuFe_age, thick = tk+1, linestyle = 2, color = 0

legend, label, box =0, charthick = tk, /bottom, /left, color = 0

if (plottops eq 1) then begin
   device, /close_file
   set_plot, entry_device
endif
endif

;^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^;                                                               
;                               5PLOT19                              ;                                                                
if (1 eq 0) then begin
loadct, 39
if(plottops eq 1) then begin
   entry_device = !d.name
   set_plot,'PS'
   device, bits_per_pixel = 8,/portrait, color = 1,filename = filename76
   tk = 2
endif else begin
   window, 19, retain = 1
   tk = 1
endelse

plot, FeH_host(heu_idx_per), EuFe_host_per, psym = sym(0), xthick=2, ythick = 2, charthick = 2, xtitle = '[Fe/H]', ytitle='[Eu/Fe] dwarf frac >= 0.95', /nodata, xrange = [-4, 1], charsize = 1.1, yrange = [-4, 4], Color=cgColor('black'), background=cgColor('white')
oplot,FeH_host(heu_idx_per), EuFe_host_per, psym = sym(0), symsize= 0.4, color = 250

oplot, zbins, median_EuFe, thick = tk+2, color = 0
oplot, zbins, percent_up_EuFe, linestyle = 2, thick = tk+1, color = 0
oplot, zbins, percent_down_EuFe, linestyle = 2, thick = tk+1, color = 0

legend, label, box =0, charthick = tk, /bottom, /left, color = 0

if (plottops eq 1) then begin
   device, /close_file
   set_plot, entry_device
endif
endif

;^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^;                                                   
;                               5PLOT20                              ;                                                  
if (1 eq 0) then begin
loadct, 39
if(plottops eq 1) then begin
   entry_device = !d.name
   set_plot,'PS'
   device, bits_per_pixel = 8,/portrait, color = 1,filename = filename77
   tk = 2
endif else begin
   window, 20, retain = 1
   tk = 1
endelse

plot, alog10(tform_sat(heu_sat_idx_per)), EuH_sat_per, psym = sym(0), xthick=2, ythick = 2, charthick = 2, xtitle = 'time', ytitle='[Eu/H] with dwarf >= 0.95', /nodata, xrange = [8.5, 10], charsize = 1.1, yrange = [-4, 4], Color=cgColor('black'), background=cgColor('white')
oplot,alog10(tform_sat(heu_sat_idx_per)), EuH_sat_per, psym = sym(0), symsize= 0.4, color = 250

oplot, logtbins, median_EuH_age_sat, thick = tk+1, color = 0
oplot, logtbins, percent_up_EuH_age_sat, thick = tk+1, linestyle = 2, color = 0
oplot, logtbins, percent_down_EuH_age_sat, thick = tk+1, linestyle = 2, color = 0

legend, label, box =0, charthick = tk, /bottom, /left, color = 0

if (plottops eq 1) then begin
   device, /close_file
   set_plot, entry_device
endif
endif

;^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^;              
;                               5PLOT21                              ;               
if (1 eq 0) then begin
loadct, 39
if(plottops eq 1) then begin
   entry_device = !d.name
   set_plot,'PS'
   device, bits_per_pixel = 8,/portrait, color = 1,filename = filename78
   tk = 2
endif else begin
   window, 21, retain = 1
   tk = 1
endelse

plot, alog10(tform_sat(heu_sat_idx_per)), EuFe_sat_per, psym = sym(0), xthick=2, ythick = 2, charthick = 2, xtitle = 'time', ytitle='[Eu/Fe] dwarf frac >= 0.95', /nodata, xrange = [8.5, 10], charsize = 1.1, yrange = [-4, 4], Color=cgColor('black'), background=cgColor('white')
oplot,alog10(tform_sat(heu_sat_idx_per)), EuFe_sat_per, psym = sym(0), symsize= 0.4, color = 250

oplot, logtbins, median_EuFe_age_sat, thick = tk+1, color = 0
oplot, logtbins, percent_up_EuFe_age_sat, thick = tk+1, linestyle = 2, color = 0
oplot, logtbins, percent_down_EuFe_age_sat, thick = tk+1, linestyle = 2, color = 0

legend, label, box =0, charthick = tk, /bottom, /left, color = 0

if (plottops eq 1) then begin
   device, /close_file
   set_plot, entry_device
endif
endif

;^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^;      
;                               5PLOT22                              ;          
if (1 eq 0) then begin
loadct, 39
if(plottops eq 1) then begin
   entry_device = !d.name
   set_plot,'PS'
   device, bits_per_pixel = 8,/portrait, color = 1,filename = filename79
   tk = 2
endif else begin
   window, 22, retain = 1
   tk = 1
endelse

plot, FeH_sat(heu_sat_idx_per), EuFe_sat_per, psym = sym(0), xthick=2, ythick = 2, charthick = 2, xtitle = '[Fe/H]', ytitle='[Eu/Fe] dwarf frac >= 0.95', /nodata, xrange = [-4, 1], charsize = 1.1, yrange = [-4, 4], Color=cgColor('black'), background=cgColor('white')
oplot,FeH_sat(heu_sat_idx_per), EuFe_sat_per, psym = sym(0), symsize= 0.4, color = 250

oplot, zbins, median_EuFe_sat, thick = tk+2, color = 0
oplot, zbins, percent_up_EuFe_sat, linestyle = 2, thick = tk+1, color = 0
oplot, zbins, percent_down_EuFe_sat, linestyle = 2, thick = tk+1, color = 0

legend, label, box =0, charthick = tk, /bottom, /left, color = 0

if (plottops eq 1) then begin
   device, /close_file
   set_plot, entry_device
endif
endif
retall
;^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^;
;^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^;
plot6:
zmin1 = -5.0
zmax1 = 5.0
nzbins = 51

dz = (zmax1 - zmin1)/(nzbins-1)
zbins1 = findgen(nzbins)*dz + zmin1 ;y-axis
index = where(FeH_host lt -2.5 and Ofe_host gt 0)
index1 = where(seu_host gt 0)
index11 = where(dwarfeu gt 0.85)
help, index

;Host
;x-axis
hist = histogram(EuH_host(index), binsize = dz, min=zmin1, max=zmax1)
hist_FeH1 = histogram(FeH_host(index1), binsize = dz, min = zmin1, max=zmax1)
hist_EuH1 = histogram(EuH_host(index1), binsize = dz, min = zmin1, max=zmax1)
hist_EuFe1 = histogram(EuFe_host(index1), binsize = dz, min = zmin1, max=zmax1)

hist_FeH11 = histogram(FeH_host(index11), binsize = dz, min = zmin1, max=zmax1)
hist_EuH11 = histogram(EuH_host(index11), binsize = dz, min = zmin1, max=zmax1)
hist_EuFe11 = histogram(EuFe_host(index11), binsize = dz, min = zmin1, max=zmax1)

;for plotting histograms
histlike, zbins1, hist, zbins1p, histp             ;dist of Eu/H

histlike, zbins1, hist_FeH1, zbins1p, hist_FeH1p   ;dist of Fe/H
histlike, zbins1, hist_EuH1, zbins1p, hist_EuH1p   ;dist of Eu/H
histlike, zbins1, hist_EuFe1, zbins1p, hist_EuFe1p ;dist of Eu/Fe

histlike, zbins1, hist_FeH11, zbins1p, hist_FeH11p   ;dist of Fe/H
histlike, zbins1, hist_EuH11, zbins1p, hist_EuH11p   ;dist of 85% Eu/H
histlike, zbins1, hist_EuFe11, zbins1p, hist_EuFe11p ;dist of 85% Eu/Fe

;Sat
index1s = where(seu_sat gt 0)
index11s = where(dwarfeu gt 0.85)

hist_FeH1s = histogram(FeH_sat(index1s), binsize = dz, min = zmin1, max=zmax1)
hist_EuH1s = histogram(EuH_sat(index1s), binsize = dz, min = zmin1, max=zmax1)
hist_EuFe1s = histogram(EuFe_sat(index1s), binsize = dz, min = zmin1, max=zmax1)

hist_FeH11s = histogram(FeH_sat(index11s), binsize = dz, min = zmin1, max=zmax1)
hist_EuH11s = histogram(EuH_sat(index11s), binsize = dz, min = zmin1, max=zmax1)
hist_EuFe11s = histogram(EuFe_sat(index11s), binsize = dz, min = zmin1, max=zmax1)

;for plotting histograms
histlike, zbins1, hist_FeH1s, zbins1p, hist_FeH1sp   ;dist of Fe/H
histlike, zbins1, hist_EuH1s, zbins1p, hist_EuH1sp   ;dist of Eu/H
histlike, zbins1, hist_EuFe1s, zbins1p, hist_EuFe1sp ;dist of Eu/Fe

histlike, zbins1, hist_FeH11s, zbins1p, hist_FeH11sp   ;dist of Fe/H
histlike, zbins1, hist_EuH11s, zbins1p, hist_EuH11sp   ;dist of 85% Eu/H
histlike, zbins1, hist_EuFe11s, zbins1p, hist_EuFe11sp ;dist of 85% Eu/Fe

;^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^;
;^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^;
;                                6PLOT1                              ;
filename80 = '80EuH__hist_idx1.eps'
if (1 eq 0) then begin
loadct, 39
if(plottops eq 1) then begin
   entry_device = !d.name
   set_plot,'PS'
   device, bits_per_pixel = 8,/portrait, color = 1,filename = filename80
   tk = 2
endif else begin
   window, 24, retain = 2
   tk = 1
endelse
plot, zbins1p, histp, xthick=tk, ythick = tk, charthick = tk, xtitle= '[Eu/H]', ytitle='distribution dN/d[Eu/H]', xrange = [-5,0],charsize = 1.5, thick = tk
label = ['[Fe/H] < -2.5, no smooth']
legend, label, box=0, /top, /left

;plot, zbins1p, hist_EuH1p, xthick=tk, ythick = tk, charthick = tk, xtitle = '[Eu/H]', ytitle='Metallicity distribution', xrange = [-5, 0],charsize = 1.5, thick = tk

if (plottops eq 1) then begin
   device, /close_file
   set_plot, entry_device
endif
endif

;^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^;
;                                6plot2                              ;
filename81 = '81FeH_Allstar_hist_idx1.eps'
if (1 eq 1) then begin
loadct, 39
if(plottops eq 1) then begin
   entry_device = !d.name
   set_plot,'PS'
   device, bits_per_pixel = 8,/portrait, color = 1,filename = filename81
   tk = 2
endif else begin
   window, 25, retain = 2
   tk = 1
endelse
plot, zbins1p, hist_FeH1p/total(hist_FeH1p)/0.1, xthick=tk, ythick = tk, charthick = tk, xtitle = '[Fe/H]', ytitle='Metallicity distribution', xrange = [-3, 2],charsize = 1.5, thick = tk

if (plottops eq 1) then begin
   device, /close_file
   set_plot, entry_device
endif
endif

;^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^;
;                                6plot3                              ;
filename82 = '82EuH_Allstar_hist_idx1.eps'
if (1 eq 1) then begin
loadct, 39
if(plottops eq 1) then begin
   entry_device = !d.name
   set_plot,'PS'
   device, bits_per_pixel = 8,/portrait, color = 1,filename = filename82
   tk = 2
endif else begin
   window, 26, retain = 2
   tk = 1
endelse
plot, zbins1p, hist_EuH1p/total(hist_EuH1p)/0.1, xthick=tk, ythick = tk, charthick = tk, xtitle = '[Eu/H]', ytitle='Metallicity distribution', xrange = [-3, 2.5],charsize = 1.5, thick = tk
;label = ['[Fe/H] < -2.5']
;legend,label, box=0, /top, /left

if (plottops eq 1) then begin
   device, /close_file
   set_plot, entry_device
endif
endif

;^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^;
;                                6plot4                              ;
filename83 = '83EuFe_Allstar_hist_idx1.eps'
if (1 eq 1) then begin
loadct, 39
if(plottops eq 1) then begin
   entry_device = !d.name
   set_plot,'PS'
   device, bits_per_pixel = 8,/portrait, color = 1,filename = filename83
   tk = 2
endif else begin
   window, 27, retain = 2
   tk = 1
endelse
plot, zbins1p, hist_EuFe1p/total(hist_EuFe1p)/0.2, xthick=tk, ythick = tk, charthick = tk, xtitle = '[Eu/Fe]', ytitle='Metallicity distribution', xrange = [-1.5, 1.5],charsize = 1.5, thick = tk
;label = ['[Fe/H] < -2.5']
;legend,label, box=0, /top, /left

if (plottops eq 1) then begin
   device, /close_file
   set_plot, entry_device
endif
endif

;^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^;
;                                6plot5                              ;
filename84 = '84FeH_dwarf85_hist_idx1.eps'
if (1 eq 1) then begin
loadct, 39
if(plottops eq 1) then begin
   entry_device = !d.name
   set_plot,'PS'
   device, bits_per_pixel = 8,/portrait, color = 1,filename = filename84
   tk = 2
endif else begin
   window, 28, retain = 2
   tk = 1
endelse
plot, zbins1p, hist_FeH11p/total(hist_FeH11p)/0.1, xthick=tk, ythick = tk, charthick = tk, xtitle = '[Fe/H]', ytitle='Metallicity distribution', xrange = [-4.5, 0.5], yrange = [0, 1], charsize = 1.5, thick = tk
oplot,zbins1p, hist_FeH1p/total(hist_FeH1p)/0.1, color = 250

if (plottops eq 1) then begin
   device, /close_file
   set_plot, entry_device
endif
endif

;^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^;
;                                6plot6                              ;
filename85 = '85EuH_dwarf85_hist_idx1.eps'
if (1 eq 1) then begin
loadct, 39
if(plottops eq 1) then begin
   entry_device = !d.name
   set_plot,'PS'
   device, bits_per_pixel = 8,/portrait, color = 1,filename = filename85
   tk = 2
endif else begin
   window, 29, retain = 2
   tk = 1
endelse
plot, zbins1p, hist_EuH11p/total(hist_EuH11p)/0.1, xthick=tk, ythick = tk, charthick = tk, xtitle = '[Eu/H]', ytitle='Metallicity distribution', xrange = [-3, 1], yrange = [0, 1], charsize = 1.5, thick = tk
oplot, zbins1p, hist_EuH1p/total(hist_EuH1p)/0.1, color = 250

if (plottops eq 1) then begin
   device, /close_file
   set_plot, entry_device
endif
endif

;^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^;
;                                6plot7                              ;
filename86 = '86EuFe_dwarf85_hist_idx1.eps'
if (1 eq 1) then begin
loadct, 39
if(plottops eq 1) then begin
   entry_device = !d.name
   set_plot,'PS'
   device, bits_per_pixel = 8,/portrait, color = 1,filename = filename86
   tk = 2
endif else begin
   window, 30, retain = 2
   tk = 1
endelse
plot, zbins1p, hist_EuFe11p/total(hist_EuFe11p)/0.21, xthick=tk, ythick = tk, charthick = tk, xtitle = '[Eu/Fe]', ytitle='Metallicity distribution', xrange = [-1, 1], yrange = [0, 1], charsize = 1.5, thick = tk
oplot, zbins1p, hist_EuFe1p/total(hist_EuFe1p)/0.21, color = 250

if (plottops eq 1) then begin
   device, /close_file
   set_plot, entry_device
endif
endif

;^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^;
;                                6plot8                              ;
filename94 = '94FeH_Allstar_sat_hist_idx1.eps'
if (1 eq 1) then begin
loadct, 39
plottops = 1
if(plottops eq 1) then begin
   entry_device = !d.name
   set_plot,'PS'
   device, bits_per_pixel = 8,/portrait, color = 1,filename = filename94
   tk = 2
endif else begin
   window, 25, retain = 2
   tk = 1
endelse
plot, zbins1p, hist_FeH1sp/total(hist_FeH1sp)/0.1, xthick=tk, ythick = tk, charthick = tk, xtitle = '[Fe/H]', ytitle='Metallicity distribution', xrange = [-3, 2],charsize = 1.5, thick = tk

if (plottops eq 1) then begin
   device, /close_file
   set_plot, entry_device
endif
endif

;^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^;                          
;                                6plot9                              ;
filename95 = '95EuH_Allstar_sat_hist_idx1.eps'
if (1 eq 1) then begin
loadct, 39
if(plottops eq 1) then begin
   entry_device = !d.name
   set_plot,'PS'
   device, bits_per_pixel = 8,/portrait, color = 1,filename = filename95
   tk = 2
endif else begin
   window, 26, retain = 2
   tk = 1
endelse
plot, zbins1p, hist_EuH1sp/total(hist_EuH1sp)/0.1, xthick=tk, ythick = tk, charthick = tk, xtitle = '[Eu/H]', ytitle='Metallicity distribution', xrange = [-3, 2.5],charsize = 1.5, thick = tk
;label = ['[Fe/H] < -2.5']                                                                      
;legend,label, box=0, /top, /left                                                               

if (plottops eq 1) then begin
   device, /close_file
   set_plot, entry_device
endif
endif

;^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^;
;                                6plot10                             ;
filename96 = '96EuFe_Allstar_sat_hist_idx1.eps'
if (1 eq 1) then begin
loadct, 39
if(plottops eq 1) then begin
   entry_device = !d.name
   set_plot,'PS'
   device, bits_per_pixel = 8,/portrait, color = 1,filename = filename96
   tk = 2
endif else begin
   window, 27, retain = 2
   tk = 1
endelse
plot, zbins1p, hist_EuFe1sp/total(hist_EuFe1sp)/0.2, xthick=tk, ythick = tk, charthick = tk, xtitle = '[Eu/Fe]', ytitle='Metallicity distribution', xrange = [-1.5, 1.5],charsize = 1.5, thick = tk

if (plottops eq 1) then begin
   device, /close_file
   set_plot, entry_device
endif
endif

;^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^;
;                                6plot11                             ;
filename97 = '97FeH_dwarf85_sat_hist_idx1.eps'
if (1 eq 1) then begin
loadct, 39
if(plottops eq 1) then begin
   entry_device = !d.name
   set_plot,'PS'
   device, bits_per_pixel = 8,/portrait, color = 1,filename = filename97
   tk = 2
endif else begin
   window, 28, retain = 2
   tk = 1
endelse
plot, zbins1p, hist_FeH11sp/total(hist_FeH11sp)/0.1, xthick=tk, ythick = tk, charthick = tk, xtitle = '[Fe/H]', ytitle='Metallicity distribution', xrange = [-4.5, 0.5], yrange = [0, 1], charsize = 1.5, thick = tk
oplot,zbins1p, hist_FeH1sp/total(hist_FeH1sp)/0.1, color = 250

if (plottops eq 1) then begin
   device, /close_file
   set_plot, entry_device
endif
endif

;^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^;
;                                6plot12                             ;
filename98 = '98EuH_dwarf85_sat_hist_idx1.eps'
if (1 eq 1) then begin
loadct, 39
if(plottops eq 1) then begin
   entry_device = !d.name
   set_plot,'PS'
   device, bits_per_pixel = 8,/portrait, color = 1,filename = filename98
   tk = 2
endif else begin
   window, 29, retain = 2
   tk = 1
endelse
plot, zbins1p, hist_EuH11sp/total(hist_EuH11sp)/0.1, xthick=tk, ythick = tk, charthick = tk, xtitle = '[Eu/H]', ytitle='Metallicity distribution', xrange = [-3, 1], yrange = [0, 1], charsize = 1.5, thick = tk
oplot, zbins1p, hist_EuH1sp/total(hist_EuH1sp)/0.1, color = 250

if (plottops eq 1) then begin
   device, /close_file
   set_plot, entry_device
endif
endif

;^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^;
;                                6plot13                             ;
filename99 = '99EuFe_dwarf85_sat_hist_idx1.eps'
if (1 eq 1) then begin
loadct, 39
if(plottops eq 1) then begin
   entry_device = !d.name
   set_plot,'PS'
   device, bits_per_pixel = 8,/portrait, color = 1,filename = filename99
   tk = 2
endif else begin
   window, 30, retain = 2
   tk = 1
endelse
plot, zbins1p, hist_EuFe11sp/total(hist_EuFe11sp)/0.21, xthick=tk, ythick = tk, charthick = tk, xtitle = '[Eu/Fe]', ytitle='Metallicity distribution', xrange = [-1, 1], yrange = [0, 1], charsize = 1.5, thick = tk
oplot, zbins1p, hist_EuFe1sp/total(hist_EuFe1sp)/0.21, color = 250

if (plottops eq 1) then begin
   device, /close_file
   set_plot, entry_device
endif
endif

;^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^;
;^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^;
;Host
index = where(FeH_host lt -2.0 and EuFe_host gt 0.5)
index1 = where(dwarfeu gt 0.0)
index1t = where(alog10(tform_host(index)) gt 9.0)
index2t = where(alog10(tform_host(index)) gt 9.5)
index3t = where(alog10(tform_host(index)) gt 10.0)

logamin = 8.0
logamax = 11
nabins = 23
da = (logamax-logamin)/(nabins-1)
abins = findgen(nabins)*da + logamin
hist_time1 = histogram(alog10(tform_host(index)), binsize = da, min = logamin, max=logamax)
hist_time9 = histogram(alog10(tform_host(index1t)), binsize = da, min = logamin, max=logamax)
hist_time95 = histogram(alog10(tform_host(index2t)), binsize = da, min = logamin, max=logamax)
hist_time10 = histogram(alog10(tform_host(index3t)), binsize = da, min = logamin, max=logamax)

histlike, abins, hist_time9, abinsp, hist_timep9
histlike, abins, hist_time95, abinsp, hist_timep95
histlike, abins, hist_time10, abinsp, hist_timep10

;Sat
index1ts = where(alog10(tform_sat(index)) gt 9.0)
index2ts = where(alog10(tform_sat(index)) gt 9.5)
index3ts = where(alog10(tform_sat(index)) gt 10.0)

hist_time9s = histogram(alog10(tform_sat(index1ts)), binsize = da, min = logamin, max=logamax)
hist_time95s = histogram(alog10(tform_sat(index2ts)), binsize = da, min = logamin, max=logamax)
hist_time10s = histogram(alog10(tform_sat(index3ts)), binsize = da, min = logamin, max=logamax)

histlike, abins, hist_time9s, abinsp, hist_timep9s
histlike, abins, hist_time95s, abinsp, hist_timep95s
histlike, abins, hist_time10s, abinsp, hist_timep10s

;^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^;
;^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^;
;                                7plot1                              ;
filename87 = '87Age_dist_hist_idx1.eps'
plot7:
if (1 eq 1) then begin
loadct, 39
if(plottops eq 1) then begin
   entry_device = !d.name
   set_plot,'PS'
   device, bits_per_pixel = 8,/portrait, color = 1,filename = filename87
   tk = 2
endif else begin
   window, 30, retain = 2
   tk = 1
endelse
plot, abinsp, hist_timep1, xthick=tk, ythick = tk, charthick = tk, xtitle = 'log (Formation time [yr])', ytitle='Age distribution', xrange = [8, 10.2], charsize = 1.1, thick = tk
;label = ['[Fe/H] < -2.0 and [Eu/Fe] > 0.5, no smooth']
;legend, label, box =0, /top, /left

if (plottops eq 1) then begin
   device, /close_file
   set_plot, entry_device
endif
endif

;^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^;
;                                7plot2                              ;
filename88 = '88Age_dist_gt9_hist_idx1.eps'
if (1 eq 1) then begin
loadct, 39
if(plottops eq 1) then begin
   entry_device = !d.name
   set_plot,'PS'
   device, bits_per_pixel = 8,/portrait, color = 1,filename = filename88
   tk = 2
endif else begin
   window, 30, retain = 2
   tk = 1
endelse
plot, abinsp, hist_timep9/total(hist_timep9)/0.1, xthick=tk, ythick = tk, charthick = tk, xtitle = 'log (Formation time [yr])', ytitle='Age distribution', xrange = [8, 10.2], yrange = [0, 6], charsize = 1.1, thick = tk
;label = ['[Fe/H] < -2.0 and Dwarf fraction > 0']
;legend, label, box =0, /top, /left

if (plottops eq 1) then begin
   device, /close_file
   set_plot, entry_device
endif
endif

;^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^;
;                                7plot3                              ;
filename89 = '89Age_dist_gt95_hist_idx1.eps'
if (1 eq 1) then begin
loadct, 39
if(plottops eq 1) then begin
   entry_device = !d.name
   set_plot,'PS'
   device, bits_per_pixel = 8,/portrait, color = 1,filename = filename89
   tk = 2
endif else begin
   window, 30, retain = 2
   tk = 1
endelse
plot, abinsp, hist_timep95/total(hist_timep95)/0.1, xthick=tk, ythick = tk, charthick = tk, xtitle = 'log (Formation time [yr])', ytitle='Age distribution', xrange = [8, 10.2], yrange = [0, 6], charsize = 1.1, thick = tk
;label = ['[Fe/H] < -2.0 and Dwarf fraction > 0']
;legend, label, box =0, /top, /left

if (plottops eq 1) then begin
   device, /close_file
   set_plot, entry_device
endif
endif

;^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^;
;                                7plot4                              ;
filename90 = '90Age_dist_gt10_hist_idx1.eps'
if (1 eq 1) then begin
loadct, 39
if(plottops eq 1) then begin
   entry_device = !d.name
   set_plot,'PS'
   device, bits_per_pixel = 8,/portrait, color = 1,filename = filename90
   tk = 2
endif else begin
   window, 30, retain = 2
   tk = 1
endelse
plot, abinsp, hist_timep10/total(hist_timep10)/0.1, xthick=tk, ythick = tk, charthick = tk, xtitle = 'log (Formation time [yr])', ytitle='Age distribution', xrange = [8, 10.2], yrange = [0, 6], charsize = 1.1, thick = tk
;label = ['[Fe/H] < -2.0 and Dwarf fraction > 0']
;legend, label, box =0, /top, /left

if (plottops eq 1) then begin
   device, /close_file
   set_plot, entry_device
endif
endif

;^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^;
;                                7plot5                              ;
filename91 = '91Age_dist_sat_gt9_hist_idx1.eps'
if (1 eq 1) then begin
loadct, 39
if(plottops eq 1) then begin
   entry_device = !d.name
   set_plot,'PS'
   device, bits_per_pixel = 8,/portrait, color = 1, filename = filename91
   tk = 2
endif else begin
   window, 30, retain = 2
   tk = 1
endelse
plot, abinsp, hist_timep9s, xthick=tk, ythick = tk, charthick = tk, xtitle = 'log (Formation time [yr])', ytitle='Age distribution', xrange = [8, 10.2], charsize = 1.1, thick = tk
label = ['[Fe/H] < -2.0 and [Eu/Fe] > 0.5, no smooth']
legend, label, box =0, /top, /left

if (plottops eq 1) then begin
   device, /close_file
   set_plot, entry_device
endif
endif

;^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^;
;                                7plot5                              ;
filename92 = '92Age_dist_sat_gt95_hist_idx1.eps'
if (1 eq 1) then begin
loadct, 39
if(plottops eq 1) then begin
   entry_device = !d.name
   set_plot,'PS'
   device, bits_per_pixel = 8,/portrait, color = 1, filename = filename92
   tk = 2
endif else begin
   window, 30, retain = 2
   tk = 1
endelse
plot, abinsp, hist_timep95s, xthick=tk, ythick = tk, charthick = tk, xtitle = 'log (Formation time [yr])', ytitle='Age distribution', xrange = [8, 10.2], charsize = 1.1, thick = tk
label = ['[Fe/H] < -2.0 and [Eu/Fe] > 0.5, no smooth']
legend, label, box =0, /top, /left

if (plottops eq 1) then begin
   device, /close_file
   set_plot, entry_device
endif
endif

;^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^;
;                                7plot6                              ;
filename93 = '93Age_dist_sat_gt10_hist_idx1.eps'
if (1 eq 1) then begin
loadct, 39
if(plottops eq 1) then begin
   entry_device = !d.name
   set_plot,'PS'
   device, bits_per_pixel = 8,/portrait, color = 1, filename = filename93
   tk = 2
endif else begin
   window, 30, retain = 2
   tk = 1
endelse
plot, abinsp, hist_timep10s, xthick=tk, ythick = tk, charthick = tk, xtitle = 'log (Formation time [yr])', ytitle='Age distribution', xrange = [8, 10.2], charsize = 1.1, thick = tk
label = ['[Fe/H] < -2.0 and [Eu/Fe] > 0.5, no smooth']
legend, label, box =0, /top, /left

if (plottops eq 1) then begin
   device, /close_file
   set_plot, entry_device
endif
endif


end 
