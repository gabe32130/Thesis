;pro eu_fe, tipfile, grpfile, eufile, oxyfile, fefile, dyrunit, ismooth, outfile, output, plottops,filename1, filename2

;dir = '/Volumes/2/simulations_disk2/mwrun_data/'
dir = '/mn/stornext/u3/shens/scratch/Eris_data/'
dir1 = '/mn/stornext/d17/extragalactic/personal/gabrierg/rprocess/'
dir2 = '/mn/stornext/u3/shens/scratch/Eris_amiga/'
dir3 = '/mn/stornext/u3/shens/scratch/rprocess/'
name = 'L90Mpc8000_hithres.'
step = [400]
tipfile = dir + name + string(step[0], format='(i5.5)')
grpfile = dir2 + name + string(step[0], format='(i5.5)')  + '.amiga.grp';Eris_amiga
;oxyfile = tipfile + '.OxMassFrac'
;fefile = tipfile + '.FeMassFrac'
oxyfile = dir3 + name + 'oxy_smooth_128'; in rprocess dir
fefile = dir3 + name + 'fe_smooth_128' ; in rprocess dir
eufile = dir1 + name + 'seu_sm_0p01_128'; step5 output

;dsecunit = 1.2236e+18
;dyrunit = dsecunit/3600./24./365.
dyrunit = 38767.4e6

outfile = dir1 + 'Eris_host_abundance_smooth_idx1_350pc.dat'
filename1 = "Age_EuFe_mergeronly_sm_idx1.eps"
filename2 = "EuFe_FeH_mergeronly_sm_idx1.eps"

ismooth = 1
output = 1
plottops = 0
ifloor = 1
icheckh = 0
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


FeH_solar = 7.45-12
OH_solar = 8.66-12
EuH_solar = 0.52 -12 

Ofe_solar = OH_solar - FeH_solar
Eufe_solar = EuH_solar - FeH_solar 
EuO_solar = EuH_solar - OH_solar 

ageofuniverse =  lookback(1000000) ; in yrs 
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
close, 1
openr, 1, eufile
ns_eu = 0l
readf, 1, ns_eu
print, "ns_eu, ns", ns_eu, ns
seu = fltarr(ns_eu)
readf, 1, seu
close, 1

index = where(seu ne 0) 
print, "fraction of eu enriched", n_elements(index)*1.0/ns

seu_SN = dblarr(ns)
seu_SN = soxy*meu_mo_solar
seu_tot = seu + seu_SN


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

if (output eq 1) then begin 
   close, 1
   openw, 1, outfile
;printf, 1, "Assuming solar abunfance [Fe/H] = 7.45, [O/H]=8.66, [Eu/H] = 0.52" 
   printf, 1, "Age [yr], n_Fe/n_H,  n_O/n_H, n_Eu/n_H (NSNS merger only), n_Eu/n_H (merger+Type II), distance [kpc]"
   printf, 1, nhost 
   for i=0l,nhost-1 do begin 
      printf, 1, age_host[i], nFenH[i], nOxynH[i], nEunH[i], nEunH_tot[i]
   endfor   
   close, 1 
endif 

; binning in age 

;plot:
logtmin = 8.5
logtmax = 10.5 
nbins = 21

binning, tform_host, FeH_host, FeH_host, logtmin, logtmax, nbins, logtbins, median_FeH, percent_up_FeH, percent_down_FeH, number, /log
binning_ave, tform_host, sfe_host, logtmin, logtmax, nbins, logtbins, ave_sfe, rms_sfe, /log
binning_ave, tform_host, smass_host, logtmin, logtmax, nbins, logtbins, ave_smass, rms_smass, /log
ave_FeH = alog10(ave_sfe/(ave_smass*fh) *(Ah/Afe))-FeH_solar

binning, tform_host, OFe_host, OFe_host, logtmin, logtmax, nbins, logtbins, median_OFe_time, percent_up_OFe_time, percent_down_OFe_time, number, /log

binning, tform_host, EuFe_host, seu_host, logtmin, logtmax, nbins, logtbins, median_EuFe_age, percent_up_EuFe_age, percent_down_EuFe_age, number1, /log, /nozero

binning, tform_host, EuFe_host_tot, seu_host_tot, logtmin, logtmax, nbins, logtbins, median_EuFe_age_tot, percent_up_EuFe_age_tot, percent_down_EuFe_age_tot, number1_tot, /log, /nozero

;binning, tform_host(index), FeH_host(index), FeH_host(index), logtmin, logtmax, nbins, logtbins, median_FeH, percent_up_FeH, percent_down_FeH, number, /log
;binning, tform_host(index), OFe_host(index), OFe_host(index), logtmin, logtmax, nbins, logtbins, median_OFe_time, percent_up_OFe_time, percent_down_OFe_time, number, /log

;binning, tform_host(index), EuFe_host(index), seu_host(index), logtmin, logtmax, nbins, logtbins, median_EuFe_age, percent_up_EuFe_age, percent_down_EuFe_age, number1, /log, /nozero

;binning, tform_host(index), EuFe_host_tot(index), seu_host_tot(index), logtmin, logtmax, nbins, logtbins, median_EuFe_age_tot, percent_up_EuFe_age_tot, percent_down_EuFe_age_tot, number1_tot, /log, /nozero

zmin = -4.0
zmax = 1.0 
nzbins = 51

binning, FeH_host, OFe_host, sOxy_host, zmin, zmax, nzbins, zbins, median_OFe, percent_up_OFe, percent_down_OFe,number_oxy,/nozero 

binning, FeH_host, EuFe_host, seu_host, zmin, zmax, nzbins, zbins, median_EuFe, percent_up_EuFe, percent_down_EuFe,number,/nozero 

binning, FeH_host, EuFe_host_tot, seu_host_tot, zmin, zmax, nzbins, zbins, median_EuFe1, percent_up_EuFe1, percent_down_EuFe1, number_tot,/nozero 


;binning, FeH_host(index), OFe_host(index), sOxy_host(index), zmin, zmax, nzbins, zbins, median_OFe, percent_up_OFe, percent_down_OFe,number_oxy, /nozero 

;binning, FeH_host(index), EuFe_host(index), seu_host(index), zmin, zmax, nzbins, zbins, median_EuFe, percent_up_EuFe, percent_down_EuFe,number, /nozero 

;binning, FeH_host(index), EuFe_host_tot(index), seu_host_tot(index),
;zmin, zmax, nzbins, zbins, median_EuFe1, percent_up_EuFe1,
;percent_down_EuFe1, number_tot, /nozero
;plot: 
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
;                                PLOT1                               ;
;retall 
plot: 

if (1 eq 0) then begin 
plottops = 0
loadct, 39 
if(plottops eq  1) then begin
   entry_device = !d.name
   set_plot,'PS'
   device, bits_per_pixel = 8,/portrait, color = 1,filename ='Age_dist_lowZ_smooth.eps'
   tk = 3
endif else begin 
   window, 3, retain = 2
   tk = 1
endelse 
;plot, zbins1p, histp, xthick=tk, ythick = tk, charthick = tk, xtitle = '[Eu/H]', ytitle='distribution dN/d[Eu/H]', xrange = [-5, 0], charsize = 1.5, thick = tk
;label = ['[Fe/H] < -2.5, no smooth']
plot, abinsp, hist_timep1, xthick=tk, ythick = tk, charthick = tk, xtitle = 'log (Formation time [yr])', ytitle='Age distribution', xrange = [8, 10.2], charsize = 1.5, thick = tk
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
;retall
;plottops = 0
loadct, 39 
if(plottops eq  1) then begin
   entry_device = !d.name
   set_plot,'PS'
   device, bits_per_pixel = 8,/portrait, color = 1,filename ='Age_FeH_smooth_128.eps'
   tk = 3
endif else begin 
   window, 0, retain = 2
   tk = 1
endelse 
plot, alog10(tform_host(0:nhost-1:100)), FeH_host(0:nhost-1:100), psym = 3, xthick=tk, ythick = tk, charthick = tk, xtitle = 'log (Formation time [yr])', ytitle='[Fe/H]', /nodata, xrange = [8, 10.5], xstyle = 1, yrange = [-4, 0.5], charsize = 1.5 
oplot,alog10(tform_host(0:nhost-1:100)), FeH_host(0:nhost-1:100),psym= 3, color = 250

;plot, alog10(tform_host(index(0:nindex-1:100))), FeH_host(index(0:nindex-1:100)), psym = 3, xthick=tk, ythick = tk, charthick = tk, xtitle = 'log (Formation time [yr])', ytitle='[Fe/H]', /nodata, xrange = [8, 10.5], xstyle = 1, yrange = [-4, 0.5], charsize = 1.5 
;oplot,alog10(tform_host(index(0:nindex-1:100))), FeH_host(index(0:nindex-1:100)),psym = 3, color = 250

oplot, logtbins, median_FeH, thick = tk
oplot, logtbins, ave_FeH, thick = tk+3, linestyle  =1  
oplot, logtbins, percent_down_FeH, thick = tk, linestyle = 2 
oplot, logtbins, percent_up_FeH, thick = tk, linestyle = 2 
if (plottops eq 1) then begin 
   device, /close_file
   set_plot, entry_device 
endif 

;--------------------------------------------------------------------;
;                                PLOT3                               ;
;retall 
loadct, 39 
if(plottops eq  1) then begin
   entry_device = !d.name
   set_plot,'PS'
   device, bits_per_pixel = 8,/portrait, color = 1,filename = filename1
   tk = 3
endif else begin 
   window, 1, retain = 2
   tk = 1
endelse 
plot, alog10(tform_host(0:nhost-1:50)), EuFe_host(0:nhost-1:50), psym = sym(1), xthick=tk, ythick = tk, charthick = tk, xtitle = 'log (Formation time [yr])', ytitle='[Eu/Fe]', /nodata, xrange = [8.6, 10.5], xstyle = 1, yrange = [-3, 2], charsize = 1.5 
oplot,alog10(tform_host(0:nhost-1:50)), EuFe_host(0:nhost-1:50),psym = sym(1), color = 70, symsize = 0.3

;plot,alog10(tform_host(index(0:nindex-1:100))), EuFe_host(index(0:nindex-1:100)), psym = sym(1), xthick=tk, ythick = tk, charthick = tk, xtitle = 'log (Formation time [yr])', ytitle='[Eu/Fe]', /nodata, xrange = [8.6, 10.5], xstyle = 1, yrange = [-3, 2], charsize = 1.5 
;oplot,alog10(tform_host(index(0:nindex-1:100))),EuFe_host(index(0:nindex-1:100)),psym = sym(1), color = 70, symsize = 0.3

oplot, logtbins, median_EuFe_age, thick = tk+2 
oplot, logtbins, percent_down_EuFe_age, thick = tk+2, linestyle = 2 
oplot, logtbins, percent_up_EuFe_age, thick = tk+2, linestyle = 2 

;label=['merger + SN II, 128 Nbr']
legend, label, box =0, charthick = tk, /bottom, /left

if (plottops eq 1) then begin 
   device, /close_file
   set_plot, entry_device 
endif 

;--------------------------------------------------------------------;
loadct, 39 
if(plottops eq  1) then begin
   entry_device = !d.name
   set_plot,'PS'
   device, bits_per_pixel = 8,/portrait, color = 1,filename =filename2
   tk = 3
endif else begin 
   window, 2, retain = 2
   tk = 1
endelse 
plot, FeH_host(0:nhost-1:50),EuFe_host(0:nhost-1:50), xthick=tk, ythick = tk, charthick = tk, xtitle = '[Fe/H]', ytitle='[Eu/Fe]', /nodata, xrange = [-4, 0.5], charsize = 1.5, yrange = [-4, 4]
oplot, FeH_host(0:nhost-1:50),EuFe_host(0:nhost-1:50), psym = sym(1), symsize= 0.3, color = 70

;plot, FeH_host(index(0:nindex-1:50)),EuFe_host(index(0:nindex-1:50)), xthick=tk, ythick = tk, charthick = tk, xtitle = '[Fe/H]', ytitle='[Eu/Fe]', /nodata, xrange = [-4, 0.5], charsize = 1.5, yrange = [-4, 4]
;oplot, FeH_host(index(0:nindex-1:50)),EuFe_host(index(0:nindex-1:50)), psym = sym(1),symsize= 0.3, color = 70

oplot, zbins, median_EuFe, thick = tk+2 
oplot, zbins, percent_down_EuFe, linestyle = 2, thick = tk
oplot, zbins, percent_up_EuFe, linestyle = 2, thick = tk

;label=['NSNS merger only, 16 Neighbours']
legend, label, box =0, charthick = tk, /bottom, /left

if (plottops eq 1) then begin 
   device, /close_file
   set_plot, entry_device 
endif 

end 
