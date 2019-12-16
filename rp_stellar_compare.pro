;dir = '/Users/shens/mwrun/'
;dir1 = '/Volumes/2/mwrun_data/'
dir1 = '/mn/stornext/u3/shens/scratch/Eris_data/'
dir2 = '/mn/stornext/d17/extragalactic/personal/gabrierg/rprocess/'

tipfile = dir1 + "L90Mpc8000_hithres.00400"
eufile = dir1 + "L90Mpc8000_hithres.sid_paint2";what is this file?
eufile_step = dir2 + 'L90Mpc8000_hithres.00100.stellar_rp';why run 00100?

lunit = 90. 
msolunit = 1.078e17
fH = 0.764 
vunit = 2269.98 
kpcunit = 90000.

FeH_solar = 7.45-12
OH_solar = 8.66-12
EuH_solar = 0.52 -12 

Ofe_solar = OH_solar - FeH_solar
Eufe_solar = EuH_solar - FeH_solar 
EuO_solar = EuH_solar - OH_solar 

dsecunit = 1.2236e+18
dyrunit = dsecunit/3600./24./365.
;ageofuniverse =  1.3739577e+10 ; in yrs 
ageofuniverse =  lookback(1000000) ; in yrs 

Ah  = 1.00794
Ao  = 15.9994
Afe = 55.845
Aeu = 151.964 
meu_mo_solar = 10^(EuO_solar)*Aeu/Ao 

read_tipsy_header, tipfile, h

ns = h.nstar
seu = dblarr(ns)
nseu = 0l
record= {sid:0L,gid:0L,stepform:0,timeform:0.0d, massform:0.0d, eufrac:0.0d}
close, 1
openr, 1, eufile
readu, 1, nseu  
print, "nseu =", nseu
sinfo = replicate(record, nseu)
readu, 1, sinfo 
close, 1 

for i=0l,nseu-1 do begin 
   seu[sinfo[i].sid] = sinfo[i].eufrac
endfor 

;index = where(sinfo.stepform eq 99)
;help, index
;print,"min max IDs", min(sinfo(index).sid),  min(sinfo(index).sid)

openr, 1, eufile_step
istart = 0l
iend = 0l
readu,1, istart
readu, 1, iend 
ns_step = iend-istart+1 
eufrac_nosm = fltarr(ns_step)
eufrac_sm = fltarr(ns_step)

for i=0l, ns_step-1 do begin 
   readu, 1, x
   eufrac_nosm[i] = x
endfor 

for i=0l, ns_step-1 do begin 
   readu, 1, x
   eufrac_sm[i] = x
endfor 

close, 1 

print, "ns_step", ns_step
print,"istart, iend", istart, iend 

print, eufrac_nosm[0:10]
print, eufrac_sm[0:10]

print, seu[istart:istart+10]
seu_sub = seu[istart:iend]
index1 = where(seu_sub ne 0)
index2 = where(eufrac_nosm ne 0)
index3 = where(eufrac_sm ne 0)

print, min(seu_sub(index1))
print, min(eufrac_nosm(index2))
print, min(eufrac_sm(index3))

;print, n_elements(where(seu_sub ne 0))
;print, n_elements(where(eufrac_nosm ne 0))
;print, n_elements(where(eufrac_sm ne 0))

loadct, 39
window, 0, retain = 2 

plot, eufrac_nosm(index2), psym = sym(1), symsize = 0.5, /ylog, yrange = [1e-15, 1e-8]
oplot, eufrac_sm(index3), psym = sym(1), symsize = 0.5, color = 250


end 
