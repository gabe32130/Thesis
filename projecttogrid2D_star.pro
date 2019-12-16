dir = '/mn/stornext/u3/shens/scratch/Eris_data/'
dir1 = '/mn/stornext/d17/extragalactic/personal/gabrierg/rprocess_halos/'             
dir2 = '/mn/stornext/d17/extragalactic/personal/gabrierg/Eris_AHF/'
name = 'L90Mpc8000_hithres.'
step = [400]
tipfile = dir + name + string(step[0], format='(i5.5)')
grpfile = dir2 + name + string(step[0], format='(i5.5)')  + '.amiga.grp' ;Eris_amiga
posCenter = dir2 + name + string(step[0], format='(i5.5)')  + '.amiga.stat'
oxyfile = dir1 + name + 'oxy_smooth_128'      ; in rprocess_halos dir
fefile = dir1 + name + 'fe_smooth_128'        ; in rprocess_halos dir
eufile = dir1 + name + 'seu_sm_0p01_128'      ; step5 output smoothing
dwarfile = dir1 + name + 'seu_dwarf_0p01_128' ; step5 output with DWARF FRACTION  

kpcunit = 90000.
msolunit=1.078e17
;kpcunit = 25000. 
;msolunit = 2.31e15
rhounit = 1.00169e-29
ageofuniverse = 1.37412e+10
tu = 38.80e9
t2 = ageofuniverse 
t1 = ageofuniverse - 1.0e9
;goto, plot

;--------------------------------------------------------------------; 
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

ns = h.nstar

print, "Read in oxygen file"
close, 1
openr, 1, oxyfile
ns1 = 0l
readf, 1, ns1
print, "ns1, ns", ns1, ns
soxy1 = fltarr(ns1)
readf, 1, soxy1
close, 1

print, "Read in iron file"
openr, 1, fefile
ns1 = 0l
readf, 1, ns1
print, "ns1, ns", ns1, ns
sfe1 = fltarr(ns1)
readf, 1, sfe1
close, 1

sfe = sfe1 + fe_floor
soxy = soxy1 + o_floor

ns = h.nstar

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

;we look at host particles
ihost = where(sgrp eq 1)

print, "number of stars in host", n_elements(ihost)
nhost =  n_elements(ihost)

;we look at non-host particles
isat = where(sgrp gt 1)
print, "number of stars in satellites", n_elements(isat)
nsat =  n_elements(isat)

print, "Read in dwarf file"
close, 1
openr, 1, dwarfile
dwarf = 0l
readf, 1, dwarf
print, "dwarf, ns", dwarf, ns
dwarfeu = fltarr(dwarf)
readf, 1, dwarfeu
close, 1
dwarfindex = where(dwarfeu ne 0)
print, "fraction of eu enrichment from dwarfs",n_elements(dwarfindex)*1.0/ns

dwarfeu_SN = dblarr(ns)
dwarfeu_SN = soxy*meu_mo_solar
dwarfeu_tot = dwarfeu + dwarfeu_SN

;--------------------------------------------------------------------;
n1=401
n2=401
d1kpc = 50.0;can change b/w 20-100
d2kpc = 50.0; 50 or 100 & 150
dd1_p = 2.0*d1kpc/(n1-1) 
dd2_p = 2.0*d2kpc/(n2-1)

d1binsp = findgen(n1)*dd1_p -d1kpc ;+ 0.5*dd1_p 
d2binsp = findgen(n2)*dd2_p -d2kpc ;+ 0.5*dd2_p 

statfile = dir2 + name +  string(step[0], format = '(i5.5)') +'.grpcenter'
velfile = dir2 + name +  string(step[0], format = '(i5.5)') +'.grpvc'
outfile1 = dir1 + name +  string(step[0], format = '(i5.5)') +'.sigma'
;outfile2 =  dir1 + name +  string(step[0], format = '(i5.5)') +'.eusigma'
;outfile3 =  dir1 + name +  string(step[0], format = '(i5.5)') +'.deusigma'

rtipsy, tipfile, h, g, d, s 
redshift = (1.0-h.time)/h.time
lunit = kpcunit/(1+redshift)
;goto, plot
 
if (1 eq 0) then begin 
   close, 1
   openr, 1, HIfile 
   ntot=0L
   readf, 1, ntot 
   fHI_tot = fltarr(ntot)
   readf, 1, fHI_tot
   close, 1 
   fHI = fHI_tot(0:h.nsph-1)
endif  

print, "Read in stat file"
grpid = 1
iproj = 1
close, 1
openr, 1, statfile 
nn = 0
readf, 1, nn
xc_all = fltarr(nn)
yc_all = fltarr(nn)
zc_all = fltarr(nn)
readf, 1, xc_all
readf, 1, yc_all
readf, 1, zc_all 
close, 1

xc = xc_all[grpid]
yc = yc_all[grpid]
zc = zc_all[grpid]

print, "Read in velocity file"
close, 1
openr, 1, velfile 
nn = 0
readf, 1, nn
vxc_all = fltarr(nn)
vyc_all = fltarr(nn)
vzc_all = fltarr(nn)
readf, 1, vxc_all
readf, 1, vyc_all
readf, 1, vzc_all 
close, 1

vxc = vxc_all[grpid]
vyc = vyc_all[grpid]
vzc = vzc_all[grpid]

d1 = d1kpc/lunit 
d2 = d2kpc/lunit 
balign = 1  

;--------------------------------------------------------------------;
for k=13,400 do begin
   dwarfiles = dir1 + name + string(k, format='(i5.5)') + 'seu_dwarf_0p01_128'
   outfile5 = dir1 + name +  string(k, format = '(i5.5)') +'.sigma5'

   print, "Read in dwarf file"
   print, dwarfiles
;endfor
;retall
   close, 1
   openr, 1, dwarfiles
   dwarfs = 0l
   readf, 1, dwarfs
;   print, "dwarf, ns", dwarfs, ns
   dwarfeus = fltarr(dwarfs)
   readf, 1, dwarfeus
   close, 1
   dwarfindexs = where(dwarfeus ne 0)
   print, "fraction of eu enrichment from dwarfs",n_elements(dwarfindexs)*1.0/ns
   
   dwarfeu_SNs = dblarr(ns)
   dwarfeu_SNs = soxy*meu_mo_solar
   dwarfeu_tots = dwarfeus + dwarfeu_SNs
   
   
   if (balign eq 1) then begin 
      g.pos[0] = g.pos[0] - xc 
      g.pos[1] = g.pos[1] - yc 
      g.pos[2] = g.pos[2] - zc 
      d.pos[0] = d.pos[0] - xc 
      d.pos[1] = d.pos[1] - yc 
      d.pos[2] = d.pos[2] - zc 
      s.pos[0] = s.pos[0] - xc 
      s.pos[1] = s.pos[1] - yc 
      s.pos[2] = s.pos[2] - zc 
      
      g.vel[0] = g.vel[0] - vxc 
      g.vel[1] = g.vel[1] - vyc 
      g.vel[2] = g.vel[2] - vzc 
      d.vel[0] = d.vel[0] - vxc 
      d.vel[1] = d.vel[1] - vyc 
      d.vel[2] = d.vel[2] - vzc 
      s.vel[0] = s.vel[0] - vxc 
      s.vel[1] = s.vel[1] - vyc 
      s.vel[2] = s.vel[2] - vzc 
      
      limit = 5.0/lunit
      print, "align to the xy plane"
   ;align_star, g, d, s, limit 
      align, g, d, s, limit, jjx, jjy, jjz
      
      xc = 0.
      yc = 0. 
      zc = 0. 
      vxc = 0. 
      vyc = 0. 
      vzc = 0.
      
   endif 
plot:

;--------------------------------------------------------------------;
   ido2d = 1
   if (ido2d eq 1) then begin 
;   d1min = xc-d1 
;   d1max = xc+d1
;   d2min = yc-d2
;   d2max = yc+d2

;    d1min = yc-d1 
;    d1max = yc+d1
;    d2min = zc-d2
;    d2max = zc+d2

      d1min = zc-d1 
      d1max = zc+d1
      d2min = xc-d2
      d2max = xc+d2

      dd1 = (d1max - d1min)/n1  ; in code unit, comoving
      dd2 = (d2max - d2min)/n2
   
      d1bins = findgen(n1)*dd1 + d1min + 0.5*dd1 
      d2bins = findgen(n2)*dd2 + d2min + 0.5*dd2 

      mass2D = fltarr(n1, n2)
      eumass2D = fltarr(n1,n2) 
      deumass2D = fltarr(n1,n2)
      deumass2Ds = fltarr(n1,n2)
      frac2D = fltarr(n1, n2)

;   index= where(s.pos[0] ge d1min and s.pos[0] lt d1max and s.pos[1] ge d2min and s.pos[1] lt d2max);xy
;   index= where(s.pos[1] ge d1min and s.pos[1] lt d1max and s.pos[2] ge d2min and s.pos[2] lt d2max);zx
      index= where(s.pos[2] ge d1min and s.pos[2] lt d1max and s.pos[0] ge d2min and s.pos[0] lt d2max) ;yz
  
      for ii=0L, n_elements(index)-1 do begin 
         if (ii/10000*10000 eq ii) then print, "ii", ii
         j = floor((s(index[ii]).pos[0] -d1min)/dd1)
         k = floor((s(index[ii]).pos[2] -d2min)/dd2)
         if (j lt n1 and k lt n2) then begin 
            mass2D(j, k) = mass2D(j, k) + s(index[ii]).mass
            eumass2D(j, k) = eumass2D(j, k) + s(index[ii]).mass*seu(index[ii])
            deumass2D(j, k) = deumass2D(j, k) + s(index[ii]).mass*seu(index[ii])*dwarfeu(index[ii])
            deumass2Ds(j, k) = deumass2Ds(j, k) + s(index[ii]).mass*seu(index[ii])*dwarfeus(index[ii])
            frac2D(j, k) = deumass2D(j, k) +  s(index[ii]).mass*dwarfeu(index[ii])
         endif 
      endfor 
      
      sigma = mass2D*msolunit/(dd1_p*dd2_p)    ; solar mass per pc^2
      eusigma = eumass2D*msolunit/(dd1_p*dd2_p) ; in solar mass per kpc^2  
      deusigma = deumass2D*msolunit/(dd1_p*dd2_p) ; in solar mass per yr per kpc^2 
      deusigmas = deumass2Ds*msolunit/(dd1_p*dd2_p)
      frac2D = frac2D/mass2D    ; mass weighted dwarf fraction  

   endif 

;--------------------------------------------------------------------;  
   if(ido2d eq 1) then begin 
      ;close, 1
      ;openw, 1, outfile1 
      ;printf, 1, n1, n2, grpid, iproj, d1kpc, d2kpc, xc, yc, zc 
      ;printf, 1, d1binsp
      ;printf, 1, d2binsp
      ;printf, 1, sigma
      ;printf, 1, eusigma
      ;printf, 1, deusigma
      ;printf, 1, frac2D
      ;close, 1
      
      close, 1
      openw, 1, outfile5
      printf, 1, n1, n2, grpid, iproj, d1kpc, d2kpc, xc, yc, zc
      printf, 1, d1binsp
      printf, 1, d2binsp
      printf, 1, deusigmas
      close, 1
 
   endif 
endfor
end 
