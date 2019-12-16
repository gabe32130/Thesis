dir0 = "/mn/stornext/d17/extragalactic/personal/gabrierg/"
dir = "/mn/stornext/d17/extragalactic/personal/gabrierg/rprocess/"
dir1 = "/mn/stornext/u3/shens/scratch/Eris_data/"
name = 'L90Mpc8000_hithres.'
outfile1 = dir0 + name + 'seu_nosm_0p01_128'
outfile2 = dir0 + name + 'seu_sm_0p01_128'

;outfile1 = dir + name + 'oxy_smooth'
;outfile2 = dir + name + 'fe_smooth'

tipfile = dir1 +  "L90Mpc8000_hithres.00400"
read_tipsy_header, tipfile, h
ns = h.nstar
nfiles = 400
seu_nosm = fltarr(ns) 
seu_sm = fltarr(ns)
filename = indgen(nfiles)+1

for j=0, nfiles-1 do begin 
   print, "step j=", j+1
   eufile = dir + name + string(filename[j], format = '(i5.5)')+".stellar_rp";_0p01_128"
 ;   eufile = dir + name + string(filename[j], format = '(i5.5)')+".smooth_oxy_fe"
  if (file_test(eufile) eq 0) then begin
     print, "Eu file doesn't exist", filename[j]
     continue
  endif 
  close, 1
  openr, 1, eufile
  istart = 0l
  iend = 0l
  readu,1, istart
  readu, 1, iend 
  ns_step = iend-istart+1 
 ; print, "ns_step", ns_step
 ; print,"istart, iend", istart, iend 

  eufrac_nosm = fltarr(ns_step)
  eufrac_sm = fltarr(ns_step)

  for i=0l, ns_step-1 do begin 
     readu, 1, x
     eufrac_nosm[i] = x
  endfor

  for i=0l, ns_step-1 do begin 
     readu, 1, x
     eufrac_sm[i] = x

     if (eufrac_sm[i] lt 0) then begin
        print, "eufrac_sm", eufrac_sm[i]
        continue
     endif
  endfor
 ; print, "non zero eu",n_elements(where(eufrac_nosm ne 0))

  close, 1 
  seu_nosm[istart:iend] = eufrac_nosm
  seu_sm[istart:iend] = eufrac_sm
  
endfor 

write:
close, 1
openw, 1, outfile1
printf, 1, ns
printf, 1, seu_nosm

close, 1

close, 1
openw, 1, outfile2
printf, 1, ns
printf, 1, seu_sm

close, 1




end 
