dir0 = "/mn/stornext/d17/extragalactic/personal/gabrierg/" ;test directory dump
dir = "/mn/stornext/d17/extragalactic/personal/gabrierg/rprocess_halos/"
dir1 = "/mn/stornext/u3/shens/scratch/Eris_data/"
name = 'L90Mpc8000_hithres.'
outfile1 = dir + name + 'seu_nosm_0p01_128'
outfile2 = dir + name + 'seu_sm_0p01_128'
outfile3 = dir + name + 'seu_dwarf_0p01_128'     ;ADDED TO TRACK THE DWARF FRACTION
outfile4 = dir + name + 'seu_dir_dwarf_0p01_128' ;ADDED TO TRACK DIRECT DWARF FRACTION

;outfile1 = dir + name + 'oxy_smooth'
;outfile2 = dir + name + 'fe_smooth'

tipfile = dir1 +  "L90Mpc8000_hithres.00400"
read_tipsy_header, tipfile, h
ns = h.nstar
nfiles = 400
seu_nosm = fltarr(ns)
seu_sm = fltarr(ns)
seu_dwarf = fltarr(ns)          ;ADDED TO TRACK THE DWARF FRACTION
dir_eudwarf = fltarr(ns)        ;ADDED TO TRACK DIRECT DWARF FRACTION
;for i = 14, nfiles-1 do begin
filename = indgen(nfiles)+1
;--------------------------------------------------------------------;
;                                BLOCK1                              ;
;creat outfiles
for j=0, nfiles-1 do begin 
   print, "step j=", j+1
   eufile = dir + name + string(filename[j], format = '(i5.5)')+".stellar_rp" ;_0p01_128"
;   print, eufile
   ;eufile = dir + name + string(filename[j], format = '(i5.5)')+".smooth_oxy_fe"
   if (file_test(eufile) eq 0) then begin
      print, "Eu file doesn't exist", filename[j]
      continue
   endif

;write info in outfiles
   close, 1
   openr, 1, eufile
   istart = 0l
   iend = 0l
   readu,1, istart
   readu, 1, iend 
   ns_step = iend-istart+1 
;  print, "ns_step", ns_step
;  print,"istart, iend", istart, iend 
  
   eufrac_nosm = fltarr(ns_step)
   eufrac_sm = fltarr(ns_step)
   eufrac_dwarf = fltarr(ns_step)       ;ADDED TO TRACK THE DWARF FRACTION
   eufrac_dwarfdir = fltarr(ns_step)    ;ADDED TO TRACK DIRECT DWARF FRACTION 


;--------------------------------------------------------------------;
;                               BLOCK1.1                             ;
   for i=0l, ns_step-1 do begin 
      readu, 1, x
      eufrac_nosm[i] = x
;     print, "x =", x
      if (eufrac_nosm[i] lt 0) then begin
         print, "eufrac_nosm", eufrac_nosm[i]
         continue
      endif
   endfor
   
   for i=0l, ns_step-1 do begin 
      readu, 1, x
      eufrac_sm[i] = x
      if (eufrac_sm[i] lt 0) then begin
         print, "eufrac_sm", eufrac_sm[i]
         continue
      endif
   endfor
   ;print, "non zero eu",n_elements(where(eufrac_nosm ne 0))
   
   for i=0l, ns_step-1 do begin ;ADDED TO TRACK THE DWARF FRACTION
      readu, 1, x
      eufrac_dwarf[i] = x
      if (eufrac_dwarf[i] lt 0) then begin
         print, "eufrac_dwarf", eufrac_dwarf[i]
         continue
      endif
   endfor
   
   for i=0l, ns_step-1 do begin ;ADDED TO TRACK DIRECT DWARF FRACTION
      readu, 1, x
      eufrac_dwarfdir[i] = x
     ;if (eufrac_dwarfdir[i] gt 1e-37) then begin
      if (eufrac_dwarfdir[i] lt 0) then begin
         print, "eufrac_dwarfdir", eufrac_dwarfdir[i]
         continue
      endif
   endfor
   
   close, 1
   seu_nosm[istart:iend] = eufrac_nosm
   seu_sm[istart:iend] = eufrac_sm
   seu_dwarf[istart:iend] = eufrac_dwarf         ;ADDED TO TRACK THE DWARF FRACTION
   dir_eudwarf[istart:iend] = eufrac_dwarfdir    ;ADDED TO TRACK DIRECT DWARF FRACTION
   
;endfor

;retall
;--------------------------------------------------------------------;
;                                BLOCK2                              ;
outfile5 = dir + name + string(filename[j], format = '(i5.5)')+"seu_dwarf_0p01_128"
write:
close, 1                        ;ADDED TO TRACK THE DWARF FRACTION
openw, 1, outfile5
printf, 1, ns
printf, 1, seu_dwarf

close, 1
endfor


;write:
;close, 1
;openw, 1, outfile1
;printf, 1, ns
;printf, 1, seu_nosm

;close, 1

;close, 1
;openw, 1, outfile2
;printf, 1, ns
;printf, 1, seu_sm

;close, 1

;close, 1                        ;ADDED TO TRACK THE DWARF FRACTION
;openw, 1, outfile3
;printf, 1, ns
;printf, 1, seu_dwarf

;close, 1

;close, 1                        ;ADDED TO TRACK DIRECT DWARF FRACTION
;openw, 1, outfile4
;printf, 1, ns
;printf, 1, dir_eudwarf

;close, 1
   
print, "Done"
end
