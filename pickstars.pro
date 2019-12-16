; ramdonly pick up stellar particles to paint r-process
; elements. Should not choose the same star to paint twice 

msolunit = 1.078e17
rho_mean = 0.042 
lunit = 90.
dsecunit = 1.22341e+18
dyrunit = dsecunit/3600./24./365.


dir = "/mn/stornext/d7/shens/Eris_data/" ;where the snapshots are
;dir1 = "/mn/stornext/d17/extragalactic/personal/gabrierg/rprocess/";where fiducial outfiles go
;dir1 = "/mn/stornext/d17/extragalactic/personal/gabrierg/rprocess_idx2p0/"
dir1 = "/mn/stornext/d17/extragalactic/personal/gabrierg/rprocess_halos/"
name = "L90Mpc8000_hithres."
;file = dir + "NSmergers_per_timestep_delay100Myr_idx1p0_yield0p01.dat"
file = dir + "NSmergers_per_timestep_delay100Myr_idx1p0_yield0p05.dat" ;fiducial model
;file = dir + "NSmergers_per_timestep_delay100Myr_idx2p0_yield0p05.dat"
;file = dir1 + "NSNSmergers_per_timestep.dat"
;goto, after
print, "star list file", file
z0file = dir+name + '00400'
rtipsy, z0file, h, g, d, s
tform = s.tform*dyrunit 

after:
nlines = file_lines(file)
steps = intarr(nlines)
redshift = fltarr(nlines)
time = fltarr(nlines)
nsrate = fltarr(nlines)

close, 1
openr, 1, file
for i=0, nlines-1 do begin 
   readf, 1, x1, x2, x3, x4
   steps[i]=floor(x1)
   redshift[i] = x2
   time[i] = x3
   nsrate[i] = x4
endfor

close, 1

print, "total amount of mergers", total(nsrate)

;retall 

npick = round(nsrate) ; 10 time more stars selected
print, total(npick)
filename = indgen(400)+1
;filename = [63]
nfiles = n_elements(filename)
count = 0
for i=0, nfiles-1 do begin 
   nmerger = npick[filename[i]-1]
   count = count + nmerger 
endfor 
sid_list = lonarr(count)

count = 0
sid_old = [-1, -1, -1]
;seed = 78934
for i=0, nfiles-1 do begin 
   tipfile = dir+name + string(filename[i], format='(i5.5)')
   outfile = dir1+name + string(filename[i], format='(i5.5)')+ '.sid_merger_0p05'
   read_tipsy_header, tipfile, h
   red = (1-h.time)/h.time
   tnow = lookback(100000)-lookback(red)
  ; print, "current time", tnow
   ns = 0L
   ns = h.nstar
   nmerger = npick[filename[i]-1]
   ;print, "merger rates for timestep", filename[i], "     is", nmerger
   if (nmerger ne 0) then begin 
      sid = lonarr(nmerger)
      count = count + nmerger
      random:
      for j=0, nmerger-1 do begin 
        loop: number = floor(randomu(seed)*ns)
        ;print, "stellar age", tnow - tform(number)
         if (tnow - tform(number) lt 1e8) then begin 
            ;print, "stellar age < 100 Myr"
            goto, loop
         endif 
         sid[j]=number
      endfor 
      
      uni_id = sid(uniq(sid, sort(sid)))
      ;print, "unique pick for timestep", filename[i], "     is", n_elements(uni_id) 
  ;    if (n_elements(uni_id) ne n_elements(sid)) then begin 
  ;       print, "elements not unique", i
  ;       goto, random 
   ;   endif 
      
      ;for j=0, nmerger-1 do begin 
       ;  test:
        ; for k=0, n_elements(sid_old)-1 do begin 
         ;   if (sid[j] eq sid_old[k]) then begin 
          ;     print, "element the same as a previous timpstep", i
           ;    sid[j]=floor(randomu(seed)*ns)
            ;   goto, test
           ; endif;
         ;endfor 
      ;endfor 
      close, 1
      openw, 1, outfile 
      printf, 1, nmerger
      for j=0, nmerger-1 do begin 
         printf, 1, sid[j]
      endfor 
      close, 1
      sid_old = [sid_old, sid]
      sid_list[count-nmerger:count-1] = sid 
   endif else begin
      close, 1
      openw, 1, outfile
      printf, 1, nmerger
      close, 1
   endelse 
endfor

x = sid_list(uniq(sid_list, sort(sid_list)))
print, "unique ids", n_elements(x)
print, "total ids", n_elements(sid_list)

end 
