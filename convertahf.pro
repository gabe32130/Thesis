;dir = "/Users/shens/mwrun/"
;dir = "/Volumes/3/dg1_new2_hi/"
;dir = "/Volumes/2/agora_1e11_large/"
;dir = '/data/shens/Eris2_RTforce/'
;root = 'L90Mpc8000_hithres_rtforce.'
;dir =  '/mn/stornext/u3/shens/extragal/Eris_AHF/'
dir = '/mn/stornext/d17/extragalactic/personal/gabrierg/Eris_AHF/'
;dir1 =  '/mn/stornext/u3/shens/scratch/Eris_data/'
;dir = '/mn/stornext/u3/shens/scratch/PONOS/low_res/'
;dir =  '/mn/stornext/u3/shens/scratch/Eris2k_data/'
root =  'L90Mpc8000_hithres.'
;root = 'SE_noBH_LR.'
;root = 'dg1_new2_hi1.'
;dir = "/Volumes/3/AHFfiles/"
;filename = [15, 19, 23, 27, 31, 34, 35, 39, 43, 45, 47, 51, 55, 59]
;filename = indgen(22)*5+200
filename = indgen(400) + 1
;filename = [28, 39, 55, 90, 328]
;filename = [305]
;filename = [54]
;filename = [64, 98, 174, 253,305]
;filename = indgen(22)*16 + 160 
;filename = [36, 60]
;filename = [140]
;filename = [16, 32, 48, 64, 80, 96, 112, 144]
;filename = [32, 36, 48, 50, 59, 64, 73, 80, 92, 96, 104, 112, 120,
;128, 141, 144, 160, 168, 176, 185, 192, 205, 208, 224, 228, 240, 256]

;filename = indgen(32)*16+16
;filename = [30, 46]
nfile = n_elements(filename)
redshift = fltarr(nfile)
;outfile = dir + "g15784_nd.zinfo.dat"
;outfile = dir + "mwmc.zinfo.dat"
print, filename
;retall
for i=87, nfile-1 do begin ;change back to i=0
  ; file = dir + 'g15784_nd.' + string(filename[i], format = '(i5.5)')
  file = dir + root + string(filename[i], format = '(i5.5)')
 ; file = dir + 'dg1_new2.' + string(filename[i], format = '(i5.5)')
   read_tipsy_header, file, h1
   redshift[i] = (1.-h1.time)/h1.time
   newversion = 0
   if (h1.nstar ne 0) then begin   
      new_grp_stat, file
      print, "done file", file, "redshift = ", redshift[i]
   endif else begin 
      print, "no star formation", i
   endelse
endfor 
close, 1

;openw, 1, outfile
;for i =0, nfile-1 do begin 
;   printf, 1, filename[i], redshift[i]
;endfor 

;close,1 

end
