PRO rtipsy, filename, $
                    header, gp, dp, sp

header = {time:0.0d0,$
          nbodies:0l,$
          ndim:0l,$
          nsph:0l,$
          ndark:0l,$
          nstar:0l,$
          dummy:0l}

openr,lun,filename,/get_lun,/xdr
readu,lun,header

;; Define data structure

gas_particle = { mass:    0.0, $
                 pos:     fltarr(3), $
                 vel:     fltarr(3), $
                 rho:     0.0, $
                 temp:    0.0, $
                 eps:     0.0, $
                 metals:  0.0, $
                 phi:     0.0 }

dark_particle = { mass:    0.0, $
                  pos:     fltarr(3), $
                  vel:     fltarr(3), $
                  eps:     0.0, $
                  phi:     0.0 }


star_particle = { mass:    0.0, $
                 pos:     fltarr(3), $
                 vel:     fltarr(3), $
                 metals:  0.0, $
                 tform:   0.0, $
                 eps:     0.0, $
                 phi:     0.0 }

if header.nsph ne 0 then begin
    gp = replicate(gas_particle,header.nsph)
    readu, lun, gp 
endif

if header.ndark ne 0 then begin
    dp = replicate(dark_particle,header.ndark)
    readu, lun, dp
endif
if header.nstar ne 0 then begin
    sp = replicate(star_particle,header.nstar)
    readu, lun, sp
endif

free_lun, lun
END
