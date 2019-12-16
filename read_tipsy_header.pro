pro read_tipsy_header,filename,tipsyheader

tipsyheader = {time:0.0d0,$
               nbodies:0l,$
               ndim:0l,$
               nsph:0l,$
               ndark:0l,$
               nstar:0l,$
               dummy:0l}

openr,lun,filename,/get_lun,/xdr
readu,lun,tipsyheader
free_lun,lun

return
end
