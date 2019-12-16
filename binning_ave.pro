pro binning_ave, x, y, xmin, xmax, n, xbins, ave_value, rms_value, log=log
  dx=(xmax-xmin)/(n-1)
  xbins = findgen(n)*dx + xmin 
  ave_value = fltarr(n)
  rms_value = fltarr(n)
  sum = fltarr(n)
  sum2 = fltarr(n)

  if(keyword_set(log)) then begin 
     hist =histogram(alog10(x), binsize = dx, min=xmin, max=xmax,reverse_indices = ridx1)
  endif else begin
     hist =histogram(x, binsize = dx, min=xmin, max=xmax,reverse_indices = ridx1)
  endelse 
 
  for j=0, n -2 do begin 
     if (ridx1[j] ne ridx1[j+1]) then begin 
        idx= ridx1[ridx1[j]:ridx1[j+1]-1]
    
        sum[j]=total(y(idx))
        ave_value[j]=sum[j]/hist[j]
        sum2[j]=total((y[idx]-ave_value[j])^2)
        rms_value[j]=sqrt(sum2[j]/(hist[j]-1))
     endif 
  endfor
end 
