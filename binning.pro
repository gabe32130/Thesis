pro binning, x, y, y_test, xmin, xmax, n, xbins, median_value, percentile_up, percentile_down, number,log=log, nozero=nozero
  dx=(xmax-xmin)/(n-1)
  xbins = findgen(n)*dx + xmin 
 
  median_value = fltarr(n)
  percentile_up = fltarr(n)
  percentile_down = fltarr(n)
  number = lonarr(n)
 
  if(keyword_set(log)) then begin 
     hist =histogram(alog10(x), binsize = dx, min=xmin, max=xmax,reverse_indices = ridx1)
  endif else begin
     hist =histogram(x, binsize = dx, min=xmin, max=xmax,reverse_indices = ridx1)
  endelse 
  
  for j=0, n -2 do begin 
     if (ridx1[j] ne ridx1[j+1]) then begin 
        idx= ridx1[ridx1[j]:ridx1[j+1]-1]
        
        if(keyword_set(nozero)) then begin 
           ii = where(y_test(idx) ne 0.0)
           if (ii[0] ne -1) then begin 
              number[j]=n_elements(ii)
              median_value[j]=median(y[idx[ii]])
              percentile_up[j] = percentile_84(y[idx(ii)])
              percentile_down[j] = percentile_16(y[idx(ii)])
            
           endif 
        endif else begin 
           median_value[j]=median(y[idx])
           percentile_up[j] = percentile_84(y[idx])
           percentile_down[j] = percentile_16(y[idx])
           
        endelse 
     endif 
  endfor 
return
end 
