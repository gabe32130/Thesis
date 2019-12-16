pro histlike, x, y, xout, yout 
xmin = min(x)
xmax = max(x)
nn = n_elements(x) 
dx = abs(x[1]-x[0])
xbins = findgen(nn)*dx + xmin 

xout = fltarr(2*nn)
yout = fltarr(2*nn)

xout[2*lindgen(nn)] = x
xout[2*lindgen(nn)+1] = x 

yout[2*lindgen(nn)] = y;  [y, 0] 
yout[2*lindgen(nn)+1] = y; [y, 0] 
yout =  shift(yout, 1)

end 
