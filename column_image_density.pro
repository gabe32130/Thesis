msolunit = 1.078e17
rho_mean = 0.0462 
lunit = 90.
kpcunit = 90000.
gmperccunit = 1.00169e-29
kmpersecunit = 2269.98
dhubble0 = 2.8944 ;code
vunit = 2269.98
omega_l = 0.76
omega_m = 0.24
m_pro = 1.6726d-24
Ao  = 15.9994
Ac = 12.0107
Afe = 55.845
Asi = 28.0855
Amg = 24.30506

kpctocm = 3.086d21
msun = 1.99d33
mpro = 1.67e-24
dir = '/mn/stornext/u3/shens/scratch/Eris_data/'
dir2 = '/mn/stornext/u3/shens/scratch/Eris_amiga/'
dir1 = '/mn/stornext/d17/extragalactic/personal/gabrierg/rprocess_halos/'
name = 'L90Mpc8000_hithres.'
filename = [400]

tipfile = dir + name + string(filename[0], format = '(i5.5)')
rvirfile = dir2 + name + string(filename[0], format = '(i5.5)')+'.rvir' ;where?
gridfile = dir1 + name + string(filename[0], format = '(i5.5)')+'.sigma' ;same as .sigstar_align_1D_bash?

;--------------------------------------------------------------------;
;item1 = ['log!l10!n (!7R!3!lSFR!n/M'+sunsymbol()+' kpc!u-2!n yr!u-1!n)']
item1 = ['log!l10!n (!7R!3!l*!n/M'+sunsymbol()+' kpc!u-2!n)']
item2 = ['log!l10!n (!7R!3!lEu!n/M'+sunsymbol()+' kpc!u-2!n)']
item3 = ['log(<f_d>)']
;tname1 = ['-3.0', '-2.0', '-1.0', '0.0', '1.0']
;tname1 = ['5.0', '6.0', '7.0', '8.0']
tname1 = ['4.5','5.5', '6.5', '7.5', '8.5', '9.5', '10.5']
;tname1 = ['6.0', '6.5', '7.0', '7.5', '8.0', '8.5']
;tname1 = ['-2.0', '-1.0', '0.0', '1.0']

ndiv = 6

read_tipsy_header, tipfile, h
z = (1.0-h.time)/h.time
print, "redshift", z
;goto, plot

nd1 = 0
nd2 = 0 
grpid = 0
iproj = 0
print, "Read in sig1 file"
close, 1
openr, 1, gridfile
readf, 1, nd1, nd2, grpid, iproj, d1l, d2l,xc, yc,zc 
d1binsp = fltarr(nd1)
d2binsp = fltarr(nd2)
sigma = fltarr(nd1, nd2)
eusigma = fltarr(nd1,nd2)
deusigma = fltarr(nd1, nd2) 
frac2D = fltarr(nd1, nd2)

readf, 1, d1binsp 
readf, 1, d2binsp 
readf, 1, sigma
readf, 1, eusigma
readf, 1, deusigma
readf, 1, frac2D
close, 1

;--------------------------------------------------------------------;
openr, 1, rvirfile
 nn = 0
 readf, 1, nn
 rvir_all = fltarr(nn)
 readf, 1, rvir_all
 close, 1
rvir = rvir_all[1]*lunit*1000./(1+z)
print, "virial radius from host", rvir

plot:
rmax = rvir
rmin = -rvir
xp = findgen(101)*((rmax-rmin)/100.) + rmin 
yp = sqrt(rvir^2. - xp^2.)

;max1 = 21.0 
;min1 = 16.0
;min1 = -3.0 
;max1 = 1.0 
;max = 18.0 
;min = 11.0
;sig_sfr = sig2/1e8
 
max = max(alog10(sigma))
;max = 1.0
min = max - 5.75
;max = 8.0
;min =5.0
max1 = max(alog10(eusigma))
min1 = max1 - 5.0 

max2 = max(alog10(deusigma))
min2 = max2 - 5.0 

max4 = max(alog10(frac2D))
min4 = max4 - 7.5

print, "min, max", min, max 

loadct, 39
image1 = bytscl(alog10(sigma), max=max, min=min, top=255)
image2 = bytscl(alog10(eusigma), max=max1, min=min1, top=255)
image3 = bytscl(alog10(deusigma), max=max2, min=min2, top=255)
image4 = bytscl(alog10(frac2D), max=max4, min=min4, top=255)
;xr = 250
;idx1 = where(d1binsp ge -xr and d1binsp le xr)
;idx2 = where(d2binsp ge -xr and d2binsp le xr)

;n1 = n_elements(idx1)
;n2 = n_elements(idx2)

;sig1_cut = sig2(idx1[0]:idx1[n1-1], idx2[0]:idx2[n2-1])
;image1_cut = bytscl(alog10(sig1_cut), max=max, min=min, top=255)
 
;--------------------------------------------------------------------;
if (1 eq 0) then begin 
write: 
   outfile = 'sigma_gas_z2_250k.dat'
   close, 1 
   openw, 1, outfile
   printf, 1, min, max
   printf, 1, n1, n2 
   printf, 1, image1_cut
   close, 1 
endif 

!p.multi = [0, 1, 1]
xnulltick = replicate(' ', 5)
ynulltick = xnulltick

xwidth = 8.5
ywidth = 11.0 
xs = 4.0
ys = 4.0
xoffset = (xwidth - xs)*0.5
yoffset = (ywidth - ys)*0.5 

xr = 50.0 ;scale
xran = [-xr, xr]
yran = [-xr, xr]

;--------------------------------------------------------------------;
plottops = 1
if (plottops eq 1) then begin 
   entry_device = !d.name
   set_plot,'PS'
;   device, bits_per_pixel = 8,/portrait, color = 1,filename ='Eris_sigma_stars_xy50.eps', xsize = xs, ysize =ys, xoffset=xoffset, yoffset=yoffset, /inches 

;   device, bits_per_pixel = 8,/portrait, color = 1,filename ='Eris_eusigma_stars_xy50.eps', xsize = xs, ysize =ys, xoffset=xoffset, yoffset=yoffset, /inches

   device, bits_per_pixel = 8,/portrait, color = 1,filename ='Eris_deusigma_stars_xz50.eps', xsize = xs, ysize =ys, xoffset=xoffset, yoffset=yoffset, /inches

;   device, bits_per_pixel = 8,/portrait, color = 1,filename ='Eris_frac_stars_xz50.eps', xsize = xs, ysize =ys, xoffset=xoffset, yoffset=yoffset, /inches

   tk = 3
   csize = 1.0
endif else begin 
   tk = 1
   window, 0, retain = 2
endelse 

xw = [0.12, 0.95]
yw = [0.12, 0.95]

position0 = [xw[0], yw[0], xw[1], yw[1]]
barpos0 = [xw[0]+0.03, yw[0]+0.07, xw[1]-0.03, yw[0]+0.10]

loadct, 0
bw = 0
cgplot, d1binsp, d2binsp, /nodata, xrange = xran, yrange = yran, charsize = csize, charthick = tk, ytitle = ' ', axiscolor = bw, xstyle = 1, ystyle = 1, xtitle = ' ', position = position0
;loadct, 39 ; xtitle = 'x [kpc]',, position = position0,
loadct, 39
cgimage, image4, /noerase, /overplot
;retall
loadct, 0
cgplot, d1binsp, d2binsp, /nodata, xrange = xran, yrange = yran, charsize = csize, charthick = tk, /noerase, axescolor = 254, xtickname = xnulltick, ytickname = ynulltick, xthick = tk, ythick = tk, xstyle = 1, ystyle = 1, position = position0

legend, item3[0], /top,/right, box=0, charthick =tk, textcolor = 255, charsize =0
;legend, itemt[0], /top,/left, box=0, charthick =tk, textcolor = 254, charsize =0
loadct, 39
colorbar, ticknames = tname1,divisions = ndiv, charsize = 0, color = 255, charthick = tk,position = barpos0;, title = item1[0]

if (plottops eq 1) then begin 
   device, /close_file
   set_plot, entry_device
endif

end
