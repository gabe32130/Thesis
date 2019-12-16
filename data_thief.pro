pro set_plot_axes,x_norm,y_norm,xmin,xmax,ymin,ymax,xlog,ylog
;
; determines the axis scaling and orientation for orthogonal axes and 
; arbitrary orientation angles
;

basic_colors,black,white,red,green,blue,yellow

;
x_norm = fltarr(3)
y_norm = fltarr(3)
;
plot,[0,1],[0,1],xsty=4,ysty=4,/nodata,/noerase
;
; get the origin first
;
print,'click at the plot origin'
cursor,a,b,/down
x_norm(0) = a & y_norm(0) = b
oplot,[a],[b],color=red,/psym
;
; x-axis end point
;
print,'click at the x-axis end point'
cursor,a,b,/down
x_norm(1) = a & y_norm(1) = b
oplot,[a],[b],color=red,/psym
;
; y-axis end point
;
print,'click at the y-axis end point'
cursor,a,b,/down
x_norm(2) = a & y_norm(2) = b
oplot,[a],[b],color=red,/psym
;
; determine scaling from x and y-axis ranges
;
print,'enter x and y-axis ranges (xmin,xmax,ymin,ymax):'
readf,0,xmin,xmax,ymin,ymax
;
; axis types (log or linear)
;
print,'x-axis type (linear=0, log=1):'
readf,0,xlog
;
print,'y-axis type:'
readf,0,ylog
;
return
end

pro compute_coords,x,y,x_norm,y_norm,xmin,xmax,ymin,ymax,xlog,ylog
;
; given the data axes defined by x_norm,y_norm,xmin,xmax,ymin,ymax,
; this procedure converts (x,y) from normal coordinates to data coordinates
;
; first find the x and y distances of the pick point in normal coordinates
;
xaxis = sqrt((x_norm(1) - x_norm(0))^2 + (y_norm(1) - y_norm(0))^2)
x_bar = ((x - x_norm(0))*(x_norm(1) - x_norm(0)) + $
     (y - y_norm(0))*(y_norm(1) - y_norm(0)))/xaxis
;
yaxis = sqrt((x_norm(2) - x_norm(0))^2 + (y_norm(2) - y_norm(0))^2)
y_bar = ((x - x_norm(0))*(x_norm(2) - x_norm(0)) + $
     (y - y_norm(0))*(y_norm(2) - y_norm(0)))/yaxis
;
; convert to data coordinates
;
if (xlog eq 1) then begin
   x = exp(x_bar*(alog(xmax) - alog(xmin))/xaxis + alog(xmin))
endif else begin
   x = x_bar*(xmax - xmin)/xaxis + xmin
endelse
;
if (ylog eq 1) then begin
   y = exp(y_bar*(alog(ymax) - alog(ymin))/yaxis + alog(ymin))
endif else begin
   y = y_bar*(ymax - ymin)/yaxis + ymin
endelse
;
return
end

pro data_thief,xx,yy
;
; allows the user to establish the coordinate axes and pick points 
; to be read into x and y data arrays
;
basic_colors,black,white,red,green,blue,yellow
;
xt = fltarr(1000)     ; temporary data arrays 
yt = fltarr(1000)
;
; set appropriate cursor
;
;device,cursor_standard=68
device,cursor_standard=34
;
; establish the plot axes
;
set_plot_axes,x0,y0,xmx,ymx,xmin,xmax,ymin,ymax
;
; determine the desired data points using the mouse pointer
;
print,'pick the data points with the left or center mouse buttons'
print,'press the right button when done'
;
; determine pick point in normal coordinates
;
cursor,x,y,/down
oplot,[x],[y],color=red,/psym
;
; convert to data coordinates
;
compute_coords,x,y,x0,y0,xmx,ymx,xmin,xmax,ymin,ymax
print,x,",",y,", $"
xt(0) = x & yt(0) = y
;
; loop until done
;
npts=1
while (!err ne 4) do begin
   cursor,x,y,/down
   oplot,[x],[y],color=red,/psym
   compute_coords,x,y,x0,y0,xmx,ymx,xmin,xmax,ymin,ymax
   print,x,",",y,", $"
   xt(npts) = x & yt(npts) = y
   npts = npts + 1
endwhile
;
; define and fill the output data array
;
npts = npts - 1
xx = fltarr(npts)
yy = fltarr(npts)
for i=0,npts-1 do begin
   xx(i) = xt(i) & yy(i) = yt(i)
endfor
;
; reset the cursor to the standard pv-wave cursor
;
device,cursor_standard=34
;
return
end

