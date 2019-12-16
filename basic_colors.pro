pro basic_colors,black,white,red,green,blue,yellow,cyan,magenta
;
; procedure to define eight basic colors. see pvwave user guide p 14-13
;

if (!d.n_colors le 256L) then begin
rd  = [0,1,1,0,0,1,0,1,0,0,0,0]
grn = [0,1,0,1,0,1,1,0,1.0,1.0,0.7,0.3]
bl  = [0,1,0,0,1,0,1,1,0.3,0.8,1.,1.]
tvlct,fix(255*rd),fix(255*grn),fix(255*bl)
black = 0
white = 1
red = 2
green = 3
blue = 4
yellow = 5
cyan = 6
magenta=7
endif else begin
black = 0
red = 255L
green = 255L*256L
blue = 255L*256L*256L
yellow = red+green
cyan = green+blue
magenta = blue+red
white = blue+red+green

endelse


return
end
