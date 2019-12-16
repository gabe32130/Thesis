pro align,gas,dark,stars,limit, jjx, jjy, jjz

J=0.
Jxs=0.
Jys=0.
Jzs=0.
Js=0.
mt=0.

r = where(sqrt(gas.pos[0]*gas.pos[0]+gas.pos[1]*gas.pos[1]+gas.pos[2]*gas.pos[2])  lt limit)
mt = total(gas[r].mass)
Jxs = total(gas[r].mass*(gas[r].pos[1]*gas[r].vel[2]-gas[r].pos[2]*gas[r].vel[1]))/mt
Jys = total(gas[r].mass*(gas[r].pos[2]*gas[r].vel[0]-gas[r].pos[0]*gas[r].vel[2]))/mt
Jzs = total(gas[r].mass*(gas[r].pos[0]*gas[r].vel[1]-gas[r].pos[1]*gas[r].vel[0]))/mt
Js=sqrt(Jxs*Jxs+Jys*Jys+Jzs*Jzs)

if Js gt 0. then begin
           jjx=Jxs/Js
           jjy=Jys/Js
           jjz=Jzs/Js
           costh=jjz
           sinth=sqrt(1.0-jjz*jjz)
if  sinth gt 0.0 then begin
              sinph=jjy/sinth
              cosph=jjx/sinth
endif
endif 
if sinth le 0.  then begin
           cosph = 1.0
           sinph = 0.0
endif

        ax=costh*cosph
        bx=costh*sinph
        cx=-sinth
        ay=-sinph
        by=cosph
        cy=0.0
        az=sinth*cosph
        bz=sinth*sinph
        cz=costh

; /**** translate star particles 

           txs=stars.pos[0]
           tys=stars.pos[1]
           tzs=stars.pos[2]
           stars.pos[0]=(ax*txs+bx*tys+cx*tzs)
           stars.pos[1]=(ay*txs+by*tys+cy*tzs)
           stars.pos[2]=(az*txs+bz*tys+cz*tzs)

           txs=stars.vel[0]
           tys=stars.vel[1]
           tzs=stars.vel[2]
           stars.vel[0]=(ax*txs+bx*tys+cx*tzs)
           stars.vel[1]=(ay*txs+by*tys+cy*tzs)
           stars.vel[2]=(az*txs+bz*tys+cz*tzs)
; /**** translate gas particles 

           txs=gas.pos[0]
           tys=gas.pos[1]
           tzs=gas.pos[2]
           gas.pos[0]=(ax*txs+bx*tys+cx*tzs)
           gas.pos[1]=(ay*txs+by*tys+cy*tzs)
           gas.pos[2]=(az*txs+bz*tys+cz*tzs)

           txs=gas.vel[0]
           tys=gas.vel[1]
           tzs=gas.vel[2]
           gas.vel[0]=(ax*txs+bx*tys+cx*tzs)
           gas.vel[1]=(ay*txs+by*tys+cy*tzs)
           gas.vel[2]=(az*txs+bz*tys+cz*tzs)
; translate dark matter particles  
           txs=dark.pos[0]
           tys=dark.pos[1]
           tzs=dark.pos[2]
           dark.pos[0]=(ax*txs+bx*tys+cx*tzs)
           dark.pos[1]=(ay*txs+by*tys+cy*tzs)
           dark.pos[2]=(az*txs+bz*tys+cz*tzs)

           txs=dark.vel[0]
           tys=dark.vel[1]
           tzs=dark.vel[2]
           dark.vel[0]=(ax*txs+bx*tys+cx*tzs)
           dark.vel[1]=(ay*txs+by*tys+cy*tzs)
           dark.vel[2]=(az*txs+bz*tys+cz*tzs)

return
end
