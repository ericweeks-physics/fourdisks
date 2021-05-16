; computeaabb		Eric R. Weeks, 7-8-14
;

function fourcross,v1,v2
    ; compute (v1 cross v2) dot z-hat
    return,v1[0,*]*v2[1,*] - v1[1,*]*v2[0,*]
end

function computeaabb,traj,unit=unit,sqrt=sqrt

v12=[traj(2,*)-traj(0,*),traj(3,*)-traj(1,*)]
v34=[traj(6,*)-traj(4,*),traj(7,*)-traj(5,*)]
v13=[traj(4,*)-traj(0,*),traj(5,*)-traj(1,*)]
v24=[traj(6,*)-traj(2,*),traj(7,*)-traj(3,*)]
v14=[traj(6,*)-traj(0,*),traj(7,*)-traj(1,*)]
v23=[traj(4,*)-traj(2,*),traj(5,*)-traj(3,*)]

if (keyword_set(unit)) then begin
	vv12=sqrt(total(v12*v12,1))
	vv13=sqrt(total(v13*v13,1))
	vv14=sqrt(total(v14*v14,1))
	vv23=sqrt(total(v23*v23,1))
	vv24=sqrt(total(v24*v24,1))
	vv34=sqrt(total(v34*v34,1))
	e=0
	v12(e,*)=v12(e,*)/vv12 & v13(e,*)=v13(e,*)/vv13
	v14(e,*)=v14(e,*)/vv14 & v23(e,*)=v23(e,*)/vv23
	v24(e,*)=v24(e,*)/vv24 & v34(e,*)=v34(e,*)/vv34
	e=1
	v12(e,*)=v12(e,*)/vv12 & v13(e,*)=v13(e,*)/vv13
	v14(e,*)=v14(e,*)/vv14 & v23(e,*)=v23(e,*)/vv23
	v24(e,*)=v24(e,*)/vv24 & v34(e,*)=v34(e,*)/vv34
endif

if (keyword_set(sqrt)) then begin
	vv12=sqrt(sqrt(total(v12*v12,1)))
	vv13=sqrt(sqrt(total(v13*v13,1)))
	vv14=sqrt(sqrt(total(v14*v14,1)))
	vv23=sqrt(sqrt(total(v23*v23,1)))
	vv24=sqrt(sqrt(total(v24*v24,1)))
	vv34=sqrt(sqrt(total(v34*v34,1)))
	e=0
	v12(e,*)=v12(e,*)/vv12 & v13(e,*)=v13(e,*)/vv13
	v14(e,*)=v14(e,*)/vv14 & v23(e,*)=v23(e,*)/vv23
	v24(e,*)=v24(e,*)/vv24 & v34(e,*)=v34(e,*)/vv34
	e=1
	v12(e,*)=v12(e,*)/vv12 & v13(e,*)=v13(e,*)/vv13
	v14(e,*)=v14(e,*)/vv14 & v23(e,*)=v23(e,*)/vv23
	v24(e,*)=v24(e,*)/vv24 & v34(e,*)=v34(e,*)/vv34
endif

xx=reform(v12(0,*)*v34(0,*)+v12(1,*)*v34(1,*))
yy=reform(v13(0,*)*v24(0,*)+v13(1,*)*v24(1,*))
zz=reform(v14(0,*)*v23(0,*)+v14(1,*)*v23(1,*))
aa=(xx+yy)/sqrt(2.0)

if (not keyword_set(unit) and not keyword_set(sqrt)) then begin
	bb=-sqrt(1.5)*zz
	; this works because:  a-hat = (+1, +1, 0) / sqrt(2)
	;                      b-hat = (+1, -1, -2) / sqrt(6)
	;                      c-hat = (+1, -1, +1) / sqrt(3)
	; (c-hat can be worked out by taking a cross b; also can
	;  see by inspection that a dot c = 0, b dot c = 0.)
	;
	; Then, note that (xx,yy,zz) dot c-hat = 0.  This can
	; be seen by writing out the expressions for xx,yy,zz;
	; for example, xx = (v12 dot v34) = (r2 - r1) dot (r4-r3)
	; where r1,r2,r3,r4 are the vector positions of the particles.
	; Bottom line:  xx,yy,zz are not independent.  In fact
	; knowing that the c component is strictly zero shows
	; that xx-yy+zz = 0, true both mathematically and also
	; true based on numerical check with simulation data.  Can't
	; get more true than that!  So, (xx,yy,zz) only has two
	; degrees of freedom, not the three that one might think.
	;
	; Therefore, rather than computing b = (xx,yy,zz)dot(b-hat),
	; I can compute (xx,yy,zz)dot(b-hat - c-hat/sqrt(2)) = 
	; (xx,yy,zz) dot (0,0,-sqrt(3/2)).
	;
	; Why bother?  Numerically, it's slightly more accurate.  It
	; has to give the same answer, but avoids depending on
	; cancellations when subtracting variables from one another.
	; See next line below for the calculation using b-hat.
	;
endif else begin
	bb=(xx-yy-2.0*zz)/sqrt(6.0)
	; this equals the above simpler expression, *iff* they
	; are regular vectors and not unit vectors
endelse

pp=fourcross(v13,v12)
qq=fourcross(v14,v12)
rr=fourcross(v13,v14)
ss=fourcross(v23,v24)

v13=0
v12=0
v14=0
v23=0
v24=0
v34=0

nx=pp+qq+rr+ss
ny=pp+qq-rr-ss
nz=pp-qq+rr-ss

nrad=sqrt(nx*nx+ny*ny+nz*nz)

na=n_elements(aa)
ix=abs(nx/nrad)
iy=abs(ny/nrad)
iz=abs(nz/nrad)
jj=fltarr(2,na)
w=where((ix ge iy) and (ix ge iz),nw)
if (nw gt 0) then begin
	jj(0,w)=iy(w)/ix(w) & jj(1,w)=iz(w)/ix(w)
endif
w=where((iy gt ix) and (iy ge iz),nw)
if (nw gt 0) then begin
	jj(0,w)=ix(w)/iy(w) & jj(1,w)=iz(w)/iy(w)
endif
w=where((iz gt ix) and (iz gt iy),nw)
if (nw gt 0) then begin
	jj(0,w)=ix(w)/iz(w) & jj(1,w)=iy(w)/iz(w)
endif
; jj has been projected onto one face of cube

anx=abs(nx) & any=abs(ny) & anz=abs(nz)
vv=(-2.0*anx+any+anz)/sqrt(6.0)
ww=(anz-any)/sqrt(2.0)

res=fltarr(18,na)
res(0,*)=aa
res(1,*)=bb
res(2,*)=atan(bb,aa)
res(3,*)=sqrt(aa*aa+bb*bb)
res(4,*)=nx
res(5,*)=ny
res(6,*)=nz
res(7,*)=nrad
res(8,*)=jj[0,*]
res(9,*)=jj[1,*]
res(10,*)=pp
res(11,*)=qq
res(12,*)=rr
res(13,*)=ss
res(14,*)=vv
res(15,*)=ww
res(16,*)=atan(vv,ww)-!DPI/6.0
res(17,*)=sqrt(vv*vv+ww*ww)

w=where(res(16,*) lt -!DPI,nw)
if (nw gt 0) then res(16,w)=res(16,w)+2.0*!DPI

return,res
end
