; fourtau4.pro
; started 7-7-14 by Eric Weeks, modified from my circtau3.pro
;
; v4: I'm calling this v4 since it was v3 for circtau, see note below!
; v5: generalize since now I have four different theta's
; v5 8-4-14:  adding in period keyword
;
;
; circtau3.pro
; started 9-6-13 by Eric Weeks, modified from Gary Hunter's
; 'ghtclifetimes2.pro'
;
; (I'm calling my version 3 since he has version 2; there is no
; circtau2.pro)

FUNCTION fourtau5,theta0,rad,drad=drad,dtheta=dtheta,debug=debug, $
   period=period

twopi = 2.0*3.14159265358

if (not keyword_set(drad)) then drad = 6.0
if (not keyword_set(dtheta)) then dtheta = (15.0)*(twopi/360.0)
;  = 15 degrees

theta=reform(theta0)+twopi/6.0
w=where(theta gt twopi*0.5,nw)
if (nw gt 0) then theta(w) = theta(w) - twopi

if (keyword_set(debug)) then begin
	aa=rad*cos(theta) & bb=rad*sin(theta)
	plot,aa,bb,/iso,psym=3
endif

nrow = n_elements(theta)
sh = bytarr(nrow)
flag=bytarr(nrow)
angl = twopi/3.0
t1 = -angl + dtheta
t2 = -dtheta
w=where((rad gt drad) and (theta gt t1) and (theta lt t2),nw)
if (nw gt 0) then sh(w) = 1
w=where((theta gt t1) and (theta lt t2),nw)
if (nw gt 0) then flag(w) = 1
if (keyword_set(debug)) then begin
	oplot,[0,cos(t1)]*20,[0,sin(t1)]*20
	oplot,[0,cos(t2)]*20,[0,sin(t2)]*20
endif

t1 = dtheta
t2 = angl - dtheta
w=where((rad gt drad) and (theta gt t1) and (theta lt t2),nw)
if (nw gt 0) then sh(w) = 2
w=where((theta gt t1) and (theta lt t2),nw)
if (nw gt 0) then flag(w) = 2
if (keyword_set(debug)) then begin
	oplot,[0,cos(t1)]*20,[0,sin(t1)]*20
	oplot,[0,cos(t2)]*20,[0,sin(t2)]*20
endif

t1 = -angl - dtheta
t2 = angl + dtheta
w=where((rad gt drad) and ((theta lt t1) or (theta gt t2)),nw)
if (nw gt 0) then sh(w) = 3
w=where(((theta lt t1) or (theta gt t2)),nw)
if (nw gt 0) then flag(w) = 3
if (keyword_set(debug)) then begin
	oplot,[0,cos(t1)]*20,[0,sin(t1)]*20
	oplot,[0,cos(t2)]*20,[0,sin(t2)]*20
endif

if (keyword_set(debug)) then begin
	th=findgen(1000)/1000.*(twopi+0.01)
	aa=drad*cos(th) & bb=drad*sin(th)
	oplot,aa,bb,psym=3
endif


u = [uniq(sh)]
nu=n_elements(u)
sh2 = sh
; u = the LAST element this is true for

ndu=0
if (nu gt 1) then begin
	for i=0,nu-1 do begin
		state0 = sh(u(i))
		if (state0 eq 0) then begin
			if (i eq 0) then begin
				; the initial state is 0, make it determinate.
				; later we're going to ignore this event anyway!
				w=where(sh gt 0)
				sh2(0:w(0)) = sh(w(0))
			endif else begin
				state2 = sh(u(i-1))
				sh2(u(i-1)+1:u(i)) = state2
			endelse
		endif
	endfor

	; patch added 7-15-14, to deal with trailing zero's
	if (sh2(nrow-1) eq 0) then begin
		sh2(u(nu-1):*)=state0
	endif

	w=where(sh2 eq 0,nw)
	if (nw gt 0) then message,'why is anything still zero?'
	w=where((sh gt 0) and (sh ne sh2),nw)
	if (nw gt 0) then message,'why did something nonzero change?'

	u2=uniq(sh2)
	du=shift(u2,-1)-u2
	ndu = n_elements(du)
	if (ndu ge 3) then res = du(0:ndu-2) else res = 0.0
endif else begin
	res = 0.0
endelse

if (ndu gt 0) then period = (1.0*nrow)/(1.0*ndu) else period=0

if ((nu gt 1) and keyword_set(debug)) then begin
	w=where(du eq 1,nw)
	if (nw gt 0) then begin
		shiftdu = shift(u2,-1)
		print,u2(w(0)),shiftdu(w(0)),du(w(0))
		u3=u2(w(0))
		print,sh2(u3-3:u3+3)
		print,sh(u3-3:u3+3)
		print,flag(u3-3:u3+3)
		print,theta(u3-3:u3+3)
		print,rad(u3-3:u3+3)
	endif
endif

return, res

end
