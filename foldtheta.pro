; foldtheta.pro		Eric R. Weeks, 7-14-14
;
; code taken out of evolvefour02.pro

function foldtheta,th,en
; th = theta
; en = energy
;
; returns [theta,energy,dtheta(dt=1)]

twopi=2.0*3.14159265

dth = reform(th)-shift(reform(th),1)
dth(0) = 0.0
dth=shift(dth,-1)
w=where(dth gt 3.14159265,nw)
if (nw gt 0) then dth(w) = dth(w)-twopi
w=where(dth lt -3.14159265,nw)
if (nw gt 0) then dth(w) = dth(w)+twopi

th2=th+(twopi/3.0) & th3=th2+(twopi/3.0)
th4=[th,th2,th3]
en4=[en,en,en]; energy
dth4=[dth,dth,dth]; dtheta
w=where(th4 gt twopi/2.0,nw)
if (nw gt 0) then th4(w)=th4(w)-twopi
th5=[th4,-th4] & en5=[en4,en4] & dth5=[dth4,dth4]
result=fltarr(3,n_elements(en5))
result(0,*)=th5
result(1,*)=en5
result(2,*)=dth5

return,result
end
