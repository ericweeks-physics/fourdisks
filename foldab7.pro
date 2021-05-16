; foldab		Eric Weeks, 7-7-14
;
; taking the code out of evolvefour02.pro
;
; version 7:  for evolvefour07.pro, fix max & min

function foldab7,ab,themin=themin,themax=themax,_extra=eee
; _extra=eee fed to hist2d

twopi=2.0*3.14159265
cs=cos(twopi/3.0)
sn=sin(twopi/3.0)

ab(0:1,0) = themin & ab(0:1,1) = themax
tmphist = hist2d(ab(0,*),ab(1,*),_extra=eee)
ab2=ab
ab2(0,*)=ab(0,*)*cs - ab(1,*)*sn
ab2(1,*)=ab(0,*)*sn + ab(1,*)*cs
ab2(0:1,0) = themin & ab2(0:1,1) = themax
tmphist2=hist2d(ab2(0,*),ab2(1,*),_extra=eee)
ab3=ab2
ab3(0,*)= ab(0,*)*cs + ab(1,*)*sn
ab3(1,*)=-ab(0,*)*sn + ab(1,*)*cs
ab3(0:1,0) = themin & ab3(0:1,1) = themax
tmphist3=hist2d(ab3(0,*),ab3(1,*),_extra=eee)
tmphist4=tmphist+tmphist2+tmphist3
tmphist4 = tmphist4+reverse(tmphist4,2)

return,tmphist4
end
