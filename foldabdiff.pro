; foldab		Eric Weeks, 7-7-14
;
; taking the code out of evolvefour02.pro
;
; foldabdiff:  7-10-2019, taken from foldab.pro

function foldabdiff,ab,abdiff,themin=themin,themax=themax,_extra=eee
; _extra=eee fed to hist2d

twopi=2.0*3.14159265
cs=cos(twopi/3.0)
sn=sin(twopi/3.0)

ab(*,0)=themin
ab(*,1)=themax

tmphist = avgbin2d(ab(0,*),ab(1,*),abdiff, $
    xmin=themin,ymin=themin,xmax=themax,ymax=themax,_extra=eee)
tmphist(*,*,0) *= tmphist(*,*,1)

ab2=ab
ab2(0,*)=ab(0,*)*cs - ab(1,*)*sn
ab2(1,*)=ab(0,*)*sn + ab(1,*)*cs
ab2(*,0)=themin
ab2(*,1)=themax
tmphist2=avgbin2d(ab2(0,*),ab2(1,*),abdiff, $
    xmin=themin,ymin=themin,xmax=themax,ymax=themax,_extra=eee)
tmphist2(*,*,0) *= tmphist2(*,*,1)

ab3=ab2
ab3(0,*)= ab(0,*)*cs + ab(1,*)*sn
ab3(1,*)=-ab(0,*)*sn + ab(1,*)*cs
ab3(*,0)=themin
ab3(*,1)=themax
tmphist3=avgbin2d(ab3(0,*),ab3(1,*),abdiff, $
    xmin=themin,ymin=themin,xmax=themax,ymax=themax,_extra=eee)
tmphist3(*,*,0) *= tmphist3(*,*,1)
tmphist4=tmphist+tmphist2+tmphist3
tmphist4 = tmphist4+reverse(tmphist4,2)

return,tmphist4
end
