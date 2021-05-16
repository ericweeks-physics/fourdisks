; fourmd2.pro -- started 4-26-2020 by Eric Weeks, from evolvefour07md.pro
;
; v3:  pre-calculate everything, then interpolate to fit dt
; v3.3:  make sure everything is double precision
; v4:  fix interpolation routine (match "redomd4.pro")
; 
; this version will be event-driven.

function whenwillhitwall,pinfo,wallpos

; see my notes from 4-25-2020 (MSWord document, black notebook)
pvec = pinfo[0:1]
vvec = pinfo[2:3]
r2 = pvec[0]*pvec[0]+pvec[1]*pvec[1]
v2 = vvec[0]*vvec[0]+vvec[1]*vvec[1]; always positive
pdv = pvec[0]*vvec[0]+pvec[1]*vvec[1]; dot product
dr2 = wallpos - r2;  this should always be positive
dt1 = (-pdv + sqrt(pdv*pdv + dr2*v2))/v2

; we're going to catch this problem elsewhere
; if (dt1 lt 0) then message,'problem with dt1'

return,dt1
end

; ============================================================

function whentwodiskswillhit,circs,wa,wb

; when will disks wa and wb collide?
dx = circs[0,wa]-circs[0,wb]
dy = circs[1,wa]-circs[1,wb]
q = [dx,dy];  relative position
wx = circs[2,wa]-circs[2,wb]
wy = circs[3,wa]-circs[3,wb]
w = [wx,wy];  relative velocity
qw = q[0]*w[0]+q[1]*w[1]; dot product
q2 = q[0]*q[0]+q[1]*q[1]; dot product
wsqr = w[0]*w[0]+w[1]*w[1]; dot product
; the "4" in the next line is because each disk has radius 1.
; When they contact, their center-center separation is 2.  The
; 2 becomes squared in the next line, so thus 4.  Thus today,
; August 22, 2020, I have verified that the disks have radius 1.
; See my calculation in my black lab notebook as of today.
disc = qw*qw + (4.0-q2)*(wsqr)
dt1 = 9e9;  no collision (will collide with wall first, presumably)
if (disc ge 0.0) then begin
	temp = sqrt(disc)
	dt1a = (-qw + temp) / wsqr
	dt1b = (-qw - temp) / wsqr
	temp = [dt1a,dt1b]
	w1=where(temp gt 0,nw1)
	if (nw1 gt 0) then begin
		if (nw1 eq 1) then begin
			dt1 = temp[w1]
		endif else begin
			dt1 = min(temp)
		endelse
	endif
endif

return,dt1[0]; make sure it's a scalar
end

; ============================================================

function updateschedule,circs,wallpos,ndisks

result=dblarr(ndisks*(ndisks+1)/2)
for i=0,ndisks-1 do result[i] = whenwillhitwall(circs[*,i],wallpos)
count=ndisks
for i=0,ndisks-2 do begin
	for j=i+1,ndisks-1 do begin
		result[count] = whentwodiskswillhit(circs,i,j)
		count++
	endfor
endfor

return,result
end


function mdcalc,circs,wallpos,sqrtwallpos,totaltime,ndisks
; updates 'circs' with MD step

result=dblarr(ndisks*4+2,20000000L)
count=0L
abmat = bytarr(2,6)
abmat[0,*] = [0,0,0,1,1,2]
abmat[1,*] = [1,2,3,2,3,3]
; to later compute distance between all pairs of particles

sched = updateschedule(circs,wallpos,ndisks)
oldsched=sched
dt1 = min(sched,mi)
dt2 = double(totaltime);  how much time left to simulate this function-call
newcircs = circs
timesum = 0.0d

while (dt1 lt dt2) do begin
	; advance disks to the collision point
	newcircs(0:1,*) += newcircs(2:3,*)*dt1
	timesum += dt1
	result(0:ndisks*4-1,count) = newcircs
	result(ndisks*4+0,count) = mi
	;if (mi ge ndisks) then begin
	;	result(ndisks*4+1,count)=abmat(0,mi-ndisks)
	;	result(ndisks*4+2,count)=abmat(1,mi-ndisks)
	;endif
	result(ndisks*4+1,count) = timesum
	count++
	; print,dt1,mi,dt2
	if (mi lt ndisks) then begin
		; disk mi is colliding with wall, adjust its velocity
		vvec = newcircs[2:3,mi]
		tanvec = [-newcircs[1,mi],newcircs[0,mi]]
		tann = sqrt(tanvec[0]*tanvec[0]+tanvec[1]*tanvec[1])
		tanvec /= tann;  now it's a unit vector
		relveltan = tanvec * (vvec[0]*tanvec[0]+vvec[1]*tanvec[1])
		relvelperp = vvec - relveltan
		newcircs[2:3,mi] -= relvelperp*2
	endif else begin
		; disks wa and wb collide, adjust their velocities
		wa = abmat[0,mi-ndisks]
		wb = abmat[1,mi-ndisks]
		wx = newcircs[2,wa]-newcircs[2,wb]
		wy = newcircs[3,wa]-newcircs[3,wb]
		wvec = [wx,wy];  relative velocity
		tanx = (newcircs[1,wa] - newcircs[1,wb])
		tany = (newcircs[0,wb] - newcircs[0,wa])
		tann = sqrt(tanx*tanx+tany*tany)
		tanvec = [tanx,tany] / tann
		relveltan = tanvec * (wvec[0]*tanvec[0]+wvec[1]*tanvec[1])
		relvelperp = wvec - relveltan
		newcircs[2:3,wa] -= relvelperp
		newcircs[2:3,wb] += relvelperp
	endelse

	; done moving disk(s)
	oldoldsched=oldsched
	oldsched=sched
	sched = updateschedule(newcircs,wallpos,ndisks)
	if (mi ge ndisks) then sched[mi] = 9e9
	; if (mi lt ndisks) then --> could potentially hit same
	; wall again if trajectory is nearly tangent to wall
	dt2 -= dt1; now we have less time we need to simulate
	dt1 = min(sched,mi); get new minimum
	if (dt1 lt 0) then message,'dt1 is negative'
	if ((count mod 100000) eq 0) then begin
		print,count,' time remaining:',dt2,dt2+timesum
	endif
endwhile

return,result(*,0L:count-1L)
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function getemaxmin,theav
	twopi = 2.0*3.14159265
	thetas = theav(0,*)
	thetas = thetas mod twopi/3.0
	w=where(thetas lt 0,nw)
	if (nw gt 0) then thetas(w)=thetas(w)+twopi/3.0
	theta0 = twopi/6.0 & dtheta0 = twopi/72.0
	w=where(thetas gt (theta0-dtheta0) and thetas lt (theta0+dtheta0))
	emin = total(theav(1,w)*theav(5,w))
	emin = emin / total(theav(5,w))
	w=where(thetas lt (dtheta0) or thetas gt (twopi/3.0-dtheta0),nw)
	if (nw eq 0) then begin
		emax = 0.0
	endif else begin
		emax = total(theav(1,w)*theav(5,w))
		emax = emax / total(theav(5,w))
	endelse
	return,[emax,emin]
end

function findenergybarrier,thetahist,doplot=doplot
	w=where(thetahist(1,*) gt 0)
	freeenergy=-alog(thetahist(1,w))
	w=where(thetahist(0,*) gt -3.0 and thetahist(0,*) lt 3.0)
	tmpmin=min(freeenergy(w),mi,max=tmpmax)

	w=where(thetahist(3,*) gt 0)
	freeenergy2=-alog(thetahist(3,w))
	w=where(thetahist(0,*) gt -3.0 and thetahist(0,*) lt 3.0)
	tmpmin2=min(freeenergy2(w),mi,max=tmpmax2)

	w=where(thetahist(5,*) gt 0)
	freeenergy3=-alog(thetahist(5,w))
	w=where(thetahist(0,*) gt -3.0 and thetahist(0,*) lt 3.0)
	tmpmin3=min(freeenergy3(w),mi,max=tmpmax3)

	if (keyword_set(doplot)) then begin
		themax=((tmpmax2-tmpmin2) > (tmpmax-tmpmin))
		themax=((tmpmax3-tmpmin3) > themax)
		plot,thetahist(0,*),freeenergy-tmpmin,xr=[-1,1]*3.14159/2,/xs, $
			yr=[0,themax],pos=[0.1,0.55,0.9,0.95]
		oplot,thetahist(0,*),freeenergy2-tmpmin2,l=1
		oplot,thetahist(0,*),freeenergy3-tmpmin3,l=2
	endif

	return,[tmpmax-tmpmin,tmpmax2-tmpmin2,tmpmax3-tmpmin3]
end


function fourmd4,epsilon,circs,maxtime, $
    noinit=noinit,seed=seed,filename=filename,trajsave=trajsave, $
    noplot=noplot,nopic=nopic,dt=dt

version=4.0
ndisks=4
if (not keyword_set(dt)) then dt=1d-2
if (not keyword_set(noinit)) then begin
	;below code for three disks:
	;inv3=1.0/sqrt(3.0)
	;circs=[[0,2*inv3],[1,-inv3],[-1,-inv3]]
	;
	;below code for four disks:
	inv4=(1.0d + epsilon)*sqrt(2.0d)
	circs = dblarr(4,4);  [x,y,vx,vy]
	circs(0:1,*)=[[0,inv4],[0,-inv4],[inv4,0],[-inv4,0]]
	flag = 1b
	while (flag eq 1b) do begin
		theta=randomu(seed,4)*3.14159265*2.0
		circs(2,*)=cos(theta)
		circs(3,*)=sin(theta)
		; random direction, unit length initial velocities

		; let's try to minimize angular momentum
		angmomentum = 0.0
		e=0 & angmomentum += circs(0,e)*circs(3,e)-circs(1,e)*circs(2,e)
		e=1 & angmomentum += circs(0,e)*circs(3,e)-circs(1,e)*circs(2,e)
		e=2 & angmomentum += circs(0,e)*circs(3,e)-circs(1,e)*circs(2,e)
		sintheta = angmomentum/inv4
		if (abs(sintheta) lt 1.0) then begin
			theta=asin(sintheta)
			flag = 0b
			e=3
			circs(2,e)=cos(theta)
			circs(3,e)=sin(theta)
			angmomentum += circs(0,e)*circs(3,e)-circs(1,e)*circs(2,e)
		endif 
	endwhile
endif
print,'angular momentum = ',angmomentum

largerad=3.0d + epsilon
sqrtwallpos = largerad-1.0d
wallpos = (sqrtwallpos)^2

flag=0b
plotflag=0b
twopi=2.0*3.14159265
mm=40.
mm2=1.6

allcircs = mdcalc(circs,wallpos,sqrtwallpos,maxtime,ndisks)
maxicount = 1000000L
dtsmall = dt*0.1d
time0 = min(allcircs(-1,*))
allcircs(-1,*) -= time0
time0 = min(allcircs(-1,*),max=timemax)
timeoffset = (time0 mod dt)
count=0L
nel = n_elements(allcircs(*,0))
x = allcircs(-1,*)
forhisto=fltarr(ndisks*2+1,maxicount)
result=fltarr(ndisks*2+1,(maxtime/dt)+210000)
while (time0 lt timemax) do begin
	time1 = time0 + maxicount * dtsmall
	print,'time = ',time0,' / ',timemax,'   epsilon = ',epsilon
	thetimes = dindgen(maxicount)*dtsmall+time0
	for i=0,ndisks-1 do begin
		y = allcircs(0+i*4,*)
		forhisto(0+i*2,*) = interpol(y,x,thetimes)
		y = allcircs(1+i*4,*)
		forhisto(1+i*2,*) = interpol(y,x,thetimes)
	endfor
	if (time1 ge timemax) then begin
		; don't bother saving this data
	endif else begin
		; so now 'forhisto' is full of data to histogram
		; next five lines: save the data every 'dt'
		;w2=where(((thetimes-timeoffset) mod dt) le dtsmall*0.999, nw2)
		w2=lindgen(maxicount/10L)*10L
		nw2=n_elements(w2)
		result(*,count:count+nw2-1) = forhisto(*,w2)
		count += nw2

		; next long bit is for histograms, accumulation of data, etc
		ab=computeaabb(forhisto)
		tmphist4 = foldab7(ab,themin=-mm,themax=mm,bin=0.05)
		; kludge to get right limits:
		ab(8:9,0) = 0.0 & ab(8:9,1)=1.0
		tmp3d=avgbin2d(ab(8,*),ab(9,*),ab(7,*),xbin=0.002,ybin=0.002)
		tmp3d(*,*,0)=tmp3d(*,*,0)*tmp3d(*,*,1)
		thetas = foldtheta(ab(2,*),forhisto(ndisks*2,*))
		; force min & max happy:
		thetas(0,0)=-(twopi/2.0) & thetas(0,1)=+(twopi/2.0)
		tmpav=avgbin(thetas(0,*),thetas(1,*),binsize=0.01,/quiet)
		tmpav(1,*)=tmpav(1,*)*tmpav(5,*)
		tmpav(2,*)=tmpav(2,*)*tmpav(5,*)
		tmpdiff=avgbin(thetas(0,*),thetas(2,*)^2,binsize=0.01,/quiet,/more)
		bladiffn = (ab(0,*)-shift(ab(0,*),1))^2+(ab(1,*)-shift(ab(1,*),1))^2
		bladiffn = sqrt(bladiffn(1:*))
		tmpabdiff = foldabdiff(ab(0:1,1:*),bladiffn, $
			themin=-mm,themax=mm,xbin=0.05,ybin=0.05)

		thetaq = foldtheta(ab(16,*),forhisto(ndisks*2,*))
		; force min & max happy:
		thetaq(0,0)=-(twopi/2.0) & thetaq(0,1)=+(twopi/2.0)
		tmpavq=avgbin(thetaq(0,*),thetaq(1,*),binsize=0.01,/quiet)
		tmpavq(1,*)=tmpavq(1,*)*tmpavq(5,*)
		tmpavq(2,*)=tmpavq(2,*)*tmpavq(5,*)
		tmpdiffq=avgbin(thetaq(0,*),thetaq(2,*)^2,binsize=0.01,/quiet,/more)

		ab=computeaabb(forhisto,/unit)
		tmphist5 = foldab7(ab,themin=-mm2,themax=mm2,bin=0.002)
		thetas2 = foldtheta(ab(2,*),forhisto(ndisks*2,*))
		; force min & max happy:
		thetas2(0,0)=-(twopi/2.0) & thetas2(0,1)=+(twopi/2.0)
		tmpavu=avgbin(thetas2(0,*),thetas2(1,*),binsize=0.01,/quiet)
		tmpavu(1,*)=tmpavu(1,*)*tmpavu(5,*)
		tmpavu(2,*)=tmpavu(2,*)*tmpavu(5,*)
		tmpdiffu=avgbin(thetas2(0,*),thetas2(2,*)^2,binsize=0.01,/quiet,/more)
		bladiffn = (ab(0,*)-shift(ab(0,*),1))^2+(ab(1,*)-shift(ab(1,*),1))^2
		bladiffn = sqrt(bladiffn(1:*))
		tmpabudiff = foldabdiff(ab(0:1,1:*),bladiffn, $
			themin=-mm2,themax=mm2,xbin=0.002,ybin=0.002)

		if (flag eq 0b) then begin
			flag = 1b
			theabhist = tmphist4
			theabdiff = tmpabdiff
			theabudiff = tmpabudiff
			theabhistu = tmphist5
			the3dhist = tmp3d
			theav = tmpav & theav(3,*) = 9e9
			nrowhist = n_elements(theav(0,*))
			thetahist = fltarr(10,nrowhist)
			thetahist(0,*)=theav(0,*)
			theavu = tmpavu & theavu(3,*) = 9e9
			theavq = tmpavq & theavq(3,*) = 9e9
		endif else begin
			theabhist = theabhist + tmphist4
			theabdiff = theabdiff + tmpabdiff
			theabudiff = theabudiff + tmpabudiff
			theav([1,2,5],*) = theav([1,2,5],*) + tmpav([1,2,5],*)
			w=where(tmpav(5,*) gt 0); FIX 6-16-14 (evolvesoft07)
			theav(3,w) = theav(3,w) < tmpav(3,w)
			theav(4,*) = theav(4,*) > tmpav(4,*)
			theabhistu = theabhistu + tmphist5
			theavu([1,2,5],*) = theavu([1,2,5],*) + tmpavu([1,2,5],*)
			w=where(tmpavu(5,*) gt 0)
			theavu(3,w) = theavu(3,w) < tmpavu(3,w)
			theavu(4,*) = theavu(4,*) > tmpavu(4,*)
			theavq([1,2,5],*) = theavq([1,2,5],*) + tmpavq([1,2,5],*)
			w=where(tmpavq(5,*) gt 0)
			theavq(3,w) = theavq(3,w) < tmpavq(3,w)
			theavq(4,*) = theavq(4,*) > tmpavq(4,*)
			the3dhist = the3dhist + tmp3d
		endelse
		thetahist(1,*)=thetahist(1,*)+theav(5,*)
		thetahist(3,*)=thetahist(3,*)+theavu(5,*)
		thetahist(5,*)=thetahist(5,*)+theavq(5,*)
		thetahist(2,*)=thetahist(2,*)+tmpdiff(1,*)*tmpdiff(5,*)
		thetahist(4,*)=thetahist(4,*)+tmpdiffu(1,*)*tmpdiffu(5,*)
		thetahist(6,*)=thetahist(6,*)+tmpdiffq(1,*)*tmpdiffq(5,*)
		thetahist(7,*)=thetahist(7,*)+tmpdiff(6,*)*tmpdiff(5,*)
		thetahist(8,*)=thetahist(8,*)+tmpdiffu(6,*)*tmpdiffu(5,*)
		thetahist(9,*)=thetahist(9,*)+tmpdiffq(6,*)*tmpdiffq(5,*)
		if (not keyword_set(nopic)) then begin
			; **** PLOTTING STUFF ****
			pic1=(small(theabdiff(*,*,0)/(theabdiff(*,*,1)+1.0),n=4))
			pic1(0,0)=0. & pic1(0,299)=0. & pic1(299,299)=0. & pic1(299,0)=0.
			pic1=bytscl(pic1)
			foo5=-alog(theabhist)
			w=where(theabhist eq 0) & foo5(w)=0
			mmm=min(foo5) & foo5(w) = mmm
			pic2=bytscl(small(foo5,n=4))
			nel=n_elements(pic1(0,*))
			tmp=(theabudiff(*,*,0)/(theabudiff(*,*,1)+1.0))
			;pic4=bytscl(small(tmp(50:1549,50:1549),n=5))
			pic4=bytscl(small(tmp,n=4))
			foo6=-alog(theabhistu(50:1549,50:1549)+1.0)
			w=where(theabhistu(50:1549,50:1549) eq 0) & foo6(w)=0
			mmm=min(foo6) & foo6(w) = mmm
			;pic6=bytscl(small(foo6,n=5))
			tmp=fltarr(1600,1600) & tmp(50:1549,50:1549)=foo6
			pic6=bytscl(small(tmp,n=4))
			pic3=bytarr(nel*2,nel*2)
			pic3(0:nel-1,0:nel-1)=pic2
			pic3(nel:*,0:nel-1)=pic6
			pic3(0:nel-1,nel:*)=pic1
			pic3(nel:*,nel:*)=pic4
			tvscl,pic3
			plot,total(foo5(*,795:805),2)/11.0-min(foo5), $
				pos=[0.05,0.05,0.3,0.25],/xs,/noerase
		endif
		if (not keyword_set(noplot)) then begin
			if (plotflag eq 0b) then begin
				dsh=findenergybarrier(thetahist,/doplot)
			endif else begin
				dsh=findenergybarrier(thetahist)
			endelse
			if (plotflag eq 2b) then begin
				foo5=-alog(theabhistu)
				w=where(theabhistu eq 0) & foo5(w)=0
				m=min(foo5) & foo5(w) = m
				tvscl,small(foo5)
			endif
			if (plotflag eq 1b) then begin
				foo5=-alog(theabhist)
				w=where(theabhist eq 0) & foo5(w)=0
				m=min(foo5) & foo5(w) = m
				tvscl,small(foo5)
			endif
			if (plotflag eq 0b) then begin
				; diffusion
				w=where(thetahist(1,*) gt 0)
				plot,thetahist(0,w),thetahist(7,w)/thetahist(1,w), $
					xr=[-1.5,1.5],pos=[0.1,0.1,0.30,0.5],/noerase
				if (n_elements(tau) gt 1) then begin
					bla=abs(thetahist(0,*) mod (3.14159*2./3.))
					w=where(bla lt 0.1 and thetahist(1,*) gt 0)
					dmin1=mean(thetahist(7,w)/thetahist(1,w))
					w=where(bla lt 0.1 and thetahist(3,*) gt 0)
					dmin2=mean(thetahist(8,w)/thetahist(3,w))
					w=where(bla lt 0.1 and thetahist(5,*) gt 0)
					dmin3=mean(thetahist(9,w)/thetahist(5,w))
					fee=3.14159/3.0
					w=where(bla gt (fee-0.1) and bla lt (fee+0.1))
					foo=thetahist(*,w)
					w=where(foo(1,*) gt 0)
					dmax1=mean(foo(7,w)/foo(1,w))
					w=where(foo(3,*) gt 0)
					dmax2=mean(foo(8,w)/foo(3,w))
					w=where(foo(5,*) gt 0)
					dmax3=mean(foo(9,w)/foo(5,w))
					hor,[dmin1,dmax1]
					;print,'Dratios: ',dmax1/dmin1,dmax2/dmin2,dmax3/dmin3
				endif
				w=where(thetahist(3,*) gt 0)
				plot,thetahist(0,w),thetahist(8,w)/thetahist(3,w), $
					xr=[-1.5,1.5],pos=[0.40,0.1,0.65,0.5],/noerase
				if (n_elements(tau) gt 1) then hor,[dmin2,dmax2]
				w=where(thetahist(5,*) gt 0)
				plot,thetahist(0,w),thetahist(9,w)/thetahist(5,w), $
					xr=[-1.5,1.5],pos=[0.75,0.1,1.00,0.5],/noerase
				if (n_elements(tau) gt 1) then hor,[dmin3,dmax3]
			endif
			plotflag = plotflag + 1b
			if (plotflag ge 2b) then plotflag = 0b
		endif else begin
			dsh=findenergybarrier(thetahist)
		endelse
		emaxmin=getemaxmin(theav)
		emaxminu=getemaxmin(theavu)
		deh = (emaxmin(0)-emaxmin(1))
		dehu = (emaxminu(0)-emaxminu(1))

		; find free energy barriers:  ab
		foo5=-alog(theabhist)
		w=where(theabhist eq 0) & foo5(w)=0
		mmm=min(foo5) & foo5(w) = mmm
		nnn=n_elements(foo5(0,*))
		ddd=shift(dist(nnn),nnn/2,nnn/2)
		w=where(ddd lt 50)
		ctr=mean(foo5(w))-mmm
		xscan=total(foo5(0:nnn/2,nnn/2-3:nnn/2+3),2)/7
		tmp=max(xscan,mi)
		abmax=min(xscan(mi:*))-mmm

		; find free energy barriers:  abu
		foo5=-alog(theabhistu)
		w=where(theabhistu eq 0) & foo5(w)=0
		mmm=min(foo5) & foo5(w) = mmm
		nnn=n_elements(foo5(0,*))
		ddd=shift(dist(nnn),nnn/2,nnn/2)
		w=where(ddd lt 100)
		uctr=mean(foo5(w))-mmm
		xscan=total(foo5(0:nnn/2,nnn/2-3:nnn/2+3),2)/7
		tmp=max(xscan,mi)
		abumax=min(xscan(mi:*))-mmm
		print,'barriers: ab abu ctr uctr',abmax,abumax,float(ctr),float(uctr)
	endelse

	time0 = time1
endwhile

result=result(*,0L:count-1L)

w=where(theav(5,*) gt 0)
theav(1,w)=theav(1,w)/theav(5,w)
theav(2,w)=theav(2,w)/theav(5,w)
w=where(theavu(5,*) gt 0)
theavu(1,w)=theavu(1,w)/theavu(5,w)
theavu(2,w)=theavu(2,w)/theavu(5,w)
w=where(theavq(5,*) gt 0)
theavq(1,w)=theavq(1,w)/theavq(5,w)
theavq(2,w)=theavq(2,w)/theavq(5,w)

tmp0=reform(the3dhist(*,*,0))
tmp1=reform(the3dhist(*,*,1))
w=where(tmp1 gt 0)
tmp0(w)=tmp0(w)/tmp1(w)
the3dhist(*,*,0)=tmp0
the3dhist(*,*,1)=tmp1

w=where(thetahist(1,*) gt 0)
thetahist(2,w)=thetahist(2,w)/thetahist(1,w)
w=where(thetahist(3,*) gt 0)
thetahist(4,w)=thetahist(4,w)/thetahist(3,w)


if (keyword_set(filename)) then begin
	if (keyword_set(trajsave)) then write_gdf,result,'traj.'+filename

	write_gdf,theabhist,'ab.'+filename
	write_gdf,theabhistu,'abu.'+filename
	write_gdf,the3dhist,'uv.'+filename
	foo5=-alog(theabhist)
	w=where(theabhist eq 0) & foo5(w)=0
	m=min(foo5) & foo5(w) = m
	write_gdf,foo5,'fab.'+filename
	write_tiff,'fab.'+filename+'.tif',bytscl(foo5)
	foo5=-alog(theabhistu)
	w=where(theabhistu eq 0) & foo5(w)=0
	m=min(foo5) & foo5(w) = m
	write_gdf,foo5,'fabu.'+filename
	write_tiff,'fabu.'+filename+'.tif',bytscl(foo5)

	ab=computethetar(result)
	tau=fourtau5(ab(0,*),ab(1,*),drad=8.0,dtheta=0.349); 20 deg
	life1=mean(tau)
	print,'lifetime = ',life1,n_elements(tau),' eps = ',epsilon
	write_gdf,tau,'tau.'+filename

	w=where(theav(5,*) gt 0,nw) & w=w(1:nw-2)
	write_gdf,theav(*,w),'energy.'+filename
	w=where(theavu(5,*) gt 0,nw) & w=w(1:nw-2)
	write_gdf,theavu(*,w),'energyu.'+filename
	w=where(theavq(5,*) gt 0,nw) & w=w(1:nw-2)
	write_gdf,theavq(*,w),'energyq.'+filename
	pic1=bytscl(theabdiff(*,*,0)/(theabdiff(*,*,1)+1.0))
	write_tiff,'abdiff.'+filename+'.tif',pic1
	write_gdf,theabdiff,'abdiff.'+filename
	pic4=bytscl(theabudiff(*,*,0)/(theabudiff(*,*,1)+1.0))
	write_tiff,'abudiff.'+filename+'.tif',pic4
	write_gdf,theabudiff,'abudiff.'+filename

	write_gdf,allcircs,'allcircs.'+filename

	thetahist(1,*)=theav(5,*)
	w=where(thetahist(1,*) gt 0)
	thetahist(1,w)=-alog(thetahist(1,w))
	w=where(thetahist(3,*) gt 0)
	thetahist(3,w)=-alog(thetahist(3,w))
	w=where(thetahist(5,*) gt 0)
	thetahist(5,w)=-alog(thetahist(5,w))
	write_gdf,thetahist,'theta.'+filename

	summary=[epsilon,dt, $ ; 0-2
	   life1,dsh[0],deh,dsh[1],dehu, $; 3-7
       version];  8
	write_gdf,summary,'sum.'+filename
	;print,'lifetime1 = ',mean(life1),n_elements(life1)
	;print,'lifetime2 = ',mean(life2),n_elements(life2)
endif

!P.MULTI=0

return,result
end
