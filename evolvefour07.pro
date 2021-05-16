; evolvefour01 -- started 7-6-14 by Eric Weeks
;
; version 02:  try to save some data as we go
; version 03:  add in unit vector code, diffusivity
; version 04:  add in 3-D phase space, upgrade lifetime calculation
; version 05:  add in period, from fourtau5
; version 06:  compute & save diffusivity as a function of (a,b)
; version 07:  add diffusivity(au,bu) and fix uv file bug
;
; taken from evolvesoft07.pro
;
; note that softsteps08, 09 is my idea for forward flux sampling
;
; evolvesoft07 -- started 8-14-13 by Eric Weeks
; try to make some subroutines
; add in ability to be 'hard' for T < 1e-15
;
; evolvesoft06 -- started 8-13-13 by Eric Weeks
; started from v04, add in some stuff from v06
; modify the step size to be less random
; evolvesoft01 -- started 4-5-12 by Eric Weeks
; taken from evolvecirc01 which was
; based heavily on Gary Hunter's matlab code (evolvebrownianthreecircmodxy.m)

function computeallenergy,circs,wallpos,sqrtwallpos,ndisks,energyptr

nel=ndisks*(ndisks+1)
nel=nel/2
; this gives n*(n-1)/2 disk-disk contacts and n disk-wall contacts
energy=fltarr(nel); by default this is zero

for j=0,ndisks-1 do begin
	radw=circs(0,j)*circs(0,j)+circs(1,j)*circs(1,j)
	if (radw gt wallpos) then begin
		dr=sqrt(radw)-sqrtwallpos
		energy(j)=dr^6
	endif
endfor
for e1 = 0,ndisks-2 do begin
	for e2 = e1+1,ndisks-1 do begin
		dx=circs(0,e1)-circs(0,e2)
		dy=circs(1,e1)-circs(1,e2)
		dr=dx*dx+dy*dy
		if (dr lt 4.0) then begin
			dr=2.0-sqrt(dr) & energy(energyptr(e1,e2))=dr^6
		endif
	endfor
endfor

print,'initial energy =',energy

return,energy
end

function browniansoftstep,circs,wallpos,sqrtwallpos,drs,s,seed,ptrs, $
		energy,temp=temp,ndisks,energyptr
; tries to update 'circs' with three Brownian steps

thernd=randomu(seed,ndisks)

for j=0,ndisks-1 do begin
	newenergy=energy
	testpos=circs(*,s(j))+drs(*,j)
	; checkwalloverlap
	radw=testpos(0)*testpos(0)+testpos(1)*testpos(1)
	if (radw gt wallpos) then begin
		dr=sqrt(radw)-sqrtwallpos
		newenergy(s(j))=dr^6
	endif
	; checkparticleoverlap -- with the other particles
	for k=0,ndisks-2 do begin
		dxy2 = circs(*,ptrs(k,s(j)))-testpos
		dxy2=dxy2(0)*dxy2(0)+dxy2(1)*dxy2(1)
		if (dxy2 lt 4.0) then begin
			dr=2.0-sqrt(dxy2);  this is the overlap
			newenergy(energyptr(s(j),ptrs(k,s(j))))=dr^6
		endif
	endfor

	denergy = total(newenergy-energy)
	if (denergy le 0.0) then begin
		circs(*,s(j))=testpos
		energy=newenergy
	endif else begin
		if (temp gt 1e-18) then begin
			boltz = exp(-denergy/temp)
			if (thernd(j) lt boltz) then begin
				circs(*,s(j))=testpos
				energy=newenergy
			endif
		endif 
		; if temp lt 1e-15, then don't move because
		; it would cause an overlap
	endelse
endfor


return,circs
end

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


function evolvefour07,epsilon,circs,steps,diffn=diffn, $
    savesteps=savesteps,temp=temp,noinit=noinit,seed=seed, $
    filename=filename,moresteps=moresteps,trajsave=trajsave, $
    noplot=noplot,nopic=nopic

version=7.1
ndisks=4
if (not keyword_set(diffn)) then diffn=1e-4
if (not keyword_set(savesteps)) then savesteps=10
if (not keyword_set(temp)) then temp=1.0
if (not keyword_set(noinit)) then begin
	;below code for three disks:
	;inv3=1.0/sqrt(3.0)
	;circs=[[0,2*inv3],[1,-inv3],[-1,-inv3]]
	;
	;below code for four disks:
	inv4=(1.0+epsilon)*sqrt(2.0)
	circs=[[0,inv4],[0,-inv4],[inv4,0],[-inv4,0]]
endif
if (not keyword_set(moresteps)) then moresteps=1L
ttlsteps=(1.0*moresteps)*(1.0*steps)

largerad=3.0+epsilon
sqrtwallpos = largerad-1.0
wallpos = (sqrtwallpos)^2
result=fltarr(ndisks*2+1,(steps/savesteps)*moresteps+1)
sdiffn=sqrt(diffn)
; ptrs: these point to the other disks
; 3 disks:  ptrs=[[1,2],[2,0],[0,1]]
ss=indgen(ndisks*ndisks) mod ndisks
ss(indgen(ndisks)*(ndisks+1))=-1;   remove self-pointing array elements
w=where(ss ge 0)
ptrs=reform(ss(w),ndisks-1,ndisks)

; ss:  these are the different possible orders for trial moves
; below line for 3 disks exactly:
; ss=[[0,1,2],[0,2,1],[1,0,2],[1,2,0],[2,0,1],[2,1,0],[2,1,0]]
ss=ecombinatoric(indgen(ndisks))
nss=n_elements(ss(0,*))

count=ndisks
energyptr = intarr(ndisks,ndisks)
for i=0,ndisks-2 do begin
	for j=i+1,ndisks-1 do begin
		energyptr(i,j)=count
		energyptr(j,i)=count
		count=count+1
	endfor
endfor

energy=computeallenergy(circs,wallpos,sqrtwallpos,ndisks,energyptr)

count=0L
icount=0L
maxicount=999999L
forhisto=fltarr(ndisks*2+1,maxicount+1L)
flag=0b
plotflag=0b
twopi=2.0*3.14159265
mm=40.
mm2=1.6


for k=0L,moresteps-1L do begin 
for i=0L,steps-1L do begin
	if ((i mod 200000) eq 0) then begin
		if ((i mod 1000000) eq 0) then begin
			message,string((i*1.0+(k*1.0)*(steps*1.0))/ttlsteps)+ $
			" eps="+string(epsilon)+ $
			"  T="+string(temp),/inf
		endif
		scount=0L
		; below lines generate some random steps, faster to generate
		; a whole bunch all at once.
		rndtheta = randomu(seed,ndisks,200000)*2.0*3.14159265358
		rndsteps = fltarr(2,ndisks,200000)
		rndsteps(0,*,*)=sdiffn*cos(rndtheta)
		rndsteps(1,*,*)=sdiffn*sin(rndtheta)
		rndtheta=0; save on memory
	endif
	rands=randomu(seed)*nss
	s=ss(*,rands)
	nowsteps=rndsteps(*,*,scount)
	circs = browniansoftstep(circs,wallpos,sqrtwallpos,nowsteps,s,seed, $
		ptrs,energy,temp=temp,ndisks,energyptr)
	scount = scount+1L
	etotal=total(energy)
	if ((i mod savesteps) eq 0) then begin
		result(0:ndisks*2-1,count)=circs
		result(ndisks*2,count)=etotal
		count++
	endif

	; next long bit is for histograms, accumulation of data, etc
	forhisto(0:ndisks*2-1,icount)=circs
	forhisto(ndisks*2,icount)=etotal
	if (icount ge maxicount) then begin
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
			pic4=bytscl(small(tmp(50:1549,50:1549),n=5))
			foo6=-alog(theabhistu(50:1549,50:1549)+1.0)
			w=where(theabhistu(50:1549,50:1549) eq 0) & foo6(w)=0
			mmm=min(foo6) & foo6(w) = mmm
			pic6=bytscl(small(foo6,n=5))
			pic3=bytarr(nel*2,nel*2)
			pic3(0:nel-1,0:nel-1)=pic2
			pic3(nel:*,0:nel-1)=pic1
			pic3(0:nel-1,nel:*)=pic6
			pic3(nel:*,nel:*)=pic4
			tvscl,pic3
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
			if (plotflag eq 3b) then begin
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
					print,'Dratios: ',dmax1/dmin1,dmax2/dmin2,dmax3/dmin3
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
			;plotflag = plotflag + 1b
			;if (plotflag ge 2b) then plotflag = 0b
		endif else begin
			dsh=findenergybarrier(thetahist)
		endelse
		emaxmin=getemaxmin(theav)
		emaxminu=getemaxmin(theavu)
		deh = (emaxmin(0)-emaxmin(1))/temp
		dehu = (emaxminu(0)-emaxminu(1))/temp

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

		; ----- check lifetime
		ab=computeaabb(result(*,0L:count-1L))
		tau=fourtau5(ab(16,*),ab(17,*),drad=8.0,dtheta=0.349,period=per);20 deg
		print,'lifetime = ',mean(tau)*savesteps,n_elements(tau),per*savesteps
		icount = -1L
	endif
	icount = icount + 1L

endfor; i loop
endfor; k loop

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

	ab=computeaabb(result(*,0L:count-1L))
	tau=fourtau5(ab(16,*),ab(17,*),drad=8.0,dtheta=0.349); 20 deg
	life1=mean(tau)*savesteps
	print,'lifetime = ',life1,n_elements(tau)
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

	thetahist(1,*)=theav(5,*)
	w=where(thetahist(1,*) gt 0)
	thetahist(1,w)=-alog(thetahist(1,w))
	w=where(thetahist(3,*) gt 0)
	thetahist(3,w)=-alog(thetahist(3,w))
	w=where(thetahist(5,*) gt 0)
	thetahist(5,w)=-alog(thetahist(5,w))
	write_gdf,thetahist,'theta.'+filename

	summary=[epsilon,temp,ttlsteps,savesteps,diffn, $ ; 0-4
	   life1,dsh[0],deh,dsh[1],dehu, $; 5-9
       version];  10
	write_gdf,summary,'sum.'+filename
	;print,'lifetime1 = ',mean(life1),n_elements(life1)
	;print,'lifetime2 = ',mean(life2),n_elements(life2)
endif

!P.MULTI=0

return,result
end
