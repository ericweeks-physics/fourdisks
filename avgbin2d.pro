;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
; NAME:
;  avgbin2d
;
; WARNING:
;  To display the output as an image use tvscl or bytscl to scale the
;  range between 0 and 255!!
;+
; PURPOSE:
;  Calculates the average value of f(x,y) versus x and y. It is the 2d
;  equivalent of avgbin.
;-
; INPUT:
;  DATA1: An array with three columns and n rows. 
;    data1[0,*] = x-values
;    data1[1,*] = y-values
;    data1[2,*] = f(x,y)-values
;  DATA2: An array with f(x,y)-values. If this is defined then data1 has 
;    only the x and y-values.
;  DATA3: An array with f(x,y)-values. If this is defined then data1=x-values
;    and data2=y-values.
;
; PARAMS:
;  XBINSIZE: Size of bins in the x direction. 
;  YBINSIZE: Size of bins in the y direction. 
;  Probably you should have XBINSIZE GE dx and YBINSIZE GE dy (where dx
;  and dy are the smallest x and y intervals in your data).
;
; OUTPUT:
;  Returns an MxN matrix where M = (xmax-xmin)/xbinsize and 
;  N = (ymax-min)/ybinsize. Each pixel represents a range in x and y
;  described by xbinsize and ybinsize. The value of result[i,j] is the
;  average of f(x,y) over the ranges in x and y spanned by that pixel.
;
;  To display the output as an image use tvscl or bytscl to scale the
;  range between 0 and 255!!
;
; EXAMPLE:
;  i=avgbin2d(f, xbinsixe=0.01, ybinsize=0.5)
;
; TO DO:
;  Figure out exactly how this works and especially how hist_nd
;  works and then improve documentation!!
;
; HISTORY:
;  10 Nov 2004 - Gianguido Cianci
;   8 Jun 2005 - Eric Weeks (allow for 1-3 input arrays)
;
; OTHER:
;  Really I did not write this ingenious little piece of idl magic. A
;  guy called JD Smith suggested this on google's group
;  comp.lan.idl-pvwave. It is much faster than gavgbin2d.pro, the
;  naive implementation. From what I understand presently, this could
;  possibly be extended to more dimensions...
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

FUNCTION avgbin2d, data1, data2, data3, xbinsize=xbin, ybinsize=ybin, $
    xmin=xmin,ymin=ymin,xmax=xmax,ymax=ymax

Compile_Opt IDL2

if (keyword_set(data3)) then begin
    n1 = n_elements(data1)
    n2 = n_elements(data2)
    n3 = n_elements(data3)
    if ((n1 ne n2) or (n2 ne n3)) then begin
        message,'Hey!  Your three data sets are not the same size.',/inf
        message,'I am going to quit.'
    endif
    dataxyz = fltarr(3,n1)
    dataxyz[0,*] = data1
    dataxyz[1,*] = data2
    dataxyz[2,*] = data3
endif else begin
    if (keyword_set(data2)) then begin
        ; assume data1 = [x,y]  and data2 = [z]
        n1 = n_elements(data1[0,*])
        n2 = n_elements(data2)
        if (n1 ne n2) then begin
            message,'Hey!  Your two data sets are not the same size.',/inf
            message,'I am going to quit.'
        endif
        dataxyz = fltarr(3,n1)
        dataxyz[0:1,*] = data1
        dataxyz[2,*] = data2
    endif else begin
        ; function like Gianguido's original code
        dataxyz = data1
    endelse
endelse


xmin2=min(dataxyz[0,*], max=xmax2)
ymin2=min(dataxyz[1,*], max=ymax2)
if (not keyword_set(xmin)) then xmin=xmin2
if (not keyword_set(ymin)) then ymin=ymin2
if (not keyword_set(xmax)) then xmax=xmax2
if (not keyword_set(ymax)) then ymax=ymax2
; print,xmin,ymin,xmax,ymax,xbin,ybin

;Determine how many bins to use (i.e. the number of pixels of the
;output picture)
IF (n_elements(xbin) GT 0) THEN BEGIN 
  nbx=ceil((xmax-xmin)/xbin)
endif else begin
   xbin=(xmax-xmin)/128.0
   nbx=ceil((xmax-xmin)/xbin)
ENDELSE

IF (n_elements(ybin) GT 0) THEN BEGIN 
   nby=ceil((ymax-ymin)/ybin)
endif else begin
   ybin=(ymax-ymin)/128.0
   nby=ceil((ymax-ymin)/ybin)
ENDELSE


im=hist_nd(dataxyz[0:1,*],NBINS=[nbx,nby],REVERSE_INDICES=ri)
F_im=make_array(/FLOAT,size(im,/DIMENSIONS))
FOR j=0,n_elements(im)-1 DO BEGIN 
   IF ri[j+1] GT ri[j] THEN BEGIN 
      F_im[j] = total(dataxyz[2,ri[ri[j]:ri[j+1]-1]])/n_elements(dataxyz[2,ri[ri[j]:ri[j+1]-1]])
   ENDIF
ENDFOR

result=fltarr(nbx,nby,2)

;this contains avg f(x,y) as a function of x and y
result[*,*,0]=F_im

;this contains number of times f exists for a given x and y
result[*,*,1]=im


return, result
END 

