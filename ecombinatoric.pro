; ecombinatoric -- 7/6/14 Eric R. Weeks
;
; sadly, this assumes intarr and that's your only option at the moment.

function ecombinatoric,list

nl=n_elements(list)

if (nl eq 2) then begin
	result=[[list[0],list[1]],[list[1],list[0]]]
endif else begin
	temp=ecombinatoric(list(1:*))
	n=n_elements(temp(0,*))
	m=n_elements(temp(*,0))
	result=intarr(nl,n*(m+1))
	count=0
	for i=0,n-1 do begin
		for j=0,m do begin
			if (j eq 0) then begin
				result(*,count)=[list[0],temp[*,i]]
			endif else if (j eq m) then begin
				result(*,count)=[temp[*,i],list[0]]
			endif else begin
				result(*,count)=[temp[0:j-1,i],list[0],temp[j:*,i]]
			endelse
			count = count + 1
		endfor
	endfor
endelse

return,result
end
