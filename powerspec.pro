pro powerspec,box1,box2,dx,nrad,k,pk
;By Saleem Zaroubi

a=size(box1,/dimension)
nbox=float(a[0])
kx=shift(dindgen(nbox)-nbox/2,nbox/2)*2.*!pi/float(nbox*dx)
ky=kx 
kz=kx
absk=fltarr(nbox,nbox,nbox)
for i=0,nbox-1 do for j=0,nbox-1 do absk(i,j,*)=sqrt(kx(i)^2+ky(j)^2+kz(*)^2)
nrad=nrad
kradmin=kx(0)
kradmax= max(kx)
;dlogk=alog10(kradmax/kradmin)/float(nrad+1)
;krad=kradmin*10.^(findgen(nrad+1)*dlogk)
krad=findgen(nrad+1)/float(nrad)*max(absk)
pk=fltarr(nrad)

Fbox1=fft(box1,/double)*sqrt((float(nbox))^3)
Fbox2=fft(box2,/double)*sqrt((float(nbox))^3)

print,total(box1*box2),total(abs(Fbox1*CONJ(Fbox2)))
for i=0,nrad-1 do begin
   nn=where(absk gt krad(i) and absk le krad(i+1))
   if total(nn) ge 0 then pk(i) = mean(fbox1(nn)*CONJ(fbox2(nn)))
endfor
k=((krad(0:nrad-1)^(3.)+krad(1:nrad)^(3.))/2.)^(1./3.)
return
end
