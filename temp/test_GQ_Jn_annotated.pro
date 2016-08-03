;fastmoments
common Jval,j
common qval,q
common mval,m
common muval, mu

GMJ = 126686535.d0 ; � 2 (km3/s2)
aJ = 71492.d0 ; eq. radius (km)
bJ = 66854.d0 ; pol. radius (km)
J2J = 14696.43d-6 ; � 0.21d-6
J3J = -0.64d-6 ;  � 0.90d-6
J4J = -587.14d-6 ; � 1.68d-6
J6J = 34.25d-6 ; � 5.22d-6
ProtJ = 9.d0*3600.d0 +55.d0*60.d0 + 29.7d0  ; rotation pd. (s)
omegaJ=2.d0*!dpi/ProtJ

qJ=omegaJ^2 * aJ^3 / GMJ
q = qJ

;now we have the value of q for Jupiter parameters

; estimate the corresponding value of m:
m = q
;estimate ell:
gues = 0.2d0
gues0=[gues*1.1d0, gues, gues*0.9d0]
l = fx_root(gues0, 'lfunc',tol=1.d-15)

print,'q, m, l: ',q,m,l

l2 = l*l

; iterate to get Maclaurin value of polar radius b :

for iter = 0,20 do begin
bmac = aJ/sqrt(1.d0+l2)

print, 'bJ, bmac = ',bJ, bmac

sJ3 = bmac*aJ^2

mmac =omegaJ^2 * sJ3 / GMJ
m=mmac

print, 'mmac = ', mmac
gues=0.2
gues0=[gues*1.1d0, gues, gues*0.9d0]

root=fx_root(gues0, 'lfunc',tol=1.d-15)
l=root
l2=l^2
l2mac = real_part(l2)
print,form='(f18.15)',l2mac
endfor

; after 21 iterations, we now
; have the value of m and ell that correspond to the prescribed q
;+++++++++++++++++++++++++++

j = 0.d0*dindgen(15)
;estimate Jn by inserting the analytic expression:

for k=0,14 do begin
n=2+2*k
j[k]=3.d0*(-1.d0)^(2.d0+k)/(n+1.d0)/(n+3.d0) * ( (l2)/(1.d0+l2) )^(1.d0+k)
print,'Jn ',n,j[k]
endfor

;j[0] = q / 2.d0

nq = 48
; the following canned routine fills in the gaussian quadrature
;  points xq and weights wq

gauss_leg_quadr,nq,xq,wq,xrange=[0,1]

xis = 0.9d0 + 0.d0 * dindgen(nq)
; initial evaluation of xi(mu)
for i=0,nq-1 do begin
gues = xis[i]
mu = xq[i]
gues0=[gues*1.1d0, gues, gues*0.9d0]
xis[i] = fx_root(gues0, 'Qfunc4',tol=1.d-15)
endfor

; the roots of Qfunc4 give the equipotential surface's radius
;  at each value of xq

xisprev = xis

;p3=plot(xq,xis)

denom=total(wq* xis^3,/double)
; denominator of Eq. (9)

denom_anal = 1.d0 / sqrt(1.d0+l2)
;  analytic value for denominator

print, denom / denom_anal - 1.d0


;  do 100 iterations by going between Eq. (9) and Eq. (14)
for iterj = 0,100 do begin

for nord=2,30,2 do begin
num = total(wq * xis^(nord+3) * legendre(xq,nord,/double) ,/double)
k = (nord-2)/2
j[k] = -(3.d0/(nord+3.d0)) * num / denom
endfor

;print,iterj
;print, j

;  solve for roots of Eq. (14)
for i=0,nq-1 do begin
gues = 0.9d0
mu = xq[i]
gues0=[gues*1.1d0, gues, gues*0.9d0]
xis[i] = fx_root(gues0, 'Qfunc4',tol=1.d-15)

endfor
;p3=plot(xq,xis-xisprev,overplot=1)

xisprev = xis

endfor

; solve Eq. (14) for the special value mu=1 (pole)
;  which is not a gaussian quadrature point but
;  gives us xi0 = b/a
mu = 1.d0
gues = 0.9d0
gues0=[gues*1.1d0, gues, gues*0.9d0]
xi0 = fx_root(gues0, 'Qfunc4',tol=1.d-15)

l2 = 1.d0/xi0^2 - 1.d0

;  exhibit the difference between ell^2 from CMS
;  and ell^2 analytic:
print,'l2 - l2mac = ',l2-l2mac

xiell = sqrt(1.d0/(1.d0+l2*xq^2))

p4=plot(xq,xis-xiell)

xiellmac = sqrt(1.d0/(1.d0+l2mac*xq^2))
p4=plot(xq,xis-xiellmac,color='red',overplot=1)

denom_anal = 1.d0 / sqrt(1.d0+l2)

print,'gquad, analytic: '
print,denom, denom_anal,denom_anal/denom - 1.d0,form='(2f25.18,e15.5)'

for nord=2,30,2 do begin
num = total(wq * xis^(nord+3) * legendre(xq,nord,/double), /double )
num_anal = ( (-1.d0)^(nord/2) * (l2/(1.d0+l2))^(nord/2) ) / $
 ((nord + 1.d0) * sqrt(1.d0 + l2) )
print,nord,num, num_anal,num_anal/num - 1.d0,form='(i4,2f25.18,e15.5)'
endfor

end