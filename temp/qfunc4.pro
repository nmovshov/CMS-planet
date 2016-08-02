function Qfunc4, xi
common Jval,j
common qval,q
common muval, mu

jsum1 = 0.d0
jsum2 = 0.d0

for k=0,14 do begin
n=2+2*k
jsum1 = jsum1 + j[k]*legendre(mu,n,/double)/xi^n
jsum2 = jsum2 + j[k] *legendre(0.d0,n,/double)
endfor

return, (1.d0/xi) * (1.d0-jsum1) + (q/3.d0)*xi*xi* $
      (1.d0 - legendre(mu,2,/double)) -q/2.d0 - 1.d0 + jsum2

end
