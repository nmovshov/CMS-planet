function lfunc, l
common mval,m
  return, (3.d0/2.d0/l^3)* $
      ( (3.d0+l^2)*atan(l) - 3.d0*l) - m
end
