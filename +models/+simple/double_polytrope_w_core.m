function cmp = double_polytrope_w_core(N, x)
%DOUBLE_POLYTROPE_W_CORE Model planet approximated by two polytrope and core.
%    DOUBLE_POLYTROPE_W_CORE(N, x) returns an N-layer CMSPlanet object
%    constructed with two barotropes.Polytrope eos objects and a constant
%    density core. First polytrope defined by constant x(1) and index x(2) is
%    applied to layers 1:tind and second polytrope defined by constant x(3) and
%    index x(4) is applied to layers tind+1:N-1, where tind=fix(N*(1-x(5))).
%    Finally, layer N is assigned barotropes.ConstDensity eos with rho0=x(6).
%    Layer radii are then equally spaced betwen 1 and x(7).

cmp = CMSPlanet(N);
eos1 = barotropes.Polytrope(x(1), x(2));
eos2 = barotropes.Polytrope(x(3), x(4));
eos3 = barotropes.ConstDensity(x(6));
tind = fix(N*(1 - x(5)));
cmp.eos = [repmat(eos1, tind, 1); repmat(eos2, N - tind - 1, 1); eos3];
cmp.cms.lambdas = linspace(1, x(7), N);

end
