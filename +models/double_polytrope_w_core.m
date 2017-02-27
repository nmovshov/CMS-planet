function cmp = double_polytrope_w_core(N, x)
%DOUBLE_POLYTROPE_W_CORE Two polytrope and core model planet.
%    DOUBLE_POLYTROPE_W_CORE(N, x) returns an N-layer CMSPlanet object with two
%    barotropes.Polytrope eos objects. First polytrope defined by constant x(1)
%    and index x(2) is applied to layers 2:tind. Second polytrope defined by
%    constant x(3) and index x(4) is applied to layers tind+1:N-1. Transition is
%    at layer tind+1, the first layer with (ai/a0)<=x(5). Layer N is assigned a
%    barotropes.ConstDensity eos wirh rho0=x(6) and given a normalized equatorial
%    radius equal to x(7). Layer 1 is assigned a zero density
%    barotropes.ConstDensity eos and is approximately half the width of the layers
%    below. The layer spacing is designed to minimize discretization error by
%    concentrating 2/3 of the available layers in the top 0.5 of the planet (see
%    Hubbard & Militzer, 2016). However if x(7)>0.2 then layers are equally spaced
%    above x(7).

cmp = CMSPlanet(N);

if x(7) > 0.2
    dl = (1 - x(7))/N;
    cmp.cms.lambdas = [1, linspace(1 - dl/2, x(7), N - 1)]';
else
    if x(7) <= 0, x(7) = eps; end
    n1 = fix(2/3*N) - 1;
    n2 = N - n1 - 1;
    dl1 = 0.5/(n1 - 1);
    dl2 = (0.5 - x(7))/n2;
    lam1 = linspace(1 - dl1/2, 0.5, n1);
    lam2 = linspace(0.5 - dl2, x(7), n2);
    cmp.cms.lambdas = [1, lam1, lam2]';
end

eos0 = barotropes.ConstDensity(0);
eos1 = barotropes.Polytrope(x(1), x(2));
eos2 = barotropes.Polytrope(x(3), x(4));
eos3 = barotropes.ConstDensity(x(6));
assert(x(5) > x(7) && x(5)<=1)
tind = find(cmp.cms.lambdas <= x(5), 1) - 1;
tind = min(max(tind, 1), N - 1);
cmp.eos = [eos0; repmat(eos1, tind - 1, 1); repmat(eos2, N - tind - 1, 1); eos3];

end
