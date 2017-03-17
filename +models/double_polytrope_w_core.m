function cmp = double_polytrope_w_core(N, x, lamstrat)
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
%
%    DOUBLE_POLYTROPE_W_CORE(N, x, lamstrat) uses the 2-element vector lamstrat to
%    specify the layer spacing strategy. Approximately lamstrat(1) of available
%    layers will be distributed in the top lamstrat(2) of the planet. For example,
%    passing [3/4, 0.2] concentrates the layers heavily in the top 20% of the
%    planet, leaving about N/4 layers to fill the bottom 80%. Passing [r, r] gives
%    approximately equal spacing throughout the planet. (Approximately because a
%    single half-width layer of zero density is always reserved for the surface.)
%    However if x(7)>0.2 or if x(7)>(1-lamstrat(2)) then we revert to equally
%    spaced layers again. Note that when x(7) is close to (1-lamstrat(2)) you
%    might be better of using equal spacing!

narginchk(2,3)
if nargin == 2, lamstrat = [2/3, 1/2]; end
validateattributes(lamstrat, {'numeric'}, {'vector', 'numel', 2, '>', 0, '<', 1})
if x(7) <= 0, x(7) = eps; end

cmp = CMSPlanet(N);

if (x(7) > 0.2) || (x(7) >= (1 - lamstrat(2)))
    dl = (1 - x(7))/N;
    cmp.cms.lambdas = [1, linspace(1 - dl/2, x(7), N - 1)]';
else
    n1 = fix(lamstrat(1)*(N - 1));
    n2 = N - n1 - 1;
    dl1 = lamstrat(2)/(n1 - 1);
    dl2 = (1 - lamstrat(2) - x(7))/(n2);
    lam1 = linspace(1 - dl1/2, (1 - lamstrat(2)), n1);
    lam2 = linspace((1 - lamstrat(2)) - dl2, x(7), n2);
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
