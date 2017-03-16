function cmp = single_polytrope(N, x, lamstrat)
%SINGLE_POLYTROPE The simplest toy model planet.
%    SINGLE_POLYTROPE(N, x) returns an N-layer CMSPlanet object with a
%    barotropes.Polytrope with constant x(1) and index x(2) assigned to layers
%    2:N. Layer 1 is assigned a zero density barotropes.ConstDensity eos and is
%    approximately half the width of the layers below. The layer spacing is
%    designed to minimize discretization error by concentrating 2/3 of the
%    available layers in the top 0.5 of the planet.
%
%    SINGLE_POLYTROPE(N, x, lamstrat) uses the 2-element vector lamstrat to
%    specify the layer spacing strategy. Approximately lamstrat(1) of available
%    layers will be distributed in the top lamstrat(2) of the planet. For example,
%    passing [3/4, 0.2] concentrates the layers heavily in the top 20% of the
%    planet, leaving about N/4 layers to fill the bottom 80%. Passing [r, r] gives
%    approximately equal spacing throughout the planet. (Approximately because a
%    single half-width layer of zero density is always reserved for the surface.)

narginchk(2,3)
if nargin == 2, lamstrat = [2/3, 1/2]; end
validateattributes(lamstrat, {'numeric'}, {'vector', 'numel', 2, '>', 0, '<', 1})

cmp = CMSPlanet(N);

n1 = fix(lamstrat(1)*(N - 1));
n2 = N - n1 - 1;
dl1 = lamstrat(2)/(n1 - 1);
dl2 = (1 - lamstrat(2))/(n2 + 1);
lam1 = linspace(1 - dl1/2, (1 - lamstrat(2)), n1);
lam2 = linspace((1 - lamstrat(2)) - dl2, dl2, n2);
cmp.cms.lambdas = [1, lam1, lam2]';

eos0 = barotropes.ConstDensity(0);
eos1 = barotropes.Polytrope(x(1), x(2));
cmp.eos = [eos0; repmat(eos1, N-1, 1)];

end
