function cmp = double_polytrope(N, x, lamstrat)
%DOUBLE_POLYTROPE Model planet approximated by two polytropes.
%    DOUBLE_POLYTROPE(N, x) returns an N-layer CMSPlanet object with two
%    barotropes.Polytrope eos objects. First polytrope defined by constant x(1)
%    and index x(2) is applied to layers 2:tind. Second polytrope defined by
%    constant x(3) and index x(4) is applied to layers tind+1:N. Transition is at
%    layer tind+1, the first layer with (ai/a0)<=x(5). Layer 1 is assigned a zero
%    density barotropes.ConstDensity eos and is approximately half the width of
%    the layers below. The layer spacing is designed to minimize discretization
%    error by concentrating 2/3 of the available layers in the top 0.5 of the
%    planet.
%
%    DOUBLE_POLYTROPE(N, x, lamstrat) uses the 2-element vector lamstrat to
%    specify the layer spacing strategy. Approximately lamstrat(1) of available
%    layers will be distributed in the top lamstrat(2) of the planet. For example,
%    passing [3/4, 0.2] concentrates the layers heavily in the top 20% of the
%    planet, leaving about N/4 layers to fill the bottom 80%. Passing [r, r] gives
%    approximately equal spacing throughout the planet. (Approximately because a
%    single half-width layer of zero density is always reserved for the surface.)

narginchk(2,3)
if nargin == 2, lamstrat = [2/3, 1/2]; end
validateattributes(lamstrat, {'numeric'}, {'vector', 'numel', 2, '>', 0, '<', 1})
assert(x(5)>0 && x(5)<1, 'The transition must be 0<R<1.')

cmp = CMSPlanet(N);

n1 = fix(lamstrat(1)*(N - 1));
n2 = N - n1 - 1;
dl1 = lamstrat(2)/(n1 - 1);
dl2 = (1 - lamstrat(2))/(n2 + 1);
lam1 = linspace(1 - dl1/2, (1 - lamstrat(2)), n1);
lam2 = linspace((1 - lamstrat(2)) - dl2, dl2, n2);
mylambdas = [1, lam1, lam2]';

% Replace lambda_tind to match the transition
[~, tind] = min(abs(mylambdas-x(5)));
assert(tind > 2, 'The transition is too close to the surface and the first polytrope has zero layers.')
mylambdas(tind) = x(5);
cmp.cms.lambdas = mylambdas;

eos0 = barotropes.ConstDensity(0);
eos1 = barotropes.Polytrope(x(1), x(2));
eos2 = barotropes.Polytrope(x(3), x(4));

cmp.eos = [eos0; repmat(eos1, tind - 2, 1); repmat(eos2, N - tind + 1, 1)];

end
