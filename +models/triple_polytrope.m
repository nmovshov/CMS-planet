function cmp = triple_polytrope(N, x, lamstrat, forcematch)
%TRIPLE_POLYTROPE Model planet approximated by three polytropes.
%    TRIPLE_POLYTROPE(N, x) returns an N-layer CMSPlanet object with three
%    barotropes.Polytrope eos objects. First polytrope defined by constant x(1)
%    and index x(2) is applied to layers 2:tind-1. Second polytrope defined by
%    constant x(3) and index x(4) is applied to layers tind:cind-1. Third
%    polytrope defined by constant x(5) and index x(6) is applied to layers
%    cind:N. Transition from first to second polytrope is at layer tind, the layer
%    with (ai/a0) nearest x(7). Transition from second to third polytrope is at
%    layer cind, the layer with (ai/a0) nearest x(8). Layer 1 is assigned a zero
%    density barotropes.ConstDensity eos and is approximately half the width of
%    the layers below. The default layer spacing concentrates 2/3 of the available
%    layers in the top 0.5 of the planet.
%
%    TRIPLE_POLYTROPE(N, x, lamstrat) uses the 2-element vector lamstrat to
%    specify the layer spacing strategy. Approximately lamstrat(1) of available
%    layers will be distributed in the top lamstrat(2) of the planet. For example,
%    passing [3/4, 0.2] concentrates the layers heavily in the top 20% of the
%    planet, leaving about N/4 layers to fill the bottom 80%. Passing [r, r] gives
%    approximately equal spacing throughout the planet. (Approximately because a
%    single half-width layer of zero density is always reserved for the surface.)
%    Use [] (empty array) for default values.
%
%    TRIPLE_POLYTROPE(..., forcematch) if forcematch=true forces the normalized
%    radii of the transition layers to exactly match x(7) and x(8). This is
%    applied after the lambda spacing strategy.

narginchk(2,4)
if ((nargin == 2) || isempty(lamstrat)), lamstrat = [2/3, 1/2]; end
if ((nargin < 4) || isempty(forcematch)), forcematch = false; end
validateattributes(lamstrat, {'numeric'}, {'vector', 'numel', 2, '>', 0, '<', 1})
validateattributes(forcematch, {'logical'}, {'scalar'})
assert(x(7) > 0 && x(7) < 1, 'First transition radius must be in (0,1).')
assert(x(8) > 0 && x(8) < 1, 'Second transition radius must be in (0,1).')
assert(x(8) < x(7), 'Second transition must come before first transition.')

cmp = CMSPlanet(N);

n1 = fix(lamstrat(1)*(N - 1));
n2 = N - n1 - 1;
dl1 = lamstrat(2)/(n1 - 1);
dl2 = (1 - lamstrat(2))/(n2 + 1);
lam1 = linspace(1 - dl1/2, (1 - lamstrat(2)), n1);
lam2 = linspace((1 - lamstrat(2)) - dl2, dl2, n2);
lambdas = [1, lam1, lam2]';

[~, tind] = min(abs(lambdas-x(7)));
assert(tind > 2,...
    'First transition too close to surface; first polytrope has zero layers.')
[~, cind] = min(abs(lambdas-x(8)));
assert(cind > tind,...
    'Transitions are too close together; second polytrope has zero layers.')

if forcematch, lambdas(tind) = x(7); end
if forcematch, lambdas(cind) = x(8); end
cmp.cms.lambdas = lambdas;

eos0 = barotropes.ConstDensity(0);
eos1 = barotropes.Polytrope(x(1), x(2));
eos2 = barotropes.Polytrope(x(3), x(4));
eos3 = barotropes.Polytrope(x(5), x(6));

cmp.eos = [eos0;...
    repmat(eos1, tind - 2, 1);...
    repmat(eos2, cind - tind, 1);...
    repmat(eos3, N - cind + 1, 1)];

end
