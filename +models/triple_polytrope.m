function cmp = triple_polytrope(N, x)
%TRIPLE_POLYTROPE Model planet approximated by three polytropes.
%    TRIPLE_POLYTROPE(N, x) returns an N-layer CMSPlanet object with three
%    barotropes.Polytrope eos objects. First polytrope defined by constant x(1)
%    and index x(2) is applied to layers 2:tind. Second polytrope defined by
%    constant x(3) and index x(4) is applied to layers tind+1:cind. Third
%    polytrope defined by constant x(5) and index x(6) is applied to layers
%    cind+1:N. Transition from first to second polytrope is at layer tind+1, the
%    first layer with (ai/a0)<=x(7). Transition from second to third polytrope is
%    at layer cind+1, the first layer with (ai/a0)<=x(8). Layer 1 is assigned a
%    zero density barotropes.ConstDensity eos and is approximately half the width
%    of the layers below. The layer spacing is designed to minimize discretization
%    error by concentrating 2/3 of the available layers in the top 0.5 of the
%    planet (see Hubbard & Militzer, 2016).

cmp = CMSPlanet(N);

n1 = fix(2/3*N) - 1;
n2 = N - n1 - 1;
dl1 = 0.5/(n1 - 1);
dl2 = 0.5/(n2 + 1);
lam1 = linspace(1 - dl1/2, 0.5, n1);
lam2 = linspace(0.5 - dl2, dl2, n2);
cmp.cms.lambdas = [1, lam1, lam2]';

eos0 = barotropes.ConstDensity(0);
eos1 = barotropes.Polytrope(x(1), x(2));
eos2 = barotropes.Polytrope(x(3), x(4));
eos3 = barotropes.Polytrope(x(5), x(6));
assert(x(7) <= 1 && x(7) >= x(8) && x(8) >= 0)
tind = find(cmp.cms.lambdas <= x(7), 1) - 1;
if isempty(tind), tind = N; end
tind = min(max(tind, 1), N);
cind = find(cmp.cms.lambdas <= x(8), 1) - 1;
if isempty(cind), cind = N; end
cind = min(max(cind, 1), N);
assert(cind > tind,...
    'Degenerate model: transition layers must be monotonic increasing.')
cmp.eos = [eos0;...
    repmat(eos1, tind - 1, 1);...
    repmat(eos2, cind - tind, 1);...
    repmat(eos3, N - cind, 1)];

end
