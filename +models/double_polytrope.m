function cmp = double_polytrope(N, x)
%DOUBLE_POLYTROPE Model planet approximated by two polytropes.
%    DOUBLE_POLYTROPE(N, x) returns an N-layer CMSPlanet object with two
%    barotropes.Polytrope eos objects. First polytrope defined by constant x(1)
%    and index x(2) is applied to layers 2:tind. Second polytrope defined by
%    constant x(3) and index x(4) is applied to layers tind+1:N. Transition is at
%    layer tind+1, the first layer with (ai/a0)<=x(5). Layer 1 is assigned a zero
%    density barotropes.ConstDensity eos and is approximately half the width of
%    the layers below. The layer spacing is designed to minimize discretization
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
assert(x(5) >=0 && x(5)<=1)
tind = find(cmp.cms.lambdas <= x(5), 1) - 1;
if isempty(tind), tind = N; end
tind = min(max(tind, 1), N);
cmp.eos = [eos0; repmat(eos1, tind - 1, 1); repmat(eos2, N - tind, 1)];

end
