function cmp = single_polytrope(N, x)
%SINGLE_POLYTROPE The simplest toy model planet.
%    SINGLE_POLYTROPE(N, x) returns an N-layer CMSPlanet object with a
%    barotropes.Polytrope with constant x(1) and index x(2) assigned to layers
%    2:N. Layer 1 is assigned a zero density barotropes.ConstDensity eos and is
%    approximately half the width of the layers below. The layer spacing is
%    designed to minimize discretization error by concentrating 2/3 of the
%    available layers in the top 0.5 of the planet (see Hubbard & Militzer,
%    2016).

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
cmp.eos = [eos0; repmat(eos1, N-1, 1)];

end
