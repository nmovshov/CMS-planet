function cmp = single_polytrope(N, x, lamstrat)
%SINGLE_POLYTROPE The simplest toy model planet.
%    SINGLE_POLYTROPE(N, x) returns an N-layer CMSPlanet object with a
%    barotropes.Polytrope with constant x(1) and index x(2) assigned to layers
%    2:N. Layer 1 is assigned a zero density barotropes.ConstDensity eos and is
%    approximately half the width of the layers below. The default layer spacing
%    is the one referenced in lambdas.best.m.
%
%    SINGLE_POLYTROPE(N, x, lamstrat) where lamstrat is a 2-element vector lets
%    you specify the layer spacing strategy. Approximately lamstrat(1) of
%    available layers will be distributed in the top lamstrat(2) of the planet.
%    For example, passing [3/4, 0.2] concentrates the layers heavily in the top
%    20% of the planet, leaving about N/4 layers to fill the bottom 80%. A single
%    half-width layer of zero density is still reserved for the surface. To use
%    the default spacing pass lambda=[].
%
%    SINGLE_POLYTROPE(N, x, lamstrat) where lamstrat is a function handle lets you
%    specify the lambda spacing completely. Pass a handle to a function that takes
%    a single scalar integer (number of layers) and returns a vector of that
%    length with values in the interval (0, 1], for normalized layer radii. For
%    example, to set layers with equally spaced radii use
%    lamstrat=@(n)linspace(1,1/n,n). Note that the final layer radii might be
%    slightly different due to placement of the transition radius. A collection of
%    pre-made distributions is available in package +lambdas.

narginchk(2,3)
if ((nargin == 2) || isempty(lamstrat)), lamstrat = @lambdas.best; end
validateattributes(N, {'numeric'}, {'positive', 'integer'}, '', 'N', 1)
validateattributes(x, {'numeric'}, {'vector', 'numel', 2, 'nonnegative'}, 2)
validateattributes(lamstrat, {'numeric','function_handle'}, {}, '', 'lamstrat', 3)
if isnumeric(lamstrat)
    validateattributes(lamstrat, {'numeric'},...
        {'vector', 'numel', 2, '>', 0, '<', 1}, '', 'lamstrat', 3)
end

cmp = CMSPlanet(N);

if (isa(lamstrat, 'function_handle'))
    lams = lamstrat(N);
    assert(isnumeric(lams) && isvector(lams) && (numel(lams) == N),...
        '@lamstrat(N) must return a vector of length N with values in (0,1].')
    assert(all(lams > 0) && all(lams <= 1),...
        '@lamstrat(N) must return a vector of length N with values in (0,1].')
    cmp.cms.lambdas = lams;
else
    n1 = fix(lamstrat(1)*(N - 1));
    n2 = N - n1 - 1;
    dl1 = lamstrat(2)/(n1 - 1);
    dl2 = (1 - lamstrat(2))/(n2 + 1);
    lam1 = linspace(1 - dl1/2, (1 - lamstrat(2)), n1);
    lam2 = linspace((1 - lamstrat(2)) - dl2, dl2, n2);
    cmp.cms.lambdas = [1, lam1, lam2]';
end

eos0 = barotropes.ConstDensity(0);
eos1 = barotropes.Polytrope(x(1), x(2));
cmp.eos = [eos0; repmat(eos1, N-1, 1)];

end
