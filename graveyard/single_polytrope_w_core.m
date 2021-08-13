function cmp = single_polytrope_w_core(N, x, lamstrat)
%SINGLE_POLYTROPE_W_CORE Polytrope-and-core planet.
%    SINGLE_POLYTROPE_W_CORE(N, x) returns an N-layer CMSPlanet object with a
%    barotropes.Polytrope with constant x(1) and index x(2) assigned to layers
%    1:N-1 and a barotropes.ConstDensity eos with rho0=x(3) assigned to layer N
%    which has a normalized equatorial radius equal to x(4). The default layer
%    spacing is the one returned by lambdas.best(N,'I',[x(4),1]).
%
%    SINGLE_POLYTROPE_W_CORE(N, x, lamstrat) lets you specify the lambda spacing
%    completely. Pass a handle to a function that takes a single scalar integer
%    (number of layers) and returns a vector of that length with values in the
%    interval [x(4),1], for normalized layer radii.

if nargin == 0
    help('generators.single_polytrope_w_core')
    return
end
narginchk(2,3)
if ((nargin < 3) || isempty(lamstrat))
    lamstrat = @(N)lambdas.best(N,'I',[x(4),1]);
end
validateattributes(N, {'numeric'}, {'positive', 'integer'}, '', 'N', 1)
validateattributes(x, {'numeric'}, {'nonnegative', 'vector', 'numel', 4}, 2)
validateattributes(lamstrat, {'function_handle'}, {}, '', 'lamstrat', 3)
assert(x(4) > 0 && x(4) < 1, 'The core must have a normalized radius in (0,1).')

cmp = CMSPlanet;

lams = lamstrat(N);
assert(isnumeric(lams) && isvector(lams) && (numel(lams) == N),...
    '@lamstrat(N) must return a vector of length N with values in (0,1].')
assert(all(lams > 0) && all(lams <= 1),...
    '@lamstrat(N) must return a vector of length N with values in (0,1].')
lams = sort(lams, 'descend');
assert(lams(end-1) > x(4),...
    'Chosen layer spacing would put one or more layers below the core.')
if abs(lams(end) - x(4)) > (2/N)
    beep
    warning('Check lambda distribution: is the core where you wanted it?')
end
lams(end) = x(4);

cmp.ai = lams;
cmp.rhoi = ones(N,1);

eos1 = barotropes.Polytrope(x(1), x(2));
eos2 = barotropes.ConstDensity(x(3));
cmp.eos = [repmat(eos1, N - 1, 1); eos2];

for k=1:N
    cmp.rhoi(k) = cmp.eos(k).density(1e5);
end

end
