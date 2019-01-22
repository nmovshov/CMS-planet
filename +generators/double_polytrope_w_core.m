function cmp = double_polytrope_w_core(N, x, lamstrat, forcematch)
%DOUBLE_POLYTROPE_W_CORE Two polytrope and core model planet.
%    DOUBLE_POLYTROPE_W_CORE(N, x) returns an N-layer CMSPlanet object with two
%    barotropes.Polytrope eos objects. First polytrope defined by constant x(1)
%    and index x(2) is applied to layers 1:tind-1. Second polytrope defined by
%    constant x(3) and index x(4) is applied to layers tind:N-1. Transition is at
%    layer tind, the layer with (ai/a0) nearest to x(5). Layer N is assigned a
%    barotropes.ConstDensity eos wirh rho0=x(6) and given a normalized equatorial
%    radius equal to x(7). The default layer spacing is the one returned by
%    lambdas.best(N,'I',[x(7),1]).
%
%    DOUBLE_POLYTROPE_W_CORE(N, x, lamstrat) lets you specify the lambda spacing
%    completely. Pass a handle to a function that takes a single scalar integer
%    (number of layers) and returns a vector of that length with values in the
%    interval [x(7),1], for normalized layer radii. For example, to set layers
%    with equally spaced radii in the envelope use
%    lamstrat=@(n)linspace(1,x(7),n).
%
%    DOUBLE_POLYTROPE_W_CORE(..., forcematch) if forcematch=true forces the
%    normalized radius of the transition layer to exactly match x(5). This is
%    applied after the lambda spacing strategy. The default is forcematch=true.

if nargin == 0
    help('generators.double_polytrope_w_core')
    return
end
narginchk(2,4)
if ((nargin < 3) || isempty(lamstrat))
    lamstrat = @(n)lambdas.best(n,'I',[x(7),1]);
end
if ((nargin < 4) || isempty(forcematch)), forcematch = true; end
validateattributes(N, {'numeric'}, {'positive', 'integer'}, '', 'N', 1)
validateattributes(x, {'numeric'}, {'nonnegative', 'vector', 'numel', 7}, 2)
validateattributes(lamstrat, {'function_handle'}, {}, '', 'lamstrat', 3)
validateattributes(forcematch, {'logical'}, {'scalar'}, '', 'forcematch', 4)
assert(x(5) > 0 && x(5) < 1, 'Transition (normalized) radius must be in (0,1).')
assert(x(7) > 0 && x(7) < 1, 'The core must have a normalized radius in (0,1).')
assert(x(7) < x(5), 'Core must be below transition radius.')

cmp = CMSPlanet;

lams = lamstrat(N);
assert(isnumeric(lams) && isvector(lams) && (numel(lams) == N),...
    '@lamstrat(N) must return a vector of length N with values in (0,1].')
assert(all(lams > 0) && all(lams <= 1),...
    '@lamstrat(N) must return a vector of length N with values in (0,1].')
lams = sort(lams, 'descend');
assert(lams(end-1) > x(7),...
    'Chosen layer spacing would put one or more layers below the core.')
if abs(lams(end) - x(7)) > (2/N)
    beep
    warning('Check lambda distribution: is the core where you wanted it?')
end
lams(end) = x(7);

[~, tind] = min(abs(lams-x(5)));
assert(tind > 1,...
    'Transition too close to surface; first polytrope has zero layers.')
assert(tind < N,...
    'Transition too close to core; second polytrope has zero layers.')
if forcematch, lams(tind) = x(5); end

cmp.ai = lams;
cmp.rhoi = ones(N,1);

eos1 = barotropes.Polytrope(x(1), x(2));
eos2 = barotropes.Polytrope(x(3), x(4));
eos3 = barotropes.ConstDensity(x(6));

cmp.eos = [repmat(eos1, tind - 1, 1); repmat(eos2, N - tind, 1); eos3];

for k=1:N
    cmp.rhoi(k) = cmp.eos(k).density(1e5);
end

end
