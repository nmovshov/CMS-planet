function cmp = double_polytrope(N, x, lamstrat, forcematch)
%DOUBLE_POLYTROPE Model planet approximated by two polytropes.
%    DOUBLE_POLYTROPE(N, x) returns an N-layer CMSPlanet object with two
%    barotropes.Polytrope eos objects. First polytrope defined by constant x(1)
%    and index x(2) is applied to layers 1:tind-1. Second polytrope defined by
%    constant x(3) and index x(4) is applied to layers tind:N. Transition from
%    first to second polytrope is at layer tind, the layer with (ai/a0) nearest
%    x(5). The default layer spacing is the one returned by lambdas.best(N).
%
%    DOUBLE_POLYTROPE(N, x, lamstrat) lets you specify the lambda spacing. Pass a
%    handle to a function that takes a single scalar integer (number of layers)
%    and returns a vector of that length with values in the interval (0, 1], for
%    normalized layer radii. For example, to set layers with equally spaced radii
%    use lamstrat=@(n)linspace(1,1/n,n). Note that the final layer radii might be
%    slightly different due to placement of transition radii. A collection of
%    pre-made distributions is available in package +lambdas.
%
%    DOUBLE_POLYTROPE(..., forcematch) if forcematch=true forces the normalized
%    radii of the transition layer to exactly match x(5). This is applied after
%    the initial lambda spacing. The default is forcematch=true.

if nargin == 0
    help('generators.double_polytrope')
    return
end
narginchk(2,4)
if ((nargin < 3) || isempty(lamstrat)), lamstrat = @lambdas.best; end
if ((nargin < 4) || isempty(forcematch)), forcematch = true; end
validateattributes(N, {'numeric'}, {'positive', 'integer'}, '', 'N', 1)
validateattributes(x, {'numeric'}, {'nonnegative', 'vector', 'numel', 5}, 2)
validateattributes(lamstrat, {'function_handle'}, {}, '', 'lamstrat', 3)
validateattributes(forcematch, {'logical'}, {'scalar'}, '', 'forcematch', 4)
assert(x(5) > 0 && x(5) < 1, 'Transition radius must be in (0,1).')

cmp = CMSPlanet;

lams = lamstrat(N);
assert(isnumeric(lams) && isvector(lams) && (numel(lams) == N),...
    '@lamstrat(N) must return a vector of length N with values in (0,1].')
assert(all(lams > 0) && all(lams <= 1),...
    '@lamstrat(N) must return a vector of length N with values in (0,1].')

[~, tind] = min(abs(lams-x(5)));
assert(tind > 1,...
    'Transition too close to surface; first polytrope has zero layers.')

if forcematch, lams(tind) = x(5); end

cmp.ai = lams;
cmp.rhoi = ones(N,1);

eos1 = barotropes.Polytrope(x(1), x(2));
eos2 = barotropes.Polytrope(x(3), x(4));

cmp.eos = [repmat(eos1, tind - 1, 1);...
           repmat(eos2, N + 1 - tind, 1)];

for k=1:N
    cmp.rhoi(k) = cmp.eos(k).density(1e5);
end

end
