function cmp = triple_polytrope(N, x, lamstrat, forcematch)
%TRIPLE_POLYTROPE Model planet approximated by three polytropes.
%    TRIPLE_POLYTROPE(N, x) returns an N-layer CMSPlanet object with three
%    barotropes.Polytrope eos objects. First polytrope defined by constant x(1)
%    and index x(2) is applied to layers 1:tind-1. Second polytrope defined by
%    constant x(3) and index x(4) is applied to layers tind:cind-1. Third
%    polytrope defined by constant x(5) and index x(6) is applied to layers
%    cind:N. Transition from first to second polytrope is at layer tind, the layer
%    with (ai/a0) nearest x(7). Transition from second to third polytrope is at
%    layer cind, the layer with (ai/a0) nearest x(8). The default layer spacing is
%    the one returned by lambdas.best(N).
%
%    TRIPLE_POLYTROPE(N, x, lamstrat) lets you specify the lambda spacing. Pass a
%    handle to a function that takes a single scalar integer (number of layers)
%    and returns a vector of that length with values in the interval (0, 1], for
%    normalized layer radii. For example, to set layers with equally spaced radii
%    use lamstrat=@(n)linspace(1,1/n,n). Note that the final layer radii might be
%    slightly different due to placement of transition radii. A collection of
%    pre-made distributions is available in package +lambdas.
%
%    TRIPLE_POLYTROPE(..., forcematch) if forcematch=true forces the normalized
%    radii of the transition layers to exactly match x(7) and x(8). This is
%    applied after the initial lambda spacing. The default is forcematch=true.

if nargin == 0
    help('generators.triple_polytrope')
    return
end
narginchk(2,4)
if ((nargin < 3) || isempty(lamstrat)), lamstrat = @lambdas.best; end
if ((nargin < 4) || isempty(forcematch)), forcematch = true; end
validateattributes(N, {'numeric'}, {'positive', 'integer'}, '', 'N', 1)
validateattributes(x, {'numeric'}, {'nonnegative', 'vector', 'numel', 8}, 2)
validateattributes(lamstrat, {'function_handle'}, {}, '', 'lamstrat', 3)
validateattributes(forcematch, {'logical'}, {'scalar'}, '', 'forcematch', 4)
assert(x(7) > 0 && x(7) < 1, 'First transition radius must be in (0,1).')
assert(x(8) > 0 && x(8) < 1, 'Second transition radius must be in (0,1).')
assert(x(8) < x(7), 'Second transition must come before first transition.')

cmp = CMSPlanet;

lams = lamstrat(N);
assert(isnumeric(lams) && isvector(lams) && (numel(lams) == N),...
    '@lamstrat(N) must return a vector of length N with values in (0,1].')
assert(all(lams > 0) && all(lams <= 1),...
    '@lamstrat(N) must return a vector of length N with values in (0,1].')

[~, tind] = min(abs(lams-x(7)));
assert(tind > 1,...
    'First transition too close to surface; first polytrope has zero layers.')
[~, cind] = min(abs(lams-x(8)));
assert(cind > tind,...
    'Transitions are too close together; second polytrope has zero layers.')

if forcematch, lams(tind) = x(7); end
if forcematch, lams(cind) = x(8); end

cmp.ai = lams;
cmp.rhoi = ones(N,1);

eos1 = barotropes.Polytrope(x(1), x(2));
eos2 = barotropes.Polytrope(x(3), x(4));
eos3 = barotropes.Polytrope(x(5), x(6));

cmp.eos = [repmat(eos1, tind - 1, 1);...
           repmat(eos2, cind - tind, 1);...
           repmat(eos3, N - cind + 1, 1)];

for k=1:N
    cmp.rhoi(k) = cmp.eos(k).density(1e5);
end

end
