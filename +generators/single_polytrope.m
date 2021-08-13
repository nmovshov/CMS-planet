function cmp = single_polytrope(N, x, zstrat)
%SINGLE_POLYTROPE The simplest toy model planet.
%    SINGLE_POLYTROPE(N, x) returns an N-layer CMSPlanet object with a
%    barotropes.Polytrope eos, with constant x(1) and index x(2). The default
%    layer spacing is one of equal radius increments between a/a0=1 and
%    a/a0=1/N.
%
%    SINGLE_POLYTROPE(N, x, zstrat) lets you specify the zvec distribution
%    strategy. Pass a handle to a function that takes a single scalar integer
%    (number of layers) and returns a vector of that length with values in the
%    interval (0, 1], for normalized equatorial layer radii. Some pre-made
%    distributions are available in +zvecs.

if nargin == 0
    help generators.single_polytrope
    return
end

narginchk(2,3)
if ((nargin < 3) || isempty(zstrat)), zstrat = @(n)linspace(1,1/n,n); end
validateattributes(N,{'numeric'},{'positive','integer'},'','N',1)
validateattributes(x,{'numeric'},{'vector','numel',2,'nonnegative'},2)
validateattributes(zstrat,{'function_handle'},{},'','zstrat',3)

% Create and validate level distribution
zvec = zstrat(N);
assert(isnumeric(zvec) && isvector(zvec) && (numel(zvec) == N),...
    '@zstrat(N) must return a vector of length N with values in (0,1].')
assert(all(zvec > 0) && all(zvec <= 1),...
    '@zstrat(N) must return a vector of length N with values in (0,1].')

% Create planet object and assign eos
cmp = CMSPlanet();
cmp.ai = zvec;
cmp.rhoi = ones(N,1);
cmp.eos = barotropes.Polytrope(x(1), x(2));

% Initialize density at 1-bar values, just because
cmp.rhoi = cmp.rhoi*cmp.eos.density(1e5);

end
