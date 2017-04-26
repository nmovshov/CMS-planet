function lams = best(N, varargin)
%BEST Return our current notion of the best lambda distribution.
%    lams = BEST(N) returns a vector of normalized radii representing our current
%    idea of the "best" distribution for minimizing discretization error in an
%    N-layer CMS model. The actual best distribution may depend on planetary
%    conditions and/or on whether we seek to converge to a prescribed barotrope or
%    to a prescribed density profile.
%
%    lams = BEST(N, Name, Value) specifies non-default values for various
%    parameters. This is mostly used during development and experimentation, after
%    which the default values will be set to the "best" value,

narginchk(1,inf)
validateattributes(N, {'numeric'}, {'positive', 'integer', 'scalar'}, '', 'N', 1)
opts = parse_others(varargin{:});

lams = lambdas.topheavy(N, opts.skew, opts.halftop);

end

function options = parse_others(varargin)
p = inputParser;
p.addParameter('H', 0.01, @(x)isscalar(x) && x >=0 && x <= 0.1)
p.addParameter('n_per_H', 30)
p.addParameter('skew', [3/4, 1/2])
p.addParameter('halftop', true, @(x)islogical(x)&&isscalar(x))

p.parse(varargin{:})
options = p.Results;

end
