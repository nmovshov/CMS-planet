function zvec = best(N, varargin)
%BEST Return our current notion of the best zvec distribution.
%    zvec = BEST(N) returns a vector of normalized radii representing our
%    current idea of the "best" distribution for minimizing discretization
%    error in an N-layer CMS model. The actual best distribution may depend on
%    planetary conditions and/or on whether we seek to converge to a prescribed
%    barotrope or to a prescribed density profile. And it will certainly need
%    to be adjusted in the presence of density discontinuities.

narginchk(1,1)
validateattributes(N, {'numeric'}, {'positive', 'integer', 'scalar'}, '', 'N', 1)

zvec = zvecs.equal_dr(N);

end
