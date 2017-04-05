function lambdas = cosine(N, alfa)
%COSINE Return a lambda distribution with cosine-like spacing.
%    lambdas = COSINE(N) returns the vector cos(pi/2*(0:N-1)/N).
%
%    lambdas = COSINE(N, alfa) returns the vector cos(pi/2*((0:N-1)/N).^alfa).

narginchk(1,2)
if (nargin < 2) || isempty(alfa), alfa = 1; end
validateattributes(alfa, {'numeric'}, {'real', 'finite', 'positive'})

lambdas = cos(pi/2*((0:N-1)/N).^alfa);

end
