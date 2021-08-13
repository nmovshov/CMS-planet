function zvec = cosine(N, alfa)
%COSINE Return a zvec distribution with cosine-like spacing.
%    zvec = COSINE(N) returns the vector cos(pi/2*(0:N-1)/N).
%
%    zvec = COSINE(N, alfa) returns the vector cos(pi/2*((0:N-1)/N).^alfa).

narginchk(1,2)
if (nargin < 2) || isempty(alfa), alfa = 1; end
validateattributes(N,{'numeric'},{'positive','integer','scalar'})
validateattributes(alfa,{'numeric'},{'real','finite','positive'})

zvec = cos(pi/2*((0:N-1)/N).^alfa);

end
