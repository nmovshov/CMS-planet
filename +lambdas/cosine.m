function lambdas = cosine(N, alfa, halftop)
%COSINE Return a lambda distribution with cosine-like spacing.
%    lambdas = COSINE(N) returns the vector cos(pi/2*(0:N-1)/N).
%
%    lambdas = COSINE(N, alfa) returns the vector cos(pi/2*((0:N-1)/N).^alfa).
%
%    lambdas = COSINE(..., halftop) where halftop==true modifies the 2nd element
%    of lambdas to make the thickness of the first layer be exactly half that of
%    the second layer. The default is halftop=false.

narginchk(1,3)
if (nargin < 2) || isempty(alfa), alfa = 1; end
if nargin < 3, halftop = false; end
validateattributes(alfa, {'numeric'}, {'real', 'finite', 'positive'})
validateattributes(halftop, {'logical'}, {'scalar'})

lambdas = cos(pi/2*((0:N-1)/N).^alfa);

if halftop
    lambdas(2) = 1 - (1 - lambdas(3))/3;
end

end
