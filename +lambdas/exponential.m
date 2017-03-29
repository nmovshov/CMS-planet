function lambdas = exponential(N, halftop)
%EXPONENTIAL Return a lambda distribution with exponentially spaced radii.
%    lambdas = EXPONENTIAL(N) returns an N-vector of exponentially spaced
%    normalized radii.
%
%    lambdas = EXPONENTIAL(N, halftop) where halftop==true modifies the 2nd
%    element of lambdas to make the thickness of the first layer be exactly half
%    that of the second layer. The default is halftop=false.

narginchk(1,2)
if nargin == 1, halftop = false; end
validateattributes(halftop, {'logical'}, {'scalar'})

lambdas = log((exp(2)-1)*(N:-1:1)/N+1)/2;

if halftop
    lambdas(2) = 1 - (1 - lambdas(3))/3;
end

end
