function lambdas = equal_dr(N, halftop)
%EQUAL_DR Return a lambda distribution with equally spaced radii.
%    lambdas = EQUAL_DR(N) returns an N-vector of equally spaced normalized radii.
%
%    lambdas = EQUAL_DR(N, halftop) where halftop==true modifies the 2nd element
%    of lambdas to make the thickness of the first layer be exactly half that of
%    the second layer. The default is halftop=true.

narginchk(1,2)

lambdas = linspace(1, 1/N, N)';

if nargin == 1, halftop = true; end
validateattributes(halftop, {'logical'}, {'scalar'})
if halftop
    lambdas(2) = 1 - (1 - lambdas(3))/3;
end

end
