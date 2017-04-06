function lambdas = equal_dr(N, halftop)
%EQUAL_DR Return a lambda distribution with equally spaced radii.
%    lambdas = EQUAL_DR(N) returns an N-vector of equally spaced normalized radii.
%
%    lambdas = EQUAL_DR(N, halftop) where halftop==true makes the thickness of the
%    first layer be exactly half that of the other layers. The default is
%    halftop=false.

narginchk(1,2)
if nargin == 1, halftop = false; end
validateattributes(halftop, {'logical'}, {'scalar'})

dl = 1/N;
if halftop
    lambdas = [1, (1 - dl/2):-dl:1/N];
else
    lambdas = 1:-dl:1/N;
end

end
