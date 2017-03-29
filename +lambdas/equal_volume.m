function lambdas = equal_volume(N, halftop)
%EQUAL_VOLUME Return a lambda distribution making layers of equal volume.
%    lambdas = EQUAL_VOLUME(N) returns an N-vector of normalized radii such that
%    the volume of the layer between lambdas(i) and lambdas(i+1) is constant.
%
%    lambdas = EQUAL_VOLUME(N, halftop) where halftop==true modifies the 2nd
%    element of lambdas to make the thickness of the first layer be exactly half
%    that of the second layer. The default is halftop=false.

narginchk(1,2)
if nargin == 1, halftop = false; end
validateattributes(halftop, {'logical'}, {'scalar'})

lambdas = ((N:-1:1)/N).^(1/3);

if halftop
    lambdas(2) = 1 - (1 - lambdas(3))/3;
end

end
