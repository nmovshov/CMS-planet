function lambdas = equal_volume(N, halftop)
%EQUAL_VOLUME Return a lambda distribution making layers of equal volume.
%    lambdas = EQUAL_VOLUME(N) returns an N-vector of normalized radii such that
%    the volume of the layer between lambdas(i) and lambdas(i+1) is constant.

lambdas = ((N:-1:1)/N).^(1/3);

if nargin == 1, halftop = false; end
if halftop
    dl = 1 - lambdas(2);
    lambdas(2:end) = lambdas(2:end) + dl/2;
end

end
