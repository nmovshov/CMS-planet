function lambdas = exponential(N, halftop)
%EXPONENTIAL Return a lambda distribution with exponentially spaced radii.
%    lambdas = EXPONENTIAL(N) returns an N-vector of exponentially spaced
%    normalized radii.

lambdas = log((exp(2)-1)*(N:-1:1)/N+1)/2;

if nargin == 1, halftop = false; end
if halftop
    dl = 1 - lambdas(2);
    lambdas(2:end) = lambdas(2:end) + dl/2;
end

end
