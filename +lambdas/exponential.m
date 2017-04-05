function lambdas = exponential(N)
%EXPONENTIAL Return a lambda distribution with exponentially spaced radii.
%    lambdas = EXPONENTIAL(N) returns an N-vector of exponentially spaced
%    normalized radii.

lambdas = log((exp(2)-1)*(N:-1:1)/N+1)/2;

end
