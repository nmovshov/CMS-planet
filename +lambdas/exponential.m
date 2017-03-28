%% EXPONENTIAL LAMBA DISTRIBUTION
% Return an exponential lambda distribution.
function lambdas = exponential(N)

lambdas = log((exp(2)-1)*(1:N)/N+1)/2;

lambdas = flip(lambdas);

end