%% EXPONENTIAL LAMBA DISTRIBUTION
% Return an exponential lambda distribution.
function lambdas = exponential(N)

lambdas = log((exp(6)-1)*(1:N)/N+1)/6;

lambdas = flip(lambdas);

end