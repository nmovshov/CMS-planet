%% EQUIVOLUME LAMBA DISTRIBUTION
% Return a lambda distribution whose spheroids have equal volume.
function lambdas = equi_volume(N)

lambdas = power((1:N)/N, 1/3);

lambdas = flip(lambdas);

end