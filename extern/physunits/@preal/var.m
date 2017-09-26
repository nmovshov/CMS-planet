function p2 = var(p1)
%PREAL/VAR Overloaded VAR function for class preal.
% Note: Unlike datafun/var preal/var simply returns the variance of ALL
% elements of the array, resulting in one scalar.

global useUnitsFlag

if ~(useUnitsFlag) % If physunits is disabled...
    p2=var(double(p1)); % ... treat as double.
    return
end

p2 = sum((p1 - mean(p1)).^2)/(numel(p1) - 1);
