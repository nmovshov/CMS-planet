function p2 = mean(p1)
%PREAL/MEAN Overloaded MEAN function for class preal.
% Note: Unlike datafun/mean preal/mean simply returns the arithmetic mean of ALL
% elements of the array, resulting in one scalar.

global useUnitsFlag

if ~(useUnitsFlag) % If physunits is disabled...
    p2=mean(double(p1)); % ... treat as double.
    return
end

p2 = sum(p1)/numel(p1);
