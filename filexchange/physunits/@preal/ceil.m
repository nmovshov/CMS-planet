function p3 = ceil(p1)
%PREAL/CEIL Overloaded CEIL function for class preal.

global useUnitsFlag

if ~(useUnitsFlag) % If physunits is disabled...
    p3=ceil(double(p1)); % ... treat as double.
    return
end

p3=preal(ones(size(p1)));
for k=1:numel(p3)
    p3(k).value=ceil(p1(k).value);
    p3(k).units=p1(k).units;
end