function p3 = fix(p1)
%PREAL/FIX Overloaded FIX function for class preal.

global useUnitsFlag

if ~(useUnitsFlag) % If physunits is disabled...
    p3=fix(double(p1)); % ... treat as double.
    return
end

p3=preal(ones(size(p1)));
for k=1:numel(p3)
    p3(k).value=fix(p1(k).value);
    p3(k).units=p1(k).units;
end