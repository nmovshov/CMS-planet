function p3 = floor(p1)
%PREAL/FLOOR Overloaded FLOOR function for class preal.

global useUnitsFlag

if ~(useUnitsFlag) % If physunits is disabled...
    p3=floor(double(p1)); % ... treat as double.
    return
end

p3=preal(ones(size(p1)));
for k=1:numel(p3)
    p3(k).value=floor(p1(k).value);
    p3(k).units=p1(k).units;
end