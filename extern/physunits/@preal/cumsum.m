function p3 = cumsum(p1)
%PREAL/CUMSUM Overloaded CUMSUM function for class preal.

global useUnitsFlag

if ~(useUnitsFlag) % If physunits is disabled...
    p3=cumsum(double(p1)); % ... treat as double.
    return
end

if isvector(p1)
    p3=preal(zeros(size(p1)));
    p3(1)=p1(1);
    for k=2:numel(p1)
        p3(k)=p3(k-1)+p1(k);
    end
else
    error('preal/cumsum is only defined for vectors')
end
