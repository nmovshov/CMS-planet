function L = ispure(d)
%ISPURE Overloaded PREAL/ISPURE function to run when physunits is disabled.

if isnumeric(d)
    L = true(size(d));
else
    L = false(size(d));
end