function L = isa(~,typeStr)
%PREAL/ISA Overloaded ISA function for class preal.
%
% ISA(p,'numeric') returns true, so that calls to validateattributes() with
% class category 'numeric' will pass a preal variable.

if any(strcmpi(typeStr,{'preal','numeric'}))
    L = true;
else
    L = false;
end
