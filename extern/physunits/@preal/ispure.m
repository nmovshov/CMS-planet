function L = ispure(p)
%PREAL/ISPURE Return true if preal is non-dimensional.

L = true(size(p));
tol=0.001;
for k=1:numel(L)
    if any(abs(p(k).units)>tol)
        L(k) = false;
    end
end