function y = midpts(x)
%MIDPTS Values between vector elements.
%   MIDPTS(x) returns (x(1:end-1) + x(2:end))/2.

y = (x(1:end-1) + x(2:end))/2;
