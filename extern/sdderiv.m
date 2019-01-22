function yp = sdderiv(x,y)
%SDDERIV Scattered data derivative.
%   SDDERIV(x,y) returns the numerical derivative of the one-dimensional function
%   known by its values y on the scattered coordinates x. The derivative is
%   estimated using the 3-point scheme, f'(x0)=(f(x0+h)-f(x0-h))/2h. SDDERIV calls
%   griddedInterpolant(x,y,'pchip') to construct an interpolant which it then uses
%   to evaluate the function at points x+h and x-h. The step size h is chosen to
%   minimize the sum of roundoff and truncation errors. A new optimal step size is
%   chosen for each coordinate in x, allowing an accurate estimate when x spans
%   many orders of magnitude.
%   
%   NOTE: SDDERIV is often significantly more accurate than the matlab built-in
%   function gradient, but it is also much slower.

if nargin == 0 && nargout == 0
    help sdderiv
    return
end

validateattributes(x,{'numeric'},{'vector'})
validateattributes(y,{'numeric'},{'vector'})
assert(isequal(size(x),size(y)),'x and y data must have same size and shape.')

[xi,ind] = unique(x);
yi = y(ind);
F = griddedInterpolant(xi,yi,'pchip');
yp = nan(size(x));
for k=1:length(yp)
    x0 = x(k);
    xc = max(abs(x0),1);   % characteristic function length scale
    ef = (eps(xc)*10);     % wild guess at function's roundoff error
    ef = ef^(1/3);         % optimal power of ef for 3pt-scheme
    h = ef*xc;             % optimal step size
    xx = [x0 - h, x0 + h];
    yy = F(xx);
    yp(k) = 0.5*(yy(2) - yy(1))/h;
end
