function p3=colon(p1,dp,p2)
%PREAL/COLON Oveloaded colon operator for class PREAL.

global useUnitsFlag

if nargin<2 || nargin>3
    error('don''t be a dick')
end

if nargin==2
    p2=dp;
    dp=preal(1,get(p1,'units'));
end

if ~(useUnitsFlag) % If physunits is disabled...
    p3=colon(double(p1),double(dp),double(p2)); % ... treat as double.
    return
end

p1=preal(p1);
p2=preal(p2);
dp=preal(dp);
tol=0.001;

if ~isscalar(p1)||~isscalar(p2)||~isscalar(dp)
    error('Input points must be scalar.')
end

if any(abs(p1.units-p2.units)>tol)
    error('Input points must be of same dimension.')
end

if any(abs(p1.units-dp.units)>tol)
    error('Input points must be of same dimension.')
end

u=get(p1,'units');
p0=preal(1,u);
p3=colon(double(p1),double(dp),double(p2))*p0;
