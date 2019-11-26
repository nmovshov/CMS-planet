function [Js, out] = cms(zvec, dvec, qrot, tol, maxiter, sskip, J0s)
%CMS Return gravity coefficients using the Concentric Maclaurin Spheroids method.
%   Js = CMS(zvec, dvec, qrot)
%   [Js, out] = CMS(zvec, dvec, qrot, tol, maxiter, sskip, J0s)
%
% Inputs
% ------
% zvec : 1d array, positive real
%     Equatorial radii of layers where density is specified. Units are unimportant
%     as values will be normalized to outer radius.
% dvec : 1d array, positive real
%     Layer densities. Units are unimportant as values will be normalized to the
%     mean (bulk) density. The density should be monotonically non-increasing with
%     zvec, but this is not enforced. (Note: these are layer densities and NOT the
%     density deltas of concentric spheroids.)
% qrot : scalar, nonnegative
%     Dimensionless rotation parameter. Recall q = w^2a0^3/GM.
% tol : scalar, positive, (tol=1e-6)
%     Convergence tolerance on fractional change in Js in successive iterations.
% maxiter : scalar, positive, integer, (maxiter=100)
%     Maximum number of iterations of CMS algorithm.
% sskip : scalar, nonnegative, integer (sskip=0)
%     Step size for skip-and-spline approximation. The shape functions (zetas)
%     will be explicitly calculated only on every sskip layer, and
%     spline-interpolated in between. This can result in significant speedup with
%     minimal loss of precision, but the sskip value must be chosen by trial and
%     error to fit the required precision and number of layers. For example, a
%     1024-layer model can benefit from almost 16x-speedup by specifying sskip=16
%     while retaining a 10^-6 relative precision on J2.
% J0s : struct
%     J-like values representing initial state. This is not just for speeding up
%     convergence. Mostly it's a mechanism to preserve state between calls.
%
% Outputs
% -------
% Js : 1-by-16 vector, real
%     Even harmonic gravity coefficients J0 to J30. Typically only J2 to J10 are
%     helpful. J0 is included as a sanity check and test of convergence.
% out : struct
%     A structure holding other quantities calculated in the course of running
%     cms. Including out.zetas and out.JLike that together define the converged
%     hydrostatic shape.
%
%CMS(zvec, dvec, qrot, tol=1e-6, maxiter=100, sskip=0, J0s=[])

%% Input parsing
if nargin == 0
    fprintf('Usage:\n  CMS(zvec, dvec, qrot, tol=1e-6, maxiter=100, sskip=0)\n')
    return
end
narginchk(3,7);
if nargin < 4 || isempty(tol), tol = 1e-6; end
if nargin < 5 || isempty(maxiter), maxiter = 100; end
if nargin < 6 || isempty(sskip), sskip = 0; end
if nargin < 7 || isempty(J0s), J0s = struct(); end
validateattributes(zvec,{'numeric'},{'finite','nonnegative','vector'},'','zvec',1)
validateattributes(dvec,{'numeric'},{'finite','nonnegative','vector'},'','dvec',2)
validateattributes(qrot,{'numeric'},{'finite','nonnegative','scalar'},'','qrot',3)
validateattributes(tol,{'numeric'},{'finite','positive','scalar'},'','tol',4)
validateattributes(maxiter,{'numeric'},{'positive','scalar','integer'},'','maxiter',5)
validateattributes(sskip,{'numeric'},{'nonnegative','scalar','integer'},'','sskip',6)
validateattributes(J0s,{'struct'},{'scalar'},'','J0s',7)
assert(length(zvec) == length(dvec),...
    'length(zvec)=%d~=%d=length(dvec)',length(zvec),length(dvec))
[zvec, I] = sort(zvec);
dvec = dvec(I);
zvec = flipud(zvec(:)); % now it's a column for sure
dvec = flipud(dvec(:)); % now it's a column for sure
if zvec(end) == 0, zvec(end) = eps; end

%% Normalize radii and density
dro = [dvec(1); diff(dvec)];
m = sum(dro.*zvec.^3);
robar = m/zvec(1)^3;
zvec = zvec/zvec(1);
dvec = dvec/robar;

%% Initialize local variables (in CMS notation)
lambdas = zvec;
deltas = [dvec(1); diff(dvec)];
nlay = length(lambdas);
nangles = 48;
kmax = 30;

% Initialize zetas as spherical
zetas = ones(nlay, nangles);

% Initialize J-like quantities for spherical planet
if isempty(fieldnames(J0s))
    Jlike = allocate_spherical_Js(nlay, kmax, lambdas, deltas);
else
    Jlike = J0s;
end

% Abscissas and weights for Gaussian quadrature
[mus, gws] = gauleg(0, 1, nangles);

% Precompute Legendre polynomials for fixed colatitudes (gauss quad)
Ps.Pnmu(kmax+1,nangles) = 0;
Ps.Pnzero(kmax+1,1) = 0;
for k=0:kmax
    Ps.Pnmu(k+1,1:nangles) = Pn(k, mus);
    Ps.Pnzero(k+1,1) = Pn(k, 0);
end

% Precompute powers of ratios of lambdas
lamratpow = nan(kmax+2,nlay,nlay);
for ii=1:nlay
    for jj=1:nlay
        for kk=1:kmax+2
            lamratpow(kk,ii,jj) = ...
                (lambdas(ii)/lambdas(jj))^(kk-1);
        end
    end
end

%% The loop (see Hubbard, 2012 and ./notes/CMS.pdf)
Js = Jlike.Jn;
for iter=1:maxiter
    % Update shape with current gravity
    if sskip == 0
        new_zetas = update_zetas(Jlike, Ps, lamratpow, qrot, zetas);
    else
        new_zetas = skipnspline_zetas(Jlike, Ps, lamratpow, qrot, zetas, zvec, sskip);
    end
    
    % Update gravity with current shape
    new_Jlike = update_Js(lambdas, deltas, new_zetas, Ps, gws);
    new_Js = new_Jlike.Jn;
    
    % Check for convergence of J0-J8 to terminate...
    dJs = abs(Js - new_Js)./abs(Js+eps);
    if all(dJs(1:5) < tol), break, end
    
    % ... or update to new values and continue
    Jlike = new_Jlike;
    Js = new_Js;
    zetas = new_zetas;
end

% It's not always a disaster if maxiter is reached, but we'd like to know
if iter == maxiter
    warning('CMS:maxiter','Shape may not be fully converged.')
end

%% Return
Js = new_Js; % may as well use the latest...
out.dJs = dJs;
out.iter = iter;
out.zetas = new_zetas;
out.lambdas = lambdas;
out.deltas = deltas;
out.JLike = new_Jlike;
out.mus = mus;
out.gws = gws;
out.Ps = Ps;

end

%% Helper functions
function Js = allocate_spherical_Js(nlay,nmom,lambdas,deltas)
Js.tilde = zeros(nlay,(nmom+1));
Js.tildeprime = zeros(nlay,(nmom+1));
Js.tildeprimeprime = zeros(nlay,1);
Js.Jn = zeros(1,nmom/2+1);
den = sum(deltas.*lambdas.^3);
Js.tilde(:,1) = -1*(deltas.*lambdas.^3)/den;
Js.tildeprime(:,1) = -1.5*(deltas.*lambdas.^3)/den;
Js.tildeprimeprime(:) = 0.5*(deltas.*lambdas.^3)/den;
Js.Jn(1) = sum(Js.tilde(:,1));
end

function y = Pn(n,x)
% Fast implementation of ordinary Legendre polynomials of low degree.
switch n
    case 0
        y = ones(size(x));
    case 1
        y = x;
    case 2
        y = 0.5*(3*x.^2 - 1);
    case 3
        y = 0.5*(5*x.^3 - 3*x);
    case 4
        y = (1/8)*(35*x.^4 - 30*x.^2 + 3);
    case 5
        y = (1/8)*(63*x.^5 - 70*x.^3 + 15*x);
    case 6
        y = (1/16)*(231*x.^6 - 315*x.^4 + 105*x.^2 - 5);
    case 7
        y = (1/16)*(429*x.^7 - 693*x.^5 + 315*x.^3 - 35*x);
    case 8
        y = (1/128)*(6435*x.^8 - 12012*x.^6 + 6930*x.^4 - 1260*x.^2 + 35);
    case 9
        y = (1/128)*(12155*x.^9 - 25740*x.^7 + 18018*x.^5 - 4620*x.^3 + 315*x);
    case 10
        y = (1/256)*(46189*x.^10 - 109395*x.^8 + 90090*x.^6 - 30030*x.^4 + 3465*x.^2 - 63);
    case 11
        y = (1/256)*(88179*x.^11 - 230945*x.^9 + 218790*x.^7 - 90090*x.^5 + 15015*x.^3 - 693*x);
    case 12
        y = (1/1024)*(676039*x.^12 - 1939938*x.^10 + 2078505*x.^8 - 1021020*x.^6 + 225225*x.^4 - 18018*x.^2 + 231);
    otherwise
        assert(isvector(x))
        Pnm = legendre(n,x);
        y = Pnm(1,:);
        if ~isrow(x), y = y'; end
end
end

function [x,w] = gauleg(x1,x2,n)
%GAULEG Calculate abscissas and weights for Gauss-Legendre n-point quadrature.
%   [x,w] = GAULEG(x1,x2,n) returns the abscissas x and weights w that can be
%   used to evaluate the definite integral, I, of a function well approximated
%   by an (2n - 1) degree polynomial in the interval [x1,x2] using the
%   Gauss-Legendre formula:
%
%       I = sum(w.*f(x))
%
%   Algorithm
%     This function is based on the C++ implementation of a routine with the
%     same name in Numerical Recipes, 3rd Edition. But in several places I opt
%     for readability over performance, on the assumption that this function is
%     most likely to be called in a setup routine rather than in an inner-loop
%     computation.
%
%   Example
%     fun = @(x)sin(x);
%     [x,w] = gauleg(0,pi,6);
%     I_adaptive = integral(fun,0,pi)
%     I_gaussleg = sum(w.*fun(x))
%
% Author: Naor Movshovitz (nmovshov at google dot com)
%         Earth and Planetary Sciences, UC Santa Cruz
%
% Reference: William H. Press, Saul A. Teukolsky, William T. Vetterling, and
% Brian P. Flannery. 2007. Numerical Recipes 3rd Edition: The Art of Scientific
% Computing (3 ed.). Cambridge University Press, New York, NY, USA.

% Input parsing and minimal assertions
narginchk(3,3)
nargoutchk(2,2)
validateattributes(x1,{'numeric'},{'scalar','finite','real'},1)
validateattributes(x2,{'numeric'},{'scalar','finite','real'},2)
validateattributes(n,{'numeric'},{'scalar','finite','integer','>=',2},3)
assert(x2 > x1, 'Interval must be positive.');

% Local variables
tol = 1e-14;
m = ceil(n/2);
xmid = (x1 + x2)/2;
dx = (x2 - x1);
x = NaN(1,n);
w = NaN(1,n);

% Main loop
for j=1:m
    % Get j-th root of Legendre polynomial Pn, along with Pn' value there.
    z = cos(pi*((j - 1) + 0.75)/(n + 0.5)); % initial guess for j-th root
    while true
        % Calculate Pn(z) and Pn-1(z) and Pn'(z)
        p = NaN(1,n+1);
        p(1) = 1;
        p(2) = z;
        for k=2:n
            pkm1 = p(k);
            pkm2 = p(k-1);
            pk = (1/k)*((2*k - 1)*z*pkm1 - (k - 1)*pkm2);
            p(k+1) = pk;
        end
        pn = p(end);
        pp = (n*p(end-1) - n*z*p(end))/(1 - z^2);
        
        % And now Newton's method (we are hopefully very near j-th root)
        oldz = z;
        z = z - pn/pp;
        if abs(z - oldz) < tol, break, end
    end
    
    % Now use j-th root to get 2 abscissas and weights
    x(j)     = xmid - z*dx/2; % Scaled abscissa left of center
    x(n+1-j) = xmid + z*dx/2; % Scaled abscissa right of center
    w(j)     = dx/((1 - z^2)*pp^2);
    w(n+1-j) = w(j);
end

% Verify and return
assert(all(isfinite(x)))
assert(all(isfinite(w)))
end

function newzetas = update_zetas(Js, Ps, lamrats, qrot, oldzetas)
% Update level surfaces using current value of Js.

% Loop over layers (outer) and colatitudes (inner)
nlay = size(lamrats,2);
nangles = size(Ps.Pnmu,2);
newzetas = NaN(nlay,nangles);
for j=1:nlay
    for alfa=1:nangles
        oldzeta = oldzetas(j,alfa);
        newzetas(j,alfa) = zeta_j_of_alfa(j, alfa, Js, Ps, lamrats, qrot, oldzeta);
    end
end
end

function newzetas = skipnspline_zetas(Js, Ps, lamrats, qrot, oldzetas, zvec, sskip)
% Update layer shapes using current value of Js.

nlay = size(Js.tilde,1);
ind = 1:sskip:nlay;
nangles = size(Ps.Pnmu,2);
Y = NaN(length(ind),nangles);

% skip...
parfor j=1:length(ind)
    for alfa=1:nangles
        oldzeta = oldzetas(ind(j),alfa);
        Y(j,alfa) = zeta_j_of_alfa(ind(j), alfa, Js, Ps, lamrats, qrot, oldzeta);
    end
end

% ... and spline
newzetas = NaN(nlay,nangles);
for alfa = 1:nangles
    newzetas(:,alfa) = spline(zvec(ind), Y(:,alfa), zvec);
end

end

function newJs = update_Js(lambdas, deltas, zetas, Ps, gws)
% Single-pass update of gravitational moments by Gaussian quad.

nlay = length(lambdas);
kmax = length(Ps.Pnzero)-1;

% Do common denominator in eqs. (48)
denom = 0;
for j=1:nlay
    fun = zetas(j,:).^3;
    I = gws*fun'; % gauss quad formula
    denom = denom + deltas(j)*lambdas(j)^3*I;
end

% Do J tilde, eq. (48a)
new_tilde = zeros(nlay,kmax+1);
for ii=1:nlay
    for kk=0:kmax
        if rem(kk, 2), continue, end
        fun = Ps.Pnmu(kk+1,:).*zetas(ii,:).^(kk+3);
        I = gws*fun'; % gauss quad formula
        new_tilde(ii,kk+1) = -(3/(kk + 3))*deltas(ii)*lambdas(ii)^3*I/denom;
    end
end

% Do J tilde prime, eqs. (48b and 48c)
new_tprime = zeros(nlay,kmax+1);
for ii=1:nlay
    for kk=0:kmax
        if rem(kk, 2), continue, end
        if kk == 2 % eq. (48c)
            fun = Ps.Pnmu(3,:).*log(zetas(ii,:));
            I = gws*fun'; % gauss quad formula
            new_tprime(ii,kk+1) = -3*deltas(ii)*lambdas(ii)^3*I/denom;
        else       % eq. (48b)
            fun = Ps.Pnmu(kk+1,:).*zetas(ii,:).^(2 - kk);
            I = gws*fun'; % gauss quad formula
            new_tprime(ii,kk+1) = -(3/(2 - kk))*deltas(ii)*lambdas(ii)^3*I/denom;
        end
    end
end

% Do J tilde double prime, eq. (48d)
new_tpprime = zeros(nlay,1);
for ii=1:nlay
    new_tpprime(ii) = 0.5*deltas(ii)*lambdas(ii)^3/denom;
end

% Return updated Js struct
newJs.tilde = new_tilde;
newJs.tildeprime = new_tprime;
newJs.tildeprimeprime = new_tpprime;
n = 0:2:kmax;
for k=1:length(n)
    newJs.Jn(k) = dot(newJs.tilde(:,n(k)+1),lambdas.^n(k));
end

end

function y = zeta_j_of_alfa(j, alfa, Js, Ps, lamrats, qrot, oldzeta)
persistent os
if isempty(os), os = optimset('TolX',1e-12); end
fun = @(x)eq52(x,j,alfa,Js,Ps,lamrats,qrot);
y = fzero(fun, oldzeta, os);
end

function y = eq52(zja, jl, alfa, Js, Ps, lamrats, qrot)
% locals
nlay = size(Js.tilde,1);
kmax = length(Ps.Pnzero)-1;
Jt = Js.tilde;
Jtp = Js.tildeprime;
Jtpp = Js.tildeprimeprime;
lamj3 = lamrats(4,jl,1);
q = qrot;
P0 = Ps.Pnzero;
Pmu = Ps.Pnmu(:,alfa);
zetpow = zja.^(0:kmax+1);
zetipow = zja.^-(0:kmax+1);

% first double sum
y1 = 0;
for ii=jl:nlay
    for kk=0:2:kmax
        y1 = y1 + lamrats(kk+1,ii,jl)*Jt(ii,kk+1)*zetipow(kk+1)*Pmu(kk+1);
    end
end

% second double sum
y2 = 0;
for ii=jl:nlay
    for kk=0:2:kmax
        y2 = y2 + lamrats(kk+1,ii,jl)*Jt(ii,kk+1)*P0(kk+1);
    end
end

% third double sum
y3 = 0;
for ii=1:jl-1
    for kk=0:2:kmax
        y3 = y3 + lamrats(kk+2,jl,ii)*Jtp(ii,kk+1)*(zetpow(kk+2))*Pmu(kk+1);
    end
end

% and forth double sum
y4 = 0;
for ii=1:jl-1
    for kk=0:2:kmax
        y4 = y4 + lamrats(kk+2,jl,ii)*Jtp(ii,kk+1)*P0(kk+1);
    end
end

% a single sum
y5 = 0;
for ii=1:jl-1
    y5 = y5 + lamrats(4,jl,ii)*Jtpp(ii,1)*zetpow(4);
end

% another single sum
y6 = 0;
for ii=1:jl-1
    y6 = y6 + lamrats(4,jl,ii)*Jtpp(ii,1);
end

% and the rotation term
y7 = -(1/3)*q*lamj3*zja^2*(1 - Pmu(3)) + (1/2)*q*lamj3;

% Now combine
y = (1/zja)*(y1 + y3 + y5) - y2 - y4 - y6 + y7;

end
