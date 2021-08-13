function [Js, out] = cms(zvec, dvec, qrot, varargin)
%CMS Concentric Maclaurin Spheroids equilibrium shape and gravity.
%   Js = CMS(zvec, dvec, qrot) returns 1-by-16 vector Js of gravity
%   coefficients J0 through J30 of a rotating fluid planet in hydrostatic
%   equilibrium. Coefficients are stored in ascending order so that Js(1) is
%   J0, Js(2) is J2, Js(3) is J4, etc. The mandatory inputs are a vector of
%   equatorial radii zvec, vector of corresponding layer densities dvec, and
%   rotation parameter qrot, assumed normalized to the outer layer's equatorial
%   radius (to zvec(1)).
%
%   [Js, out] = CMS(zvec, dvec, qrot, 'NAME1',VALUE1, 'NAME2',VALUE2,...)
%   accepts additional parameters as NAME/VALUE pairs, and also returns an
%   output struct holding diagnostic values and additional derived quantities,
%   including the full hydrostatic spheroid shapes.
%
% Inputs, required
% ----------------
% zvec : 1d array, positive real
%     Equatorial radii of constant density layers, indexed from the outside in,
%     i.e., zvec(1)=a0 is the outer radius of the outermost layer, zvec(2) is
%     the inner radius of the outermost layer as well as the outer radius of
%     the next layer, etc. The innermost layer extends all the way to the
%     center, so that zvec(end) is the outer radius of a central spheroid
%     layer. Units of zvec are unimportant as values will be normalized to
%     outer radius.
% dvec : 1d array, positive real
%     Layer densities. The layer lying between zvec(i) and zvec(i+1) has
%     constant density dvec(i). Units are unimportant as values will be
%     normalized to the mean (bulk) density. The density should be
%     monotonically non-increasing with zvec, but this is not enforced. (Note:
%     these are layer densities and NOT the density "deltas" of concentric
%     spheroids.)
% qrot : scalar, nonnegative
%     Dimensionless rotation parameter. Recall q = w^2a0^3/GM.
%
% Inputs, NAME/VALUE pairs
% tol : scalar, positive, (tol=1e-6)
%     Convergence tolerance on fractional change in Js in successive iterations.
% maxiter : scalar, positive, integer, (maxiter=100)
%     Maximum number of iterations of CMS algorithm.
% xlayers : scalar or vector, nonnegative, integer (xlayers=-1)
%     Layers whose shape will be explicitly calculated. The shape functions
%     (zetas) will be explicitly calculated for these layers, and
%     spline-interpolated in between. This can result in significant speedup
%     with minimal loss of precision, if the xlayers are chosen by trial and
%     error to fit the required precision and the spacing of density layers. A
%     scalar value is interpreted as a number of xlayers to be uniformaly
%     distributed among the density layers. For example, a smooth-density
%     1024-layer model can benefit from almost 16x-speedup by specifying
%     xlayers=64 while retaining a 10^-6 relative precision on J2. A vector
%     value is interpreted as indices of layers to be used as xlayers. (A
%     negative value is a shortcut to flag a full calculation instead of
%     skip-n-spline.)
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
% Algorithm
% ---------
% Concentric Maclaurin Spheroids from Hubbard 2012, 2013 ApJ/ApJL.

%% Input parsing
% Zero inputs case, usage only
if nargin == 0
    print_usage()
    return
end
narginchk(3,inf);

% Mandatory inputs
validateattributes(zvec,{'numeric'},{'finite','nonnegative','vector'},'','zvec',1)
validateattributes(dvec,{'numeric'},{'finite','nonnegative','vector'},'','dvec',2)
validateattributes(qrot,{'numeric'},{'finite','nonnegative','scalar'},'','qrot',3)
assert(length(zvec) == length(dvec),...
    'length(zvec)=%d~=%d=length(dvec)',length(zvec),length(dvec))
[zvec, I] = sort(zvec);
dvec = dvec(I);
zvec = flipud(zvec(:)); % now it's a column for sure
dvec = flipud(dvec(:)); % now it's a column for sure
if zvec(end) == 0, zvec(end) = eps; end

% Optional arguments
opts = parsem(varargin{:});

% Normalize radii and density
dro = [dvec(1); diff(dvec)];
m = sum(dro.*zvec.^3);
robar = m/zvec(1)^3;
zvec = zvec/zvec(1);
dvec = dvec/robar;

%% Define and initialize local variables (in CMS notation)
lambdas = zvec;
deltas = [dvec(1); diff(dvec)];
nlay = length(lambdas);
nangles = opts.nangles;
kmax = opts.kmax;

% Define down-sampled variabels (for skip-n-spline)
if isscalar(opts.xlayers)
    if opts.xlayers > 0
        sskip = max(fix(nlay/opts.xlayers), 1);
        xind = 1:sskip:nlay;
    else
        xind = 1:nlay;
    end
else
    warning('CMS:xind','Experimental feature, use with care.')
    xind = opts.xlayers;
end
xlambdas = lambdas(xind);
xdvec = dvec(xind);
xdeltas = [xdvec(1); diff(xdvec)];
nxlay = length(xlambdas);

% Initialize spherical layer shapes
zetas = NaN(nlay,nangles);
xzetas = ones(nxlay, nangles);

% Initialize J-like quantities for spherical planet
if isempty(fieldnames(opts.J0s))
    Jlike = allocate_spherical_Js(nxlay, kmax, xlambdas, xdeltas);
else
    Jlike = opts.J0s;
end

% Abscissas and weights for Gaussian quadrature
[mus, gws] = gauleg(0, 1, nangles);

% Precompute Legendre polynomials for fixed colatitudes (for gauss quad)
Ps.Pnmu(kmax+1,nangles) = 0;
Ps.Pnzero(kmax+1,1) = 0;
for k=0:kmax
    Ps.Pnmu(k+1,1:nangles) = Pn(k, mus);
    Ps.Pnzero(k+1,1) = Pn(k, 0);
end

% Precompute powers of ratios of lambdas (only for explicit layers)
if opts.prerat
    lamratpow = nan(kmax+2,nxlay,nxlay);
    for ii=1:nxlay
        for jj=1:nxlay
            for kk=1:kmax+2
                lamratpow(kk,ii,jj) = ...
                    (xlambdas(ii)/xlambdas(jj))^(kk-1);
            end
        end
    end
else
    lamratpow = @(kk,ii,jj)(xlambdas(ii)./xlambdas(jj)).^(kk-1);
end

%% The loop (see Hubbard, 2012 and ./notes/CMS.pdf)
Js = Jlike.Jn;
for iter=1:opts.maxiter
    % Update shape with current gravity
    new_xzetas = update_zetas(Jlike, Ps, lamratpow, qrot, xzetas);
    for alfa = 1:nangles
        zetas(:,alfa) = spline(xlambdas, new_xzetas(:,alfa), lambdas);
    end
    
    % Update gravity with current shape
    new_Jlike = update_Js(lambdas, deltas, zetas, xind, Ps, gws);
    new_Js = new_Jlike.Jn;
    
    % Check for convergence of J0-J8 to terminate...
    dJs = abs(Js - new_Js)./abs(Js+eps);
    dJs(abs(new_Js) < 1e-12) = 0; % a hack for nonrotating planets
    if all(dJs(1:5) < opts.tol), break, end
    
    % ... or update to new values and continue
    Jlike = new_Jlike;
    Js = new_Js;
    xzetas = new_xzetas;
end

% It's not always a disaster if maxiter is reached, but we'd like to know
if iter == opts.maxiter
    warning('CMS:maxiter','Shape may not be fully converged.')
end

%% Return
Js = new_Js; % may as well use the latest...
out.dJs = dJs;
out.iter = iter;
out.zetas = zetas;
out.lambdas = lambdas;
out.deltas = deltas;
out.JLike = new_Jlike;
out.mus = mus;
out.gws = gws;
out.Ps = Ps;
out.xind = xind;

end

%% Helper functions
function print_usage()
fprintf('Usage:\n\tcms(zvec,dvec,qrot,''name'',value)\n')
fprintf('Name-Value pair arguments:\n')
fprintf('tol - Convergence tolerance for gravity coefficients [ positive real {1e-6} ]\n')
fprintf('maxiter - Number of iterations allowed for relaxation to equilibrium shape [ positive integer {60} ]\n')
fprintf('xlayers - Solve shape functions on xlayers and spline the rest [ integer scalar or vector {-1} ]\n')
fprintf('prerat - Precalculate powers of ratios of lambdas (trades memory for speed) [ {true} | false ]\n')
fprintf('J0s - J-like values representing initial state [ scalar struct {[]} ]\n')
end

function options = parsem(varargin)
p = inputParser;
p.FunctionName = 'cms.m';

p.addParameter('tol',1e-6,@(x)isscalar(x)&&isreal(x)&&x>0)
p.addParameter('maxiter',60,@(x)isscalar(x)&&isreal(x)&&x>0&&mod(x,1)==0)
p.addParameter('xlayers',-1,@(x)validateattributes(x,{'numeric'},{'vector','integer'}))
p.addParameter('prerat',true,@(x)isscalar(x)&&islogical(x))
p.addParameter('J0s',struct(),@(x)isscalar(x)&&isstruct(x))

% undocumented or obsolete options
p.addParameter('nangles',48,@(x)isscalar(x)&&(x>0)&&(mod(x,1)==0)) % colatitudes defining level surface
p.addParameter('kmax',30,@(x)isscalar(x)&&(x>6)&&(mod(x,2)==0)) % degree to cut mulitpole expansion
p.addParameter('TolX',1e-12,@(x)isscalar(x)&&(x>0)) % termination tolerance for root finding

% Parse name-value pairs and return
p.parse(varargin{:})
options = p.Results;
end

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
nlay = size(oldzetas, 1);
nangles = size(oldzetas, 2);
newzetas = NaN(nlay,nangles);
for j=1:nlay
    for alfa=1:nangles
        oldzeta = oldzetas(j,alfa);
        newzetas(j,alfa) = zeta_j_of_alfa(j, alfa, Js, Ps, lamrats, qrot, oldzeta);
    end
end
end

function newJs = update_Js(lambdas, deltas, zetas, xind, Ps, gws)
% Single-pass update of gravitational moments by Gaussian quad.

nlay = length(lambdas);
kmax = length(Ps.Pnzero)-1;
xlambdas = lambdas(xind);
dvec = cumsum(deltas);
xdvec = dvec(xind);
xdeltas = [xdvec(1); diff(xdvec)];
xzetas = zetas(xind, :);
nxlay = length(xlambdas);

% Do common denominator in eqs. (48) (USING FULL DENSITY PROFILE)
denom = 0;
for j=1:nlay
    fun = zetas(j,:).^3;
    I = gws*fun'; % gauss quad formula
    denom = denom + deltas(j)*lambdas(j)^3*I;
end

% Do J tilde, eq. (48a)
new_tilde = zeros(nxlay,kmax+1);
for ii=1:nxlay
    for kk=0:kmax
        if rem(kk, 2), continue, end
        fun = Ps.Pnmu(kk+1,:).*xzetas(ii,:).^(kk+3);
        I = gws*fun'; % gauss quad formula
        new_tilde(ii,kk+1) = -(3/(kk + 3))*xdeltas(ii)*xlambdas(ii)^3*I/denom;
    end
end

% Do J tilde prime, eqs. (48b and 48c)
new_tprime = zeros(nxlay,kmax+1);
for ii=1:nxlay
    for kk=0:kmax
        if rem(kk, 2), continue, end
        if kk == 2 % eq. (48c)
            fun = Ps.Pnmu(3,:).*log(xzetas(ii,:));
            I = gws*fun'; % gauss quad formula
            new_tprime(ii,kk+1) = -3*xdeltas(ii)*xlambdas(ii)^3*I/denom;
        else       % eq. (48b)
            fun = Ps.Pnmu(kk+1,:).*xzetas(ii,:).^(2 - kk);
            I = gws*fun'; % gauss quad formula
            new_tprime(ii,kk+1) = -(3/(2 - kk))*xdeltas(ii)*xlambdas(ii)^3*I/denom;
        end
    end
end

% Do J tilde double prime, eq. (48d)
new_tpprime = zeros(nxlay,1);
for ii=1:nxlay
    new_tpprime(ii) = 0.5*xdeltas(ii)*xlambdas(ii)^3/denom;
end

% And finally, the external Js deserve full grid resolution
full_tilde = zeros(nlay,kmax+1);
for ii=1:nlay
    for kk=0:kmax
        if rem(kk, 2), continue, end
        fun = Ps.Pnmu(kk+1,:).*zetas(ii,:).^(kk+3);
        I = gws*fun'; % gauss quad formula
        full_tilde(ii,kk+1) = -(3/(kk + 3))*deltas(ii)*lambdas(ii)^3*I/denom;
    end
end

% Return updated Js struct
newJs.tilde = new_tilde;
newJs.tildeprime = new_tprime;
newJs.tildeprimeprime = new_tpprime;
newJs.fulltilde = full_tilde;
n = 0:2:kmax;
for k=1:length(n)
    newJs.Jn(k) = dot(full_tilde(:,n(k)+1),lambdas.^n(k));
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
