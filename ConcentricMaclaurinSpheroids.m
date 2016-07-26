classdef ConcentricMaclaurinSpheroids < handle
    %CONCENTRICMACLAURINSPHEROIDS Implementation of CMS shape model.
    %   This class implements the iterative relaxation of concentric Maclaurin
    %   spheroids from starting (dimensionless) radii and densities to a
    %   self-consistent hydrostatic equilibrium shape, as explained in Hubbard
    %   (2013).
    
    %% Properties
    properties (GetAccess = public, SetAccess = public)
        opts    % holds CMS project-wide options
        lambdas % normalized layer equatorial radii
        deltas  % normalized density steps
        mus     % colatitude cosines
        zetas   % normalized and scaled level-surface radii
        Js      % rescaled dimensionless gravity moments
        Pn      % coefficients of Legendre polynomials
    end
    
    %% The constructor
    methods
        function obj = ConcentricMaclaurinSpheroids(opts)
            if nargin == 0
                opts = cmsset();
            end
            
            % Pre-allocation and simple assignments
            obj.lambdas = linspace(1, opts.rcore, opts.nlayers)';
            obj.deltas = zeros(opts.nlayers, 1);
            obj.deltas(1) = 1;
            obj.mus = linspace(0,1,opts.nangles);
            obj.zetas = ones(opts.nlayers, opts.nangles);
            obj.Js.tilde = zeros(opts.nlayers, (opts.kmax+1));
            obj.Js.tilde_prime = zeros(opts.nlayers, (opts.kmax+1));
            obj.Js.pprime = zeros(opts.nlayers, 1);
            obj.opts = opts;
            
            % Approximate degree 0 Js
            den = sum(obj.deltas.*obj.lambdas.^3);
            obj.Js.tilde(:,1) = -(obj.deltas.*obj.lambdas.^3)/den;
            obj.Js.tilde_prime(:,1) = -1.5*(obj.deltas.*obj.lambdas.^3)/den;
            obj.Js.pprime(:) = 0.5*obj.deltas/den;
            
            %TODO: setup better deltas
            %TODO: setup better mus?
            
            % Precompute Legendre polynomial coefficients up to degree kmax
            
        end
    end
    
    %% Ordinary methods
    methods (Access = public)
        function update_zetas(obj)
            % Update level surfaces using current value of Js.
            pbar = (obj.opts.verbosity > 1);
            if pbar
                progressbar('updating zetas')
            end
            nbz = numel(obj.zetas);
            for ii=1:size(obj.zetas, 1)
                for alfa=1:size(obj.zetas, 2)
                    obj.zetas(ii,alfa) = obj.zeta_j_of_mu(ii, obj.mus(alfa));
                    if pbar
                        progressbar(((ii -1)*obj.opts.nangles + alfa)/nbz);
                    end
                end
            end
        end
        
        function dJ = update_Js(obj)
            % Single-pass update of gravitational moments (dispatcher).
            
            t_J_pass = tic;
            
            % Dispatch based on chosen integration method.
            switch lower(obj.opts.J_integration_method)
                case 'adaptive'
                    dJ = obj.update_Js_adaptive();
                otherwise
                    error('Unknown integration method: %s.',obj.opts.J_integration_method)
            end
            
            toc(t_J_pass)            
        end
        
    end % public methods
    
    methods (Access = public) % to become private
        function dJ = update_Js_adaptive(obj)
            % Single-pass update of gravitational moments using built-in integral.
            
            pbar = (obj.opts.verbosity > 1);
            if pbar, progressbar('updating Js'), end
            nbJ = numel(obj.Js.tilde) + numel(obj.Js.tilde_prime) + ...
                  numel(obj.Js.pprime);
            
            % Do common denominator in eqs. (40-43)
            denom = 0;
            for j=1:obj.opts.nlayers
                fun = @(mu)obj.zeta_j_of_mus(j, mu).^3;
                I = integral(fun, 0, 1, 'RelTol', obj.opts.IntTol);
                denom = denom + obj.deltas(j)*obj.lambdas(j)^3*I;
            end
            
            % Do J tilde, eq. (40)
            new_tilde = zeros(size(obj.Js.tilde));
            for ii=1:obj.opts.nlayers
                for kk=0:obj.opts.kmax
                    if rem(kk, 2), continue, end
                    fun = @(mu)Pn(kk,mu).*obj.zeta_j_of_mus(ii, mu).^(kk+3);
                    I = integral(fun, 0, 1, 'RelTol', obj.opts.IntTol);
                    new_tilde(ii,kk+1) = -(3/(2*kk + 3))*obj.deltas(ii)*obj.lambdas(ii)^3*I/denom;
                end
            end
            
            % Do J tilde prime, eqs. (41,42)
            new_tprime = zeros(size(obj.Js.tilde_prime));
            for ii=1:obj.opts.nlayers
                for kk=0:obj.opts.kmax
                    if rem(kk, 2), continue, end
                    if kk == 2
                        % eq. (42)
                        fun = @(mu)Pn(2,mu).*log(obj.zeta_j_of_mus(ii, mu));
                        I = integral(fun, 0, 1, 'RelTol', obj.opts.IntTol);
                        new_tprime(ii,kk+1) = -3*obj.deltas(ii)*obj.lambdas(ii)^3*I/denom;
                    else
                        % eq. (41)
                        fun = @(mu)Pn(kk,mu).*obj.zeta_j_of_mus(ii, mu).^(2 - kk);
                        I = integral(fun, 0, 1, 'RelTol', obj.opts.IntTol);
                        new_tprime(ii,kk+1) = -(3/(2 - kk))*obj.deltas(ii)*obj.lambdas(ii)^3*I/denom;
                    end
                end
            end
            
            % Do J double prime, eq. (27)
            new_pprime = zeros(size(obj.Js.pprime));
            denom = 0;
            for j=1:obj.opts.nlayers
                fun = @(mu)obj.zeta_j_of_mus(j, mu).^3;
                I = obj.lambdas(j)*integral(fun, 0, 1, 'RelTol', obj.opts.IntTol);
                denom = denom + obj.deltas(j)*I;
            end
            denom = 2*denom;
            for ii=1:obj.opts.nlayers
                new_pprime(ii) = obj.deltas(ii)/denom;
            end
            
            % Return max change in Js and update Js in obj.
            dJ = max(abs(new_tilde(:) - obj.Js.tilde(:)));
            dJ = max(dJ, max(abs(new_tprime(:) - obj.Js.tilde_prime(:))));
            dJ = max(dJ, max(abs(new_pprime(:) - obj.Js.pprime(:))));
            
            obj.Js.tilde = new_tilde;
            obj.Js.tilde_prime = new_tprime;
            obj.Js.pprime = new_pprime;
        end
        
        function y = zeta_j_of_mu(obj,jlayer,mu)
            % Find lvl surface of jth layer at colat mu.
            assert(jlayer > 0 && jlayer <= obj.opts.nlayers)
            assert(mu >= 0 && mu <=1)
            if jlayer == 1
                fun = @(x)eq50(x,mu,obj.Js.tilde,obj.lambdas,obj.opts.qrot);
            else
                fun = @(x)eq51(x,jlayer,mu,obj.Js.tilde,obj.Js.tilde_prime,obj.Js.pprime,obj.lambdas,obj.opts.qrot);
            end
           %y = fzero(fun, [0.6, 1.02]);
            y = fzero(fun, 1);
        end
        
        function y = zeta_j_of_mus(obj,jlayer,mus)
            % ML quad requires y=f(x) take and return vector...
            y = nan(size(mus));
            for alfa=1:length(mus)
                if jlayer == 1
                    fun = @(x)eq50(x,mus(alfa),obj.Js.tilde,obj.lambdas,obj.opts.qrot);
                else
                    fun = @(x)eq51(x,jlayer,mus(alfa),obj.Js.tilde,obj.Js.tilde_prime,obj.Js.pprime,obj.lambdas,obj.opts.qrot);
                end
                %y = fzero(fun, [0.6, 1.02]);
                y(alfa) = fzero(fun, 1);
            end
        end
    end % private methods
        
    %% Static methods
    methods (Static)
        
    end % End of static methods
    
end % End of classdef

%% Class-related functions
function y = eq50(zeta0,mu,Jt,lambda,qrot)
% Equation 50 in Hubbard (2013)

% summation limits
nlayers = length(lambda);
kmax = size(Jt, 2) - 1;

% Legendre polynomials at 0 and mu
P0(kmax+1) = 0;
Pmu(kmax+1) = 0;
for k=0:kmax
    P0(k+1) = Pn(k, 0);
    Pmu(k+1) = Pn(k, mu);
end

% Double sum in eq. (47)
x1 = 0;
for ii=1:nlayers
    for kk=2:kmax % (note ind shift, start ind, odd J=0)
        x1 = x1 + Jt(ii,kk+1)*lambda(ii)^kk*P0(kk+1);
    end
end

% Double sum in eq. (50)
x2 = 0;
for ii=1:nlayers
    for kk=2:kmax % (note ind shift, start ind, odd J=0)
        x2 = x2 + Jt(ii,kk+1)*lambda(ii)^kk*zeta0^(-kk)*Pmu(kk+1);
    end
end

% And combine
U0 = 1 + 0.5*qrot - x1;
U = (1/zeta0)*(1 - x2) + 1/3*qrot*zeta0^2*(1 - Pmu(3));
y = U - U0;
end

function y = eq51(zeta_j,jj,mu,Jt,Jtp,Jpp,lambda,qrot)
% Equation 51 in Hubbard (2013)

% summation limits
nlayers = length(lambda);
kmax = size(Jt, 2) - 1;

% Legendre polynomials at 0 and mu
P0(kmax+1) = 0;
Pmu(kmax+1) = 0;
for k=0:kmax
    P0(k+1) = Pn(k, 0);
    Pmu(k+1) = Pn(k, mu);
end

% Double sum, row 1
x1 = 0;
for ii=jj:nlayers
    for kk=0:kmax % (note ind shift, start ind, odd J=0)
        x1 = x1 + Jt(ii,kk+1)*(lambda(ii)/lambda(jj))^kk*zeta_j^(-kk)*Pmu(kk+1);
    end
end

% Double sum, row 2
x2 = 0;
for ii=1:jj-1
    for kk=0:kmax
        x2 = x2 + Jtp(ii,kk+1)*(lambda(jj)/lambda(ii))^(kk+1)*zeta_j^(kk+1)*Pmu(kk+1);
    end
end

% Single sum, row 2
x3 = 0;
for ii=1:jj-1
    x3 = x3 + Jpp(ii)*lambda(jj)^3*zeta_j^3;
end

% Double sum, row 3
x4 = 0;
for ii=jj:nlayers
    for kk=0:kmax
        x4 = x4 + Jt(ii,kk+1)*(lambda(ii)/lambda(jj))^kk*P0(kk+1);
    end
end

% Double sum, row 4
x5 = 0;
for ii=1:jj-1
    for kk=0:kmax
        x5 = x5 + Jtp(ii,kk+1)*(lambda(jj)/lambda(ii))^(kk+1)*P0(kk+1);
    end
end

% Single sum, row 4
x6 = 0;
for ii=1:jj-1
    x6 = x6 + Jpp(ii)*lambda(jj)^3;
end

% And combine
y = -(1/zeta_j)*(x1 + x2 + x3) + ...
    (1/3)*qrot*lambda(jj)^3*zeta_j^2*(1 - Pmu(3)) + ...
    (x4 + x5 + x6) - 0.5*qrot*lambda(jj)^3;

end

function y = Pn_builtin(n,x)
% Ordinary fully normalized Legendre polynomial using built-in legendre (slow).
assert(isvector(x))
Pnm = legendre(n,x);
y = Pnm(1,:);
if ~isrow(x), y = y'; end
end

function y = Pn(n,x)
% Ordinary fully normalized Legendre polynomial of order n at x in [0,1].
y = ones(size(x));

end
