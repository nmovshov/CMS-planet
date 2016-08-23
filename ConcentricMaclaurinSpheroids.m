classdef ConcentricMaclaurinSpheroids < handle
    %CONCENTRICMACLAURINSPHEROIDS Implementation of CMS shape model.
    %   This class implements the iterative relaxation of concentric Maclaurin
    %   spheroids from starting (dimensionless) radii and densities to a
    %   self-consistent hydrostatic equilibrium shape, as explained in Hubbard
    %   (2013).
    
    %% Properties
    properties (Access = public)
        opts    % holds CMS user configurable options
        lambdas % normalized layer equatorial radii
        deltas  % normalized density steps
    end
    properties (SetAccess = private)
        mus     % colatitude cosines
        zetas   % normalized and scaled level-surface radii
        Js      % rescaled, dimensionless, layer gravity moments
        bs      % normalized polar radii
        fs      % layer flattening (oblateness)
        ars     % layer aspect ratio
    end
    properties (Access = private)
        Pnmu    % values of Legendre polynomials at fixed colatitudes
        Pnzero  % values of Legendre polynomials at equator
        gws     % weight factors for Gauss integration (correspond to mus)
        cooked = false  % flag indicating successful convergence
    end
    
    %% The constructor
    methods
        function obj = ConcentricMaclaurinSpheroids(varargin)
            %CONCENTRICMACLAURINSPHEROIDS Class constructor.
            
            if nargin == 0
                % Use all default options
                op = cmsset();
            elseif (nargin == 1) && (isstruct(varargin{1}))
                % Use options struct returned by cmsset()
                op = varargin{1};
            else
                % Use Name/Value pairs
                op = cmsset(varargin{:});
            end
            
            % Pre-allocation and default assignments
            obj.lambdas = linspace(1, op.rcore, op.nlayers)';
            if (op.nlayers == 1), obj.lambdas = 1; end % N=1 special case
            obj.deltas = zeros(op.nlayers, 1);
            obj.deltas(1) = 1;
            obj.mus = linspace(0,1,op.nangles); % may be modified by gauss
            obj.zetas = ones(op.nlayers, op.nangles);
            obj.Js.tilde = zeros(op.nlayers, (op.kmax+1));
            obj.Js.tilde_prime = zeros(op.nlayers, (op.kmax+1));
            obj.Js.pprime = zeros(op.nlayers, 1);
            obj.Js.Jn = NaN(1, op.kmax+1);
            obj.opts = op;
            
            % Approximate degree 0 Js
            den = sum(obj.deltas.*obj.lambdas.^3);
            obj.Js.tilde(:,1) = -(obj.deltas.*obj.lambdas.^3)/den;
            obj.Js.tilde_prime(:,1) = -1.5*(obj.deltas.*obj.lambdas.^3)/den;
            obj.Js.pprime(:) = 0.5*obj.deltas/den;
            obj.Js.Jn(1) = sum(obj.Js.tilde(:,1));
            
            %TODO: setup better deltas
            
            % Get mus and weights for Gaussian quadrature
            if strcmpi(obj.opts.J_integration_method, 'gauss')
                [obj.mus, obj.gws] = gauleg(0, 1, obj.opts.nangles);
            end
            
            % Precompute Legendre polynomials for fixed colatitudes (gauss quad)
            for k = 0:obj.opts.kmax
                obj.Pnmu(k+1,:) = Pn(k, obj.mus);
                obj.Pnzero(k+1,1) = Pn(k, 0);
            end
        end % Constructor
    end % Constructor block
    
    %% Public methods
    methods (Access = public)
        function relax(obj)
            % Iterate calculation of gravitational moments until converged.
            
            % Optional communication
            t_rlx = tic;
            verb = obj.opts.verbosity;
            if (verb > 0)
                fprintf('Relaxing CMS to self-consistent level surfaces...\n\n')
            end
            
            % Main loop
            dJ = Inf;
            iter = 1;
            while (dJ > obj.opts.dJtol) && (iter <= obj.opts.MaxIter)
                t_pass = tic;
                if (verb > 0)
                    fprintf('Pass %d (of %d max)...\n', iter, obj.opts.MaxIter)
                end
                obj.update_zetas;
                dJ = obj.update_Js();
                if (verb > 0)
                    fprintf('Pass %d (of %d max)...done. (%g sec.)\n',...
                        iter, obj.opts.MaxIter, toc(t_pass))
                    if (verb > 1)
                        fprintf('dJ = %g; required tolerance = %g.\n\n',...
                            dJ, obj.opts.dJtol)
                    else
                        fprintf('\n')
                    end
                end
                iter = iter + 1;
            end
            if (iter < obj.opts.MaxIter), obj.cooked = true; end % converged!
            
            % Optional communication
            if (verb > 0)
                msg = 'Relaxing CMS to self-consistent level surfaces...done.';
                fprintf([msg, '\n'])
                fprintf('Total elapsed time %g sec.\n', toc(t_rlx))
            end
        end
        
        function update_zetas(obj)
            % Update level surfaces using current value of Js.

            % Optional communication
            verb = obj.opts.verbosity;
            if (verb > 0)
                t_z_pass = tic;
                fprintf('Updating zetas....')
            end
            
            % Loop over layers (outer) and colatitudes (inner)
            nlayers = size(obj.zetas, 1);
            nangles = size(obj.zetas, 2);
            for ii=1:nlayers
                for alfa=1:nangles
                    obj.zetas(ii,alfa) = obj.zeta_j_of_alfa(ii,alfa);
                    %obj.zetas(ii,alfa) = obj.zeta_j_of_mu(ii, obj.mus(alfa));
                end
            end
            
            % Optional communication
            if (verb > 1)
                t_z_pass = toc(t_z_pass);
                fprintf('done. (%g sec.)\n', t_z_pass)
            elseif (verb > 0)
                fprintf('done.\n')
            end
        end
        
        function dJ = update_Js(obj)
            % Single-pass update of gravitational moments (dispatcher).
            
            % Optional communication
            verb = obj.opts.verbosity;
            if (verb > 0)
                t_J_pass = tic;
                fprintf('Updating Js....')
            end
            
            % Dispatch based on chosen integration method.
            switch lower(obj.opts.J_integration_method)
                case 'gauss'
                    dJ = obj.update_Js_gauss();
                case 'adaptive'
                    dJ = obj.update_Js_adaptive();
                otherwise
                    error('Unknown integration method: %s.',...
                           obj.opts.J_integration_method)
            end
            obj.Js.Jn = obj.Jn();
            
            % Optional communication
            if (verb > 1)
                t_J_pass = toc(t_J_pass);
                fprintf('done. (%g sec).\n', t_J_pass)
            elseif (verb > 0)
                fprintf('done.\n')
            end
        end
        
        function TF = validate(obj)
            % Run a few simple sanity checks, warn if fail.
            
            TF = false;
            warning('off','backtrace')
            
            % Don't bother if obj is not even cooked
            if (~obj.cooked)
                warning(['Object may not be fully converged, try running',...
                    ' %s.relax()'],inputname(1))
                return
            end
            
            % Check for properly normalized zetas
            for j=1:obj.opts.nlayers
                er = abs(1 - obj.zeta_j_of_mu(j,0));
                if er > 1e-14
                    warning('Layer %d normalized radius deviates by %g',j,er)
                    return
                end
            end
            
            % Check for properly normalized Js
            er = obj.Jn(0) + 1;
            if abs(er) > 1e-15
                warning('J0 + 1 = %g (should be near zero)',er)
            end
            
            % Check for weak concentricity
            ind = find(diff(obj.bs) > 0);
            if any(ind)
                warning('Concentricity may be broken')
            end
            for k=1:length(ind)
                warning('Polar radius b_%d > b_%d',ind(k) + 1, ind(k))
            end
            
            % Check for strong concentricity
            ind = find((obj.bs(1:end-1) - obj.lambdas(2:end)) < 0);
            if any(ind)
                warning('Strong concentricity may be broken')
            end
            for k=1:length(ind)
                warning('b_%d < a_%d',ind(k), ind(k) + 1)
            end
            
            % TODO: more checks
            
            % If you can read this you passed.
            TF = true;
            warning('on','backtrace')
        end
        
        function J = Jn(obj,n)
            % Convenience method: return external, not re-scaled J multipoles
            if exist('n','var')
                validateattributes(n,{'numeric'},...
                                     {'nonnegative','even','<=',obj.opts.kmax})
            else
                n = 0:2:obj.opts.kmax;
            end
            for k=1:length(n)
                J(k) = dot(obj.Js.tilde(:,n(k)+1),obj.lambdas.^n(k));%#ok<AGROW>
            end
        end
        
        function ah = plot(obj)
            % Visualize a CMS object, return axes handle.
            
            mu = [1, fliplr(obj.mus), 0];
            th = acos(mu);            % 0 to pi/2
            th = [th, fliplr(pi-th)]; % 0 to pi
            th = [th, (pi + th)];     % 0 to 2pi
            
            try % requires R2016a or above
                % Prepare axes
                fh = figure;
                pax = polaraxes;
                pax.ThetaZeroLocation = 'top';
                pax.ThetaDir = 'clockwise';
                pax.ThetaAxisUnits = 'rad';
                cmap = parula(obj.opts.nlayers);
                hold(pax, 'on')
                
                % Add distinct outer surface
                xi0 = [obj.bs(1), fliplr(obj.zetas(1,:)), 1];
                xi0 = [xi0, fliplr(xi0)];
                xi0 = [xi0, fliplr(xi0)];
                polarplot(pax, th, xi0, '-k', 'linewidth', 2)
                
                % Return handle if requested
                if (nargout == 1), ah = pax; end
                
            catch
                warning('Fallback for R2016a and earlier not yet implemented')
                % Return handle if requested
                if (nargout == 1), ah = pax; end
            end
            
        end
    end % public methods
    
    %% Private methods
    methods (Access = private) % to become private
        function dJ = update_Js_gauss(obj)
            % Single-pass update of gravitational moments by Gaussian quad.
                        
            % Do common denominator in eqs. (40-43)
            denom = 0;
            for j=1:obj.opts.nlayers
                fun = obj.zetas(j,:).^3;
                I = obj.gws*fun'; % gauss quad formula
                denom = denom + obj.deltas(j)*obj.lambdas(j)^3*I;
            end
            
            % Do J tilde, eq. (40)
            new_tilde = zeros(size(obj.Js.tilde));
            for ii=1:obj.opts.nlayers
                for kk=0:obj.opts.kmax
                    if rem(kk, 2), continue, end
                    fun = obj.Pnmu(kk+1,:).*obj.zetas(ii,:).^(kk+3);
                    I = obj.gws*fun'; % gauss quad formula
                    new_tilde(ii,kk+1) = -(3/(kk + 3))*obj.deltas(ii)*obj.lambdas(ii)^3*I/denom;
                end
            end
            
            % Do J tilde prime, eqs. (41,42)
            new_tprime = zeros(size(obj.Js.tilde_prime));
            for ii=1:obj.opts.nlayers
                for kk=0:obj.opts.kmax
                    if rem(kk, 2), continue, end
                    if kk == 2
                        % eq. (42)
                        fun = obj.Pnmu(3,:).*log(obj.zetas(ii,:));
                        I = obj.gws*fun'; % gauss quad formula
                        new_tprime(ii,kk+1) = -3*obj.deltas(ii)*obj.lambdas(ii)^3*I/denom;
                    else
                        % eq. (41)
                        fun = obj.Pnmu(kk+1,:).*obj.zetas(ii,:).^(2 - kk);
                        I = obj.gws*fun'; % gauss quad formula
                        new_tprime(ii,kk+1) = -(3/(2 - kk))*obj.deltas(ii)*obj.lambdas(ii)^3*I/denom;
                    end
                end
            end
            
            % Do J double prime, eq. (27)
            new_pprime = zeros(size(obj.Js.pprime));
            denom = 0;
            for j=1:obj.opts.nlayers
                fun = obj.zetas(j,:).^3;
                I = obj.lambdas(j)*obj.gws*fun';
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
        
        function dJ = update_Js_adaptive(obj)
            % Single-pass update of gravitational moments using built-in integral.
            
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
                    new_tilde(ii,kk+1) = -(3/(kk + 3))*obj.deltas(ii)*obj.lambdas(ii)^3*I/denom;
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
        
        function y = zeta_j_of_alfa(obj,jlayer,alfa)
            % Find lvl surface of jth layer at colat mu(alfa).
            assert(jlayer > 0 && jlayer <= obj.opts.nlayers)
            assert(alfa > 0 && alfa <= obj.opts.nangles)
            if jlayer == 1
                fun = @(x)obj.eq50_alfa(x,alfa);
            else
                fun = @(x)obj.eq51_alfa(x,jlayer,alfa);
            end
           %y = fzero(fun, [0.6, 1.02]);
            y = fzero(fun, 1);
        end
        
        function y = eq50_alfa(obj,zeta0,alfa)
            % Equation 50 in Hubbard (2013) for fixed colatitudes.
            
            % Local variables
            nlayers = obj.opts.nlayers;
            kmax = obj.opts.kmax;
            Jt = obj.Js.tilde;
            lambda = obj.lambdas;
            qrot = obj.opts.qrot;
            P0 = obj.Pnzero;
            Pmu = obj.Pnmu(:,alfa);
            
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

        function y = eq51_alfa(obj,zeta_j,jj,alfa)
            % Equation 51 in Hubbard (2013) for fixed colatitudes.
            
            % Local variables
            nlayers = obj.opts.nlayers;
            kmax = obj.opts.kmax;
            Jt = obj.Js.tilde;
            Jtp = obj.Js.tilde_prime;
            Jpp = obj.Js.pprime;
            lambda = obj.lambdas;
            qrot = obj.opts.qrot;
            P0 = obj.Pnzero;
            Pmu = obj.Pnmu(:,alfa);
            
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
        
        function y = zeta_j_of_mus(obj,jlayer,mus)
            % ML quad requires y=f(x) take and return vector...
            assert(jlayer > 0 && jlayer <= obj.opts.nlayers)
            assert(all(mus >= 0) && all(mus <=1))
            y = nan(size(mus));
            for alfa=1:numel(mus)
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
    
    %% Access methods
    methods
        function set.opts(obj,val)
            % Certain parameters are not alowed to change after construction:
            forbiddenFields = {'kmax','nangles','nlayers'};
            if ~isempty(obj.opts) % for the call in the constructor
                for k=1:length(forbiddenFields)
                    if val.(forbiddenFields{k}) ~= obj.opts.(forbiddenFields{k})
                        msg = ['Changing %s in an existing obj makes no ',...
                               'sense; create a new CMS object instead.'];
                        error(msg,forbiddenFields{k})
                    end
                end
            end
            
            % For permitted opts, filter through cmsset again.
            obj.opts = cmsset(val);
            obj.cooked = false; %#ok<MCSUP>
        end
        
        function set.lambdas(obj,val)
            obj.lambdas = val;
            obj.cooked = false; %#ok<MCSUP>
        end
        
        function set.deltas(obj,val)
            obj.deltas = val;
            obj.cooked = false; %#ok<MCSUP>
        end
        
        function val = get.bs(obj)
            val = NaN(size(obj.lambdas));
            if obj.cooked
                for j=1:obj.opts.nlayers
                    val(j) = obj.zeta_j_of_mu(j,1)*obj.lambdas(j);
                end
            end
        end
        
        function val = get.fs(obj)
            val = NaN(size(obj.lambdas));
            if obj.cooked
                for j=1:obj.opts.nlayers
                    a = 1;
                    b = obj.zeta_j_of_mu(j,1);
                    val(j) = (a - b)/a;
                end
            end
        end
        
        function val = get.ars(obj)
            val = NaN(size(obj.lambdas));
            if obj.cooked
                for j=1:obj.opts.nlayers
                    a = 1;
                    b = obj.zeta_j_of_mu(j,1);
                    val(j) = b/a;
                end
            end
        end
    end % access methods

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

%% Input parsing and minimal assertions
narginchk(3,3)
nargoutchk(2,2)
validateattributes(x1,{'numeric'},{'scalar','finite','real'},1)
validateattributes(x2,{'numeric'},{'scalar','finite','real'},2)
validateattributes(n,{'numeric'},{'scalar','finite','integer','>=',2},3)
assert(x2 > x1, 'Interval must be positive.');

%% Local variables
tol = 1e-14;
m = ceil(n/2);
xmid = (x1 + x2)/2;
dx = (x2 - x1);
x = NaN(1,n);
w = NaN(1,n);

%% Main loop
for j=1:m
    % Get j-th root of Legenre polynomial Pn, along with Pn' value there.
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

%% Verify and return
assert(all(isfinite(x)))
assert(all(isfinite(w)))
end
