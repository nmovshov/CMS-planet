classdef CMSPlanet < handle
    %CMSPLANET Interior model of rotating fluid planet.
    %   This class implements a model of a rotating fluid planet using the
    %   Concentric Maclaurin Spheroids technique. A CMSPlanet object is defined
    %   by a given mass, equatorial radius, rotation period, and barotrope.
    
    %% Properties
    properties (Access = public)
        name    % model name
        desc    % model description (one line)
        a0      % equatorial radius
        M       % total mass
    end
    properties (Dependent)
        ai      % layer equatorial radii
        bi      % layer polar radii
        si      % layer mean radii
        rhoi    % layer densities
        Pi      % layer surface pressure (top of layer)
    end
    properties
        cms     % a CMS object
        eos     % a barotrope object
    end
    properties (Dependent)
        opts    % struct with opts, both CMS and CMSPlanet
        qrot    % rotation parameter 
        nlayers % layers of constant density 
        Js      % external gravity moments (even, 0:2:kmax degree)
        NMoI    % normalized moment of inertia
        s0      % mean radius
        b0      % polar radius
        f0      % flattening, a.k.a, oblateness: (a - b)/a
        rho0    % reference density (uses equatorial radius)
        rho_s   % mean density (uses mean radius)
        M_calc  % mass from current state of cms
        beta    % so called mass renormalization factor (M/M_calc)
        P_c     % central pressure
        P_mid   % layer internal pressure (avg. of surface pressures)
    end
    properties (Access = private)
        
    end
    
    %% The constructor
    methods
        function obj = CMSPlanet(varargin)
            %CMSPLANET Class constructor.
            
            % Empty, integer, or name/value arguments.
            warning off CMS:obsolete
            if nargin == 0
                op = cmsset();
            else
                if isnumeric(varargin{1})
                    varargin = {'nlayers',varargin{1},varargin{2:end}};
                end
                op = cmsset(varargin{:});
            end
            warning on CMS:obsolete
            
            % Call InitPlanet to finish the construction
            obj.InitPlanet(op);
        end
    end % End of constructor block
    
    %% Public methods
    methods (Access = public)
        function ET = relax_to_HE(obj)
            ET = obj.cms.relax();
        end
        
        function ET = relax_to_barotrope(obj)
            % Iterate relaxation to HE and density updates until converged.
            
            if isempty(obj.eos)
                warning('Set valid barotrope first (obj.eos = <barotrope>)')
                return
            end
            
            t_rlx = tic;
            
            % Optional communication
            verb = obj.opts.verbosity;
            fprintf('Relaxing CMP to desired barotrope...\n\n')
            if (verb > 2)
                try
                    sbj = ['CMP.relax_to_barotrope() started on ',...
                        getenv('computername')];
                    sendmail(obj.opts.email,sbj)
                catch
                end
            end
            
            % Main loop
            dBar = Inf;
            iter = 1;
            while (abs(dBar) > obj.opts.dBtol) && (iter <= obj.opts.MaxIterBar)
                t_pass = tic;
                fprintf('Baropass %d (of max %d)...\n',...
                    iter, obj.opts.MaxIterBar)
                if (verb > 0), fprintf('\n'), end
                
                obj.relax_to_HE;
                
                if (verb > 0), fprintf('\n'), end
                dM = (obj.M_calc - obj.M)/obj.M;
                dro = obj.update_densities;
                dBar = var(dro);
                
                if (verb > 0), fprintf('\n'), end
                fprintf('Baropass %d (of max %d)...done. (%g sec.)\n',...
                    iter, obj.opts.MaxIterBar, toc(t_pass))
                if (verb > 0)
                    fprintf(['|drho| < %g; var(drho) = %g; dM = %g;'...
                        ' (required tolerance = %g).\n\n'],...
                        max(abs(double(dro))), var(double(dro)), double(dM),...
                        obj.opts.dBtol)
                else
                    fprintf('\n')
                end
                if (verb > 2)
                    try
                        sbj = ['CMP.relax_to_barotrope() on ',...
                            getenv('computername')];
                        msg{1} = sprintf(...
                            'Baropass %d (of max %d)...done. (%g sec.)',...
                            iter, obj.opts.MaxIterBar, toc(t_pass));
                        msg{2} = sprintf(...
                            'dBar = %g; required tolerance = %g.',...
                            double(dBar), obj.opts.dBtol);
                        sendmail(obj.opts.email,sbj,msg)
                    catch
                    end
                end
                iter = iter + 1;
            end
            
            % Flags and maybe warnings
            if (dBar > obj.opts.dBtol)
                msg = ['Planet may not have fully relaxed to desired eos.\n',...
                    'Try increasing the number of layers '...
                    'and/or convergence tolerance (%s.opts.dBtol) ',...
                    'and/or iteration limit (%s.opts.MaxIterBar).'];
                warning off backtrace
                warning(msg, inputname(1), inputname(1))
                warning on backtrace
            end
            
            ET = toc(t_rlx);
            
            % Optional communication
            fprintf('Relaxing CMP to desired barotrope...done.\n')
            try
                fprintf('Total elapsed time %s\n',lower(seconds2human(ET)))
            catch
                fprintf('Total elapsed time %g sec.\n', ET)
            end
            if (verb > 2)
                try
                    sbj = ['CMP.relax_to_barotrope() finished on ',...
                        getenv('computername')];
                    sendmail(obj.opts.email,sbj)
                catch
                end
            end
        end
        
        function dro = update_densities(obj)
            t_rho = tic;
            verb = obj.opts.verbosity;
            if (verb > 0)
                fprintf('  Updating layer densities...')
            end
            P = obj.Pi;
            newro = obj.eos.density((P(1:end-1) + P(2:end))/2);
            newro(end+1) = obj.eos.density((P(end) + obj.P_c)/2);
            dro = ((newro - obj.rhoi)./obj.rhoi);
            if (verb > 0)
                fprintf('done. (%g sec.)\n', toc(t_rho))
            end
            obj.rhoi = newro;
        end
        
        function ah = plot_equipotential_surfaces(obj)
            % Visualize a CMSPlanet object by plotting equipotential contours.
            
            % Require R2016a to use the amazing polarplot features
            if verLessThan('matlab','9')
                error('CMS plotting requires R2016a or later')
            end
            
            ah = obj.cms.plot();
            
            % Add a colorbar
            if ~isempty(obj.rhoi)
                ch = colorbar;
                ch.Label.String =...
                    sprintf('\\times %.0f kg/m^3', double(max(obj.rhoi)));
                ch.Label.FontSize = 10;
            end
            
            % Indicate equatorial radius
            if ~isempty(obj.a0)
                ah.RTickLabel{end} =...
                    sprintf('1.0\\times{}%g km',double(obj.a0)/1e3);
            end
        end
        
        function ah = plot_barotrope(obj)
            % Plot P(rho) of current model and of input barotrope.
            
            % Prepare the canvas
            fh = figure;
            set(fh, 'defaultTextInterpreter', 'latex')
            set(fh, 'defaultLegendInterpreter', 'latex')
            ah = axes;
            hold(ah, 'on')
            
            % Prepare the data
            x_cms = double(obj.rhoi);
            y_cms = double(obj.P_mid);
            if ~isempty(obj.eos)
                x_bar = logspace(-4, 4); % expected range in SI units
                y_bar = double(obj.eos.pressure(x_bar));
            end
            
            % Plot the lines (pressure in GPa)
            lh(1) = stairs(x_cms, y_cms/1e9);
            if ~isempty(obj.eos)
                lh(2) = line(x_bar, y_bar/1e9);
            end
            
            % Style and annotate lines
            lh(1).LineWidth = 2;
            lh(1).Color = [0.31, 0.31, 0.31];
            if isempty(obj.name)
                lh(1).DisplayName = 'CMS model';
            else
                lh(1).DisplayName = obj.name;
            end
            if ~isempty(obj.eos)
                lh(2).Color = 'r';
                if isempty(obj.eos.name)
                    lh(2).DisplayName = 'input barotrope';
                else
                    lh(2).DisplayName = obj.eos.name;
                end
            end
            
            % Style and annotate axes
            ah.Box = 'on';
            ah.XScale = 'log';
            ah.YScale = 'log';
            xlim([min(x_cms),max(x_cms)])
            xlabel('$\rho$ [kg/m$^3$]')
            ylabel('$P$ [GPa]')
            
            % Legend
            gh = legend('show','location','nw');
            gh.FontSize = 11;
            
        end
    end % End of public methods block
    
    %% Private methods
    methods (Access = private)
        function InitPlanet(obj,op)
           % (Re)Initialize a CMSPlanet object.
           
           obj.cms = ConcentricMaclaurinSpheroids(op);
        end
    end % End of private methods block
    
    %% Access methods
    methods
        function val = get.opts(obj)
            val = obj.cms.opts;
        end
        
        function set.opts(obj,val)
            obj.cms.opts = val;
        end
        
        function val = get.nlayers(obj)
            val = obj.cms.nlayers;
        end
        
        function set.nlayers(~,~)
            msg = ['Changing number of layers in an existing CMSPlanet ',...
                'makes no sense; create a new object instead.'];
            error(msg)
        end
        
        function val = get.qrot(obj)
            val = obj.cms.qrot;
        end
        
        function set.qrot(obj,val)
            try
                (val == 1); %#ok<EQEFF>
                obj.cms.qrot = double(val);
            catch ME
                if (strfind(ME.message,'Attempt to compare')) == 1
                    error('Check dimensions of qrot')
                else
                    rethrow(ME)
                end
            end
        end
        
        function val = get.Js(obj)
            val = obj.cms.Jn;
        end
        
        function set.M(obj,val)
            try
                assert(isnumeric(val))
                assert(isscalar(val))
                assert(double(val) > 0)
                obj.M = val;
            catch
                error('Mass must be a positive scalar.')
            end
        end
        
        function set.a0(obj,val)
            try
                assert(isnumeric(val))
                assert(isscalar(val))
                assert(double(val) > 0)
                obj.a0 = val;
            catch
                error('Radius must be a positive scalar.')
            end
        end
        
        function val = get.ai(obj)
            val = [];
            if ~isempty(obj.a0)
                val = obj.a0*obj.cms.lambdas;
            end
        end
        
        function set.ai(obj,val)
            validateattributes(val,{'numeric'},{'real','finite','vector'},...
                '','ai')
            assert(numel(val) == obj.nlayers,...
                'length(ai) = %g ~= nlayers = %g',...
                numel(val),obj.nlayers)
            assert(all(double(val) > 0), 'layer radii must be positive')
            obj.a0 = val(1);
            obj.cms.lambdas = val/obj.a0;
        end
        
        function val = get.rhoi(obj)
            val = [];
            if ~isempty(obj.rho0)
                val = cumsum(obj.cms.deltas);
                val = val*obj.rho0;
            end
        end
        
        function set.rhoi(obj,val)
            validateattributes(val,{'numeric'},{'real','finite','vector'},...
                '','rhoi')
            assert(numel(val) == obj.nlayers,...
                'length(rhoi) = %g ~= nlayers = %g',...
                numel(val),obj.nlayers)
            assert(all(double(val) >= 0), 'layer densities must be nonnegative')
            assert(~isempty(obj.a0) && ~isempty(obj.M),...
                'cannot set layer densities before setting mass and radius')
            obj.cms.deltas = [val(1); diff(val(:))]/obj.rho0;
        end
        
        function set.eos(obj,val)
            if ~isa(val,'barotropes.Barotrope')
                error('eos must be a valid instance of class Barotrope')
            end
            obj.eos = val;
        end
        
        function val = get.s0(obj)
            val = obj.a0*obj.cms.ss(1);
        end
        
        function val = get.b0(obj)
            val = obj.a0*obj.cms.bs(1);
        end
        
        function val = get.f0(obj)
            val = (obj.a0 - obj.b0)/obj.a0;
        end
        
        function val = get.rho0(obj)
            if isempty(obj.M) || isempty(obj.a0)
                val = [];
            else
                val = obj.M/(4*pi/3*obj.a0^3);
            end
        end

        function val = get.rho_s(obj)
            if isempty(obj.M) || isempty(obj.s0)
                val = [];
            else
                val = obj.M/(4*pi/3*obj.s0^3);
            end
        end

        function val = get.M_calc(obj)
            if isempty(obj.rhoi), val = []; return, end
            drho = [obj.rhoi(1); diff(obj.rhoi)];
            val = (4*pi/3)*(obj.a0^3)*sum(drho.*obj.cms.Vs);
        end
        
        function val = get.beta(obj)
            val = obj.M/obj.M_calc;
        end
        
        function val = get.bi(obj)
            val = [];
            if ~isempty(obj.a0)
                val = obj.a0*obj.cms.bs;
            end
        end
        
        function val = get.si(obj)
            val = [];
            if ~isempty(obj.a0)
                val = obj.a0*obj.cms.ss;
            end
        end
        
        function val = get.Pi(obj)
            %#ok<*PROP> 
            if isempty(obj.rho0), val = []; return, end
            try
                si = setUnits; % if you have physunits
            catch
                si = setFUnits; % if you don't have physunits
            end
            G = double(si.gravity);
            if isequal(obj.opts.equipotential_squeeze, 'polar')
                U = double(obj.cms.equiUpu*G*obj.M/obj.a0);
            else
                U = double(mean(obj.cms.Upu, 2)*G*obj.M/obj.a0);
            end
            rho = double(obj.rhoi);
            val = zeros(obj.nlayers, 1);
            val(2:end) = cumsum(rho(1:end-1).*diff(U));
            val = val*si.Pa;
        end
        
        function val = get.P_c(obj)
            if isempty(obj.rho0), val = []; return, end
            try
                si = setUnits; % if you have physunits
            catch
                si = setFUnits; % if you don't have physunits
            end
            G = si.gravity;
            if isequal(obj.opts.equipotential_squeeze, 'polar')
                U = obj.cms.equiUpu*G*obj.M/obj.a0;
            else
                U = mean(obj.cms.Upu, 2)*G*obj.M/obj.a0;
            end
            U_center = -G*obj.M/obj.a0*...
                sum(obj.cms.Js.tilde_prime(:,1).*obj.cms.lambdas.^-1);
            P = obj.Pi; % matlab answers 307052
            rho = obj.rhoi;
            val = P(end) + rho(end)*(U_center - U(end));
        end
        
        function val = get.P_mid(obj)
            P = obj.Pi;
            val = (P(1:end-1) + P(2:end))/2;
            val(end+1) = (P(end) + obj.P_c)/2;
        end
        
        function val = get.NMoI(obj)
            val = obj.cms.NMoI;
        end
        
        function set.P_c(~,~)
            error('Pi is a calculated property and cannot be assigned')
        end
        
        function set.Pi(~,~)
            error('Pi is a calculated property and cannot be assigned')
        end
    end % End of access methods block
    
    %% Static methods
    methods (Static)
        
    end % End of static methods block
end % End of classdef
