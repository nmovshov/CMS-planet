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
        Mi      % layer masses
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
        J2      % convenience name for Js(2)
        J4      % convenience name for Js(3)
        J6      % convenience name for Js(4)
        NMoI    % normalized moment of inertia
        s0      % mean radius
        b0      % polar radius
        f0      % flattening, a.k.a, oblateness: (a - b)/a
        rho0    % reference density (uses equatorial radius)
        rho_s   % mean density (uses mean radius)
        M_calc  % mass from current state of cms
        P_c     % central pressure
        P_mid   % layer internal pressure (avg. of surface pressures)
    end
    properties (SetAccess = private)
        betanorm % mass renormalization factor (M/M_calc)
    end
    
    %% The constructor
    methods
        function obj = CMSPlanet(nlay, varargin)
            %CMSPLANET Class constructor.
            
            % The constructor only dispatches to InitCMS().
            if nargin == 0
                error(['Required argument missing: specify number of layers',...
                    ' as first input argument.'])
            end
            validateattributes(nlay,...
                {'numeric'},{'positive','integer','scalar'},...
                '','nlayers',1)
            warning off CMS:obsolete
            op = cmsset(varargin{:});
            obj.opts = op; % (calls set.opts which calls cmsset(op) again)
            warning on CMS:obsolete
            
            % Call InitPlanet to finish the construction
            obj.InitPlanet(nlay, op);
        end
    end % End of constructor block
    
    %% Public methods
    methods (Access = public)
        function ET = relax_to_HE(obj)
            ET = obj.cms.relax();
        end
        
        function ET = relax_to_barotrope(obj)
            % Relax shape function, Js, and density simultaneously.
            
            if isempty(obj.eos)
                warning('Set valid barotrope first (obj.eos = <barotrope>)')
                return
            end
            
            if isempty(obj.rho0)
                warning('Make sure mass (%s.M) and radius (%s.a0) are set',...
                    inputname(1), inputname(1))
                return
            end
            
            t_rlx = tic;
            
            % Optional communication
            verb = obj.opts.verbosity;
            if (verb > 0)
                fprintf('Relaxing CMP to desired barotrope...\n\n')
            end
            if (verb > 3)
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
                if (verb > 0)
                    fprintf('Baropass %d (of max %d)...\n',...
                        iter, obj.opts.MaxIterBar)
                end
                obj.cms.update_zetas;
                if isequal(obj.opts.equipotential_squeeze, 'polar')
                    obj.cms.update_polar_radii;
                end
                dJ = obj.cms.update_Js;
                if (verb > 1), fprintf('  '), end
                dro = obj.update_densities;
                
                % The convergence tolerance is the max of dJs and drhos
                vdro = var(dro(~isnan(dro)));
                dBar = max([vdro, dJ]);
                
                if (verb > 0)
                    fprintf('Baropass %d (of max %d)...done. (%g sec)\n',...
                        iter, obj.opts.MaxIterBar, toc(t_pass))
                end
                if (verb > 1)
                    fprintf(['var(drho) = %g; dJ = %g;'...
                        ' (required tolerance = %g).\n\n'],...
                        double(vdro), double(dJ),obj.opts.dBtol)
                elseif (verb > 0)
                    fprintf('\n')
                end
                if (verb > 3)
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
            
            % Renormalize densities; save betanorm constant
            obj.betanorm = obj.M/obj.M_calc;
            obj.rhoi = obj.rhoi*obj.betanorm;
            
            % Update polar radii if we haven't already
            if ~isequal(obj.opts.equipotential_squeeze, 'polar')
                obj.cms.update_polar_radii;
                if (verb > 1), fprintf('\n'), end
            end
            
            % Flags and maybe warnings
            if (dBar > obj.opts.dBtol)
                msg = ['Planet may not have fully relaxed to desired eos.\n',...
                    'Try increasing the number of layers '...
                    'and/or convergence tolerance (%s.opts.dBtol) ',...
                    'and/or iteration limit (%s.opts.MaxIterBar).'];
                warning off backtrace
                warning('CMS:noconverge', msg, inputname(1), inputname(1))
                warning on backtrace
                if (verb > 1), fprintf('\n'), end
            end
            
            ET = toc(t_rlx);
            
            % Optional communication
            if (verb > 0)
                fprintf('Relaxing CMP to desired barotrope...done.\n')
                try
                    fprintf('Total elapsed time %s\n',lower(seconds2human(ET)))
                catch
                    fprintf('Total elapsed time %g sec.\n', ET)
                end
            end
            if (verb > 3)
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
            if (verb > 1)
                fprintf('  Updating layer densities...')
            end
            P = obj.P_mid;
            if isscalar(obj.eos)
                newro = obj.eos.density(P);
            else
                newro = repmat(obj.rho0, obj.nlayers, 1);
                for k=1:length(newro)
                    newro(k) = obj.eos(k).density(P(k));
                end
            end
            dro = ((newro - obj.rhoi)./obj.rhoi);
            if (verb > 2)
                fprintf('done. (%g sec)\n', toc(t_rho))
            elseif (verb > 1)
                fprintf('done.\n')
            end
            obj.rhoi = newro;
        end
        
        function ah = plot_equipotential_surfaces(obj)
            % Visualize a CMSPlanet object by plotting equipotential contours.
            
            % Require R2016a to use the amazing polarplot features
            if verLessThan('matlab','9')
                error('CMS equipotential plots require R2016a or later')
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
        
        function ah = plot_barotrope(obj, varargin)
            % Plot P(rho) of current model and optionally of input barotrope.
            
            % Don't bother if there is no pressure
            if isempty(obj.Pi)
                warning('Uninitialized object.')
                return
            end
            
            % Input parsing
            p = inputParser;
            p.FunctionName = mfilename;
            p.addParameter('axes', [],...
                @(x)isscalar(x) && isgraphics(x,'axes') && isvalid(x));
            p.addParameter('showinput', false,...
                @(x)isscalar(x) && islogical(x));
            p.addParameter('includecore', false,...
                @(x)isscalar(x) && islogical(x));
            p.addParameter('betanormalize', false,...
                @(x)isscalar(x) && islogical(x));
            p.parse(varargin{:})
            pr = p.Results;
            
            % Prepare the canvas
            if isempty(pr.axes)
                fh = figure;
                set(fh, 'defaultTextInterpreter', 'latex')
                set(fh, 'defaultLegendInterpreter', 'latex')
                ah = axes;
                hold(ah, 'on')
            else
                ah = pr.axes;
                hold(ah, 'on')
            end
            
            % Prepare the data
            x_cms = double(obj.rhoi);
            y_cms = double(obj.P_mid);
            if pr.betanormalize
                x_cms = double(obj.beta*x_cms);
                y_cms = double(obj.beta*y_cms);
            end
            if ~pr.includecore && ~isempty(obj.eos) && ...
                    isa(obj.eos(end), 'barotropes.ConstDensity')
                x_cms(end) = [];
                y_cms(end) = [];
            end
            
            if pr.showinput && ~isempty(obj.eos) && (range(x_cms) > 0)
                x_bar = linspace(min(x_cms), max(x_cms));
                if isscalar(obj.eos)
                    y_bar = double(obj.eos.pressure(x_bar));
                else
                    v = 1:length(unique(x_cms));
                    ind = interp1(unique(x_cms), v, x_bar, 'nearest', 'extrap');
                    y_bar = nan(size(x_bar));
                    for k=1:length(x_bar)
                        y_bar(k) = double(obj.eos(ind(k)).pressure(x_bar(k)));
                    end
                end
            end
            
            % Plot the lines (pressure in GPa)
            lh(1) = stairs(x_cms, y_cms/1e9);
            if exist('y_bar','var') && any(isfinite(y_bar))
                lh(2) = line(x_bar, y_bar/1e9);
            end
            
            % Style and annotate lines
            lh(1).LineWidth = 2;
            if isempty(pr.axes)
                lh(1).Color = [0.31, 0.31, 0.31];
            end
            if isempty(obj.name)
                lh(1).DisplayName = 'CMS model';
            else
                lh(1).DisplayName = obj.name;
            end
            if length(lh) == 2
                lh(2).Color = 'r';
                lh(2).DisplayName = 'input barotrope';
            end
            
            % Style and annotate axes
            if isempty(pr.axes)
                ah.Box = 'on';
                if (range(x_cms) > 0)
                    xlim([min(x_cms),max(x_cms)])
                end
                xlabel('$\rho$ [kg/m$^3$]')
                ylabel('$P$ [GPa]')
            else
                xlim('auto')
            end
            
            % Legend
            legend(ah, 'off')
            gh = legend(ah, 'show','location','nw');
            gh.FontSize = 11;
            
        end
        
        function T = report_card(obj, obs)
            % REPORT_CARD Table summary of model's vital statistics.
            
            % Minimal checks
            narginchk(1,2);
            
            % Basic table
            if isempty(obj.M_calc)
                objM = NaN;
            else
                objM = obj.M_calc;
            end
            vitals = {'Mass [kg]'; 'J2'; 'J4'; 'J6'; 'J8'; 'J10'; 'NMoI'};
            CMP1 = [objM; obj.J2; obj.J4; obj.J6;...
                obj.Js(5); obj.Js(6); obj.NMoI];
            CMP1 = double(CMP1);
            T = table(CMP1, 'RowNames', vitals);
            if ~isempty(obj.name)
                vname = matlab.lang.makeValidName(obj.name);
                T.Properties.VariableNames{'CMP1'} = vname;
            end
            if nargin == 1, return, end
            
            % Optionally compare with another CMSPlanet
            if isa(obs, 'CMSPlanet')
                if isempty(obs.M_calc)
                    obsM = NaN;
                else
                    obsM = obs.M_calc;
                end
                CMP2 = [obsM; obs.J2; obs.J4; obs.J6;...
                    obs.Js(5); obs.Js(6); obs.NMoI];
                CMP2 = double(CMP2);
                T = [T table(CMP2)];
                if ~isempty(obs.name)
                    vname = matlab.lang.makeValidName(obs.name);
                    try
                        T.Properties.VariableNames{'CMP2'} = vname;
                    catch
                        T.Properties.VariableNames{'CMP2'} = ['x_',vname];
                    end
                end
                DIFF = (CMP1 - CMP2)./CMP1;
                T = [T table(DIFF, 'VariableNames', {'fractional_diff'})];
                return
            end
            
            % Or compare with observables set
            try
                if ~isfield(obs, 'M'), obs.M = NaN; end
                if ~isfield(obs, 'NMoI'), obs.NMoI = NaN; end
                if ~isfield(obs, 'J8'), obs.J8 = NaN; end
                if ~isfield(obs, 'J10'), obs.J10 = NaN; end
                OBS = [obs.M; obs.J2; obs.J4; obs.J6;...
                    obs.J8; obs.J10; obs.NMoI];
                OBS = double(OBS);
                T = [T table(OBS)];
                if isfield(obs, 'name') && ~isempty(obs.name)
                    vname = matlab.lang.makeValidName(obs.name);
                    try
                    T.Properties.VariableNames{'OBS'} = vname;
                    catch
                        T.Properties.VariableNames{'OBS'} = ['x_',vname];
                    end
                end
                DIFF = CMP1 - OBS;
                T = [T table(DIFF, 'VariableNames', {'diff'})];
                if ~isfield(obs, 'dM'), obs.dM = NaN; end
                if ~isfield(obs, 'dNMoI'), obs.dNMoI = NaN; end
                if ~isfield(obs, 'dJ8'), obs.dJ8 = NaN; end
                if ~isfield(obs, 'dJ10'), obs.dJ10 = NaN; end
                dees = [obs.dM; obs.dJ2; obs.dJ4; obs.dJ6;...
                    obs.dJ8; obs.dJ10; obs.dNMoI];
                WE = T.diff./dees;
                T = [T table(WE, 'VariableNames', {'weighted_error'})];
                match = ((T.weighted_error < 1) & (T.weighted_error > -1)) | ...
                    isnan(T.weighted_error);
                T = [T table(match)];
            catch ME
                if any(strcmp(ME.identifier, {'MATLAB:nonStrucReference',...
                        'MATLAB:nonExistentField',...
                        'MATLAB:structRefFromNonStruct'}))
                    msg = ['To compare model to observations supply a',...
                        'struct argument with fields:', ...
                        '\n\tM\n\tdM\n\tJ2\n\tdJ2\n\tJ4\n\tdJ4\n\tJ6\n\tdJ6',...
                        '\nand optionally',...
                        '\n\tname'];
                    error(sprintf(msg)) %#ok<SPERR>
                else
                    rethrow(ME)
                end
            end
        end
        
        function beta = renormalize_density(obj)
            beta = obj.M/obj.M_calc;
            obj.betanorm = beta;
            obj.rhoi = obj.rhoi*obj.betanorm;
        end
    end % End of public methods block
    
    %% Private methods
    methods (Access = private)
        function InitPlanet(obj,nlay,op)
           % (Re)Initialize a CMSPlanet object.
           
           obj.cms = ConcentricMaclaurinSpheroids(nlay, op);
        end
        
        function ET = relax_to_barotrope_old(obj)
            % Iterate relaxation to HE and density updates until converged.
            
            if isempty(obj.eos)
                warning('Set valid barotrope first (obj.eos = <barotrope>)')
                return
            end
            
            t_rlx = tic;
            
            % Optional communication
            verb = obj.opts.verbosity;
            fprintf('Relaxing CMP to desired barotrope...\n\n')
            if (verb > 3)
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
                if (verb > 1), fprintf('\n'), end
                
                obj.relax_to_HE;
                
                if (verb > 1), fprintf('\n'), end
                dM = (obj.M_calc - obj.M)/obj.M;
                dro = obj.update_densities;
                dBar = var(dro(~isnan(dro)));
                
                if (verb > 1), fprintf('\n'), end
                fprintf('Baropass %d (of max %d)...done. (%g sec.)\n',...
                    iter, obj.opts.MaxIterBar, toc(t_pass))
                if (verb > 1)
                    fprintf(['|drho| < %g; var(drho) = %g; dM = %g;'...
                        ' (required tolerance = %g).\n\n'],...
                        max(abs(double(dro))), double(dBar), double(dM),...
                        obj.opts.dBtol)
                else
                    fprintf('\n')
                end
                if (verb > 3)
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
                warning('CMS:noconverge', msg, inputname(1), inputname(1))
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
            if (verb > 3)
                try
                    sbj = ['CMP.relax_to_barotrope() finished on ',...
                        getenv('computername')];
                    sendmail(obj.opts.email,sbj)
                catch
                end
            end
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
        
        function val = get.J2(obj)
            val = obj.Js(2);
        end
        
        function val = get.J4(obj)
            val = obj.Js(3);
        end
        
        function val = get.J6(obj)
            val = obj.Js(4);
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
        
        function set.name(obj,val)
            if ~isempty(val)
                validateattributes(val, {'char'}, {'row'})
            end
            obj.name = val;
        end
        
        function set.desc(obj,val)
            if ~isempty(val)
                validateattributes(val, {'char'}, {'row'})
            end
            obj.desc = val;
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
                'length(rhoi) == %g ~= nlayers == %g',...
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
            n = obj.nlayers; %#ok<MCSUP>
            assert(isscalar(val) || (numel(val) == n),...
                'length(eos) == %g ~= nlayers == %g',...
                numel(val),n)
            obj.eos = val(:);
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
        
        function val = get.Mi(obj)
            if isempty(obj.rhoi), val = []; return, end
            dvs = [(obj.cms.Vs(1:end-1) - obj.cms.Vs(2:end)); obj.cms.Vs(end)];
            val = (4*pi/3)*(obj.a0^3)*(dvs.*obj.rhoi);
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
                % are we using doubles...
                obj.rho0 > 1; %#ok<VUNUS>
                si = setFUnits;
                G = si.gravity;
            catch
                % ... or preals
                si = setUnits;
                G = si.gravity;
            end
            if strcmp(obj.opts.equipotential_squeeze, 'mean')
                U = mean(obj.cms.Upu, 2)*G*obj.M/obj.a0;
            else
                U = obj.cms.equiUpu*G*obj.M/obj.a0;
            end
            rho = obj.rhoi;
            val = zeros(obj.nlayers, 1)*rho(1)*U(1);
            val(2:end) = cumsum(rho(1:end-1).*diff(U));
        end
        
        function val = get.P_c(obj)
            if isempty(obj.rho0), val = []; return, end
            try
                % are we using doubles...
                obj.rho0 > 1; %#ok<VUNUS>
                si = setFUnits;
                G = si.gravity;
            catch
                % ... or preals
                si = setUnits;
                G = si.gravity;
            end
            if strcmp(obj.opts.equipotential_squeeze, 'mean')
                U = mean(obj.cms.Upu, 2)*G*obj.M/obj.a0;
            else
                U = obj.cms.equiUpu*G*obj.M/obj.a0;
            end
            U_center = -G*obj.M/obj.a0*...
                sum(obj.cms.Js.tilde_prime(:,1).*obj.cms.lambdas.^-1);
            P = obj.Pi; % matlab answers 307052
            rho = obj.rhoi;
            val = P(end) + rho(end)*(U_center - U(end));
        end
        
        function val = get.P_mid(obj)
            if isempty(obj.rho0), val = []; return, end
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
        
        function set.Mi(~,~)
            error('Mi is a calculated property and cannot be assigned')
        end
    end % End of access methods block
    
    %% Static methods
    methods (Static)
        
    end % End of static methods block
end % End of classdef
