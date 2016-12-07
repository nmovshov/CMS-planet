classdef eos_table < barotropes.Tabular
    %EOS_TABLE Example of a tabular barotrope that constructs from file.
    %   Use this class as a template for creating subclasses derived from
    %   barotropes.Tabular that initialize by reading a table from a file. All
    %   the required methods are already implemented in the base class Tabular.
    %   the derived subclasses only need to implement a constructor that takes a
    %   file name, reads the data in the appropriate format, potentially post
    %   processes missing values etc., and then assigns a vector of P values in
    %   Pascals and a vector of rho values in kg/m^3 to the P_vals and rho_vals
    %   properties derived from Tabular.

    %% Properties
    % Type 'properties(barotropes.Tabular) to see derived properties. You only
    % need to add properties that are specific to the custom eos. To save meta
    % data about the table (e.g. origin, composition, entropy...) you can add
    % any desired fields to the Tabular.meta property.
    %properties
    %end
    
    %% The constructor
    % This function reads in data from file and assigns to P_vals and rho_vals.
    % In this example the file contains a single header line followed by a
    % 3-column array so using dlmread is easiest. The columns are log(P),
    % log(rho), and log(T) in cgs units. We don't need T.
    methods
        function obj = eos_table(filename)
            if nargin == 0, return, end % matlab likes us to allow empty calls
            
            % Add input checks here if you feel like it
            narginchk(1,1)
            
            % Read in raw data
            raw = dlmread(filename,'',1,0);
            
            % Some post processing
            P = 10.^raw(:,1)*0.1; % 1st col is log(P) in dyne/cm^2
            rho = 10.^raw(:,2)*1000; % 2nd col is log(rho) in g/cm^3
            
            % Assign to base class derived properties
            obj.P_vals = P;
            obj.rho_vals = rho;
            
            % That's it but you can assign some meta data if desired
            obj.meta.original_units = 'cgs';
        end
    end
end
