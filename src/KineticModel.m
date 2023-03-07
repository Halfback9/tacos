classdef KineticModel < handle
    %
    % Dependency Optimisation Toolbox (fmincon)
    %
    % KineticModel Create a model for kinetics
    %   KineticModel(connections) Builds a matrix which represents the
    %   kinetic transfer of concentration between different compartments.
    %   
    % KineticModel Properties:
    %
    %
    %
    %   connections                         -   Matrix of connections between
    %                                           compartments
    %   connectionLocationsRows             -   Which row a connection was
    %                                           found in
    %   connectionLocationsColumns          -   Which column a connection
    %                                           was found in
    %
    %
    %
    %   nCompartments                       -   Number of compartments
    %   intialCompartmentConcentrations     -   Compartment concentrations
    %                                           at time zero
    %
    %
    %
    %   KMatrix                             -   Population transfer matrix
    %
    %
    % 
    %   nParameters                         -   Number of parameters
    %   parameterValues                     -   Current values of
    %                                           parameters
    %   parameterLabels                     -   Label of each parameters
    %
    %
    %   
    %   instrumentResponseFWHM              -   Full width half maximum of
    %                                           a Gaussian instrument
    %                                           response
    %
    %
    %
    %   timeZero                            -   Where the kinetics start
    %
    %
    %
    % KineticModel Methods:
    %
    %   generateConcentrations(times)       -   Use the kinetic model which
    %                                           has been set up to evaluate
    %                                           the concentrations of the
    %                                           compartments at each time
    %                                           point dicated by the
    %                                           argument vector "times".
    %
    %   updateParameters()                  -   Update the rate constant
    %                                           the IRF width and time zero
    %                                           point by parsing an ordered
    %                                           vector. This is for the
    %                                           fitting algorithm.
   

    properties (GetAccess = public, SetAccess = public)
        
        % Input
        connections(:, :) {mustBeNumeric}
        connectionLocationsRows (:, 1) {mustBeNumeric}
        connectionLocationsColumns (:, 1) {mustBeNumeric}

        % Compartments
        nCompartments (1,1) {mustBeInteger}
        intialCompartmentConcentrations (:, 1) {mustBeNumeric} % Column vector of initial concentrations

        % Population Transfer
        KMatrix (:,:) {mustBeNumeric}
        
        % Parameters
        nParameters (1, 1) {mustBeInteger}
        parameterValues(:, 1) {mustBeNumeric}
        parameterLabels (:, 1) %{mustBeText}

        % IRF
        instrumentResponseFWHM (1,1) {mustBeNumeric}
        
        % Times
        timeZero (1,1) {mustBeNumeric}
        

    end

    methods (Access = public)

        function obj = KineticModel(connections)
            % KineticModel Constructor
            
            obj.connections = connections;
            obj.nCompartments = length(connections);

            % Make default compartment concentrations
            obj.intialCompartmentConcentrations = [1; zeros(obj.nCompartments-1, 1)];

            obj.processConnectionsMatrix(connections);

        end

        function Concs = generateConcentrations(obj, times) 

            arguments
                obj KineticModel
                times (1, :) {mustBeNumeric}
            end

            times = times(:)';
            
            % Do an eigenvalue decomposition of the population transfer
            % matrix
            [V,Lambda]  = eig(obj.KMatrix);
            Lambda      = diag(Lambda);
            a           = V\obj.intialCompartmentConcentrations;

            % Exponential function convolved with a Gaussian IRF of given
            % FWHM (w) and t0 and using the eigenvalues from the
            % decomposition to evaluate the eigenvalues as "parallel
            % decays"
            stddev           = obj.instrumentResponseFWHM/(2*sqrt(2*log(2))); % Convert FWHM to standard deviation, since it gives simpler equations
            ExpVec_GC   = 0.5.*exp(Lambda.*(times-obj.timeZero + 0.5*Lambda*stddev.^2)).*(1+erf((times-obj.timeZero+Lambda.*stddev.^2)./(stddev.*sqrt(2))));
            
            
            % Use the eigenvectors of the K Matrix to convert the parallel
            % decays into the coupled decays.
            Concs       = (V*(a.*ExpVec_GC))';

        end

        function updateParameters(obj, parameters)
            obj.parameterValues = parameters(1:obj.nParameters, 1);
            obj.instrumentResponseFWHM = parameters(obj.nParameters + 1);
            obj.timeZero = parameters(obj.nParameters + 2);
        end

    end % methods (Access = public)

    methods (Access = private)

        function processConnectionsMatrix(obj, connections)
            [obj.connectionLocationsRows, obj.connectionLocationsColumns] = find(connections ~= 0);

            obj.nParameters = length(obj.connectionLocationsRows);
            obj.parameterValues= zeros(obj.nParameters, 1);
            obj.parameterLabels = strings(obj.nParameters, 1);

            for i = 1:obj.nParameters
                if obj.connectionLocationsRows(i) == obj.connectionLocationsColumns(i)
                    obj.parameterLabels(i) = strcat("k", num2str(obj.connectionLocationsColumns(i)), ">", "0");
                else
                    obj.parameterLabels(i) = strcat("k", num2str(obj.connectionLocationsColumns(i)), ">", num2str(obj.connectionLocationsRows(i)));
                end
            end
            obj.instrumentResponseFWHM = 0.1;
            obj.timeZero = 0;
            
        end

        function buildKMatrix(obj)
            temporaryKMatrix = zeros(obj.nCompartments); % Set initial matrix size

            % Populate with individual transfer rates
            for i = 1:obj.nParameters 
                temporaryKMatrix(obj.connectionLocationsRows(i), obj.connectionLocationsColumns(i)) = obj.parameterValues(i);
            end

            

            % Re express the diagonal elements as the negatve of the sum of
            % the rates in the corresponding column (i.e. the diagnonal is
            % the sum of all rateprocesses pulling intensity from that
            % compartment
            for i = 1:obj.nCompartments
                temporaryKMatrix(i, i) = -1*sum(temporaryKMatrix(:, i));
            end

            obj.KMatrix = temporaryKMatrix;

        end

    end
    methods % Setters

        function set.intialCompartmentConcentrations(obj, intialCompartmentConcentrations)
            % SET intialCompartmentConcentrations

            if sum(intialCompartmentConcentrations, "all") ~= 1
                warning('Sum of compartment concentrations is not 1');
            end

            obj.intialCompartmentConcentrations = intialCompartmentConcentrations(:)'; % Ensure it's a column vector
        
        end

        function set.parameterValues(obj, Parameters)
            % SET parameterValues

            if all(isfinite(Parameters))
                obj.parameterValues = Parameters(:, 1);
            else
                error('Non finite value given as parameter.');
            end
            
            obj.buildKMatrix(); % New parameter values so need to rebuild the K Matrix

        end

        function set.KMatrix(obj, KMatrix)
            % SET KMatrix

            if all(isfinite(KMatrix))
                obj.KMatrix = KMatrix;
            else
                error('Non finite value given as parameter.');
            end

        end

        function set.timeZero(obj, timeZero)
            % SET timeZero

            if isfinite(timeZero)
                obj.timeZero = timeZero;
            else
                error('The given value for time zero is not finite');
            end

        end

    end % Setters
 
end % classdef