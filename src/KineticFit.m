classdef KineticFit < handle
    % KineticFit Unified approach to fitting time resolved spectra
    %   Provides a unified interface for independent or simultaneous
    %   fitting of fluorescence lifetime traces and broadband transient
    %   spectra. 
    %   Depends upon the kineticModel class to implement 

    properties

        % General fitting parameters
        kineticsRateConstantLabels (:, 1)
        kineticsRateConstantStartingValues (:, 1) {mustBeNumeric}
        kineticsRateConstantLowerBounds (:, 1) {mustBeNumeric}
        kineticsRateConstantUpperBounds (:, 1) {mustBeNumeric}
        kineticsRateConstantFittedValues (:, 1) {mustBeNumeric}
        kineticsRateConstantFittedVariances (:, 1) {mustBeNumeric}


        % TA Fit Parameters
        TAInstrumentResponseLabel {mustBeText} = "TA IRF FWHM"; 
        TAInstrumentResponseFWHMStartingValue (1, 1) {mustBeNumeric}
        TAInstrumentResponseLowerBound (1, 1) {mustBeNumeric}
        TAInstrumentResponseUpperBound (1, 1) {mustBeNumeric}
        TAInstrumentResponseFWHMFittedValue (1, 1) {mustBeNumeric}
        TAInstrumentResponseFWHMFittedVariance (1, 1) {mustBeNumeric}
            
        TATimeZeroLabel {mustBeText} = "TA Time Zero";
        TATimeZeroStartingValue (1, 1) {mustBeNumeric}
        TATimeZeroLowerBound (1, 1) {mustBeNumeric}
        TATimeZeroUpperBound (1, 1) {mustBeNumeric}
        TATimeZeroFittedValue (1, 1) {mustBeNumeric}
        TATimeZeroFittedVariance (1, 1) {mustBeNumeric}


        % Photoluminescence Lifetime Fit Parameters
        PLInstrumentResponseLabel {mustBeText} = "PL IRF FWHM"; 
        PLInstrumentResponseFWHMStartingValue (1, 1) {mustBeNumeric}
        PLInstrumentResponseLowerBound (1, 1) {mustBeNumeric}
        PLInstrumentResponseUpperBound (1, 1) {mustBeNumeric}
        PLInstrumentResponseFWHMFittedValue (1, 1) {mustBeNumeric}
        PLInstrumentResponseFWHMFittedVariance (1, 1) {mustBeNumeric}
        
        PLTimeZeroLabel {mustBeText} = "PL Time Zero";
        PLTimeZeroStartingValue (1, 1) {mustBeNumeric}
        PLTimeZeroLowerBound (1, 1) {mustBeNumeric}
        PLTimeZeroUpperBound (1, 1) {mustBeNumeric}
        PLTimeZeroFittedValue (1, 1) {mustBeNumeric}
        PLTimeZeroFittedVariance (1, 1) {mustBeNumeric}

        PLCompartmentIntensityStartingValues (:, 1) {mustBeNumeric}
        PLCompartmentIntensityLowerBounds (:, 1) {mustBeNumeric}
        PLCompartmentIntensityUpperBounds (:, 1) {mustBeNumeric}
        PLCompartmentIntensityFittedValues (:, 1) {mustBeNumeric}
        PLCompartmentIntensityFittedVariances (:, 1) {mustBeNumeric}

        % TA Data
        TAMeasuredMean (:, :) {mustBeNumeric}
        TAMeasuredVariance (:, :) {mustBeNumeric}
        TAFittedMean (:, :) {mustBeNumeric}
        TATimes (:, 1) {mustBeNumeric}
        TAFittedSpectra (:,:) {mustBeNumeric}
        TAFittedConcentrationProfiles (:, :) {mustBeNumeric}

        % Photoluminescence Data
        PLMeasuredData (:, 1) {mustBeNumeric}
        PLMeasuredVariance (:, 1)
        PLTimes (:, 1) {mustBeNumeric}
        PLFittedData (:, 1) {mustBeNumeric}
        PLFittedConcentrationProfiles (:, :) {mustBeNumeric}


        % Model
        model KineticModel

        % Fitting Options
        includeTwoPhotonAbsorption

    end % properties

    methods (Access = public)

        function obj = KineticFit(myModel)
            arguments
                myModel KineticModel
            end
            obj.model = myModel;

            % Generate the default conditions for the fit.
            obj.kineticsRateConstantStartingValues = obj.model.parameterValues;
            obj.kineticsRateConstantLowerBounds = -inf*ones(mymodel.nParameters, 1);
            obj.kineticsRateConstantUpperBounds = inf*ones(mymodel.nParameters, 1);

            obj.TAInstrumentResponseFWHMStartingValue = 0.1;
            obj.TAInstrumentResponseLowerBound = -inf;
            obj.TAInstrumentResponseUpperBound = inf;

            obj.TATimeZeroStartingValue = 0.1;
            obj.TATimeZeroLowerBound = -inf;
            obj.TATimeZeroUpperBound = inf;

            obj.TATimeZeroStartingValue = 0.1;
            obj.TATimeZeroLowerBound = -inf;
            obj.TATimeZeroUpperBound = inf;

            obj.PLInstrumentResponseFWHMStartingValue = 22; % 22 ps standard TCSPSC IRF width
            obj.PLInstrumentResponseLowerBound = -inf;
            obj.PLInstrumentResponseUpperBound = inf;
            
            obj.PLTimeZeroStartingValue = 0;
            obj.PLTimeZeroLowerBound = -inf;
            obj.PLTimeZeroUpperBound = inf;

            obj.PLCompartmentIntensityStartingValues = ones(mymodel.nCompartments, 1);
            obj.PLCompartmentIntensityLowerBounds = -inf*ones(mymodel.nCompartments, 1);
            obj.PLCompartmentIntensityUpperBounds = inf*ones(mymodel.nCompartments, 1);
            
            obj.includeTwoPhotonAbsorption = false;
             
        end

        function [simulatedData, simulatedConcentrationProfiles, simulatedSpectra] = simulateTA(obj)

            [TAParameters, ~, ~] = obj.buildTAParameters();

            obj.model.updateParameters(TAParameters);

            if obj.includeTwoPhotonAbsorption
                twopa = exp(-0.5*((obj.TATimes' - obj.model.timeZero)/obj.model.instrumentResponseFWHM).^2);
                simulatedConcentrationProfiles = [twopa, obj.model.generateConcentrations(obj.TATimes)];
            else
                simulatedConcentrationProfiles = obj.model.generateConcentrations(obj.TATimes);
            end

            simulatedSpectra = simulatedConcentrationProfiles\obj.TAMeasuredMean;

            simulatedData = simulatedConcentrationProfiles * simulatedSpectra;
        end

        function fitTransientAbsorption(obj)

            [TAParameters, TALowerBounds, TAUpperBounds] = obj.buildTAParameters();

            optimfun = @(Parameters) obj.TAFitFunction(Parameters, obj.model, obj.TAMeasuredMean);

            options = optimoptions("fmincon", 'Display', 'iter', 'Algorithm', 'interior-point');

            [TAParameterFittedValues, ~, ~, ~, ~, ~, hessian] = fmincon(optimfun, TAParameters, [], [], [], [], TALowerBounds, TAUpperBounds, [], options);
            
            obj.model.updateParameters(TAParameterFittedValues);

            TAParameterFittedVariance = diag(inv(hessian));
            
            obj.kineticsRateConstantFittedValues = TAParameterFittedValues(1:obj.model.nParameters, 1);
            obj.kineticsRateConstantFittedVariances = TAParameterFittedVariance(1:obj.model.nParameters, 1);
            
            obj.TAInstrumentResponseFWHMFittedValue = TAParameterFittedValues(obj.model.nParameters + 1, 1);
            obj.TAInstrumentResponseFWHMFittedVariance = TAParameterFittedVariance(obj.model.nParameters + 1, 1);

            obj.TATimeZeroFittedValue = TAParameterFittedValues(obj.model.nParameters + 1, 1);
            obj.TATimeZeroFittedVariance = TAParameterFittedVariance(obj.model.nParameters + 1, 1);
            
            obj.TAFittedConcentrationProfiles = obj.model.generateConcentrations(obj.TATimes);

            obj.TAFittedSpectra = obj.TAFittedConcentrationProfiles\obj.TAMeasuredMean;

            obj.TAFittedMean = obj.model.generateConcentrations(obj.TATimes) * obj.TAFittedSpectra;

        end

        function fitFluorescenceLifetime(obj)

            [PLParameters, PLLowerBounds, PLUpperBounds] = obj.buildPLParameters();

            optimfun = @(Parameters) obj.PLFitFunction(Parameters, obj.model, obj.PLMeasuredData);

            options = optimoptions("fmincon", 'Display', 'iter', 'Algorithm', 'interior-point');
            

            [PLFittedParameterValues, ~, ~, ~, ~, ~, hessian] = fmincon(optimfun, PLParameters, [], [], [], [], PLLowerBounds, PLUpperBounds, [], options);
            
            obj.model.updateParameters(PLFittedParameterValues(1:obj.model.nParameters + 2, 1));

            PLfittedParameterVariance = diag(inv(hessian));

            obj.kineticsRateConstantFittedValues = PLFittedParameterValues(1:obj.model.nParameters, 1);
            obj.kineticsRateConstantFittedVariances = PLfittedParameterVariance(1:obj.model.nParameters, 1);
            
            obj.PLInstrumentResponseFWHMFittedValue = PLFittedParameterValues(obj.model.nParameters + 1, 1);
            obj.PLInstrumentResponseFWHMFittedVariance = PLfittedParameterVariance(obj.model.nParameters + 1, 1);

            obj.PLTimeZeroFittedValue = PLFittedParameterValues(obj.model.nParameters + 1, 1);
            obj.PLTimeZeroFittedVariance = PLfittedParameterVariance(obj.model.nParameters + 1, 1);

            obj.PLFittedConcentrationProfiles = obj.model.generateConcentrations(obj.PLTimes);
           
            obj.PLCompartmentIntensityFittedValues = PLFittedParameterValues(obj.model.nParameters + 3:end, 1);
            obj.PLCompartmentIntensityFittedVariances = PLfittedParameterVariance(obj.model.nParameters + 3:end, 1);

            obj.PLFittedData = sum(obj.PLFittedConcentrationProfiles.*PLFittedParameterValues(obj.model.nParameters + 3:end, 1)', 2);

        end

    end

    methods (Access = private) 

        function diff = TAFitFunction(obj, parameters, model, data)
            
            model.updateParameters(parameters);

            if obj.includeTwoPhotonAbsorption
                twopa = exp(-0.5*((obj.transientAbsorptionTimes' - model.timeZero)/model.instrumentResponseFWHM).^2);
                concentrationProfiles = [twopa, model.generateConcentrations(obj.transientAbsorptionTimes)];
            else
                concentrationProfiles = model.generateConcentrations(obj.transientAbsorptionTimes);
            end

            extractedSpectra = concentrationProfiles\data;
        
            predictedData = concentrationProfiles*extractedSpectra;
        
            diff = sum((predictedData - data).^2, "all");
        
        end % transientAbsorptionFitFunction

        function diff = PLFitFunction(obj, parameters, model, data)
            % TO BE TESTED
            
            model.updateParameters(parameters(1:model.nParameters + 2, 1)); 

            Concs = model.generateConcentrations(obj.emissionLifetimeTimes);

            predictedData = sum(Concs.*parameters(model.nParameters + 3:end, 1)', 2);
            
            diff = sum((predictedData - data').^2, "all");

        
        end % fluorescenceLifetimeFitFunction

        function [TAParameters, TALowerBounds, TAUpperBounds] = buildTAParameters(obj)
            
            TAParameters = [obj.kineticsRateConstantStartingValues; obj.TAInstrumentResponseFWHMStartingValue; obj.TATimeZeroStartingValue];
            TALowerBounds = [obj.kineticsRateConstantLowerBounds; obj.TAInstrumentResponseLowerBound; obj.TATimeZeroLowerBound];
            TAUpperBounds = [obj.kineticsRateConstantUpperBounds; obj.TAInstrumentResponseUpperBound; obj.TATimeZeroUpperBound];

        end

        function [PLParameters, PLLowerBounds, PLUpperBounds] = buildPLParameters(obj)
            
            PLParameters = [obj.kineticsRateConstantStartingValues; obj.PLInstrumentResponseFWHMStartingValue; obj.PLTimeZeroStartingValue; obj.PLCompartmentIntensityStartingValues];
            PLLowerBounds = [obj.kineticsRateConstantLowerBounds; obj.PLInstrumentResponseLowerBound; obj.PLTimeZeroLowerBound; obj.PLCompartmentIntensityLowerBounds];
            PLUpperBounds = [obj.kineticsRateConstantUpperBounds; obj.PLInstrumentResponseUpperBound; obj.PLTimeZeroUpperBound; obj.PLCompartmentIntensityUpperBounds];

        end

    end % methods (Access = private) 

    methods % Setters

        function set.TATimes(obj, Times)

            if issorted(Times) && all(isfinite(Times)) && length(Times) == length(unique(Times))
                obj.TATimes = Times(:)'; % Ensure it's set as row vector
            else
                error('A problem was found in the vector of times, either there are non finite numbers, the list of times is does continuously increase or decrease, or there are duplicate time values');
            end

        end % set.TATimes

    end % methods Setters

end % KineticFit