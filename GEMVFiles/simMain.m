function[V2XData,vehiclesLatLon,vehiclesHeight,buildingsLatLon,foliageLatLon,...
    numTimesteps,numVehiclesPerTimestep,vehicleMidpointsLatLon,RSUsLatLon]...
    = simMain(vehiclesFile,bBoxVehicles,RSUFile,staticFile,foliageFile,...
    toleranceBuildings,commRange,numRowsPerVehicle,numCommPairs,...
    numCommPairsV2I,vehDimensionParams,lengthThreshold,maxNLOSbRange,...
    maxNLOSvRange,antennaHeight,polarization,vehReflRelPerm,...
    buildingReflRelPerm,freq,txPower,antennaGain,PLENLOSb,smallScale,...
    minDensityRange,NLOSvModelType,addNLOSvLoss,useReflDiffr,verbose,...
    V2XNames,commPairArray,commPairArrayV2I)
% SIMMAIN   Sets up simulation environment.
%   Loads the object outlines (vehicles, buildings, foliage).
%   Initiates output variables.
%   Runs simulation for each time step.
% Input:
%   see simSettings.m
% Output:
% (in V2XData struct)
%   largeScalePwrCell:  received power for each communication pair based on
%                       the large-scale path loss model
%   smallScaleVarCell:  small-scale variations for each comm. pair
%   communicationPairsCell:  list of comm. pairs that are within
%                       the maximum comm. range for the link type
%                       they form. Max comm. range for link type is
%                       defined in simSettings.m
%   communicationPairsCellAll: list of all generated communication pairs.
%                       NB: the list can, but may not contain all the
%                       communication pairs in the system (determined by
%                       variable numCommPairs in simSettings.m (if set to
%                       Inf, commPairsAll contains all pairs within comm.
%                       range in the system)
%   commPairTypeCell:   comm. pair link type (1-LOS; 2-NLOSb; 3-NLOSv)
%   effectivePairRangeCell: effective communication range for a pair.
%                       It is set to Inf if comm. pair if the vehicles are
%                       not within the maximum comm. range for the link
%                       type
%   numNeighborsPerVehPerIntervalCell: number of neighbors per vehicle and
%                       per timestep
% (in V2XData struct)
%
%   vehiclesLatLon:     latitude and longitude coordinates of vehicle
%                       outlines
%   vehiclesHeight:     height of vehicles
%   buildingsLatLon:    latitude and longitude coordinates of buildings
%   foliageLatLon:      latitude and longitude coordinates of foliage
%   numCommPairsPerTimestep: number of comm. pair per timestep
%   numTimesteps:       number of simulation timesteps
%   numVehiclesPerTimestep: number of vehicles per simulation timestep
%   vehicleMidpointsLatLon: latitude and longitude coordinates of vehicle
%                       	midpoints
%
% Copyright (c) 2014-2015, Mate Boban

% Load vehicle outlines
if ~isempty(vehiclesFile)
    [SUMOFile,vehicles,vehiclesLatLon,numTimesteps,numVehicles,...
        numVehiclesPerTimestep,vehiclesHeight,vehicleMidpoints,...
        vehicleMidpointsLatLon] = ...
        loadFunctions.loadVehicles(vehiclesFile,bBoxVehicles,...
        numRowsPerVehicle,lengthThreshold,vehDimensionParams,verbose);
else
    error('Vehicle array cannot be empty.');
end

%% Load RSUs
if ~isempty(RSUFile)
    RSUsLatLon = load(RSUFile);
    RSUs = RSUsLatLon;
    % Convert Lat/Lon to UTM
    [xxRSU,yyRSU,~] = ...
        externalCode.deg2utm.deg2utm(RSUsLatLon(:,2),RSUsLatLon(:,3));
    RSUs(:,2) = yyRSU;
    RSUs(:,3) = xxRSU;
else
    RSUsLatLon = [];
    RSUs = [];
end

%% Load building and foliage outlines
[buildings,buildingsLatLon,boundingBoxesBuildings,objectCellsBuildings,...
    BigBoxesBuildings,foliage,foliageLatLon,boundingBoxesFoliage,...
    objectCellsFoliage,BigBoxesFoliage] =...
    loadFunctions.loadBuildingsFoliage(staticFile,foliageFile,...
    toleranceBuildings,verbose);

%% Load communication pairs (if any)
commPairArray = load(commPairArray);
commPairArrayV2I = load(commPairArrayV2I);

%% Run simulation for each timestep
if round(numTimesteps)~=numTimesteps
    error('Something is wrong with the input file containing vehicle polygons');
end

% Initialize cells and arrays used to collect info for all timesteps
for ii=1:length(V2XNames)
    V2XData.(V2XNames{ii}).largeScalePwrCell = cell(numTimesteps,1);
    V2XData.(V2XNames{ii}).smallScaleVarCell = cell(numTimesteps,1);
    V2XData.(V2XNames{ii}).communicationPairsCell = cell(numTimesteps,1);
    V2XData.(V2XNames{ii}).communicationPairsCellAll = cell(numTimesteps,1);
    V2XData.(V2XNames{ii}).commPairTypeCell = cell(numTimesteps,1);
    V2XData.(V2XNames{ii}).effectivePairRangeCell = cell(numTimesteps,1);
    V2XData.(V2XNames{ii}).numNeighborsPerVehPerIntervalCell = cell(numTimesteps,1);
    V2XData.(V2XNames{ii}).communicationPairsLatLon = zeros(0);
    V2XData.(V2XNames{ii}).numCommPairsPerTimestep = zeros(numTimesteps,1);
end
% Counter for rows in vehicle array in case of SUMOFile (i.e., when
% there are non-fixed number of vehicles per timestep)
vehicleRowCounter = 0;
for kk = 1:numTimesteps
    % Get the current vehicle data (outlines and height)
    if ~SUMOFile
        currVehicles = vehicles((kk-1)*numVehicles*numRowsPerVehicle+...
            1:kk*numVehicles*numRowsPerVehicle,:);
        currVehiclesHeight = vehiclesHeight((kk-1)*numVehicles+...
            1:kk*numVehicles);
        currVehicleMidpoints = vehicleMidpoints((kk-1)*numVehicles+...
            1:kk*numVehicles,:);
        % Set the number of vehicles for current timestep
        numVehiclesPerTimestep=numVehicles;
    else
        currVehicles = vehicles(vehicleRowCounter+1:vehicleRowCounter...
            +numVehiclesPerTimestep(kk)*numRowsPerVehicle,:);
        currVehiclesHeight = vehiclesHeight(vehicleRowCounter/numRowsPerVehicle+...
            1:vehicleRowCounter/numRowsPerVehicle+numVehiclesPerTimestep(kk));
        currVehicleMidpoints = vehicleMidpoints(vehicleRowCounter/numRowsPerVehicle+...
            1:vehicleRowCounter/numRowsPerVehicle+numVehiclesPerTimestep(kk),:);
        vehicleRowCounter = vehicleRowCounter+size(currVehicles,1);
    end
    % Simulate current timestep
    [V2XDataTimestep]...
        = simOneTimestep(numCommPairs,numCommPairsV2I,currVehicles,...
        RSUs,BigBoxesBuildings,BigBoxesFoliage,objectCellsBuildings,...
        objectCellsFoliage,commRange,numRowsPerVehicle,...
        boundingBoxesBuildings,boundingBoxesFoliage,currVehiclesHeight,...
        maxNLOSbRange,maxNLOSvRange,antennaHeight,polarization,...
        vehReflRelPerm,buildingReflRelPerm,freq,txPower,antennaGain,...
        PLENLOSb,smallScale,minDensityRange,NLOSvModelType,addNLOSvLoss,...
        verbose,useReflDiffr,commPairArray,commPairArrayV2I,...
        currVehicleMidpoints,V2XNames);
    
    for ii=1:length(V2XNames)
        % Store the simulation results from current timestep in cells
        if ~isempty(V2XDataTimestep.(V2XNames{ii}))
            V2XData.(V2XNames{ii}).largeScalePwrCell{kk,1} = ...
                V2XDataTimestep.(V2XNames{ii}).largeScalePwr;
            V2XData.(V2XNames{ii}).smallScaleVarCell{kk,1} = ...
                V2XDataTimestep.(V2XNames{ii}).smallScaleVar;
            V2XData.(V2XNames{ii}).commPairTypeCell{kk,1} = ...
                V2XDataTimestep.(V2XNames{ii}).commPairType;
            V2XData.(V2XNames{ii}).communicationPairsCellAll{kk,1} = ...
                V2XDataTimestep.(V2XNames{ii}).communicationPairs;
            V2XData.(V2XNames{ii}).communicationPairsCell{kk,1} = ...
                V2XDataTimestep.(V2XNames{ii}).communicationPairs(...
                V2XDataTimestep.(V2XNames{ii}).effectivePairRange~=Inf,:);
            V2XData.(V2XNames{ii}).effectivePairRangeCell{kk,1} = ...
                V2XDataTimestep.(V2XNames{ii}).effectivePairRange;
            V2XData.(V2XNames{ii}).numCommPairsPerTimestep(kk) = ...
                length(V2XDataTimestep.(V2XNames{ii}).largeScalePwr);
            
            %% Analyze the neighborhood of each vehicle/RSU
            if strcmpi(V2XNames{ii},'v2v')
                % V2V pairs: get both columns to calculate neighbors
                % (other vehicles only)
                inputNodes = ...
                    V2XData.(V2XNames{ii}).communicationPairsCell{kk,1}(:);
            else
                % V2I pairs; get RSUs only (first column); neighbors
                % (other vehicles only)
                inputNodes = ...
                    V2XData.(V2XNames{ii}).communicationPairsCell{kk,1}(:,1);
            end
            uniqueNeighborsCommPairs = unique(inputNodes);
            numNeighborsPerVehPerInterval = zeros(length(uniqueNeighborsCommPairs),1);
            for jj=1:length(uniqueNeighborsCommPairs)
                numNeighborsPerVehPerInterval(jj) = ...
                    sum(inputNodes==uniqueNeighborsCommPairs(jj));
            end
            % Store vehicle/RSU IDs (from object Midpoints) in first column of the
            % cell, and number of vehicles in second
            V2XData.(V2XNames{ii}).numNeighborsPerVehPerIntervalCell{kk,1} =...
                uniqueNeighborsCommPairs;
            V2XData.(V2XNames{ii}).numNeighborsPerVehPerIntervalCell{kk,2} =...
                numNeighborsPerVehPerInterval;
        else
            V2XData.(V2XNames{ii}) = [];
        end
    end
    fprintf('Time-step: %i\n',kk);
end
end