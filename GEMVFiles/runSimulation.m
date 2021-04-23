% RUNSIMULATION Runs the simulation. 
%   Loads the settings from simSettings.m. 
%   Saves the simulation output to .mat and .csv files. 
%   Generates Google Earth Visualization.
%   For details on variables, see simSettings.m and simMain.m
%
% Usage: runSimulation
%
% Copyright (c) 2014-2015, Mate Boban

clear all;
clear global;

fprintf('GEMV^2: Geometry-based Efficient propagation Model for V2V communication\n');
fprintf('Copyright (c) 2014, Mate Boban\n');

% Load the simulation settings (contained and explained in simSettings.m)
simSettings;

%% Run simMain (for explanation of returned variables, see simMain.m)    
[V2XData,vehiclesLatLon,vehiclesHeight,buildingsLatLon,foliageLatLon,...
    numTimesteps,numVehiclesPerTimestep,vehicleMidpointsLatLon,RSUsLatLon] = ...
    simMain(vehiclesFile,bBoxVehicles,RSUFile,staticFile,...
    foliageFile,toleranceBuildings,commRange,numRowsPerVehicle,...
    numCommPairs,numCommPairsV2I,vehDimensionParams,lengthThreshold,...
    maxNLOSbRange,maxNLOSvRange,antennaHeight,polarization,vehReflRelPerm,...
    buildingReflRelPerm,freq,txPower,antennaGain,PLENLOSb,smallScale,minDensityRange,...
    NLOSvModelType,addNLOSvLoss,useReflDiffr,verbose,V2XNames,...
    commPairArray,commPairArrayV2I);

%% Save output in a matlab .mat and comma-separated .csv file.
% Description of output variables stored in V2XData struct:
%   
%   largeScalePwr:              received power for each comm. pair based on
%                               the large-scale path loss model
%   smScaleVar:                 small-scale variation for each comm. pair
%   commPairs:                  list of comm. pairs that are within
%                               the maximum comm. range for the link type
%                               they form. Max comm. range for link type is
%                               defined in simSettings.m
%   commPairsAll:               list of all generated communication pairs.
%                               NB: the list can, but may not contain all
%                               the communication pairs in the system
%                               (determined by variable numCommPairs in
%                               simSettings.m (if set to Inf, commPairsAll
%                               contains all pairs within comm. range in
%                               the system)
%   commPairType:               comm. pair link type (1-LOS; 2-NLOSb; 3-NLOSv)
%   effectivePairRange:         effective communication range for a pair.
%                               It is set to Inf if comm. pair if the
%                               vehicles are not within the maximum comm.
%                               range for the link type
%   numNeighborsPerVehicle:     Number of neighboring vehicles for which
%                               the received power is above rec. threshold
%   numCommPairsPerTimestep:    Number of feasible communication pairs for
%                               each timestep (i.e., number of pairs for
%                               which effectivePairRange~=Inf). This
%                               variable is useful for getting the
%                               per-timestep information from e.g.,
%                               largeScalePwr or smScaleVar

% Save output to .mat file
fileName = ['outputSim/',date,'_outputSim.mat'];
save(fileName,'V2XData');
timeNow = now;

for ii=1:length(V2XNames)
    if ~isempty(V2XData.(V2XNames{ii}))
        % Turn cells to arrays
        largeScalePwrCelltoArray = cell2mat(V2XData.(V2XNames{ii}).largeScalePwrCell);
        smallScaleVarCelltoArray = cell2mat(V2XData.(V2XNames{ii}).smallScaleVarCell);
        communicationPairsCelltoArray = cell2mat(V2XData.(V2XNames{ii}).communicationPairsCell);
        communicationPairsCellAlltoArray = cell2mat(V2XData.(V2XNames{ii}).communicationPairsCellAll);
        commPairTypeCelltoArray = cell2mat(V2XData.(V2XNames{ii}).commPairTypeCell);
        effectivePairRangeCelltoArray = cell2mat(V2XData.(V2XNames{ii}).effectivePairRangeCell);
        numNeighborsPerVehPerIntervalCelltoArray = cell2mat(V2XData.(V2XNames{ii}).numNeighborsPerVehPerIntervalCell);
        numCommPairsPerTimestep = V2XData.(V2XNames{ii}).numCommPairsPerTimestep; 
        
        % Save variables to .csv files
        fileNamelargeScalePwrCell = ['outputSim/',date,'_largeScalePwr_',V2XNames{ii},'.csv'];
        fileNamesmallScaleVarCell = ['outputSim/',date,'_smScaleVar_',V2XNames{ii},'.csv'];
        fileNamecommunicationPairsCell = ['outputSim/',date,'_commPairs_',V2XNames{ii},'.csv'];
        fileNamecommunicationPairsCellAll = ['outputSim/',date,'_commPairsAll_',V2XNames{ii},'.csv'];
        fileNamecommPairTypeCell = ['outputSim/',date,'_commPairType_',V2XNames{ii},'.csv'];
        fileNameeffectivePairRangeCell = ['outputSim/',date,'_effectivePairRange_',V2XNames{ii},'.csv'];
        fileNamenumNeighborsPerNodePerIntervalCell = ['outputSim/',date,'_numNeighborsPerNode_',V2XNames{ii},'.csv'];
        fileNamenumCommPairsPerTimestepCell = ['outputSim/',date,'_numCommPairsPerTimestep_',V2XNames{ii},'.csv'];
        
        % Write the output to .csv files
        dlmwrite(fileNamelargeScalePwrCell,largeScalePwrCelltoArray);
        dlmwrite(fileNamesmallScaleVarCell,smallScaleVarCelltoArray);
        dlmwrite(fileNamecommunicationPairsCell,communicationPairsCelltoArray);
        dlmwrite(fileNamecommunicationPairsCellAll,communicationPairsCellAlltoArray);
        dlmwrite(fileNameeffectivePairRangeCell,effectivePairRangeCelltoArray);
        dlmwrite(fileNamenumNeighborsPerNodePerIntervalCell,numNeighborsPerVehPerIntervalCelltoArray);
        dlmwrite(fileNamenumCommPairsPerTimestepCell,numCommPairsPerTimestep);
    end

    %% Google Earth Visualization
    if GEVisualize
        if ~isempty(V2XData.(V2XNames{ii}))
            if strcmpi(V2XNames{ii},'v2v')
                %objectMidpointsLatLon = vehicleMidpointsLatLon;
                %objectHeight = vehiclesHeight;
                RSUMidpointsLatLon = [];
                RSUHeight = [];
            elseif strcmpi(V2XNames{ii},'v2i')
                RSUMidpointsLatLon = RSUsLatLon(:,[2,3]);
                RSUHeight = RSUsLatLon(:,4);
            else
                disp('Do not know how to plot links other than V2V and V2I. Skipping plotting...');
                continue
            end            
            % Run Google Earth output function
                plotFunctions.GEOutput(plotPoly,plotRxPwr,plotNeighbors,...
                    vehicleMidpointsLatLon,vehiclesLatLon,vehiclesHeight,...
                    RSUMidpointsLatLon,RSUHeight,numRowsPerVehicle,...
                    buildingsLatLon,foliageLatLon,verbose,timeNow,...
                    V2XData.(V2XNames{ii}).numNeighborsPerVehPerIntervalCell,...
                    V2XData.(V2XNames{ii}).numCommPairsPerTimestep,...
                    V2XData.(V2XNames{ii}).communicationPairsCell,0,...
                    V2XData.(V2XNames{ii}).largeScalePwrCell,...
                    V2XData.(V2XNames{ii}).smallScaleVarCell,...
                    V2XNames{ii},numTimesteps,numVehiclesPerTimestep);      
        end
    end
end