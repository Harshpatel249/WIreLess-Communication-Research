function[V2XDataTimestep]...    
    = simOneTimestep(numCommPairs,numCommPairsV2I,vehicles,RSUs,RTreeB,RTreeF,...
    objectCellsBuildings,objectCellsFoliage,commRange,numRowsPerVehicle,...
    boundingBoxesBuildings,boundingBoxesFoliage,vehiclesHeight,...
    MaxNLOSbRange,MaxNLOSvRange,antennaHeight,polarization,vehReflRelPerm,...
    buildingReflRelPerm,freq,txPower,antennaGain,...
    PLENLOSb,smallScale,minDensityRange,NLOSvModelType,addNLOSvLoss,verbose,...
    useReflDiffr,commPairArray,commPairArrayV2I,vehicleMidpoints,V2XNames)
% SIMONETIMESTEP Executes one simulation timestep. 
%   Builds vehicle, buildings, and foliage R-tree
%   Generates random communication pairs (if not provided)
%   Determines the link type (LOS, NLOSv, NLOSb) between each comm. pair
%   Calculates large- and small-scale signal variation for each comm. pair
%   Calculates received power for each comm. pair
%   (optionally) Calculates 1-interaction reflections and diffractions
% Input:
%   See simSettings.m and simMain.m
%
% Output (in V2XDataTimestep struct, separate for each link type: V2V and V2I)
%   largeScalePwr:      large-scale received power (including path-loss and
%                       shadowing) for each communication pair
%   smallScaleVar:      small-scale signal variations for each
%                       communication pair
%   communicationPairs: array containing pairwise IDs of vehicles
%                       participating in communication pairs
%   commPairType:       type of link between a communication pair: LOS = 1,
%                       NLOSb/NLOSf = 2, NLOSv = 3 
%   effectivePairRange: for each communication pair, effectivePairRange is 
%                       the maximum allowed distance from any object to
%                       both Tx and Rx (i.e., it defines an ellipse with Tx
%                       and Rx as foci). Used to search for objects around
%                       Tx and Rx that could have impact on signal
%                       variation (see paper for details).
%   objectCellsVehicles:cell containing outline of all vehicles
%
% Copyright (c) 2014-2015, Mate Boban

ticForTimestep=tic;
% Speed of light (m/s) constant
c = 299792458;
% Build the vehicle R-tree
[bBoxesVehicles,objectCellsVehicles,RTreeV,~] = ...
    RTree.prepareData(vehicles,verbose);

% NB: in comments, the term "link" and "communication pair" are used
% interchangeably 
 
%% Generate and select communication pairs (if not provided)
if isempty(commPairArray) && isempty(commPairArrayV2I)
    % Get all communicating vehicle pairs within commRange
    if exist('rangesearch') == 5 || exist('rangesearch') == 2
        % Rangesearch command is available in Statistics Toolbox (at
        % least from Matlab version 2011b)
        if numCommPairs>0
            tic
            communicationPairs = rangesearch...
                (vehicleMidpoints,vehicleMidpoints,commRange);
            if verbose
                fprintf(['Matlab rangesearch command takes %f seconds '...
                    'to map the V2V communication pairs.\n'],toc);
            end
        else
            communicationPairs = [];
        end
        % Find V2I (RSU-vehicle) comm. pairs
        if ~isempty(RSUs) && numCommPairsV2I>0
            tic
            communicationPairsV2I = rangesearch...
                (RSUs(:,[2 3]),vehicleMidpoints,commRange);
            if verbose
                fprintf(['Matlab rangesearch command takes %f seconds '...
                    'to map the V2I communication pairs.\n'],toc);
            end
        else
            communicationPairsV2I = [];
        end
    else
        % Use the much slower function for getting all communicating
        % vehicle pairs within commRange
        if numCommPairs>0
            communicationPairs = commPairs.mapCommPairs...
                (vehicleMidpoints,commRange,verbose);
        else
            communicationPairs = [];
        end
        if ~isempty(RSUs) && numCommPairsV2I>0
            communicationPairsV2I = commPairs.mapCommPairs...
                (vehicleMidpoints,commRange,verbose,RSUs(:,[2 3]));
        else
            communicationPairsV2I = [];
        end
    end
    % Get a given number of random communication pairs
    if numCommPairs>0
        communicationPairs = commPairs.getRandCommPairs...
            (communicationPairs,numCommPairs);
    end
    if ~isempty(RSUs) && numCommPairsV2I>0
        communicationPairsV2I = commPairs.getRandCommPairs...
            (communicationPairsV2I,numCommPairsV2I,'V2I');
        if exist('rangesearch') == 5 || exist('rangesearch') == 2
            communicationPairsV2I = fliplr(communicationPairsV2I);
        end
    end
    % NB: Depending on maximum communication ranges, there might be no
    % communication pairs for the current timestep
    if isempty(communicationPairs)
        fprintf(['There are no V2V communication pairs given the current\n'...
            'maximum communication ranges! Moving to the next timestep'...
            ' (if any)...\n']);
        V2XDataTimestep.V2V = [];
        
    end
    if isempty(communicationPairsV2I)
        if ~isempty(RSUs)
            fprintf(['There are no V2I communication pairs given the current\n'...
                'maximum communication ranges! Moving to the next timestep'...
                ' (if any)...\n']);
        end
        V2XDataTimestep.V2I = [];
    end
    if isempty(communicationPairs) && isempty(communicationPairsV2I)
        return
    end
else
    if size(commPairArray,2)==2 && max(max(commPairArray))<=...
            size(vehicleMidpoints,1)
        communicationPairs = commPairArray;
    else
        error(['Something is wrong with the supplied communication '...
            'pairs file! The structure of the file is: two columns,'...
            'each row represents vehicle ID for a communication pair.']);
    end
end

%% Calculate the distance between randomly selected communication pairs.
if ~isempty(communicationPairs)
    commPairVeh1 = vehicleMidpoints(communicationPairs(:,1),:);
    commPairVeh2 = vehicleMidpoints(communicationPairs(:,2),:);
end
if ~isempty(communicationPairsV2I)
    commPairV2I1 = RSUs(communicationPairsV2I(:,1),[2 3]);
    commPairV2I2 = vehicleMidpoints(communicationPairsV2I(:,2),:);
end

if ~isempty(communicationPairs)
    commPairDist = sqrt((commPairVeh1(:,1)-commPairVeh2(:,1)).^2+...
        (commPairVeh1(:,2)-commPairVeh2(:,2)).^2);
end
if ~isempty(communicationPairsV2I)
    commPairDistV2I = sqrt((commPairV2I1(:,1)-commPairV2I2(:,1)).^2+...
        (commPairV2I1(:,2)-commPairV2I2(:,2)).^2);
end

% Initialize RSU midpoints and heights variables
RSUHeights = [];
RSUMidpoints = [];

%% New code to perform both V2V and V2I calculations
for ii=1:length(V2XNames)
    if strcmpi(V2XNames{ii},'v2v') 
        if isempty(communicationPairs)
            continue
        else
            % V2V links;
            combinedMidpoints = vehicleMidpoints;
            objectHeights = vehiclesHeight;
        end
    elseif strcmpi(V2XNames{ii},'v2i')
        if isempty(communicationPairsV2I)
            continue
        else
            % V2I links;
            communicationPairs = communicationPairsV2I;
            commPairDist = commPairDistV2I;            
            combinedMidpoints = {RSUs(:,[2 3]),vehicleMidpoints};
            RSUHeights = RSUs(:,4);
            RSUMidpoints = RSUs(:,[2 3]);
        end
    else
        error('Unknown link type');
    end
%% Get links obstructed by buildings    
if ~isempty(RTreeB)
    obstructedCommPairsB = polygonManipulation.getObstructingObjectsV2X...
        (communicationPairs,objectCellsBuildings,combinedMidpoints,RTreeB,-1,1,verbose);
else
    obstructedCommPairsB = zeros(size(communicationPairs,1),1);
end

%% Get links obstructed by foliage (among those not obstructed by buildings)
if ~isempty(RTreeF)
    obstructedCommPairsFList = polygonManipulation.getObstructingObjectsV2X...
        (communicationPairs(~obstructedCommPairsB,:),objectCellsFoliage,...
        combinedMidpoints,RTreeF,-1,0,verbose);
    obstructedCommPairsF2 = ~cellfun(@isempty,obstructedCommPairsFList);
    obstructedCommPairsF = zeros(size(communicationPairs,1),1);
    obstructedCommPairsF(~obstructedCommPairsB) = obstructedCommPairsF2;
else
    obstructedCommPairsF = zeros(size(communicationPairs,1),1);
end
    
%% Get links obstructed by vehicles (among those not obstructed by buildings or foliage)
obstructedCommPairsV = polygonManipulation.getObstructingObjectsV2X...
    (communicationPairs(~obstructedCommPairsB & ~obstructedCommPairsF,:),...
    objectCellsVehicles,combinedMidpoints,RTreeV,numRowsPerVehicle,0,verbose);    
    % Convert array to cell; merge NLOSb and NLOSv. -1 represents comm.
    % pairs that are obstructed by buildings 
    obstructedCommPairs = num2cell(obstructedCommPairsB.*-1);    
    obstructedCommPairs(~(obstructedCommPairsB|obstructedCommPairsF),1)...
        = obstructedCommPairsV;

%% Categorize links/comm. pairs 
%   LOS (line of sight)
%   NLOSb (obstructed by buildings)
%   NLOSf (obstructed by foliage)
%   NLOSv (obstructed by vehicles) 
LOSPairs = cellfun(@isempty,obstructedCommPairs);
NLOSbPairs = logical(obstructedCommPairsB);
NLOSfPairs = logical(obstructedCommPairsF);
NLOSvPairs = ~(NLOSfPairs|NLOSbPairs|LOSPairs);

%% Calculate effective range for each link/comm. pair 
% effectivePairRange is the maximum distance from any object to both Tx and
% Rx in a comm. pair (i.e., effectivePairRange =
% dist(Tx,point)+dist(point,Rx)). Defines an ellipse with Tx and Rx as
% foci. Used to search for objects around Tx and Rx that could have impact
% on signal variation (see paper for details).
effectivePairRange = ones(size(communicationPairs,1),1);
% Effective range for LOS pairs
LOSPairsRange = max(commPairDist(LOSPairs),commRange);
% Effective range for NLOSb pairs
NLOSbPairsRange = max(commPairDist(NLOSbPairs | NLOSfPairs),MaxNLOSbRange);
% Effective range for NLOSv pairs
NLOSvPairsRange = max(commPairDist(NLOSvPairs),MaxNLOSvRange);
% If distance for a NLOSb/NLOSf/NLOSv comm. pair is above respective
% maximum range, assume the pair cannot communicate and set the range to
% Inf. Used for filtering "unfeasible" links
NLOSbPairsRange(NLOSbPairsRange>MaxNLOSbRange)=Inf;
NLOSvPairsRange(NLOSvPairsRange>MaxNLOSvRange)=Inf;
% Set the effective range
effectivePairRange(LOSPairs) = LOSPairsRange;
effectivePairRange(NLOSbPairs | NLOSfPairs) = NLOSbPairsRange;
effectivePairRange(NLOSvPairs) = NLOSvPairsRange; 

%% Calculate the "real" LOS distance: distance between Tx and Rx that
% accounts for the height of the antennas (i.e., not only location of
% vehicle midpoints).
if strcmpi(V2XNames{ii},'v2v')
    realLOSDists = sqrt(commPairDist.^2+(vehiclesHeight(communicationPairs...
        (:,1))+antennaHeight-(vehiclesHeight(communicationPairs(:,2))+...
        antennaHeight)).^2);
else
    % Distance from RSUs to vehicles
    realLOSDists = sqrt(commPairDist.^2+(RSUs(communicationPairs(:,1),4)...
        -(vehiclesHeight(communicationPairs(:,2))+...
        antennaHeight)).^2);
end
 
%% Get only those comm. pairs that are designated as "feasible" 
feasibleCommPairs = communicationPairs(effectivePairRange~=Inf,:);
feasibleEffRange = effectivePairRange(effectivePairRange~=Inf);

% Get the IDs of vehicles and RSUs that are in feasibleCommPairs
IDRSUMidpoints = [];
if strcmpi(V2XNames{ii},'v2v')
    % Get vehicle IDs
    IDVehMidpoints = ismember(1:size(vehicleMidpoints,1),...
        unique(feasibleCommPairs));
else
    % Get RSU IDs
    IDRSUMidpoints = ismember(1:size(RSUs,1),...
        unique(feasibleCommPairs(:,1)));
    % Get vehicle IDs
    IDVehMidpoints = ismember(1:size(vehicleMidpoints,1),...
        unique(feasibleCommPairs));
end

% If no RSUs, pass an empty array
if strcmpi(V2XNames{ii},'v2v') || isempty(RSUs)
    RSUMidpoints=[];
%elseif ~isempty(RSUs)
else
    RSUMidpoints=RSUs(:,[2,3]);
%else
%    RSUMidpoints=[];
end

%% Get objects inside comm. pair ellipse
% Get buildings inside comm. pair ellipse
if ~isempty(boundingBoxesBuildings) && ~isempty(feasibleCommPairs)
    [jointBuildings,jointDistances] = ...
        polygonManipulation.getObjectsInsideEllipse...
        (feasibleCommPairs,vehicleMidpoints,RSUMidpoints,commRange,feasibleEffRange,...
        IDVehMidpoints,IDRSUMidpoints,boundingBoxesBuildings,verbose);
else
    jointBuildings = cell(1,size(feasibleCommPairs,1));
    jointDistances = cell(1,size(feasibleCommPairs,1));
end

% Get vehicles inside comm. pair ellipse
if ~isempty(bBoxesVehicles)
    if ~isempty(feasibleCommPairs)
    [jointVehicles,jointVehicleDistances] = ...
        polygonManipulation.getObjectsInsideEllipse...
        (feasibleCommPairs,vehicleMidpoints,RSUMidpoints,commRange,feasibleEffRange,...
        IDVehMidpoints,IDRSUMidpoints,bBoxesVehicles,verbose);
    else
        jointVehicles = cell(1,size(feasibleCommPairs,1));
        jointVehicleDistances = cell(1,size(feasibleCommPairs,1));
    end
else
    error(['There needs to be at least one vehicle in the system for '...
        'the simulation to make sense!']);
end

% Get foliage inside comm. pair ellipse
if ~isempty(boundingBoxesFoliage) && ~isempty(feasibleCommPairs)
    [jointFoliage,~] = ...
        polygonManipulation.getObjectsInsideEllipse...
        (feasibleCommPairs,vehicleMidpoints,RSUMidpoints,commRange,feasibleEffRange,...
        IDVehMidpoints,IDRSUMidpoints,boundingBoxesFoliage,verbose);
else
    jointFoliage = cell(1,size(feasibleCommPairs,1));    
end

% Exclude the vehicles in the comm. pair from the jointVehicles and
% jointVehicleDistances 
if ~isempty(feasibleCommPairs) 
    fh = str2func('commPairs.excludeCommPairJointVeh');
    [jointVehicles,jointVehicleDistances] = cellfun(fh,jointVehicles,...
        jointVehicleDistances,num2cell(feasibleCommPairs,2),'uni',false);
end

%% Set the comm. pair/link type
% LOS = 1, NLOSb/NLOSf = 2, NLOSv = 3
commPairType = ones(size(communicationPairs,1),1);
commPairType(LOSPairs) = 1;
commPairType(NLOSbPairs | NLOSfPairs) = 2;
commPairType(NLOSvPairs) = 3;

% Get the link type of all feasibleCommPairs and their distances
feasibleCommPairType = commPairType(effectivePairRange~=Inf,:);
feasibleCommPairDists = realLOSDists(effectivePairRange~=Inf);
% Distance for NLOSb comm. pairs (used for reflections and diffractions) 
feasibleCommPairDistsNLOSb = feasibleCommPairDists(feasibleCommPairType==2);

% Get the correct antenna gain
if strcmpi(V2XNames{ii},'v2v')
    Gt = antennaGain.veh;
    Gr = antennaGain.veh;    
elseif strcmpi(V2XNames{ii},'v2i')
    Gt = antennaGain.veh;
    Gr = antennaGain.RSU;
end

%% Calculate power for NLOSf links
if sum(NLOSfPairs(effectivePairRange~=Inf))>0
    % Get obstructing foliage IDs for NLOSf comm. pairs
    obstructedCommPairsFListTemp = cell(size(communicationPairs,1),1);
    obstructedCommPairsFListTemp(~obstructedCommPairsB) = obstructedCommPairsFList;
    obstructedCommPairsFListObsOnly2 = obstructedCommPairsFListTemp...
        (NLOSfPairs==1 & effectivePairRange~=Inf);
    % Calculate power for NLOSf links
    [powerNLOSf,~] = LOSNLOS.NLOSf(objectCellsFoliage,...
        NLOSfPairs(effectivePairRange~=Inf),feasibleCommPairs,...
        obstructedCommPairsFListObsOnly2,feasibleCommPairDists,...
        vehicleMidpoints,txPower,c,freq,Gt,Gr,verbose);
else
    powerNLOSf = [];    
end
    
%% Calculate reflections and diffractions
if ~useReflDiffr
    finalReflRays=[];
    reflEfields=[];
    reflDist=[];
    diffrEfields=[];
    finalDiffrRays=[];
    diffrDist=[];
else
    % NB: the code below calculates V2X diffractions and reflections (i.e.,
    % V2V and V2I); if you want to calculate reflection for I2V, then RSU
    % power from RSUs(:,5) needs to be plugged in instead of txPower in
    % reflect and diffract functions
    [finalReflRays,feasibleCommPairsNLOSb,jointBuildingsNLOSb,...
        jointDistancesNLOSb,reflEfields,reflDist] = reflections.reflect...
        (feasibleCommPairs,feasibleCommPairType,antennaHeight,jointBuildings,...
        jointDistances,jointVehicles,objectCellsBuildings,objectCellsVehicles,...
        objectCellsFoliage,numRowsPerVehicle,vehicleMidpoints,...
        vehiclesHeight,RTreeB,RTreeV,RTreeF,buildingReflRelPerm,vehReflRelPerm,...
        polarization,txPower,Gt,verbose);
    if ~isempty(finalReflRays)
        [diffrEfields,finalDiffrRays,diffrDist] = diffractions.diffract...
            (finalReflRays,feasibleCommPairsNLOSb,feasibleCommPairDistsNLOSb,...
            feasibleEffRange,feasibleCommPairType,jointBuildingsNLOSb,...
            jointDistancesNLOSb,objectCellsBuildings,objectCellsVehicles,...
            vehicleMidpoints,vehiclesHeight,antennaHeight,RTreeB,RTreeV,...
            txPower,Gt,Gr,c,freq,verbose);
    else
        diffrEfields=[];
        finalDiffrRays=[];
        diffrDist=[];
    end
end

%% Calculate rec. power for LOS and NLOSv links
[powerNLOSv,~] = LOSNLOS.LOSNLOSv(objectCellsVehicles,LOSPairs,...
    NLOSvPairs,communicationPairs,obstructedCommPairs,effectivePairRange,...
    realLOSDists,vehiclesHeight,vehicleMidpoints,RSUHeights,RSUMidpoints,...
    antennaHeight,txPower,c,freq,Gt,Gr,polarization,NLOSvModelType,...
    addNLOSvLoss,V2XNames{ii},verbose);

%% Calculate small-scale signal variations
smallScaleVar=zeros(0);
if sum(effectivePairRange<Inf)>0
smallScaleVar = smallScaleVariations.smallScaleVariation(feasibleCommPairs,...
    jointVehicles,jointBuildings,jointFoliage,feasibleEffRange,...
    effectivePairRange,minDensityRange,objectCellsBuildings,objectCellsFoliage,...
    LOSPairs,NLOSvPairs,NLOSfPairs,NLOSbPairs,smallScale,V2XNames{ii});
end
    
%% Combine large-scale signal variation results for all comm. pair types in one array
largeScalePwr = powerCalculations.largeScaleVariations...
    (communicationPairs,LOSPairs,NLOSvPairs,NLOSbPairs,NLOSfPairs,...
    effectivePairRange,useReflDiffr,powerNLOSv,powerNLOSf,finalReflRays,...
    reflEfields,realLOSDists,reflDist,finalDiffrRays,diffrEfields,...
    diffrDist,c,freq,txPower,Gr,Gt,PLENLOSb,V2XNames{ii});

%% Store variables
V2XDataTimestep.(V2XNames{ii}).largeScalePwr = largeScalePwr;
V2XDataTimestep.(V2XNames{ii}).smallScaleVar = smallScaleVar;
V2XDataTimestep.(V2XNames{ii}).commPairType = commPairType;
V2XDataTimestep.(V2XNames{ii}).communicationPairs = communicationPairs;
V2XDataTimestep.(V2XNames{ii}).effectivePairRange = effectivePairRange;
V2XDataTimestep.(V2XNames{ii}).objectCellsVehicles = objectCellsVehicles;

% In case of V2I links, add the RSU Tx power difference to I2V "direction"
% of the link (the default power is calculated for V2I "direction");
if strcmpi(V2XNames{ii},'v2i')
    RSUPower = RSUs(communicationPairs(effectivePairRange~=Inf,1),5);
    V2XDataTimestep.(V2XNames{ii}).largeScalePwrI2V = ...
        largeScalePwr - txPower + RSUPower;
end
end
