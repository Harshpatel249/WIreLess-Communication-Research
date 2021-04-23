% Simulation settings
%
% Contains the most relevant simulation parameters. Loaded each time
% simulation is run.
%
% Copyright (c) 2014, Mate Boban

%% Different "pre-installed" scenarios
%   1 - Porto entire city
%   2 - Porto quarter of the city
%   3 - Cologne
%   4 - Newcastle
%   5 - Bologna (scenario with both V2V and V2I)
scenario = 1;

%% Names of link types (supported: V2V and V2I)
V2XNames = {'V2V'};

%% Parameters related to vehicle, foliage, and building outlines
% Load the file containing vehicle outlines/polygons in a three column
% array [vehicleID|Latitude|Longitude] or a five column array
% [vehicleID|midpointLat|midpointLong|vehicleType|bearing] provided by
% SUMO. In the latter case, vehicle polygons will be generated by
% generateVehiclePolygons.m using specified vehicle dimension statistics.
% Vehicles in the entire city of Porto
%   vehiclesFile = 'inputPolygon/vehiclesPorto';
% Also in the package:
% Vehicles in one quarter of Porto (bottom left quarter)
%   vehiclesFile = 'inputPolygon/vehiclesPortoQuarter'; 
% SUMO simulation in the city of Cologne, Germany
%   vehiclesFile = 'inputMobilitySUMO/CologneSUMOMobility.xml'; 
% SUMO simulation in the city of Newcastle, UK
%   vehiclesFile = 'inputMobilitySUMO/NewcastleSUMOMobility.xml';
% SUMO simulation in the city of Bologna, Italy (suitable for V2I
% simulation combined with Bologna RSUs) 
%   vehiclesFile = 'bologna/bolognaVehicles25sec.xml';

%% Load the file containing outlines of static objects (buildings, foliage)
% Buildings in entire city of Porto
%   staticFile = 'inputPolygon/buildingsPorto'; 
% Also in the package:
% Buildings in one quarter of Porto (bottom left quarter)
%   staticFile = 'inputPolygon/buildingsPortoQuarter'; 
% Buildings in part of Cologne (OpenStreetMap .osm)
%   staticFile = 'inputPolygon/CologneCityPart.osm'; 
% Buildings in part of Newcastle, UK (OpenStreetMap .osm)
%   staticFile = 'inputPolygon/Newcastle.osm';
% Buildings in part of Bologna, Italy (OpenStreetMap .osm)
%   staticFile = 'inputPolygon/bolognaLeftHalfRSU1-3-10.osm';
%   staticFile = 'inputPolygon/bolognaRightHalfRSU2-4-5-6-7-8-9.osm';

%% Load the file containing Locations, height, and power of RSUs in a five
% column array: [RSU_ID|Latitude|Longitude|Height(m)|TxPower(dBm)]
% If you need V2V only, use empty string
RSUFile = '';

% Bounding box (bottom-left, top-right) for removing vehicles outside of
% the desired area (applies for SUMO-generated mobility files). 
% If bounding box is not needed, set it to [-Inf,-Inf,Inf,Inf]; 
% For the included Cologne SUMO mobility dataset, to match the included
% CologneCityPart.osm static file, set the bBoxVehicles as follows: 
% bBoxVehicles = [6.93,50.913,7,50.96];
% bBoxVehicles = [11.3259 44.4829 11.3629 44.51];
bBoxVehicles = [-Inf,-Inf,Inf,Inf];
 
switch scenario
    case 1
        % Vehicles in the entire city of Porto
        vehiclesFile = 'inputMobilitySUMO/c6r.xml';
        % Buildings in entire city of Porto
         staticFile = 'inputPolygon/map.osm';  
         RSUFile = 'inputRSU/RSULocations-1-3-10.txt';
         V2XNames = {'V2I'};
    case 2
        % Vehicles in one quarter of Porto (bottom left quarter)
        vehiclesFile = 'inputPolygon/vehiclesPortoQuarter';
        % Buildings in entire city of Porto
        staticFile = 'inputPolygon/buildingsPortoQuarter';
    case 3
        % SUMO simulation in the city of Cologne, Germany
        vehiclesFile = 'inputMobilitySUMO/CologneSUMOMobility.xml';
        bBoxVehicles = [6.93,50.913,7,50.96];        
        % Buildings in part of Cologne (OpenStreetMap .osm)
        staticFile = 'inputPolygon/CologneCityPart.osm';
    case 4
        % SUMO simulation in the city of Newcastle, UK
        vehiclesFile = 'inputMobilitySUMO/NewcastleSUMOMobility.xml';
        % Buildings in part of Newcastle, UK (OpenStreetMap .osm)
        staticFile = 'inputPolygon/newCastle.osm';
    case 5
        % SUMO simulation in the city of Bologna, Italy (suitable for V2I
        % simulation combined with Bologna RSUs)
        % vehiclesFile = 'inputMobilitySUMO/bolognaVehicles50sec.xml';
        vehiclesFile = 'inputMobilitySUMO/bolognaVehicles25sec.xml';
        bBoxVehicles = [11.3259 44.4829 11.3629 44.51];
        % Buildings in part of Bologna, Italy (OpenStreetMap .osm)
        staticFile = 'inputPolygon/bolognaRightHalfRSU2-4-5-6-7-8-9.osm';
        % Locations and characteristics of 15 RSUs in the city of Bologna, Italy.
        % Details available in Gozalvez et al.:
        % http://ieeexplore.ieee.org/xpls/abs_all.jsp?arnumber=6194400        
        RSUFile = 'inputRSU/RSULocations-2-4-5-6-7-8-9.txt';  
        
        % OR       
        % (make sure to delete vehicle pre-processed file when shifting
        % between two parts of Bologna!)
        % bBoxVehicles = [-Inf -Inf 11.3259 44.51];
        % RSUFile = 'inputRSU/RSULocations-1-3-10.txt';
        % staticFile = 'inputPolygon/bolognaLeftHalfRSU1-3-10.osm';
        
        % Names of link types
        V2XNames = {'V2V', 'V2I'};
end

% Load the file containing outlines of foliage. If foliage is contained in
% staticFile above, use 'inputPolygon/emptyStaticFile' as foliageFile.
foliageFile = 'inputPolygon/emptyStaticFile';

% Number of rows per vehicle outline (by default, there are 5 sides per
% vehicle outline defined by 5 points, with the first point repeating to
% close the polygon (i.e., 6 rows in total).
numRowsPerVehicle = 6;

% Tolerance (in meters) for simplifying building outlines. Useful if
% building outlines are overly complex. For more info, see:
% http://www.mathworks.com/matlabcentral/fileexchange/21132-line-simplification
toleranceBuildings = 0;

%{ 
% Which model is used for NLOSv links (models contained in LOSNLOSv.m):
%   1: Simple model: add attenuation based on the number of obstructing
%   vehicles, irrespective of their dimensions. NB: requires addNLOSvLoss
%   variable to be set up.
%   2: Bullington method for knife-edge diffraction ("Radio propagation for
%   vehicular communications", K. Bullington, IEEE Transactions on
%   Vehicular Technology, Vol. 26, Issue 4): Top and side diffraction are
%   accounted for.
%   3: Multiple knife-edge based on ITU-R method (Propagation by
%   diffraction, Recommendation P.526, Feb. 2007, available at
%   http://www.itu.int/rec/R-REC-P.526-13-201311-I/en). Top and side
%   diffraction are accounted for.
%
% Details on how the NLOSv models are implemented are available in: "Impact
% of vehicles as obstacles in vehicular ad hoc networks", Boban et al.,
% IEEE Journal on Selected Areas in Communications, Vol. 29, Issue 1.
%}
NLOSvModelType = 3;

% Additional loss (in dB) due to obstructing vehicles (min, med, max).
% Applies for NLOSvModelType = 1 only (for details, see LOSNLOSv.m).
addNLOSvLoss = [2 6 10];

%% Parameters related to communicating vehicle pairs
% Array containing communicating vehicle pairs. NB1: vehicles
% in pairs are fetched based on the row ID (e.g., if [4,10] is a pair, 4th
% and 10th vehicle in the vehicle list are selected. NB2: if same pairs are
% used for all time-steps, then perform "repmat" on the array.
commPairArray = 'inputPolygon/emptyCommPairFile';
commPairArrayV2I = 'inputPolygon/emptyCommPairFile';

% In case commPairArray is empty, numCommPairs defines how many
% communication pairs are generated. If this number is larger than the
% total number of pairs in the system, then it is replaced by the latter.
% Conversely, if you want to analyze all possible communication pairs,
% simply set numCommPairs to Inf;
numCommPairs = 50;
numCommPairsV2I = Inf;

%% Dimensions of vehicles (in meters)
% Mean and standard deviation for normally distributed dimensions
carMeanHeight = 1.5;
carStdDev =.1;
truckMeanHeight = 3; 
truckStdDev = .1;

carMeanWidth = 1.75;
carStdDevWidth = .1;
truckMeanWidth = 2;
truckStdDevWidth = .1;

carMeanLength = 5;
carStdDevLength = .4;
truckMeanLength = 9;
truckStdDevLength = .6;

% Putting vehicle dimensions in an array for more compact argument passing 
vehDimensionParams = ...
    [carMeanHeight     carStdDev        truckMeanHeight     truckStdDev;...
     carMeanWidth      carStdDevWidth   truckMeanWidth      truckStdDevWidth;...
     carMeanLength     carStdDevLength  truckMeanLength     truckStdDevLength];

% For vehicles not generated by the simulation environment, any vehicle
% with length larger than lengthThreshold (in meters) is assumed to be a
% tall vehicle (e.g., a truck)
lengthThreshold = 8;

%% Max. communication range for LOS, NLOSv, and NLOSb (commRange = LOSRange)
commRange = 500; 
maxNLOSvRange = 400;
maxNLOSbRange = 300;

%% Propagation-related parameters
% Antenna height (in meters, added atop vehicle height); 
antennaHeight = .1;

% Antenna polarization: 0 - vertical; 1 - horizontal
polarization = 0;

% Operating frequency (in Hz)
freq = 5.89*10^9;

% Transmit power (in dBm) for vehicles. RSU input file contains Tx power
% for each of the RSUs. NB: if you want to use different Tx power per
% vehicle, RSU input file can serve as a starting point
txPower = 12;

% Vehicle and RSU antenna gains (in dBi)
antennaGain.veh = 1;
antennaGain.RSU = 1;

% Path Loss Exponent (PLE) for NLOSb/f links
PLENLOSb.V2V.NLOSb = 2.9;
PLENLOSb.V2V.NLOSf = 2.7;
PLENLOSb.V2I.NLOSb = 2.9;
PLENLOSb.V2I.NLOSf = 2.7;

%% Small-scale signal variation parameters
% (taken from the measurements described in the paper "Geometry-Based
% Vehicle-to-Vehicle Channel Modeling for Large-Scale Simulation". Details
% on the use of parameters below are described in the same paper). Other
% environments can have different min and max small scale fading.
%% V2V 
% Minimum std. dev. of small scale variation for LOS and NLOSv/b links
smallScale.V2V.minFadingLOS = 3.3;
smallScale.V2V.minFadingNLOSv = 3.8;
smallScale.V2V.minFadingNLOSb = 4.1;
smallScale.V2V.minFadingNLOSf = smallScale.V2V.minFadingNLOSb;
% Maximum std. dev. of small scale variation for LOS and NLOSv/b links
smallScale.V2V.maxFadingLOS = 5.2;
smallScale.V2V.maxFadingNLOSv = 5.3;
smallScale.V2V.maxFadingNLOSb = 6.8;
smallScale.V2V.maxFadingNLOSf = smallScale.V2V.maxFadingNLOSb;

%% V2I; calculated based on Bologna V2I measurements by Gozalvez et al.: 
% http://ieeexplore.ieee.org/xpls/abs_all.jsp?arnumber=6194400

% Minimum std. dev. of small scale variation for LOS and NLOSv/f/b links
smallScale.V2I.minFadingLOS = 1;
smallScale.V2I.minFadingNLOSv = 1;
smallScale.V2I.minFadingNLOSf = 1.5;
smallScale.V2I.minFadingNLOSb = 1.9;
% Maximum std. dev. of small scale variation for LOS and NLOSv/f/b links;
smallScale.V2I.maxFadingLOS = 3.4;
smallScale.V2I.maxFadingNLOSv = 4.2;
smallScale.V2I.maxFadingNLOSf = 3.3;
smallScale.V2I.maxFadingNLOSb = 4.7;

% Minimum radius (in meters) based on which the small scale signal
% variation model determines maximum vehicle density. The specific value is
% a guesstimate.
minDensityRange = 250;


%% Experiment-related parameters
% Determines if the simulation should perform comparison against
% measurement data (experiments = 1) or not (experiments = 0). 
% NB: measurement data needs to be separated into LOS, NLOSv, and NLOSb.
experiments = 0;

%% Experiment type: 0-V2V; 1-V2I
experimentType = 1;
if experimentType~=0 && experimentType~=1
    error('Wrong experiment type');
end

% Location of the measurement trace (NB: it should be used with appropriate
% vehicle and building/foliage outlines)
%experimentFile = 'inputExperiments/cc_high_ceed_sum_noHeader.txt';
%experimentFile = 'inputExperiments/bolognsRSU3-3-P20-h6.5_MB';
experimentFile = 'inputExperiments/bolognaRSU9-P20-h6.5.txt';

% Actual height (in meters) of the transmitting (Tx) and receiving (Rx)
% vehicle used in the experiments in Porto (located in
% inputExperiments/measurementsPorto). -Inf in the third column is used for
% distinguishing the real-world TxRxHeight in the simulation
TxRxHeight = [6.5, 1.5 -Inf];


%% Reflection- and diffraction-related parameters
% Determines if reflections and diffractions are used (useReflDiffr=1) or
% not (useReflDiffr=0) If not, only log-distance path loss is used for
% NLOSb links.
useReflDiffr=0;

% If reflections are used 
if useReflDiffr
    % Minimum fading for NLOSb data is set to 0 in case reflections and
    % diffractions are used. Detailed explanation is available in the paper
    minFadingNLOSb = 0;
end

% Relative permittivity of materials
% (http://en.wikipedia.org/wiki/Relative_permittivity).
% Buildings (concrete)
buildingReflRelPerm = 4.5;
% Vehicles (rough approximation: mix of metal, glass, rubber)
vehReflRelPerm = 6; 

%% Google Earth Visualization 
% NB: in case of large networks, generating the visualization can take very
% long time. Use accordingly.
GEVisualize = 1;
    % IF GEVisualize = 1, then parameters below indicate what is plotted.
    % Plot static objects and vehicle polygons
    plotPoly = 0;
    % Plot the received power for each of the communications pairs
    plotRxPwr = 1;
    % Plot the number of neighbors per vehicle
    plotNeighbors = 1;
    
%% Output-related parameters
% Verbose textual output in console (verbose = 1)
verbose = 1;