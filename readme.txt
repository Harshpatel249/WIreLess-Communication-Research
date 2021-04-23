# Generating results  
* Download: GEMV2 version 1.2 from [GEMV^2](http://vehicle2x.net/download/)

* Run: GEMV2PackageV1.2/runSimulation.m (After replacing the downloaded with file from the repository)

* Results available in following folders:  
KML files: GEMV2PackageV1.2/OutputKML  
csv files: GEMV2PackageV1.2/OutputSim  

----------------------------------------------------------------

# Edited matlab files:

* GEMV2PackageV1.2/simSettings.m

* Edited the file to point to the required simulation files

----------------------------------------------------------------

# Simulation files:

* Mobility file (Vehicle Mobility simulation from SUMO):  
Folder: GEMV2PackageV1.2/inputMobilitySUMO/SUMOMobility.xml

How to create:
1. Create the .sumocfg for the required area of simulation using the SUMO-OSM-Web-Wizard. [SUMO](https://sumo.dlr.de/docs/Downloads.php)
2. Convert the .sumocfg file to the mobility .xml file using the following command:
```
sumo.exe -c <InputFilename>.sumocfg --fcd-output 
<OutputMobilityFilename>.xml --fcd-output.geo true --netstate-dump <OutputMobilityFilename>.net.xml 
```
3. The output xml file contains thousands of time-steps, however we need to use around 50 time-steps so we edit the xml file accordingly. You can open
the file specified in the folder section to know more about it.

4. The final step is to copy the .xml file to GEMV2PackageV1.2/inputMobilitySUMO/ and edit the simSettings.m accordingly (vehiclesFile variable).

* Input-Polygon file (The actual map):  
Folder: GEMV2PackageV1.2/inputPolygon/map.osm

How to create:
1. Go to openstreetmap.org and select the area you want to simulate for and generate the .osm file.
2. Copy the .osm file to the specified folder (GEMV2PackageV1.2/inputPolygon/) and edit the simSettings.m accordingly (staticFile variable).

* Input-RSU file (The co-ordinates for Road-Side Units):
Folder: GEMV2PackageV1.2/inputRSU/RSULocations-1-3-10.txt

How to create:
1. Open Google Maps, go to the area you are simulating for and select any random point to place and RSU. Copy the co-ordinates of the point and place
in the specified file(GEMV2PackageV1.2/inputRSU/RSULocations-1-3-10.txt) and edit the simSettings.m accordingly (RSUFile variable).

--------------------------------------------------------------------

# Extra notes:

* simSettings.m has to be set such that it simulates the scenario that we want (it has default loaded scenarios).
    For this, we have to set the variable 'scenario = 5', if we want to simulate the 5th scenario in the switch case.
    So, now we have to go ahead and set the variables for the 5th scenario to point to our simulation files as mentioned in above steps for each file.

* The output can be visualized using the google earth pro (free to download) and load the .kml files into it.