# Footprint-Weighed-Flux-Maps
This code is used to create Footprint-Weighted Flux Maps around eddy-covariance sites.
Data input: Standard Ameriflux csv file or Standard Berkeley Biomet Lab mat files

## Initial Steps

•	Clone the “Footprint-Weighted-Flux-Maps” folder from Github in your local machine.

•	Select a folder in your local machine for L4 processing output and note it for the next step.

•	If you have not done so already, download the Ameriflux csv file for your site of interest. For speed and efficiency, you may upload a subset of the original file with the times of interest. Uploading the whole file should be fine as well but initial L4_Processing may run a bit slower.

## L4_Processing
•	Open the file L4_Process_from_AFcsv.m in the L4_processing folder to create an L4 file. For footprint processing there are multiple variables that are calculated in what I called “L4 files” These L4 files are Matlab structures with the Ameriflux csv data plus additional variables needed to run footprint including Aerodynamic canopy height and PBL height (PBL is optional).

•	Enter options
* Modify the Save Directory (from previous step)
* Select the site acronym for the site (3 letter code from Ameriflux recommended) 
* Update the filename and file path for the csv file
* Adjust parameters in BADM. Canopy height is an initial estimate. Canopy height will be calculated using the aerodynamic canopy height approach of Pennypacker and Baldocchi later on
* Adjust the variable names in data if needed. Some variable may have slightly different names (for example AF_RH_1_1_1 vs AF_RH_1_1_2). Non-existent variables can be populated with nan { e.g. data.mbar= nan(size(AF.Mdate)) }
* If you have water level above the surface and have data available update this in the water level section, otherwise disregard this step
* Run Aerodynamic canopy height. You can adjust the days to average if your dataset is shorter. 30 days is recommended for datasets longer than 1 year.
* Add PBL data if available. If not an estimate of 750 m will be used for all the footprints. Non-published tests inidicate the sensitivity of footprint models to PBL height is very low

## Calculating Footprints
•	Go to the CalculateFootprint folder, Open FP_opts_continuous.m, select the site, and update the source of the L4 file from the previous step.

•	Section How: If you only need one model, you can select only the best model (the K&M model) by adjusting opt.models=[0 0 1];

•	In section “When?” select daytime or nighttime footprints. Select the start and the end times of the run in opt.start and opt.end2

•	Select footprint contours to plot, you can choose more than 2, for example [50,60,70,80]

•	Input fluxes. Adjust the fluxes for footprint-weighted flux by choosing those fluxes of interests that are available in your dataset for which you wish to create footprint-weighted flux maps. 

•	Saving. Select a folder to Save the footprint output. For consistency this folder should be of the type '..\Footprint_Output\continuous\'. Continuous indicates that the data is being aggregated on continuous datasets. 

•	Within the continuous folder a folder named ”kml” should be created. Here all the resulting contours in kml files will be created.

•	Once all the options are ready run the code Footprint_Run_Continuous.m (click “run” or press F5)

•	To plot the resulting footprint contours in kml files that can be then seen in Google Earth, the GoogleEarth toolbox is added at the end of the Footprint_Run_Continuous.m

## To create Footprint-Weighted Flux Maps

•	Go to the folder FW_FluxMaps

•	Open the file Surface_Flux_Map and update the DataDir with the output from the previous step.

•	Modify the file to load with the output from the previous step.

•	Make any necessary changes to the fluxes to be plotted. By default, CH4, CO2 and H2O are provided.

•	Run the code. You should be able to see your initial footprint-weighted flux maps.
