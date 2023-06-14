******CRA Analysis Package (V1) ************
This package does the CRA Analysis 
Tested with R version - 3.4.3
Input file format tested: netcdf (Both observation and model forecast have to be on the same grids)

File system of CRA_package_V1.tar.gz:
obs_20210517.nc          : Obs test data
model_20210517.nc        : Model test data
Isolate_Objects.R        : Function to identify and pair objects from both model and observation based on rainfall threshold. This will provide each objects in separate NetCDF file. 
CRA_Err_Decomp_Features.R: Main function for CRA analysis
Func_Run_CRA_Analysis.sh : This is a shell script which creates a R script for creating the identified paired object files and does CRA Analysis. It also proved a pdf which shows the paired objects and snapshots of the movement of the model raster over CRA area
Run_CRA_Plot.sh : Creates figure including the necessary statistic tables
###################
User inputs for Func_Run_CRA_Analysis.sh

date=????????    #### date in the format of yyyymmdd
RUNPATH=         #### Path for rundirectory
ActINPUT_OBS=    #### path for observation file
ActINPUT_MODEL=  #### path for model file
OUT_PATH=        #### output path 
Threshold=       #### threshold rainfall for CRA Analysis 
case=            #### Experiment name
OUT_PREFIX=      #### output file prefix
Ngrids=6         #### No of grids to be shifted which is required for the RMSE minimization; 
			at present it is kept as 6 which means the model raster will be moving 
			6 grids right, left, down, and up over the CRA area. 


###################
User inputs for Run_CRA_Plot.sh

date=????????    #### date in the format of yyyymmdd
RUNPATH=         #### Path for rundirectory
ActINPUT_OBS=    #### path for observation file
ActINPUT_MODEL=  #### path for model file
OUT_PATH=        #### output path 
Threshold=       #### threshold rainfall for CRA Analysis 
case=            #### Experiment name
OUT_PREFIX=      #### output file prefix
Ngrids=6         #### No of grids to be shifted which is required for the RMSE minimization; 
ObjNum=          #### serial number of the feature analysed 
Rainlimt=        #### rainfall limit for colorbar
RainlimtScatter  #### rainfall limit for scatter plot 
shpfile=         #### path of shp file 
figname=         #### name of the figure to be saved 
INPUT_OBS=       #### observation file which includes only the rainfall object
INPUT_MODEL=     #### model file  which includes only the rainfall object

lat1,lat2,lon1 and lon2 are the lat-lon limits for the plotting purpose
###################

