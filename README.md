# coastsat_matlab
Contains my Matlab code for working with coastlines provided by Sean Vitousek

This repo resides in ../src/coastsat_matlab

Right now, this code is intended to work on the provisional data provided by Sean Vitousek named `transects_with_shorelines.mat`  

### Data files
`transects_with_shorelines.mat` - Too big to put in the repo. Contact Sean for a copy.  
`littoral_cells.csv` - List of first transect number for each littoral cell  

### Key m-files  
`plot_coastlines.m` - First bit of code to wrangle the transect format and do some analysis  
`littoral_cell_indices.m` - Extract starting transect number for each littoral cell. Write a list in `littoral_cells.csv`  

### Helper scripts
`load_colors.m` - Load some pretty colors grabbed from ColorBrewer  
`lsfit.m` - Ancient least-squares fit routine  
`transects2utm.m` - Code to convert the along-transect `Y` value to UTM  
`utm2deg.m` - Convert UTM to lat/lon  
`wgs2utm.m` - Convert lat/lon to UTM  

