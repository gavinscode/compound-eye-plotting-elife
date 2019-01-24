# compound-eye-plotting-elife
Code for calculating and plotting visual parameters from compound eye volumes, as published in Taylor et al. 2019 eLife
https://doi.org/10.1101/380527 

To perform the calculations described in the publication, enter filepaths in ConfigForBeeCEComp.m and then run the script CompAnalysisForBeeCE.m. The config file already contains metadata for the bees analysed in the publication, and the data used for analysis are labelled image stacks from Morphosource and facet measurements in csv files from DataDryad. New data can also be used by modifying the Config file appropriately. The analysis scripts are not well documented, but options exist to modify analysis parameters and the results can be saved for plotting.

The results calculated in the publications can be plotted either by running the calculation scripts (and saving the results), or downloading the calculated data from DataDryad. Run PlottingCodeForBeeCE.m to plot matching (but unformatted) plots to those published in the paper using a data file. A variety of options, either at the start of the script or throughout, can be used to modify the plotting results. 

All other scripts are dependencies. Several scripts in this repository were obtained from the Mathworks File Exchange or the Mathworks website. These files are dependencies for the scripts I wrote, although some have been modified. The original authors of those files may have given the scripts different licences from the MIT license that covers this repository, and those licenses may take precedence over the license of this repository. 

The files written for this project are with the name ‘Gavin Taylor’ at the start of each script and are covered by the MIT licence. 
Contact gavin.taylor.01@gmail.com for assistance.

Resources:
DataDryad: https://dx.doi.org/10.5061/dryad.23rj4pm
MorphoSource: https://www.morphosource.org/Detail/ProjectDetail/Show/project_id/646


