# compound-eye-plotting-elife
Code for plotting visual parameters calculated from compound eye volumes, as published in Taylor et al. 2019 eLife

To perform the calculations described in the publication, enter filepaths in ConfigForBeeCEComp.m and then run the script CompAnalysisForBeeCE.m 
The config file already contains metadata for the bees analysed in the publicaiton, and the data used for analysis are labelled image stacks availible from Morphosource and facet measruments in csv files availible from Datadryad. New data can also be used by modifying the Config file appropriately. The analysis scripts are not well documented, but options exist to modify analysis parameters and the results can be saved for plotting.

The results clacualted in the publications, please download the data from (Data_for_plots.zip) from https://dx.doi.org/10.5061/dryad.23rj4pm
Run PlottingCodeForBeeCE.m to plot matching (but unformatted) plots to those published in the paper using the data file. A variety of options, either at the start of the script or through, can be used to modify the plotting results. 

All other scripts are dependencies.

Several scripts in this repository were obtained from the Mathworks File Exchange or the Mathworks website. These files are dependencies for the scripts I wrote, although some have been modified. The original authors of those files may have given the scripts different licences from the MIT license that covers this repository, and those licenses may take precedence over the license of this repository. 

The scripts I wrote for this project are indicated at the start of each file and are covered by the MIT licence. Please contact me if you want access to the Matlab scripts used to calculate the visual paramteres from the volume information.

Contact gavin.taylor.01@gmail.com for assistance.
