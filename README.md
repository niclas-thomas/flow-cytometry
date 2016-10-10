# flowCytometry

Suite of functions and pipelines to perform automated gating and clustering of flow cytometry data.

## Running flowCytometry

1. Download R [here](http://star-www.st-andrews.ac.uk/cran/) and RStudio [here](https://www.rstudio.com/products/rstudio/download3/) and install them both in that order.
2. Download ```flowCytometry``` and unzip the folder.
3. Open RStudio.
3. Open and run ```src/setup.R``` to install the required packages.
4. Check that the gating strategy in ```gatingtemplate_bcell.csv``` in the folder ```gatingTemplates``` is suitable for your flow experiment.
5. Manually change the variables ```myDataDir``` and ```myCodeDir``` in ```pipeline/autoGating.R``` to reflect your particular system and run this script.
6. Manually change the variables ```myDataDir``` and ```myCodeDir``` in ```pipeline/autoClustering.R``` to reflect your particular system and run this script.
7. Manually change the variables ```myDataDir``` and ```myCodeDir``` in ```pipeline/pca.R``` to reflect your particular system and run this script.

## Notes

Use ```plot.gates``` and ```plot.strategy``` in ```pipeline/autoGating.R``` to define whether you want to plot gating strategy and sample gated plots. Note that plotting gated plots is quite
RAM intensive, so on lower spec machines it may not be possible to use this feature. If so, ensure that ```plot.gates==FALSE```.

Steps 6, 7 and 8 will create new folders in ```myDataDir``` to save the results to. Once the the script at each step has finished, output will be saved in subfolders within ```myDataDir```. Figures will be saved to the ```figures``` folder.