# flowCytometry

Suite of functions and pipelines to perform automated gating and clustering of flow cytometry data.

## Running flowCytometry

1. Manually change the variables ```myDataDir``` and ```myGitDir``` for your particular system.
2. Check that the gating template in the folder gatingTemplates is suitable for your flow experiment.
3. Run pipeline scripts in the order ```autoGating.R```, ```autoClustering.R```, ```pca.R```.