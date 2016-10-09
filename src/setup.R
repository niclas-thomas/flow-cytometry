list.of.packages <- c("ggplot2", "openCyto","flowMeans")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)