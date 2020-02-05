

# Nicholas Ducharme-Barth
# 05/02/2020
# calculate extra metrics
# Mean absolute percentage error, mean percentage error, mean forecast error, avg(est/pred), bias coefficient
# only calculate metrics for models/regions/replicates where other metrics were calculated

# set project directory
	project.dir = "C:/Users/nicholasd/HOME/SPC/SPC_SAM/Geostats/spatial.sampling.sim.manuscript/"
	setwd(project.dir)

# bring in metric.df
	load("Index/ResultsDF/metric.df.RData")

