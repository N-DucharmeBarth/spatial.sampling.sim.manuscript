

# Nicholas Ducharme-Barth
# 03/02/2020
# plot the results

# set working directory
	project.dir = "C:/Users/nicholasd/HOME/SPC/SPC_SAM/Geostats/spatial.sampling.sim.manuscript/"
	setwd(project.dir)

# source
	source("R/Fn/fn.plot.sim.results.summary.r")

# load
	load("SimData/simple.true.index.RData")
	load("Index/ResultsDF/ts.df.RData")
	load("Index/ResultsDF/metric.df.RData")
	load("Index/ResultsDF/nominal.df.RData")


# iterate across combinations and plot metrics
	for(q in c("Q","noQ"))
	{
		for(r in c("all",1:8))
		{
			plot.bias(metric.df,r=r,q=q,Save=TRUE,Save.Dir="Plots/Metric/")
			plot.rmse(metric.df,r=r,q=q,Save=TRUE,Save.Dir="Plots/Metric/")
			plot.mae(metric.df,r=r,q=q,Save=TRUE,Save.Dir="Plots/Metric/")
			plot.cover(metric.df,r=r,q=q,Save=TRUE,Save.Dir="Plots/Metric/")
		}
	}

# iterate across combinations and plot time series: 1 & 2
	for(q in c("Q","noQ"))
	{
		for(r in c("all",1:8))
		{
			for(s in c("Contraction", "Expansion", "Fixed", "Preferential", "Random", "Rotating"))
			{
				for(m in c("Enviro","NoEnviro","EnviroSVC","NoEnviroSVC"))
				{
					plot.ts(ts.df,nominal.df,simple.true.index,plot.type=1,r=r,s=NA,q=q,m=m,Save=TRUE,Save.Dir="Plots/Index/")
					plot.ts(ts.df,nominal.df,simple.true.index,plot.type=2,r=NA,s=s,q=q,m=m,Save=TRUE,Save.Dir="Plots/Index/")
				}
			}
		}
	}
		


