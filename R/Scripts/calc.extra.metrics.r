

# Nicholas Ducharme-Barth
# 05/02/2020
# calculate extra metrics
# Mean absolute percentage error, mean percentage error, mean forecast error, avg(est/pred), bias coefficient
# only calculate metrics for models/regions/replicates where other metrics were calculated

# load packages
	library(data.table)

# set project directory
	project.dir = "C:/Users/nicholasd/HOME/SPC/SPC_SAM/Geostats/spatial.sampling.sim.manuscript/"
	setwd(project.dir)

# bring in metric.df
	load("Index/ResultsDF/metric.df.RData")
	load("Index/ResultsDF/ts.df.RData")
	load("SimData/simple.true.index.RData")
	ts.dt = as.data.table(ts.df)
	setkey(ts.dt,Scenario,Catchability,Replicate,Model,Region)

# check equivalence
	# plot(ts.dt[Scenario=="Rotating"&Catchability=="Q"&Replicate==61&Model=="Enviro"&Region==3]$Index)
	# load("Index/Simple120/Rotating_Q/vast/61.vast_list.RData")
	# lines(vast_list$vast_output$Enviro$idx[,5]/mean(vast_list$vast_output$Enviro$idx[,5]))
	# mean(ts.dt[Scenario=="Rotating"&Catchability=="Q"&Replicate==61&Model=="Enviro"&Region==3]$Index == vast_list$vast_output$Enviro$idx[,5]/mean(vast_list$vast_output$Enviro$idx[,5]))

# identify models/regions/replicates to calculate new metrics
	unique.dt = as.data.table(metric.df)
	unique.dt = unique(unique.dt[!is.na(mgc),.(Scenario,Catchability,Replicate,Model,Region,mgc)])

# define new metric storage structure
	new.metrics = matrix(NA,nrow=nrow(unique.dt),ncol=4)
	colnames(new.metrics) = c("MAPE", "MPE", "bias.simple", "bias.coefficient")

# define bias function
	bias.coefficient = function(est,true)
	# pass two indices
	# http://kourentzes.com/forecasting/2014/12/17/measuring-the-behaviour-of-experts-on-demand-forecasting-a-complex-task/
	{
		est = scale(est)
		true = scale(true)
		error = est - true
		complex.error.vec = as.complex(error)
		z.vec = sqrt(complex.error.vec)
		MRE = mean(z.vec,na.rm=TRUE)
		a.real = Re(MRE)
		b.im = Im(MRE)
		r.polar = sqrt(a.real^2 + b.im^2)
		gamma.polar = if(a.real>0){atan(b.im/a.real)}else if(a.real == 0 & b.im>0){pi/2}else{pi/4}
		bias.k = 1 - 4*gamma.polar/pi
		return(bias.k)
	}

# iterate across unique scenario/catchability/model/region/replicate
	for(i in 1:nrow(new.metrics))
	{	
		# A  = proc.time()
		# extract estimated index
			est = ts.dt[.(unique.dt$Scenario[i],unique.dt$Catchability[i],unique.dt$Replicate[i],unique.dt$Model[i],unique.dt$Region[i])]$Index
		# get true index
			true = simple.true.index[,unique.dt$Region[i]]
		# calc new metrics
		# positive sign means est>true and negative means true>est
			new.metrics[i,"MAPE"] = 100*mean(abs((est-true)/true))
			new.metrics[i,"MPE"] = 100*mean((est-true)/true)
			new.metrics[i,"bias.simple"] = mean((est/true))
			new.metrics[i,"bias.coefficient"] = -bias.coefficient(est,true)
		# clean-up
			rm(list=c("est","true"))
		# B = proc.time()
		# B - A
	}

# append new metrics to metric.df


# save
