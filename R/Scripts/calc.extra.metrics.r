

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
	new.metrics = matrix(NA,nrow=nrow(unique.dt),ncol=5)
	colnames(new.metrics) = c("MAPE", "MPE", "bias.simple", "bias.coefficient","cover.50")

# define bias function
	bias.coefficient = function(est,true)
	# pass two indices
	# http://kourentzes.com/forecasting/2014/12/17/measuring-the-behaviour-of-experts-on-demand-forecasting-a-complex-task/
	# https://github.com/trnnick/TStools/blob/master/R/mre.R
	{
		e = est - true
		e = as.complex(e)
		e = sqrt(e)
		mre = mean(e)
		gamma = Arg(mre)
	  	bias = 1 - 4*gamma/pi
		return(bias)
	}



# iterate across unique scenario/catchability/model/region/replicate
	A  = proc.time()
	for(i in 1:nrow(new.metrics))
	{	
		# extract estimated index
			est = ts.dt[.(unique.dt$Scenario[i],unique.dt$Catchability[i],unique.dt$Replicate[i],unique.dt$Model[i],unique.dt$Region[i])]$Index
			est.se = ts.dt[.(unique.dt$Scenario[i],unique.dt$Catchability[i],unique.dt$Replicate[i],unique.dt$Model[i],unique.dt$Region[i])]$SE

		# get true index
			true = simple.true.index[,unique.dt$Region[i]]
		# calc new metrics
			new.metrics[i,"MAPE"] = 100*median(abs((est-true)/true))
			new.metrics[i,"MPE"] = 100*median((est-true)/true)
			new.metrics[i,"bias.simple"] = median((est/true))
			new.metrics[i,"bias.coefficient"] = bias.coefficient(est,true) # positive sign means est>true and negative sign means est<true
			if(mean(is.na(est.se))==0)
			{
				COVER  = rep( 0, len = length( true ) )
				UCI = est + 0.67499*est.se 
				LCI = est - 0.67499*est.se 
				for ( j in 1:length( COVER ) ) {
					if ( ( LCI[j]<=true[j] ) & ( true[j]<=UCI[j]) ){ COVER[j] = 1} 	
				}
				
				new.metrics[i,"cover.50"] = ( sum( COVER )/ length( COVER ) ) * 100
				rm(list=c("COVER","UCI","LCI"))
			}

		# clean-up
			rm(list=c("est","true","est.se",))
	}
	B = proc.time()
	B - A

# append new metrics to metric.df
	unique.dt = rbind(unique.dt,unique.dt,unique.dt,unique.dt)
	unique.dt$Metric = c(rep("MAPE",nrow(new.metrics)),rep("MPE",nrow(new.metrics)),rep("bias.simple",nrow(new.metrics)),rep("bias.coefficient",nrow(new.metrics)))
	unique.dt$Value = c(new.metrics[,"MAPE"],new.metrics[,"MPE"],new.metrics[,"bias.simple"],new.metrics[,"bias.coefficient"])
	unique.dt = unique.dt[,.(Scenario,Catchability,Replicate,Model,Region,Metric,Value,mgc)]
	# remove previously calculated "new metrics"
	"%ni%" = Negate("%in%")
	metric.df = as.data.table(metric.df)
	metric.df = metric.df[Metric %ni% c("MAPE","MPE","bias.simple","bias.coefficient","cover.50")]
	metric.df = rbind(metric.df,unique.dt)

# remove duplicate rows
	metric.df = unique(metric.df)

# save
	save(metric.df,file="Index/ResultsDF/metric.df.RData")

