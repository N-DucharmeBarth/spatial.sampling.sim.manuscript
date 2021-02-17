

# Nicholas Ducharme-Barth
# 05/02/2020
# calculate extra coverage metrics
# cover 60, 70, 80, 90, 95

# load packages
	library(data.table)

# set project directory
	project.dir = "C:/Users/nicholasd/HOME/SPC/SPC_SAM/Geostats/spatial.sampling.sim.manuscript/"
	setwd(project.dir)

# bring in metric.df
	load("Index/ResultsDF/metric.MRwR.df.RData")
	load("Index/ResultsDF/ts.df.MRwR.RData")
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
	colnames(new.metrics) = c("cover.60","cover.70","cover.80","cover.90","cover.95")

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
			if(mean(is.na(est.se))==0)
			{
				COVER  = rep( 0, len = length( true ) )
				UCI = est + qnorm(0.5*(1-0.6),lower.tail=FALSE)*est.se 
				LCI = est - qnorm(0.5*(1-0.6),lower.tail=FALSE)*est.se 
				for ( j in 1:length( COVER ) ) {
					if ( ( LCI[j]<=true[j] ) & ( true[j]<=UCI[j]) ){ COVER[j] = 1} 	
				}
				
				new.metrics[i,"cover.60"] = ( sum( COVER )/ length( COVER ) ) * 100
				rm(list=c("COVER","UCI","LCI"))
			}
			if(mean(is.na(est.se))==0)
			{
				COVER  = rep( 0, len = length( true ) )
				UCI = est + qnorm(0.5*(1-0.7),lower.tail=FALSE)*est.se 
				LCI = est - qnorm(0.5*(1-0.7),lower.tail=FALSE)*est.se 
				for ( j in 1:length( COVER ) ) {
					if ( ( LCI[j]<=true[j] ) & ( true[j]<=UCI[j]) ){ COVER[j] = 1} 	
				}
				
				new.metrics[i,"cover.70"] = ( sum( COVER )/ length( COVER ) ) * 100
				rm(list=c("COVER","UCI","LCI"))
			}
			if(mean(is.na(est.se))==0)
			{
				COVER  = rep( 0, len = length( true ) )
				UCI = est + qnorm(0.5*(1-0.8),lower.tail=FALSE)*est.se 
				LCI = est - qnorm(0.5*(1-0.8),lower.tail=FALSE)*est.se 
				for ( j in 1:length( COVER ) ) {
					if ( ( LCI[j]<=true[j] ) & ( true[j]<=UCI[j]) ){ COVER[j] = 1} 	
				}
				
				new.metrics[i,"cover.80"] = ( sum( COVER )/ length( COVER ) ) * 100
				rm(list=c("COVER","UCI","LCI"))
			}
			if(mean(is.na(est.se))==0)
			{
				COVER  = rep( 0, len = length( true ) )
				UCI = est + qnorm(0.5*(1-0.9),lower.tail=FALSE)*est.se 
				LCI = est - qnorm(0.5*(1-0.9),lower.tail=FALSE)*est.se 
				for ( j in 1:length( COVER ) ) {
					if ( ( LCI[j]<=true[j] ) & ( true[j]<=UCI[j]) ){ COVER[j] = 1} 	
				}
				
				new.metrics[i,"cover.90"] = ( sum( COVER )/ length( COVER ) ) * 100
				rm(list=c("COVER","UCI","LCI"))
			}
			if(mean(is.na(est.se))==0)
			{
				COVER  = rep( 0, len = length( true ) )
				UCI = est + qnorm(0.5*(1-0.95),lower.tail=FALSE)*est.se 
				LCI = est - qnorm(0.5*(1-0.95),lower.tail=FALSE)*est.se 
				for ( j in 1:length( COVER ) ) {
					if ( ( LCI[j]<=true[j] ) & ( true[j]<=UCI[j]) ){ COVER[j] = 1} 	
				}
				
				new.metrics[i,"cover.95"] = ( sum( COVER )/ length( COVER ) ) * 100
				rm(list=c("COVER","UCI","LCI"))
			}

		# clean-up
			rm(list=c("est","true","est.se"))
	}
	B = proc.time()
	B - A

# append new metrics to metric.df
	unique.dt = rbind(unique.dt,unique.dt,unique.dt,unique.dt,unique.dt)
	unique.dt$Metric = c(rep("cover.60",nrow(new.metrics)),rep("cover.70",nrow(new.metrics)),rep("cover.80",nrow(new.metrics)),rep("cover.90",nrow(new.metrics)),rep("cover.95",nrow(new.metrics)))
	unique.dt$Value = c(new.metrics[,"cover.60"],new.metrics[,"cover.70"],new.metrics[,"cover.80"],new.metrics[,"cover.90"],new.metrics[,"cover.95"])
	unique.dt = unique.dt[,.(Scenario,Catchability,Replicate,Model,Region,Metric,Value,mgc)]
	# remove previously calculated "new metrics"
	"%ni%" = Negate("%in%")
	metric.df = as.data.table(metric.df)
	metric.df = metric.df[Metric %ni% c("cover.60","cover.70","cover.80","cover.90","cover.95")]
	metric.df = as.data.frame(rbind(metric.df,unique.dt))

# remove duplicate rows
	metric.df = unique(metric.df)

# save
	save(metric.df,file="Index/ResultsDF/metric.MRwR.df.RData")

