

# Nicholas Ducharme-Barth
# 30/12/2019
# Sample the baseline abundance data under 6 different spatial sampling scenarios
# 1) Preferential sampling - fishery dependent
# 2) Fleet expansion - fishery dependent
# 3) Fleet contraction - fishery dependent
# 4) Seasonal closure (fixed spatial area) - fishery dependent
# 5) Seasonal closure (rotating spatial area) - fishery dependent
# 6) Random sampling - fishery independent

# Under the fishery dependent scenarios, samples will also be preferentially allocated to areas of higher underlying biomass
# All scenarios will assume the same level of observation error, a CV of 15%
# Lastly, as a baseline, the random sampling will be run with zero observation error. 

# Each spatial sampling scenario will be defined as a function with the following arguments
# 1) data.dt
# 2) ts.sampled (vector indicating which timesteps to sample between,inclusive)
# 3) n.samp (total # of samples)
# 4) cv (observation error, as a proportion)
# 5) random.seed

##########################################################
# set working directory
	setwd("C:/Users/nicholasd/HOME/SPC/SPC_SAM/Geostats/spatial.sampling.sim.manuscript/")	


##########################################################
# Define external variables for function development
	load("Background_Data/data.dt.RData")
	data.dt = data.table::as.data.table(data.dt)
	# ts.sampled = c(1,40)
	# n.samp = 20000
	# cv = 0.15
	# random.seed = 123

##########################################################
# Define sampling functions
rand.samp = function(data.dt, ts.sampled = c(1,40), n.samp = 20000, cv = 0.15, random.seed = 123)
{
	# set seed
		set.seed(random.seed)

	# subset to desired ts
		valid.ts = ts.sampled[1]:ts.sampled[2]	
		tmp.dt = data.dt[ts %in% valid.ts]

	# define sampling probability
		prob.vector = rep(1,nrow(tmp.dt))

	# sample rows
		tmp.dt = tmp.dt[sample(1:nrow(tmp.dt),size=n.samp,replace=FALSE,prob=prob.vector),][order(ts,lon,lat)]

	# add observation error
		tmp.dt$response = tmp.dt$skj.noise.patchy + rnorm(length(tmp.dt$skj.noise.patchy),0,cv*tmp.dt$skj.noise.patchy)
		tmp.dt$response = ifelse(tmp.dt$response<0,0,tmp.dt$response)

	# trim columns
		tmp.dt = tmp.dt[,.(ts,yy,qq,yyqq,lon,lat,response,id.data,id.5x5,id.raster)]

	return(tmp.dt)
}

fixed.seasonal.closure.samp = function(data.dt, ts.sampled = c(1,40), n.samp = 20000, cv = 0.15, random.seed = 123)
# no samples in the 3rd quarter of the year within 20S to 20N
{
	# set seed
		set.seed(random.seed)

	# subset to desired ts
		valid.ts = ts.sampled[1]:ts.sampled[2]	
		tmp.dt = data.dt[ts %in% valid.ts]

	# define number of samples taken per ts (total across ts must sum to n.samp)
		samp.per.ts = table(sample(valid.ts,n.samp,replace=TRUE))

	# define list to store samples from each ts
		ts.samp.list = as.list(rep(NA,length(valid.ts)))

	# iterate across ts and sample tmp.dt preferentially with respect to non-patchy abundance field (want samples to include zeros)
	# store the results in ts.samp.list
	# preferential sampling weight is sqrt(abundance) aka proportional to the square-root of abundance
		for(i in 1:length(valid.ts))
		{
			# define sampling probability
				ts.dt = tmp.dt[ts == valid.ts[i]]
				prob.vector = sqrt(ts.dt$skj.noise)

			# mask cells due to the closure
				closure.mask = which(ts.dt$qq == 3 & ts.dt$lat >=-20 & ts.dt$lat<= 20)
				if(length(closure.mask)>0)
				{
					prob.vector[closure.mask] = 0
				}

			# sample rows
				ts.dt = ts.dt[sample(1:nrow(ts.dt),size=samp.per.ts[i],replace=FALSE,prob=prob.vector),][order(ts,lon,lat)]

			# add to ts.samp.list
				ts.samp.list[[i]] = ts.dt
		}

	# combine ts samples
		tmp.dt = do.call(rbind,ts.samp.list)

	# add observation error
		tmp.dt$response = tmp.dt$skj.noise.patchy + rnorm(length(tmp.dt$skj.noise.patchy),0,cv*tmp.dt$skj.noise.patchy)
		tmp.dt$response = ifelse(tmp.dt$response<0,0,tmp.dt$response)

	# trim columns
		tmp.dt = tmp.dt[,.(ts,yy,qq,yyqq,lon,lat,response,id.data,id.5x5,id.raster)]

	return(tmp.dt)

}


rotate.seasonal.closure.samp = function(data.dt, ts.sampled = c(1,40), n.samp = 20000, cv = 0.15, random.seed = 123)
# sampling closure rotates
# Q1 the NE quadrant is closed (>155 & >15)
# Q2 the SE quadrant is closed (>155 & <15)
# Q3 the SW quadrant is closed (<155 & <15)
# Q4 the NW quadrant is closed (<155 & >15)
{
	# set seed
		set.seed(random.seed)

	# subset to desired ts
		valid.ts = ts.sampled[1]:ts.sampled[2]	
		tmp.dt = data.dt[ts %in% valid.ts]

	# define number of samples taken per ts (total across ts must sum to n.samp)
		samp.per.ts = table(sample(valid.ts,n.samp,replace=TRUE))

	# define list to store samples from each ts
		ts.samp.list = as.list(rep(NA,length(valid.ts)))

	# iterate across ts and sample tmp.dt preferentially with respect to non-patchy abundance field (want samples to include zeros)
	# store the results in ts.samp.list
	# preferential sampling weight is sqrt(abundance) aka proportional to the square-root of abundance
		for(i in 1:length(valid.ts))
		{
			# define sampling probability
				ts.dt = tmp.dt[ts == valid.ts[i]]
				prob.vector = sqrt(ts.dt$skj.noise)

			# mask cells due to the closure
				if(unique(ts.dt$qq)==1)
				{
					closure.mask = which(ts.dt$lon >155 & ts.dt$lat> 15)
					prob.vector[closure.mask] = 0
				}
				if(unique(ts.dt$qq)==2)
				{
					closure.mask = which(ts.dt$lon >155 & ts.dt$lat<= 15)
					prob.vector[closure.mask] = 0
				}
				if(unique(ts.dt$qq)==3)
				{
					closure.mask = which(ts.dt$lon <=155 & ts.dt$lat<= 15)
					prob.vector[closure.mask] = 0
				}
				if(unique(ts.dt$qq)==4)
				{
					closure.mask = which(ts.dt$lon <=155 & ts.dt$lat> 15)
					prob.vector[closure.mask] = 0
				}

			# sample rows
				ts.dt = ts.dt[sample(1:nrow(ts.dt),size=samp.per.ts[i],replace=FALSE,prob=prob.vector),][order(ts,lon,lat)]

			# add to ts.samp.list
				ts.samp.list[[i]] = ts.dt
		}

	# combine ts samples
		tmp.dt = do.call(rbind,ts.samp.list)

	# add observation error
		tmp.dt$response = tmp.dt$skj.noise.patchy + rnorm(length(tmp.dt$skj.noise.patchy),0,cv*tmp.dt$skj.noise.patchy)
		tmp.dt$response = ifelse(tmp.dt$response<0,0,tmp.dt$response)

	# trim columns
		tmp.dt = tmp.dt[,.(ts,yy,qq,yyqq,lon,lat,response,id.data,id.5x5,id.raster)]

	return(tmp.dt)

}


pref.samp = function(data.dt, ts.sampled = c(1,40), n.samp = 20000, cv = 0.15, random.seed = 123)
# sample preferentially with respect to percieved abundance aka fisheries dependent
{
	# set seed
		set.seed(random.seed)

	# subset to desired ts
		valid.ts = ts.sampled[1]:ts.sampled[2]	
		tmp.dt = data.dt[ts %in% valid.ts]

	# define number of samples taken per ts (total across ts must sum to n.samp)
		samp.per.ts = table(sample(valid.ts,n.samp,replace=TRUE))

	# define list to store samples from each ts
		ts.samp.list = as.list(rep(NA,length(valid.ts)))

	# iterate across ts and sample tmp.dt preferentially with respect to non-patchy abundance field (want samples to include zeros)
	# store the results in ts.samp.list
	# preferential sampling weight is sqrt(abundance) aka proportional to the square-root of abundance
		for(i in 1:length(valid.ts))
		{
			# define sampling probability
				ts.dt = tmp.dt[ts == valid.ts[i]]
				prob.vector = sqrt(ts.dt$skj.noise)

			# sample rows
				ts.dt = ts.dt[sample(1:nrow(ts.dt),size=samp.per.ts[i],replace=FALSE,prob=prob.vector),][order(ts,lon,lat)]

			# add to ts.samp.list
				ts.samp.list[[i]] = ts.dt
		}

	# combine ts samples
		tmp.dt = do.call(rbind,ts.samp.list)

	# add observation error
		tmp.dt$response = tmp.dt$skj.noise.patchy + rnorm(length(tmp.dt$skj.noise.patchy),0,cv*tmp.dt$skj.noise.patchy)
		tmp.dt$response = ifelse(tmp.dt$response<0,0,tmp.dt$response)

	# trim columns
		tmp.dt = tmp.dt[,.(ts,yy,qq,yyqq,lon,lat,response,id.data,id.5x5,id.raster)]

	return(tmp.dt)

}

expand.samp = function(data.dt, ts.sampled = c(1,40), n.samp = 20000, cv = 0.15, random.seed = 123)
# sample preferentially with respect to percieved abundance aka fisheries dependent
# and allow for expansion of the fleet over the middle 75% of ts.sampled
{
	# set seed
		set.seed(random.seed)

	# subset to desired ts
		valid.ts = ts.sampled[1]:ts.sampled[2]	
		tmp.dt = data.dt[ts %in% valid.ts]

	# define number of samples taken per ts (total across ts must sum to n.samp)
		samp.per.ts = table(sample(valid.ts,n.samp,replace=TRUE))

	# define maximum distance per ts
	# allow for constrained brownian excursion between min distance and max distance during the middle 75% ts
	# brownian excursion code taken from https://stats.stackexchange.com/questions/163047/simulating-a-brownian-excursion-using-a-brownian-bridge
	# specifically whuber's answer
		min.dist = 1000
		max.dist = 10000
		r.vec = rep(NA,length(valid.ts))
		ts.quantiles =  floor(quantile(valid.ts,probs = seq(0, 1, 0.125)))
		r.vec[which(valid.ts == ts.quantiles[1]):which(valid.ts == ts.quantiles[2])] = min.dist
		r.vec[which(valid.ts == ts.quantiles[8]):which(valid.ts == ts.quantiles[9])] = max.dist
		bridge.index = which(is.na(r.vec))

		bridge.n = length(bridge.index)
		bridge.times = seq(0, 1, length.out=bridge.n)
		bridge.target = max.dist - min.dist # Constraint at time=1
		bridge.dW = rnorm(bridge.n,0,0.5*min.dist)
		bridge.W = cumsum(bridge.dW)
		bridge.B = bridge.W + bridge.times * (bridge.target - bridge.W[bridge.n])   # The Brownian bridge from (0,0) to (1,target)
		bridge.range = 2*(bridge.target - 0)
		bridge.B = (bridge.B - 0) %% bridge.range
		bridge.B = pmin(bridge.B, bridge.range-bridge.B) + 0
		r.vec[bridge.index] = min.dist + bridge.B

	# define list to store samples from each ts
		ts.samp.list = as.list(rep(NA,length(valid.ts)))

	# iterate across ts and sample tmp.dt preferentially with respect to non-patchy abundance field (want samples to include zeros)
	# store the results in ts.samp.list
	# preferential sampling weight is sqrt(abundance) aka proportional to the square-root of abundance
		for(i in 1:length(valid.ts))
		{
			# define sampling probability
				ts.dt = tmp.dt[ts == valid.ts[i]]
				prob.vector = sqrt(ts.dt$skj.noise)
				prob.vector[which(ts.dt$dist.jp.km>r.vec[i])] = 0

			# sample rows
				ts.dt = ts.dt[sample(1:nrow(ts.dt),size=samp.per.ts[i],replace=TRUE,prob=prob.vector),][order(ts,lon,lat)]

			# add to ts.samp.list
				ts.samp.list[[i]] = ts.dt
		}

	# combine ts samples
		tmp.dt = do.call(rbind,ts.samp.list)

	# add observation error
		tmp.dt$response = tmp.dt$skj.noise.patchy + rnorm(length(tmp.dt$skj.noise.patchy),0,cv*tmp.dt$skj.noise.patchy)
		tmp.dt$response = ifelse(tmp.dt$response<0,0,tmp.dt$response)

	# trim columns
		tmp.dt = tmp.dt[,.(ts,yy,qq,yyqq,lon,lat,response,id.data,id.5x5,id.raster)]

	return(tmp.dt)

}

contract.samp = function(data.dt, ts.sampled = c(1,40), n.samp = 20000, cv = 0.15, random.seed = 123)
# sample preferentially with respect to percieved abundance aka fisheries dependent
# and allow for contraction of the fleet over the middle 75% of ts.sampled
{
	# set seed
		set.seed(random.seed)

	# subset to desired ts
		valid.ts = ts.sampled[1]:ts.sampled[2]	
		tmp.dt = data.dt[ts %in% valid.ts]

	# define number of samples taken per ts (total across ts must sum to n.samp)
		samp.per.ts = table(sample(valid.ts,n.samp,replace=TRUE))

	# define maximum distance per ts
	# allow for constrained brownian excursion between min distance and max distance during the middle 75% ts
	# brownian excursion code taken from https://stats.stackexchange.com/questions/163047/simulating-a-brownian-excursion-using-a-brownian-bridge
	# specifically whuber's answer
		min.dist = 1000
		max.dist = 10000
		r.vec = rep(NA,length(valid.ts))
		ts.quantiles =  floor(quantile(valid.ts,probs = seq(0, 1, 0.125)))
		r.vec[which(valid.ts == ts.quantiles[1]):which(valid.ts == ts.quantiles[2])] = max.dist
		r.vec[which(valid.ts == ts.quantiles[8]):which(valid.ts == ts.quantiles[9])] = min.dist
		bridge.index = which(is.na(r.vec))

		bridge.n = length(bridge.index)
		bridge.times = seq(0, 1, length.out=bridge.n)
		bridge.target = max.dist - min.dist # Constraint at time=1
		bridge.dW = rnorm(bridge.n,0,0.5*min.dist)
		bridge.W = cumsum(bridge.dW)
		bridge.B = bridge.W + bridge.times * (bridge.target - bridge.W[bridge.n])   # The Brownian bridge from (0,0) to (1,target)
		bridge.range = 2*(bridge.target - 0)
		bridge.B = (bridge.B - 0) %% bridge.range
		bridge.B = pmin(bridge.B, bridge.range-bridge.B) + 0
		r.vec[bridge.index] = rev(min.dist + bridge.B)

	# define list to store samples from each ts
		ts.samp.list = as.list(rep(NA,length(valid.ts)))

	# iterate across ts and sample tmp.dt preferentially with respect to non-patchy abundance field (want samples to include zeros)
	# store the results in ts.samp.list
	# preferential sampling weight is sqrt(abundance) aka proportional to the square-root of abundance
		for(i in 1:length(valid.ts))
		{
			# define sampling probability
				ts.dt = tmp.dt[ts == valid.ts[i]]
				prob.vector = sqrt(ts.dt$skj.noise)
				prob.vector[which(ts.dt$dist.jp.km>r.vec[i])] = 0

			# sample rows
				ts.dt = ts.dt[sample(1:nrow(ts.dt),size=samp.per.ts[i],replace=TRUE,prob=prob.vector),][order(ts,lon,lat)]

			# add to ts.samp.list
				ts.samp.list[[i]] = ts.dt
		}

	# combine ts samples
		tmp.dt = do.call(rbind,ts.samp.list)

	# add observation error
		tmp.dt$response = tmp.dt$skj.noise.patchy + rnorm(length(tmp.dt$skj.noise.patchy),0,cv*tmp.dt$skj.noise.patchy)
		tmp.dt$response = ifelse(tmp.dt$response<0,0,tmp.dt$response)

	# trim columns
		tmp.dt = tmp.dt[,.(ts,yy,qq,yyqq,lon,lat,response,id.data,id.5x5,id.raster)]

	return(tmp.dt)

}


##########################################################
# Sample from functions and save in appropriate folder

for(i in 1:100)
{
		save.id = i
		if(i<100){save.id = paste0("0",save.id)}
		if(i<10){save.id = paste0("0",save.id)}

	# RandomZero
		samp.dt = rand.samp(data.dt, ts.sampled = c(1,120), n.samp = 60000, cv = 0, random.seed = i)
		save(samp.dt,file=paste0("SimData/Simple120/RandomZero/samp.dt.",save.id,".RData"))
		rm(list=c("samp.dt"))
	# Random
		samp.dt = rand.samp(data.dt, ts.sampled = c(1,120), n.samp = 60000, cv = 0.15, random.seed = i)
		save(samp.dt,file=paste0("SimData/Simple120/Random/samp.dt.",save.id,".RData"))
		rm(list=c("samp.dt"))
	# Fixed
		samp.dt = fixed.seasonal.closure.samp(data.dt, ts.sampled = c(1,120), n.samp = 60000, cv = 0.15, random.seed = i)
		save(samp.dt,file=paste0("SimData/Simple120/Fixed/samp.dt.",save.id,".RData"))
		rm(list=c("samp.dt"))
	# Rotating	
		samp.dt = rotate.seasonal.closure.samp(data.dt, ts.sampled = c(1,120), n.samp = 60000, cv = 0.15, random.seed = i)
		save(samp.dt,file=paste0("SimData/Simple120/Rotating/samp.dt.",save.id,".RData"))
		rm(list=c("samp.dt"))
	# Preferential	
		samp.dt = pref.samp(data.dt, ts.sampled = c(1,120), n.samp = 60000, cv = 0.15, random.seed = i)
		save(samp.dt,file=paste0("SimData/Simple120/Preferential/samp.dt.",save.id,".RData"))
		rm(list=c("samp.dt"))
	# Expansion	
		samp.dt = expand.samp(data.dt, ts.sampled = c(1,120), n.samp = 60000, cv = 0.15, random.seed = i)
		save(samp.dt,file=paste0("SimData/Simple120/Expansion/samp.dt.",save.id,".RData"))
		rm(list=c("samp.dt"))
	# Contraction	
		samp.dt = contract.samp(data.dt, ts.sampled = c(1,120), n.samp = 60000, cv = 0.15, random.seed = i)
		save(samp.dt,file=paste0("SimData/Simple120/Contraction/samp.dt.",save.id,".RData"))
		rm(list=c("samp.dt"))
}


