

# Nicholas Ducharme-Barth
# 28/10/2020
# Simulate data for effort distribution titration


##########################################################
# set working directory
	setwd("C:/Users/nicholasd/HOME/SPC/SPC_SAM/Geostats/spatial.sampling.sim.manuscript/")	


##########################################################
# Define external variables for function development
	load("Background_Data/data.dt.RData")
	data.dt = data.table::as.data.table(data.dt)

##########################################################
# Define function
contract.titration.samp = function(data.dt, ts.sampled = c(1,40), n.samp = 20000, cv = 0.15, random.seed = 123,var.tune=0.5,min.dist=1000)
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
		max.dist = 10000
		r.vec = rep(NA,length(valid.ts))
		ts.quantiles =  floor(quantile(valid.ts,probs = seq(0, 1, 0.125)))
		r.vec[which(valid.ts == ts.quantiles[1]):which(valid.ts == ts.quantiles[2])] = max.dist
		r.vec[which(valid.ts == ts.quantiles[8]):which(valid.ts == ts.quantiles[9])] = min.dist
		bridge.index = which(is.na(r.vec))

		bridge.n = length(bridge.index)
		bridge.times = seq(0, 1, length.out=bridge.n)
		bridge.target = max.dist - min.dist # Constraint at time=1
		bridge.dW = rnorm(bridge.n,0,var.tune*min.dist)
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



n.dists = 8
min.dist.vec = seq(from=1000,to=9999,length.out=n.dists)
var.tune.vec = 100/min.dist.vec


##########################################################
# Sample from functions and save in appropriate folder

iter = 1
for(i in 1:10)
{
	for(j in 1:length(min.dist.vec))
	{
		save.id = iter
		if(iter<100){save.id = paste0("0",save.id)}
		if(iter<10){save.id = paste0("0",save.id)}

		# Contraction	
		samp.dt = contract.titration.samp(data.dt, ts.sampled = c(1,120), n.samp = 60000, cv = 0.15, random.seed = i,var.tune=var.tune.vec[j],min.dist=min.dist.vec[j])
		save(samp.dt,file=paste0("SimData/Titration/samp.dt.",save.id,".RData"))
		rm(list=c("samp.dt"))

		iter = iter +1
	}
}
