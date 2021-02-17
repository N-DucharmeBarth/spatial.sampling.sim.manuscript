

# Nicholas Ducharme-Barth
# 07/01/2021
# Simulate data for preferential sampling sensitivity with replacement


##########################################################
# set working directory
	setwd("C:/Users/nicholasd/HOME/SPC/SPC_SAM/Geostats/spatial.sampling.sim.manuscript/")	


##########################################################
# Define external variables for function development
	load("Background_Data/data.dt.RData")
	data.dt = data.table::as.data.table(data.dt)

##########################################################
# Define function

pref.samp.titration = function(data.dt, ts.sampled = c(1,40), n.samp = 20000, cv = 0.15, random.seed = 123, prob.weight.exponent=0.5)
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
	# preferential sampling weight is sqrt(abundance) aka proportional to a power (prob.weight.exponent) of abundance
		for(i in 1:length(valid.ts))
		{
			# define sampling probability
				ts.dt = tmp.dt[ts == valid.ts[i]]
				prob.vector = (ts.dt$skj.noise)^prob.weight.exponent

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


weight.vec = c(0,0.25,0.5,1,1.5,2,4,8) # seq(from=0,to=5,by=0.5)

##########################################################
# Sample from functions and save in appropriate folder

iter = 1
for(i in 1:10)
{
	for(j in 1:length(weight.vec))
	{
		save.id = iter
		if(iter<100){save.id = paste0("0",save.id)}
		if(iter<10){save.id = paste0("0",save.id)}

		# Contraction	
		samp.dt = pref.samp.titration(data.dt, ts.sampled = c(1,120), n.samp = 60000, cv = 0.15, random.seed = i, prob.weight.exponent=weight.vec[j])
		save(samp.dt,file=paste0("SimData/Titration_Preferential_wReplacement/samp.dt.",save.id,".RData"))
		rm(list=c("samp.dt"))

		iter = iter +1
	}
}

##########################################################
# Plot examples
	r = colorRampPalette(c("gray90","#80deea","#29b6f6","#1e88e5","#283593"))(32)

	# Data:
	 example.list = list('0'= pref.samp.titration(data.dt, ts.sampled = c(1,120), n.samp = 60000, cv = 0.15, random.seed = 1, prob.weight.exponent=0),
	 	                 '0.25'= pref.samp.titration(data.dt, ts.sampled = c(1,120), n.samp = 60000, cv = 0.15, random.seed = 1, prob.weight.exponent=0.25),
	 	                 '0.5'= pref.samp.titration(data.dt, ts.sampled = c(1,120), n.samp = 60000, cv = 0.15, random.seed = 1, prob.weight.exponent=0.5),
	 	                 '1'= pref.samp.titration(data.dt, ts.sampled = c(1,120), n.samp = 60000, cv = 0.15, random.seed = 1, prob.weight.exponent=1),
	 	                 '1.5'= pref.samp.titration(data.dt, ts.sampled = c(1,120), n.samp = 60000, cv = 0.15, random.seed = 1, prob.weight.exponent=1.5),
	 	                 '2'= pref.samp.titration(data.dt, ts.sampled = c(1,120), n.samp = 60000, cv = 0.15, random.seed = 1, prob.weight.exponent=2),
	 	                 '4'= pref.samp.titration(data.dt, ts.sampled = c(1,120), n.samp = 60000, cv = 0.15, random.seed = 1, prob.weight.exponent=4),
	 	                 '8'= pref.samp.titration(data.dt, ts.sampled = c(1,120), n.samp = 60000, cv = 0.15, random.seed = 1, prob.weight.exponent=8))
	samp.dt = data.table::rbindlist(example.list,idcol="Exponent")


	#Plot:
	p = ggplot2::ggplot(samp.dt, ggplot2::aes(lon,lat)) + ggplot2::ylim(-20,50) + ggplot2::facet_wrap(~Exponent) + 
		 ggplot2::xlab("Longitude") + ggplot2::ylab("Latitude") + ggplot2::ggtitle("Preferential sampling: Probability exponent") +
	     ggplot2::stat_bin_hex(bins=35) + ggplot2::scale_fill_gradientn("Number of samples",colors=r,trans="sqrt") +
	     ggthemes::theme_few(base_size = 20) + ggplot2::coord_fixed()
    ggplot2::ggsave(filename=paste0("titration.preferential.effort.distribution.wReplacement.png"), plot = p, device = "png", path = "Plots/",
  			scale = 1.25, width = 16, height = 9, units = c("in"),
  			dpi = 300, limitsize = TRUE)


