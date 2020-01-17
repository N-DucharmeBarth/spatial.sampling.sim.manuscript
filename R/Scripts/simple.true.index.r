

# Nicholas Ducharme-Barth
# 17/01/2020
# Calculate true index in each region for the simple simulations

setwd("C:/Users/nicholasd/HOME/SPC/SPC_SAM/Geostats/spatial.sampling.sim.manuscript/")	

# load packages & background data
	library(ndd.vast.utils)
	data(skj.alt2019.shp)
	load("Background_Data/data.dt.RData")

# define data
	data.dt = data.table::as.data.table(data.dt)

# create storage structure
	true.index = matrix(NA,nrow=40,ncol=9)
	colnames(true.index) = c("all",1:8)
	rownames(true.index) = 1:40

# loop over strata
	strata.sp = skj.alt2019.shp
	n.substrata = length(strata.sp)
	full.reg = strata.sp[1]
	for(i in 2:length(strata.sp)){full.reg = rgeos::gUnion(full.reg,strata.sp[i])}
	strata.sp = rbind(full.reg,strata.sp)

	for(i in 1:ncol(true.index))
	{
		extrap.points = as.data.frame(data.dt)
		sp::coordinates(extrap.points) = c("lon", "lat")
		sp::proj4string(extrap.points) = sp::proj4string(strata.sp)
		extrap.points$valid = sp::over(extrap.points,strata.sp[i])
		true.dt = data.table::as.data.table(extrap.points@data)
		true.index[,i] = true.dt[ts %in% 1:40 & !is.na(valid),.(Index=sum(skj.noise.patchy)),by=ts]$Index
	}

# normalize by the mean
	simple.true.index = apply(true.index,2,function(x)x/mean(x))

# save
	save(simple.true.index,file="SimData/simple.true.index.RData")