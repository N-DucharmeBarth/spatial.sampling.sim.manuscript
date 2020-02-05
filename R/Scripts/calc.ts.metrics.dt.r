

# Nicholas Ducharme-Barth
# 06/02/2020
# Iterate over indices and calculates RMSE & bias coefficient at each time step
# match this up with the proportion of the area sampled within each region

# load packages
	library(sp)
	library(data.table)
	library(ndd.vast.utils)
	data(skj.alt2019.shp)
	data(pacific.coast)

# set project directory
	project.dir = "C:/Users/nicholasd/HOME/SPC/SPC_SAM/Geostats/spatial.sampling.sim.manuscript/"
	setwd(project.dir)

# load true
	load("SimData/simple.true.index.RData")

# calculate the number of 1x1 cells in each spatial region & not on land
	coords = SpatialPoints(expand.grid(lon=105:215,lat=-25:55),proj4string=CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
	cells.per.region = summary(as.factor(as.character(over(coords,skj.alt2019.shp)*ifelse(is.na(over(coords,pacific.coast)),1,NA))))[1:8]
	cells.per.region = c(sum(cells.per.region),cells.per.region)
	names(cells.per.region)[1] = "all"

# bring in ts.df
	load("Index/ResultsDF/ts.df.RData")
	ts.dt = as.data.table(ts.df)
	setkey(ts.dt,Scenario,Catchability,Replicate,Model,Region)

# identify unique time series
	unique.dt = ts.dt
	unique.dt = unique(unique.dt[!is.na(mgc),.(Scenario,Catchability,Replicate,Model,Region,mgc)])

# create storage structure
	ts.metrics.dt = ts.dt[!is.na(mgc)]
	ts.metrics.dt$sample.prop = as.numeric(NA)
	ts.metrics.dt$true = as.numeric(NA)
	ts.metrics.dt$rmse = as.numeric(NA)
	ts.metrics.dt$pe = as.numeric(NA)
	ts.metrics.dt$ape = as.numeric(NA)

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

# iterate across unique time series
	first.index = seq(from=1,by=120,length.out=nrow(unique.dt))
	last.index = seq(from=120,by=120,length.out=nrow(unique.dt))

	A = proc.time()
	for(i in 1:nrow(unique.dt))
	{
		# define vars
			s = as.character(unique.dt$Scenario[i])
			q = as.character(unique.dt$Catchability[i])
			rep = unique.dt$Replicate[i]
			m = as.character(unique.dt$Model[i])
			area = as.character(unique.dt$Region[i])
			save.id = rep
			if(rep<100){save.id = paste0("0",save.id)}
			if(rep<10){save.id = paste0("0",save.id)}

		# load samp.dt
			load(paste0("SimData/Simple120/",s,"/samp.dt.",save.id,".RData"))

		# calc sample.prop
			coords =  SpatialPoints(as.matrix(samp.dt[,.(lon,lat)]),proj4string=CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
			samp.dt$region = as.character(over(coords,skj.alt2019.shp))
			samp.dt$id.1x1 = as.numeric(as.factor(paste0(samp.dt$lon,"_",samp.dt$lat)))
			if(area == "all")
			{
				samp.dt = samp.dt[!is.na(region)]
			} else {
				samp.dt = samp.dt[!is.na(region) & region == area]
			}
			sample.prop = rep(0,120)
			tmp.dt = samp.dt[,.(uCell = length(unique(id.1x1))),by=ts]
			sample.prop[tmp.dt$ts] = tmp.dt$uCell/cells.per.region[area]

		# get true
			true = simple.true.index[,area]

		# get est
			est = ts.dt[.(s,q,rep,m,area)]$Index

		# calc metrics
			rmse = sqrt(( est - true )^2)
			pe = 100*(est-true)/true
			ape = 100*abs((est-true))/true
		
		# put in ts.metrics.dt
			ts.metrics.dt$sample.prop[first.index[i]:last.index[i]] = sample.prop
			ts.metrics.dt$true[first.index[i]:last.index[i]] = true
			ts.metrics.dt$rmse[first.index[i]:last.index[i]] = rmse
			ts.metrics.dt$pe[first.index[i]:last.index[i]] = pe
			ts.metrics.dt$ape[first.index[i]:last.index[i]] = ape

		# clean-up
			rm(list=c("s","q","rep","m","area","save.id","coords","samp.dt","sample.prop","tmp.dt","true","est","rmse","pe","ape"))
	}
	B = proc.time()
	B - A

# save
	save(ts.metrics.dt,file="Index/ResultsDF/ts.metrics.dt.RData")

