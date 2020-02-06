

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

# get the spatial knot associated with each 1x1 cell
			strata.limits = data.frame('STRATA'="All_areas")
			Region = "user"
			input.grid.res = 1
			strata.sp = skj.alt2019.shp
			Data_Geostat = expand.grid(Lon=100:210,Lat=-20:50)
			Data_Geostat$id.1x1 = paste0(Data_Geostat$Lat,".",Data_Geostat$Lon)
			grid_size_km = (110 * input.grid.res)^2 # the distance between grid cells for the 2D AR1 grid if Method == "Grid"
			grid_bounds = c(floor(min(Data_Geostat[,c("Lat")])/5)*5,ceiling(max(Data_Geostat[,c("Lat")])/5)*5,floor(min(Data_Geostat[,c("Lon")])/5)*5,ceiling(max(Data_Geostat[,c("Lon")])/5)*5)
			input.grid = as.matrix(expand.grid(Lat = seq(from=grid_bounds[1],to=grid_bounds[2],by=input.grid.res),
					Lon = seq(from=grid_bounds[3],to=grid_bounds[4],by=input.grid.res), Area_km2 = grid_size_km))
			crs.en = paste0("+proj=tpeqd +lat_1=",mean(grid_bounds[1:2])," +lon_1=",round(grid_bounds[3] + (1/3)*abs(diff(grid_bounds[3:4])))," +lat_2=",mean(grid_bounds[1:2])," +lon_2=",round(grid_bounds[3] + (2/3)*abs(diff(grid_bounds[3:4])))," +datum=WGS84 +ellps=WGS84 +units=km +no_defs")
			crs.ll = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
			Extrapolation_List = make_extrapolation_info.ndd(Region = Region,strata.limits = strata.limits,input.grid = input.grid, crs.en = crs.en,crs.ll = crs.ll,strata.sp)
			 extrap.df = as.data.frame(Extrapolation_List$Data_Extrap)
			sp::coordinates(extrap.df) = ~E_km + N_km
			sp::proj4string(extrap.df) = crs.en
			pacific.coast = sp::spTransform(pacific.coast,crs.en)
			over.index.land = which(is.na(sp::over(extrap.df,pacific.coast)))
			# find overlap with data hull
			smooth.hull = smooth.hull.sp(Data_Geostat[,c("Lon","Lat")],crs.ll=crs.ll,buffer.ll=2.5,d.scalar = 0.15)
			# plot(smooth.hull)
			smooth.hull.trans = sp::spTransform(smooth.hull,crs.en)
			# identify number of extrapolation cells (within data hull and not on land) corresponding to each knot
				over.index.region = which(!is.na(sp::over(extrap.df,smooth.hull.trans)))
				over.index = intersect(over.index.land,over.index.region)
			over.index = over.index.land
			Extrapolation_List$a_el = data.frame(All_areas = Extrapolation_List$a_el[over.index,])
			Extrapolation_List$Data_Extrap = Extrapolation_List$Data_Extrap[over.index,]
			Extrapolation_List$Area_km2_x = Extrapolation_List$Area_km2_x[over.index]
			seed = 123 
			Spatial_List = make_spatial_info.ndd(n_x=150, Lon_i=Data_Geostat[,'Lon'], Lat_i=Data_Geostat[,'Lat'], Extrapolation_List = Extrapolation_List, knot_method="grid", Method="Mesh",
												  grid_size_km=grid_size_km, grid_size_LL=input.grid.res, fine_scale=FALSE, Network_sz_LL=NULL,
												  iter.max=1000, randomseed=seed, nstart=100, DirPath="/tmp", Save_Results=FALSE,crs.en = crs.en,crs.ll = crs.ll)

			Data_Geostat = cbind(Data_Geostat, knot_i = Spatial_List$knot_i)


# get the number of knots in each region
	data(pacific.coast)
	coords = SpatialPoints(cbind(Data_Geostat$Lon,Data_Geostat$Lat),proj4string=CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
	Data_Geostat$area = over(coords,skj.alt2019.shp)*ifelse(is.na(over(coords,pacific.coast)),1,NA)
	Data_Geostat.dt = as.data.table(Data_Geostat)
	knots.per.region = Data_Geostat.dt[!is.na(area),.(N=length(unique(knot_i))),by=area][order(area)]$N
	knots.per.region = c(length(unique(Data_Geostat$knot_i)),knots.per.region)
	names(knots.per.region) = c("all",1:8)

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
			if(area == "all")
			{
				tmp.dg = Data_Geostat.dt[!is.na(area)]
			} else {
				tmp.dg = Data_Geostat.dt[area == as.character(unique.dt$Region[i])]
			}

			samp.dt$id.1x1 = paste0(samp.dt$lat,".",samp.dt$lon)
			samp.dt$knot_i = tmp.dg$knot_i[match(samp.dt$id.1x1,tmp.dg$id.1x1)]

			sample.prop = rep(0,120)
			tmp.dt = samp.dt[,.(uknots = length(unique(knot_i))),by=ts]
			sample.prop[tmp.dt$ts] = tmp.dt$uknots/knots.per.region[area]

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
			rm(list=c("s","q","rep","m","area","save.id","tmp.dg","samp.dt","sample.prop","tmp.dt","true","est","rmse","pe","ape"))
	}
	B = proc.time()
	B - A

# save
	save(ts.metrics.dt,file="Index/ResultsDF/ts.metrics.dt.RData")

