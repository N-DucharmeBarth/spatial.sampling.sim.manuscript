

# Nicholas Ducharme-Barth
# 24/12/2019
# Create baseline abundance data files
# 1) based on seapodym output
# 2) introduce variability
# 3) introduce zeros

# load packages
	library(ndd.vast.utils)
	data(skj.alt2019.shp)

# bring in seapodym data
		setwd("C:/Users/nicholasd/HOME/SPC/SPC_SAM/Geostats/spatial.sampling.sim.manuscript/")	
		seapodym.raw = data.table::fread("SEAPODYM/SKJ_adultbio_1979_to_2010_allPac_INTERIM.csv")
	
	# trim columns and rename
		seapodym.raw = seapodym.raw[,c(9,8,10,12,11)]
		colnames(seapodym.raw) = c("yy","qq","lon","lat","skj")
		seapodym.raw$qq = as.numeric(substr(seapodym.raw$qq,2,2))
		seapodym.raw$yyqq = seapodym.raw$yy + c(0,0.25,0.5,0.75)[seapodym.raw$qq]
		seapodym.raw$ts = as.numeric(as.factor(seapodym.raw$yyqq))
		data.dt = seapodym.raw[lon <=210 & lat >= -20 & lat <= 50,.(ts,yy,qq,yyqq,lon,lat,skj)]

# add variability
	# set seed
		set.seed(123)

	# add noise
		cv = 0
		data.dt$skj.noise = data.dt$skj + rnorm(length(data.dt$skj),0,cv*data.dt$skj)
		data.dt$skj.noise = ifelse(data.dt$skj.noise<0,0,data.dt$skj.noise)

# add zeros
	prop.zero = 0.1
	zero.cells = as.vector(rmultinom(1,round(prop.zero*nrow(data.dt)),1/sqrt(data.dt$skj.noise)))

	zero.cells = ifelse(zero.cells>0,0,1)
	data.dt$skj.noise.patchy = zero.cells * data.dt$skj.noise

# add column for data.id
	data.dt$id.data = 1:nrow(data.dt)

# add column for 5x5 cell id
	data.dt$id.5x5 = as.numeric(as.factor(paste0(floor(data.dt$lon/5)*5,"_",floor(data.dt$lat/5)*5)))

# add column for raster cell id
	dummy.ras = raster::raster(xmn=99.5, xmx=210.5, ymn=-20.5, ymx=50.5, crs="+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0", resolution=c(1,1))
	raster::values(dummy.ras) = 1:raster::ncell(dummy.ras)  # rep(NA,raster::ncell(dummy.ras))
	data.dt$id.raster = raster::extract(dummy.ras,as.matrix(data.dt[,.(lon,lat)]))

# add column for distance from Japan (35.689722, 139.692222), Tokyo used as the distance point
	data.dt$dist.jp.km = geosphere::distHaversine(c(139.692222,35.689722),as.matrix(data.dt[,.(lon,lat)]))/1000

		# convert to rasters & plot
			# raster::values(dummy.ras) = rep(NA,raster::ncell(dummy.ras))

			# data.rasters = as.list(rep(NA,max(data.dt$ts)))
			# for(i in 1:1)
			# # for(i in 1:max(data.dt$ts))
			# {
			# 	tmp.dt = data.dt[ts == i]
			# 	data.rasters[[i]] = raster::rasterize(x=tmp.dt[,.(lon,lat)],y=dummy.ras,field=tmp.dt$skj.noise.patchy,fun="sum",na.rm=TRUE)
			# 	rm(list=c("tmp.dt"))
			# }
			# raster::plot(data.rasters[[1]],col=c("gray50",viridis::viridis(500)))

# test estimability
	# Data_Geostat = data.dt[ts %in% 1:40]
	# 		Data_Geostat = Data_Geostat[,.(skj.noise.patchy,ts,lon,lat)]
	# 		colnames(Data_Geostat) = c("Response_variable","Year","Lon","Lat")
	# 		Data_Geostat$Spp = as.factor("skj")
	# 		Data_Geostat$AreaSwept_km2 = 110^2
	# 		Data_Geostat$Vessel = "missing"
	# # subset to 10k obs for testing
	# 		Data_Geostat = as.data.frame(Data_Geostat[sample(1:nrow(Data_Geostat),20000),])
	# vast_output = fit.vast(Data_Geostat,RunDir=paste0(getwd(),"/preliminary_ignore/model_runs"),SaveDir=paste0(getwd(),"/preliminary_ignore/model_runs/data.dt"),save.output=FALSE,Q_ik = NULL,vf.re = FALSE,FieldConfig=c(Omega1 = 1, Epsilon1 = 1, Omega2 = 1, Epsilon2 = 1),RhoConfig=c(Beta1 = 0, Beta2 = 0, Epsilon1 = 0, Epsilon2 = 0),ObsModel_ez = c(1,3),fine_scale=FALSE,input.grid.res=1,knot_method = "grid",n_x=150,Version="VAST_v8_3_0",Method="Mesh",strata.sp=skj.alt2019.shp)
	# vast_output$Opt$max_gradient
	# vast_output$fit.time/60

	# # plot indices
	# 	strata.sp = skj.alt2019.shp
	# 		# define union of model region
	# 			n.substrata = length(strata.sp)
	# 			full.reg = strata.sp[1]
	# 			for(i in 2:length(strata.sp)){full.reg = rgeos::gUnion(full.reg,strata.sp[i])}
	# 			strata.sp = rbind(full.reg,strata.sp)

	#   			par(mfrow=c(3,3))
	#   			rmse.vec = rep(NA,9)
	# 			for(i in 1:9)
	# 			{
	# 				extrap.points = as.data.frame(data.dt)
	# 				sp::coordinates(extrap.points) = c("lon", "lat")
	# 	  			sp::proj4string(extrap.points) = sp::proj4string(strata.sp)
	# 	  			extrap.points$valid = sp::over(extrap.points,strata.sp[i])
	# 	  			true.dt = data.table::as.data.table(extrap.points@data)
	# 				true.index = true.dt[ts %in% 1:40 & !is.na(valid),.(Index=mean(skj.noise.patchy)),by=ts]

	# 				plot(scale(true.index$Index),pch=16,ylim=c(-2,2))
	# 				lines(scale(vast_output$Report$Index_cyl[1,,i+1]),col="hotpink",lwd=2,lty=2)
	# 				rmse.vec[i] = sqrt(mean((scale(true.index$Index)-scale(vast_output$Report$Index_cyl[1,,i+1]))^2))
	# 				rm(list=c("extrap.points","true.dt","true.index"))

	# 			}
	# 	rmse.vec
	# 	mean(rmse.vec)

# save
	save(dummy.ras,file="Background_Data/dummy.ras.RData")
	save(data.dt,file="Background_Data/data.dt.RData")

