

# Nicholas Ducharme-Barth
# 24/12/2019
# Create baseline abundance data files
# 1) based on seapodym output
# 2) introduce variability
# 3) introduce zeros

# load packages
	library(ndd.vast.utils)

# bring in seapodym data
		setwd("C:/Users/nicholasd/HOME/SPC/SPC_SAM/Geostats/spatial.sampling.sim.manuscript/")	
		seapodym.raw = data.table::fread("SEAPODYM/SKJ_adultbio_1979_to_2010_allPac_INTERIM.csv")
	
	# trim columns and rename
		seapodym.raw = seapodym.raw[,c(9,8,10,12,11)]
		colnames(seapodym.raw) = c("yy","qq","lon","lat","skj")
		seapodym.raw$qq = as.numeric(substr(seapodym.raw$qq,2,2))
		seapodym.raw$yyqq = seapodym.raw$yy + c(0,0.25,0.5,0.75)[seapodym.raw$qq]
		seapodym.raw$ts = as.numeric(as.factor(seapodym.raw$yyqq))
		data.dt = seapodym.raw[,.(ts,yy,qq,yyqq,lon,lat,skj)]

# add variability
	# set seed
		set.seed(123)

	# add noise
		cv = 0.15
		data.dt$skj.noise = data.dt$skj + rnorm(length(data.dt$skj),0,cv*data.dt$skj)
		data.dt$skj.noise = ifelse(data.dt$skj.noise<0,0,data.dt$skj.noise)

# add zeros
	prop.zero = 0.15
	zero.cells = as.vector(rmultinom(1,round(prop.zero*nrow(data.dt)),1/log(data.dt$skj.noise)))

	zero.cells = ifelse(zero.cells>0,0,1)
	data.dt$skj.noise.patchy = zero.cells * data.dt$skj.noise

		# convert to rasters & plot
			dummy.ras = raster::raster(xmn=99.5, xmx=290.5, ymn=-56.5, ymx=55.5, crs="+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0", resolution=c(1,1))
			raster::values(dummy.ras) = rep(NA,raster::ncell(dummy.ras))

			data.rasters = as.list(rep(NA,max(data.dt$ts)))
			for(i in 1:1)
			# for(i in 1:max(data.dt$ts))
			{
				tmp.dt = data.dt[ts == i]
				data.rasters[[i]] = raster::rasterize(x=tmp.dt[,.(lon,lat)],y=dummy.ras,field=tmp.dt$skj.noise.patchy,fun="sum",na.rm=TRUE)
				rm(list=c("tmp.dt"))
			}
			raster::plot(data.rasters[[1]],col=viridis::viridis(500))

# test estimability
	Data_Geostat = data.dt[ts %in% 1:10 & lon <=210 & lat >= -20 & lat <= 50]
			Data_Geostat = Data_Geostat[,.(skj.noise.patchy,ts,lon,lat)]
			colnames(Data_Geostat) = c("Response_variable","Year","Lon","Lat")
			Data_Geostat$Spp = as.factor("skj")
			Data_Geostat$AreaSwept_km2 = 110^2
			Data_Geostat$Vessel = "missing"
	vast_output = fit.vast(Data_Geostat,RunDir=paste0(getwd(),"/preliminary_ignore/model_runs"),SaveDir=paste0(getwd(),"/preliminary_ignore/model_runs/data.dt"),save.output=FALSE,Q_ik = NULL,vf.re = FALSE,FieldConfig=c(Omega1 = 1, Epsilon1 = 1, Omega2 = 1, Epsilon2 = 1),RhoConfig=c(Beta1 = 0, Beta2 = 0, Epsilon1 = 0, Epsilon2 = 0),ObsModel_ez = c(1,3),fine_scale=FALSE,input.grid.res=1,knot_method = "grid",n_x=100,Version="VAST_v8_3_0",Method="Mesh")


# save
	save(data.dt,file="Background_Data/data.dt.RData")

