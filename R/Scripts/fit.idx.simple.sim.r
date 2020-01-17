

# Nicholas Ducharme-Barth
# 17/01/2020
# Estimate the index for the simple sim using the nominal, glm, and VAST formulations
# 1) Random
# 2) RandomZero
# 3) Preferential
# 4) Fixed
# 5) Rotating
# 6) Expansion
# 7) Contraction



fit.simple.sim = function(s)
# s %in% "Fixed", "Rotating", "Expansion", "Contraction", "Preferential", "Random", & "RandomZero"
{
	setwd("C:/Users/nicholasd/HOME/SPC/SPC_SAM/Geostats/spatial.sampling.sim.manuscript/")	

	# source
		source("R/Fn/fn.fit.dglm.simple.r")
		source("R/Fn/fn.fit.hybrid.dglm.simple.r")
		source("R/Fn/fn.nominal.simple.r")
		load("SimData/simple.true.index.RData")

	# load packages
		library(ndd.vast.utils)
		library(splines)
		data(skj.alt2019.shp)

		load("Background_Data/sst.storage.df.RData")
		enviro.dt = data.table::as.data.table(sst.storage.df)
		enviro.dt$yy = as.numeric(substr(enviro.dt$time,1,4))
		enviro.dt$qq = as.numeric(substr(enviro.dt$time,5,6))
		enviro.dt = enviro.dt[x> 85 & x<225 & y>-35 & y < 65 & yy>=1979 & yy <= 1988]
		enviro.dt$ts = as.numeric(as.factor(as.character(enviro.dt$time)))
		enviro.df = data.frame(Year = enviro.dt$ts,Lon = enviro.dt$x, Lat = enviro.dt$y, sst = scale(enviro.dt$sst))
		enviro = list(formula = "~bs(sst)", covariate_data = enviro.df)

	for(i in 1:100)
	{
		effort.scenario = s
		save.id = i
		if(i<100){save.id = paste0("0",save.id)}
		if(i<10){save.id = paste0("0",save.id)}
		# bring in data
			load(paste0("SimData/Simple/",effort.scenario,"/samp.dt.",save.id,".RData"))
		# format data
			data = as.data.frame(samp.dt[,c("response","ts","lon","lat")])
			colnames(data) = c("Response_variable","Year","Lon","Lat")
		# format Data_Geostat
			Data_Geostat = as.data.frame(samp.dt[,c("response","ts","lon","lat")])
			colnames(Data_Geostat) = c("Response_variable","Year","Lon","Lat")
			Data_Geostat$Spp = as.factor("skj")
			Data_Geostat$AreaSwept_km2 = 110^2
			Data_Geostat$Vessel = "missing"
			Data_Geostat = as.data.frame(Data_Geostat)
		# get indices & save
			sd.dglm = paste0("Index/Simple/",effort.scenario,"/dglm/")
			idx_dglm = try(fit.dglm.simple(data,  agg.cell = 5, formula.stem = " ~ Year + agg.cell", data.weighting = FALSE,n.yr.rng = length(range(data$Year)[1]:range(data$Year)[2]),scale=TRUE, skj.alt2019.shp),silent=TRUE)
			if(!dir.exists(sd.dglm)){dir.create(sd.dglm, recursive = TRUE)}
			save(idx_dglm,file=paste0(sd.dglm,save.id,".idx_dglm.RData"))

			sd.hybrid = paste0("Index/Simple/",effort.scenario,"/hybrid/")
			idx_hybrid = try(fit.hybrid.dglm.simple(data, n_x = 10, formula.stem = " ~ Year * knot", data.weighting = FALSE,n.yr.rng = length(range(data$Year)[1]:range(data$Year)[2]), seed = 123,scale=TRUE, skj.alt2019.shp),silent=TRUE)
			if(!dir.exists(sd.hybrid)){dir.create(sd.hybrid, recursive = TRUE)}
			save(idx_hybrid,file=paste0(sd.hybrid,save.id,".idx_hybrid.RData"))

			sd.nominal = paste0("Index/Simple/",effort.scenario,"/nominal/")
			idx_nominal = try(nominal.simple(data,n.yr.rng = length(range(data$Year)[1]:range(data$Year)[2]),scale=TRUE, skj.alt2019.shp),silent=TRUE)
			if(!dir.exists(sd.nominal)){dir.create(sd.nominal, recursive = TRUE)}
			save(idx_nominal,file=paste0(sd.nominal,save.id,".idx_nominal.RData"))

			sd.vast = paste0("Index/Simple/",effort.scenario,"/vast/")
			idx_vast = try(fit.vast(Data_Geostat,RunDir=paste0(getwd(),"/VAST/"),SaveDir=paste0(getwd(),"/VAST/"),save.output=FALSE,Q_ik = NULL,vf.re = FALSE,FieldConfig=c(Omega1 = 1, Epsilon1 = 1, Omega2 = 1, Epsilon2 = 1),RhoConfig=c(Beta1 = 0, Beta2 = 0, Epsilon1 = 0, Epsilon2 = 0),ObsModel_ez = c(1,3),fine_scale=FALSE,input.grid.res=1,crop.extrap.by.data=TRUE,knot_method = "grid",n_x=150,Version="VAST_v8_3_0",Method="Mesh",ADREPORT=TRUE,normalize_idx=TRUE,strata.sp=skj.alt2019.shp,enviro=enviro),silent=TRUE)
			if(!dir.exists(sd.vast)){dir.create(sd.vast, recursive = TRUE)}
			save(idx_vast,file=paste0(sd.vast,save.id,".idx_vast.RData"))

			sd.vastNoEnviro = paste0("Index/Simple/",effort.scenario,"/vastNoEnviro/")
			idx_vastNoEnviro = try(fit.vast(Data_Geostat,RunDir=paste0(getwd(),"/VAST/"),SaveDir=paste0(getwd(),"/VAST/"),save.output=FALSE,Q_ik = NULL,vf.re = FALSE,FieldConfig=c(Omega1 = 1, Epsilon1 = 1, Omega2 = 1, Epsilon2 = 1),RhoConfig=c(Beta1 = 0, Beta2 = 0, Epsilon1 = 0, Epsilon2 = 0),ObsModel_ez = c(1,3),fine_scale=FALSE,input.grid.res=1,crop.extrap.by.data=TRUE,knot_method = "grid",n_x=150,Version="VAST_v8_3_0",Method="Mesh",ADREPORT=TRUE,normalize_idx=TRUE,strata.sp=skj.alt2019.shp),silent=TRUE)
			if(!dir.exists(sd.vastNoEnviro)){dir.create(sd.vastNoEnviro, recursive = TRUE)}
			save(idx_vastNoEnviro,file=paste0(sd.vastNoEnviro,save.id,".idx_vastNoEnviro.RData"))

		# clean-up
			rm(list=c("effort.scenario","save.id","samp.dt","Data_Geostat","data","idx_dglm","idx_hybrid","idx_nominal","sd.dglm","sd.hybrid","sd.nominal","sd.vast","sd.vastNoEnviro"))
	}
}

