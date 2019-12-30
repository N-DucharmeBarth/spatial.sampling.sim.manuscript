

# Nicholas Ducharme-Barth
# 30/12/2019
# Estimate the index using the VAST model from the simple, simulated examples w/ an environmental covariate
# 1) Random
# 2) RandomZero
# 3) Preferential
# 4) Fixed
# 5) Rotating
# 6) Expansion
# 7) Contraction

	setwd("C:/Users/nicholasd/HOME/SPC/SPC_SAM/Geostats/spatial.sampling.sim.manuscript/")	

# load packages
	library(ndd.vast.utils)
	library(splines)
	data(skj.alt2019.shp)

# prepare enviro
	load("Background_Data/sst.storage.df.RData")
	enviro.dt = data.table::as.data.table(sst.storage.df)
	enviro.dt$yy = as.numeric(substr(enviro.dt$time,1,4))
	enviro.dt$qq = as.numeric(substr(enviro.dt$time,5,6))
	enviro.dt = enviro.dt[x> 85 & x<225 & y>-35 & y < 65 & yy>=1979 & yy <= 1988]
	enviro.dt$ts = as.numeric(as.factor(as.character(enviro.dt$time)))
	enviro.df = data.frame(Year = enviro.dt$ts,Lon = enviro.dt$x, Lat = enviro.dt$y, sst = scale(enviro.dt$sst))
	enviro = list(formula = "~bs(sst,df=8)", covariate_data = enviro.df)

# Estimate indices
	# Random
		for(i in 1:100)
		{
			# define save id
				save.id = i
				if(i<100){save.id = paste0("0",save.id)}
				if(i<10){save.id = paste0("0",save.id)}
			# bring in data
				load(paste0("SimData/Simple/Random/samp.dt.",save.id,".RData"))
				samp.dt = data.table::as.data.table(samp.dt)
			# format Data_Geostat
				Data_Geostat = samp.dt[,.(response,ts,lon,lat)]
				colnames(Data_Geostat) = c("Response_variable","Year","Lon","Lat")
				Data_Geostat$Spp = as.factor("skj")
				Data_Geostat$AreaSwept_km2 = 110^2
				Data_Geostat$Vessel = "missing"
				Data_Geostat = as.data.frame(Data_Geostat)
			# estimate index
				SaveDir = "VAST/Simple/Random.Enviro/"
				if(!dir.exists(SaveDir)){dir.create(SaveDir, recursive = TRUE)}
				vast_output = try(fit.vast(Data_Geostat,RunDir=paste0(getwd(),"/VAST"),SaveDir=paste0(getwd(),"/Simple/Random.Enviro/"),save.output=FALSE,Q_ik = NULL,vf.re = FALSE,FieldConfig=c(Omega1 = 1, Epsilon1 = 1, Omega2 = 1, Epsilon2 = 1),RhoConfig=c(Beta1 = 0, Beta2 = 0, Epsilon1 = 0, Epsilon2 = 0),ObsModel_ez = c(1,3),fine_scale=FALSE,input.grid.res=1,knot_method = "grid",n_x=150,Version="VAST_v8_3_0",Method="Mesh",strata.sp=skj.alt2019.shp,enviro=enviro),silent=TRUE)
			# save
				save(vast_output,file=paste0(SaveDir,save.id,".vast_output.RData"))
			# clean-up
				rm(list=c("save.id","samp.dt","Data_Geostat","vast_output","SaveDir"))
		}

	# RandomZero
		for(i in 1:100)
		{
			# define save id
				save.id = i
				if(i<100){save.id = paste0("0",save.id)}
				if(i<10){save.id = paste0("0",save.id)}
			# bring in data
				load(paste0("SimData/Simple/RandomZero/samp.dt.",save.id,".RData"))
				samp.dt = data.table::as.data.table(samp.dt)
			# format Data_Geostat
				Data_Geostat = samp.dt[,.(response,ts,lon,lat)]
				colnames(Data_Geostat) = c("Response_variable","Year","Lon","Lat")
				Data_Geostat$Spp = as.factor("skj")
				Data_Geostat$AreaSwept_km2 = 110^2
				Data_Geostat$Vessel = "missing"
				Data_Geostat = as.data.frame(Data_Geostat)
			# estimate index
				SaveDir = "VAST/Simple/RandomZero.Enviro/"
				if(!dir.exists(SaveDir)){dir.create(SaveDir, recursive = TRUE)}
				vast_output = try(fit.vast(Data_Geostat,RunDir=paste0(getwd(),"/VAST"),SaveDir=paste0(getwd(),"/Simple/RandomZero.Enviro/"),save.output=FALSE,Q_ik = NULL,vf.re = FALSE,FieldConfig=c(Omega1 = 1, Epsilon1 = 1, Omega2 = 1, Epsilon2 = 1),RhoConfig=c(Beta1 = 0, Beta2 = 0, Epsilon1 = 0, Epsilon2 = 0),ObsModel_ez = c(1,3),fine_scale=FALSE,input.grid.res=1,knot_method = "grid",n_x=150,Version="VAST_v8_3_0",Method="Mesh",strata.sp=skj.alt2019.shp,enviro=enviro),silent=TRUE)
			# save
				save(vast_output,file=paste0(SaveDir,save.id,".vast_output.RData"))
			# clean-up
				rm(list=c("save.id","samp.dt","Data_Geostat","vast_output","SaveDir"))
		}		

	# Preferential
		for(i in 1:100)
		{
			# define save id
				save.id = i
				if(i<100){save.id = paste0("0",save.id)}
				if(i<10){save.id = paste0("0",save.id)}
			# bring in data
				load(paste0("SimData/Simple/Preferential/samp.dt.",save.id,".RData"))
				samp.dt = data.table::as.data.table(samp.dt)
			# format Data_Geostat
				Data_Geostat = samp.dt[,.(response,ts,lon,lat)]
				colnames(Data_Geostat) = c("Response_variable","Year","Lon","Lat")
				Data_Geostat$Spp = as.factor("skj")
				Data_Geostat$AreaSwept_km2 = 110^2
				Data_Geostat$Vessel = "missing"
				Data_Geostat = as.data.frame(Data_Geostat)
			# estimate index
				SaveDir = "VAST/Simple/Preferential.Enviro/"
				if(!dir.exists(SaveDir)){dir.create(SaveDir, recursive = TRUE)}
				vast_output = try(fit.vast(Data_Geostat,RunDir=paste0(getwd(),"/VAST"),SaveDir=paste0(getwd(),"/Simple/Preferential.Enviro/"),save.output=FALSE,Q_ik = NULL,vf.re = FALSE,FieldConfig=c(Omega1 = 1, Epsilon1 = 1, Omega2 = 1, Epsilon2 = 1),RhoConfig=c(Beta1 = 0, Beta2 = 0, Epsilon1 = 0, Epsilon2 = 0),ObsModel_ez = c(1,3),fine_scale=FALSE,input.grid.res=1,knot_method = "grid",n_x=150,Version="VAST_v8_3_0",Method="Mesh",strata.sp=skj.alt2019.shp,enviro=enviro),silent=TRUE)
			# save
				save(vast_output,file=paste0(SaveDir,save.id,".vast_output.RData"))
			# clean-up
				rm(list=c("save.id","samp.dt","Data_Geostat","vast_output","SaveDir"))
		}			

	# Fixed	
		for(i in 1:100)
		{
			# define save id
				save.id = i
				if(i<100){save.id = paste0("0",save.id)}
				if(i<10){save.id = paste0("0",save.id)}
			# bring in data
				load(paste0("SimData/Simple/Fixed/samp.dt.",save.id,".RData"))
				samp.dt = data.table::as.data.table(samp.dt)
			# format Data_Geostat
				Data_Geostat = samp.dt[,.(response,ts,lon,lat)]
				colnames(Data_Geostat) = c("Response_variable","Year","Lon","Lat")
				Data_Geostat$Spp = as.factor("skj")
				Data_Geostat$AreaSwept_km2 = 110^2
				Data_Geostat$Vessel = "missing"
				Data_Geostat = as.data.frame(Data_Geostat)
			# estimate index
				SaveDir = "VAST/Simple/Fixed.Enviro/"
				if(!dir.exists(SaveDir)){dir.create(SaveDir, recursive = TRUE)}
				vast_output = try(fit.vast(Data_Geostat,RunDir=paste0(getwd(),"/VAST"),SaveDir=paste0(getwd(),"/Simple/Fixed.Enviro/"),save.output=FALSE,Q_ik = NULL,vf.re = FALSE,FieldConfig=c(Omega1 = 1, Epsilon1 = 1, Omega2 = 1, Epsilon2 = 1),RhoConfig=c(Beta1 = 0, Beta2 = 0, Epsilon1 = 0, Epsilon2 = 0),ObsModel_ez = c(1,3),fine_scale=FALSE,input.grid.res=1,knot_method = "grid",n_x=150,Version="VAST_v8_3_0",Method="Mesh",strata.sp=skj.alt2019.shp,enviro=enviro),silent=TRUE)
			# save
				save(vast_output,file=paste0(SaveDir,save.id,".vast_output.RData"))
			# clean-up
				rm(list=c("save.id","samp.dt","Data_Geostat","vast_output","SaveDir"))
		}	

	# Rotating	
		for(i in 1:100)
		{
			# define save id
				save.id = i
				if(i<100){save.id = paste0("0",save.id)}
				if(i<10){save.id = paste0("0",save.id)}
			# bring in data
				load(paste0("SimData/Simple/Rotating/samp.dt.",save.id,".RData"))
				samp.dt = data.table::as.data.table(samp.dt)
			# format Data_Geostat
				Data_Geostat = samp.dt[,.(response,ts,lon,lat)]
				colnames(Data_Geostat) = c("Response_variable","Year","Lon","Lat")
				Data_Geostat$Spp = as.factor("skj")
				Data_Geostat$AreaSwept_km2 = 110^2
				Data_Geostat$Vessel = "missing"
				Data_Geostat = as.data.frame(Data_Geostat)
			# estimate index
				SaveDir = "VAST/Simple/Rotating.Enviro/"
				if(!dir.exists(SaveDir)){dir.create(SaveDir, recursive = TRUE)}
				vast_output = try(fit.vast(Data_Geostat,RunDir=paste0(getwd(),"/VAST"),SaveDir=paste0(getwd(),"/Simple/Rotating.Enviro/"),save.output=FALSE,Q_ik = NULL,vf.re = FALSE,FieldConfig=c(Omega1 = 1, Epsilon1 = 1, Omega2 = 1, Epsilon2 = 1),RhoConfig=c(Beta1 = 0, Beta2 = 0, Epsilon1 = 0, Epsilon2 = 0),ObsModel_ez = c(1,3),fine_scale=FALSE,input.grid.res=1,knot_method = "grid",n_x=150,Version="VAST_v8_3_0",Method="Mesh",strata.sp=skj.alt2019.shp,enviro=enviro),silent=TRUE)
			# save
				save(vast_output,file=paste0(SaveDir,save.id,".vast_output.RData"))
			# clean-up
				rm(list=c("save.id","samp.dt","Data_Geostat","vast_output","SaveDir"))
		}

	# Expansion	
		for(i in 1:100)
		{
			# define save id
				save.id = i
				if(i<100){save.id = paste0("0",save.id)}
				if(i<10){save.id = paste0("0",save.id)}
			# bring in data
				load(paste0("SimData/Simple/Expansion/samp.dt.",save.id,".RData"))
				samp.dt = data.table::as.data.table(samp.dt)
			# format Data_Geostat
				Data_Geostat = samp.dt[,.(response,ts,lon,lat)]
				colnames(Data_Geostat) = c("Response_variable","Year","Lon","Lat")
				Data_Geostat$Spp = as.factor("skj")
				Data_Geostat$AreaSwept_km2 = 110^2
				Data_Geostat$Vessel = "missing"
				Data_Geostat = as.data.frame(Data_Geostat)
			# estimate index
				SaveDir = "VAST/Simple/Expansion.Enviro/"
				if(!dir.exists(SaveDir)){dir.create(SaveDir, recursive = TRUE)}
				vast_output = try(fit.vast(Data_Geostat,RunDir=paste0(getwd(),"/VAST"),SaveDir=paste0(getwd(),"/Simple/Expansion.Enviro/"),save.output=FALSE,Q_ik = NULL,vf.re = FALSE,FieldConfig=c(Omega1 = 1, Epsilon1 = 1, Omega2 = 1, Epsilon2 = 1),RhoConfig=c(Beta1 = 0, Beta2 = 0, Epsilon1 = 0, Epsilon2 = 0),ObsModel_ez = c(1,3),fine_scale=FALSE,input.grid.res=1,knot_method = "grid",n_x=150,Version="VAST_v8_3_0",Method="Mesh",strata.sp=skj.alt2019.shp,enviro=enviro),silent=TRUE)
			# save
				save(vast_output,file=paste0(SaveDir,save.id,".vast_output.RData"))
			# clean-up
				rm(list=c("save.id","samp.dt","Data_Geostat","vast_output","SaveDir"))
		}

	# Contraction	
		for(i in 1:100)
		{
			# define save id
				save.id = i
				if(i<100){save.id = paste0("0",save.id)}
				if(i<10){save.id = paste0("0",save.id)}
			# bring in data
				load(paste0("SimData/Simple/Contraction/samp.dt.",save.id,".RData"))
				samp.dt = data.table::as.data.table(samp.dt)
			# format Data_Geostat
				Data_Geostat = samp.dt[,.(response,ts,lon,lat)]
				colnames(Data_Geostat) = c("Response_variable","Year","Lon","Lat")
				Data_Geostat$Spp = as.factor("skj")
				Data_Geostat$AreaSwept_km2 = 110^2
				Data_Geostat$Vessel = "missing"
				Data_Geostat = as.data.frame(Data_Geostat)
			# estimate index
				SaveDir = "VAST/Simple/Contraction.Enviro/"
				if(!dir.exists(SaveDir)){dir.create(SaveDir, recursive = TRUE)}
				vast_output = try(fit.vast(Data_Geostat,RunDir=paste0(getwd(),"/VAST"),SaveDir=paste0(getwd(),"/Simple/Contraction.Enviro/"),save.output=FALSE,Q_ik = NULL,vf.re = FALSE,FieldConfig=c(Omega1 = 1, Epsilon1 = 1, Omega2 = 1, Epsilon2 = 1),RhoConfig=c(Beta1 = 0, Beta2 = 0, Epsilon1 = 0, Epsilon2 = 0),ObsModel_ez = c(1,3),fine_scale=FALSE,input.grid.res=1,knot_method = "grid",n_x=150,Version="VAST_v8_3_0",Method="Mesh",strata.sp=skj.alt2019.shp,enviro=enviro),silent=TRUE)
			# save
				save(vast_output,file=paste0(SaveDir,save.id,".vast_output.RData"))
			# clean-up
				rm(list=c("save.id","samp.dt","Data_Geostat","vast_output","SaveDir"))
		}

