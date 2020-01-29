# clear workspace
	rm(list=ls())

# load packages
	library(ndd.vast.utils)
	library(rgdal)
	library(splines)
	data(skj.alt2019.shp)

# bring in environmental data
	load("sst.storage.df.RData")
	enviro.dt = data.table::as.data.table(sst.storage.df)
	enviro.dt$yy = as.numeric(substr(enviro.dt$time,1,4))
	enviro.dt$qq = as.numeric(substr(enviro.dt$time,5,6))
	enviro.dt = enviro.dt[x> 85 & x<225 & y>-35 & y < 65 & yy>=1979 & yy <= 2008]
	enviro.dt$ts = as.numeric(as.factor(as.character(enviro.dt$time)))
	enviro.df = data.frame(Year = enviro.dt$ts,Lon = enviro.dt$x, Lat = enviro.dt$y, sst = scale(enviro.dt$sst))

# add in the SVC
	load("nino.df.RData")
	nino.df = subset(nino.df,yr>=1979 & yr <= 2008)
	nino.df$ts = as.numeric(as.factor(as.character(nino.df$yrqtr)))
	nino.df = nino.df[,c("ts","index.scale")]
	colnames(nino.df) = c("Year","SVC")
	enviro.df$SVC = nino.df$SVC[match(enviro.df$Year,nino.df$Year)]

	enviro = list(formula = "~bs(sst)", covariate_data = enviro.df[,c("Year","Lon","Lat","sst")])
	enviro.svc = list(formula = "~bs(sst) + SVC", covariate_data = enviro.df[,c("Year","Lon","Lat","sst","SVC")])
	svc = list(formula = "~SVC", covariate_data = enviro.df[,c("Year","Lon","Lat","SVC")])

# bring in data
# add Catchability
	load("data.dt.RData")
	load(list.files()[grep("samp.dt",list.files())])
	source("fn.add.catchability.r")
	samp.dt$True_Abundance = data.dt$skj.noise.patchy[samp.dt$id.data]
	samp.dt = as.data.frame(samp.dt[,c("True_Abundance","ts","lon","lat")])
	colnames(samp.dt) = c("True_Abundance","Year","Lon","Lat")
	Data_Geostat = add.catchability(samp.dt,seed = as.numeric(strsplit(list.files()[grep("samp.dt",list.files())],"[.]")[[1]][3]),n.vessel.target = 90,new.entry.target=30,cv = 0.15,plot=FALSE)

# format Data_Geostat
	Data_Geostat$Spp = as.factor("skj")
	Data_Geostat$AreaSwept_km2 = 110^2
	Data_Geostat$Vessel = as.character(Data_Geostat$Vessel)
	Data_Geostat$Gear_Config = as.character(Data_Geostat$Gear_Config)
	Data_Geostat = as.data.frame(Data_Geostat)

# define Q_ik
	Q_ik = model.matrix(as.formula(paste0("Spp~Gear_Config+Class:bs(Poles)-1")),data=Data_Geostat)
	Q_ik = Q_ik[,-1]

# define the Xconfig_zcp object
	Xconfig_zcp.enviro = array( 0, dim = c( 2, 1, 4 ) )
	Xconfig_zcp.enviro[,,1:3] = 1
	Xconfig_zcp.enviro[,,4] = 2

	Xconfig_zcp.svc = array( 0, dim = c( 2, 1, 1 ) )
	Xconfig_zcp.svc[,,1] = 2

# fit indices
	A = proc.time()
	idx_vast.Enviro = try(fit.vast(Data_Geostat,RunDir=paste0(getwd()),SaveDir=paste0(getwd()),save.output=FALSE,Q_ik = Q_ik,vf.re = TRUE,FieldConfig=c(Omega1 = 1, Epsilon1 = 1, Omega2 = 1, Epsilon2 = 1),RhoConfig=c(Beta1 = 0, Beta2 = 0, Epsilon1 = 0, Epsilon2 = 0),ObsModel_ez = c(1,3),fine_scale=FALSE,input.grid.res=1,crop.extrap.by.data=TRUE,knot_method = "grid",n_x=150,Version="VAST_v8_3_0",Method="Mesh",ADREPORT=TRUE,normalize_idx=FALSE,Xconfig_zcp=NULL,slim.output = TRUE,strata.sp=skj.alt2019.shp,enviro=enviro),silent=TRUE)
	B = proc.time()
	idx_vast.NoEnviro = try(fit.vast(Data_Geostat,RunDir=paste0(getwd()),SaveDir=paste0(getwd()),save.output=FALSE,Q_ik = Q_ik,vf.re = TRUE,FieldConfig=c(Omega1 = 1, Epsilon1 = 1, Omega2 = 1, Epsilon2 = 1),RhoConfig=c(Beta1 = 0, Beta2 = 0, Epsilon1 = 0, Epsilon2 = 0),ObsModel_ez = c(1,3),fine_scale=FALSE,input.grid.res=1,crop.extrap.by.data=TRUE,knot_method = "grid",n_x=150,Version="VAST_v8_3_0",Method="Mesh",ADREPORT=TRUE,normalize_idx=FALSE,Xconfig_zcp=NULL,slim.output=TRUE,strata.sp=skj.alt2019.shp),silent=TRUE)
	C = proc.time()
	idx_vast.EnviroSVC = try(fit.vast(Data_Geostat,RunDir=paste0(getwd()),SaveDir=paste0(getwd()),save.output=FALSE,Q_ik = Q_ik,vf.re = TRUE,FieldConfig=c(Omega1 = 1, Epsilon1 = 1, Omega2 = 1, Epsilon2 = 1),RhoConfig=c(Beta1 = 0, Beta2 = 0, Epsilon1 = 0, Epsilon2 = 0),ObsModel_ez = c(1,3),fine_scale=FALSE,input.grid.res=1,crop.extrap.by.data=TRUE,knot_method = "grid",n_x=150,Version="VAST_v8_3_0",Method="Mesh",ADREPORT=TRUE,normalize_idx=FALSE,Xconfig_zcp=Xconfig_zcp.enviro,slim.output = TRUE,strata.sp=skj.alt2019.shp,enviro=enviro.svc),silent=TRUE)
	D = proc.time()
	idx_vast.NoEnviroSVC = try(fit.vast(Data_Geostat,RunDir=paste0(getwd()),SaveDir=paste0(getwd()),save.output=FALSE,Q_ik = Q_ik,vf.re = TRUE,FieldConfig=c(Omega1 = 1, Epsilon1 = 1, Omega2 = 1, Epsilon2 = 1),RhoConfig=c(Beta1 = 0, Beta2 = 0, Epsilon1 = 0, Epsilon2 = 0),ObsModel_ez = c(1,3),fine_scale=FALSE,input.grid.res=1,crop.extrap.by.data=TRUE,knot_method = "grid",n_x=150,Version="VAST_v8_3_0",Method="Mesh",ADREPORT=TRUE,normalize_idx=FALSE,Xconfig_zcp=Xconfig_zcp.svc,slim.output=TRUE,strata.sp=skj.alt2019.shp,enviro=svc),silent=TRUE)
	E = proc.time()

# load true.data and calc metrics
	load("simple.true.index.RData")
	simple.true.index = apply(simple.true.index,2,function(x)x/mean(x))

	estimate.idx.metrics = function(true,est,est.se)
	{
		# standardize first
			true = apply(true,2,function(x)x/mean(x))
			for(j in 1:ncol(est))
			{
				j.mean = mean(est[,j])
				est[,j] = est[,j]/j.mean
				est.se[,j] = est.se[,j]/j.mean
				rm(list=c("j.mean"))
			}

		# calc bias
			bias.vec = rep(NA,ncol(true))
			for(j in 1:ncol(true))
			{
				bias.vec[j] = summary( lm( est[,j] ~ true[,j] ) )$coef['true[, j]', 1]
			}

		# calc MAE
			mae.vec = rep(NA,ncol(true))
			for(j in 1:ncol(true))
			{
				mae.vec[j] = mean( abs ( est[,j] - true[,j] ) )
			}

		# calc coverage
			cover.vec = rep(NA,ncol(true))
			for(j in 1:ncol(true))
			{
				COVER  = rep( 0, len = nrow( true ) )
				UCI = est[,j] + 2*est.se[,j] 
				LCI = est[,j] - 2*est.se[,j] 
				for ( i in 1:length( COVER ) ) {
					if ( ( LCI[i]<=true[i,j] ) & ( true[i,j]<=UCI[i]) ){ COVER[i] = 1} 	
				}
				cover.vec[j] = ( sum( COVER )/ length( COVER ) ) * 100
				rm(list=c("COVER","UCI","LCI"))
			}

		return(data.frame(bias=bias.vec,mae=mae.vec,cover=cover.vec))		
	}

	NoEnviro.metrics = try(estimate.idx.metrics(true = simple.true.index,est = idx_vast.NoEnviro$idx[,2:10],est.se = idx_vast.NoEnviro$idx.se[,2:10]),silent=TRUE)
	Enviro.metrics = try(estimate.idx.metrics(true = simple.true.index,est = idx_vast.Enviro$idx[,2:10],est.se = idx_vast.Enviro$idx.se[,2:10]),silent=TRUE)
	NoEnviroSVC.metrics = try(estimate.idx.metrics(true = simple.true.index,est = idx_vast.NoEnviroSVC$idx[,2:10],est.se = idx_vast.NoEnviroSVC$idx.se[,2:10]),silent=TRUE)
	EnviroSVC.metrics = try(estimate.idx.metrics(true = simple.true.index,est = idx_vast.EnviroSVC$idx[,2:10],est.se = idx_vast.EnviroSVC$idx.se[,2:10]),silent=TRUE)

# bundle output
	vast_output = list(NoEnviro = idx_vast.NoEnviro,Enviro = idx_vast.Enviro,NoEnviroSVC = idx_vast.NoEnviroSVC,EnviroSVC = idx_vast.EnviroSVC)
	vast_metrics = list(NoEnviro = NoEnviro.metrics,Enviro = Enviro.metrics, NoEnviroSVC = NoEnviroSVC.metrics,EnviroSVC = EnviroSVC.metrics)
	vast_diagnostics = list(NoEnviro = try(idx_vast.NoEnviro[c(9,11,5,12,8,15,16,3)],silent=TRUE),Enviro = try(idx_vast.Enviro[c(9,11,5,12,8,15,16,3)],silent=TRUE), NoEnviroSVC = try(idx_vast.NoEnviroSVC[c(9,11,5,12,8,15,16,3)],silent=TRUE),EnviroSVC = try(idx_vast.EnviroSVC[c(9,11,5,12,8,15,16,3)],silent=TRUE))
	vast_list = list(vast_output = vast_output,vast_metrics = vast_metrics,vast_diagnostics = vast_diagnostics)
# save
	# vast_list = Data_Geostat
	save(vast_list,file="vast_list.RData")
	

# clear workspace
	rm(list=ls())

