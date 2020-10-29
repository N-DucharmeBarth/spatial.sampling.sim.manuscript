# clear workspace
	rm(list=ls())

# load packages
	library(ndd.vast.utils)
	library(rgdal)
	library(splines)
	data(skj.alt2019.shp)

# bring in data
	# format Data_Geostat
	load(list.files()[grep("samp.dt",list.files())])
	Data_Geostat.noQ = as.data.frame(samp.dt[,c("response","ts","lon","lat")])
	colnames(Data_Geostat.noQ) = c("Response_variable","Year","Lon","Lat")
	Data_Geostat.noQ$Spp = as.factor("skj")
	Data_Geostat.noQ$AreaSwept_km2 = 110^2
	Data_Geostat.noQ$Vessel = "missing"
	Data_Geostat.noQ = as.data.frame(Data_Geostat.noQ)

# add Catchability
	load("data.dt.RData")
	source("fn.add.catchability.r")
	samp.dt$True_Abundance = data.dt$skj.noise.patchy[samp.dt$id.data]
	samp.dt = as.data.frame(samp.dt[,c("True_Abundance","ts","lon","lat")])
	colnames(samp.dt) = c("True_Abundance","Year","Lon","Lat")
	Data_Geostat.Q = add.catchability(samp.dt,seed = as.numeric(strsplit(list.files()[grep("samp.dt",list.files())],"[.]")[[1]][3]),n.vessel.target = 90,new.entry.target=30,cv = 0.15,plot=FALSE)

# format Data_Geostat
	Data_Geostat.Q$Spp = as.factor("skj")
	Data_Geostat.Q$AreaSwept_km2 = 110^2
	Data_Geostat.Q$Vessel = as.character(Data_Geostat.Q$Vessel)
	Data_Geostat.Q$Gear_Config = as.character(Data_Geostat.Q$Gear_Config)
	Data_Geostat.Q = as.data.frame(Data_Geostat.Q)

# define Q_ik
	Q_ik = model.matrix(as.formula(paste0("Spp~Gear_Config+Class:bs(Poles)-1")),data=Data_Geostat.Q)
	Q_ik = Q_ik[,-1]

# fit indices
	A = proc.time()
	idx_vast.NoEnviro.NoQ = try(fit.vast(Data_Geostat.noQ,RunDir=paste0(getwd()),SaveDir=paste0(getwd()),save.output=FALSE,Q_ik = NULL,vf.re = FALSE,FieldConfig=c(Omega1 = 1, Epsilon1 = 1, Omega2 = 1, Epsilon2 = 1),RhoConfig=c(Beta1 = 0, Beta2 = 0, Epsilon1 = 0, Epsilon2 = 0),ObsModel_ez = c(1,3),fine_scale=FALSE,input.grid.res=1,crop.extrap.by.data=TRUE,knot_method = "grid",n_x=150,Version="VAST_v8_3_0",Method="Mesh",ADREPORT=TRUE,normalize_idx=FALSE,Xconfig_zcp=NULL,slim.output=TRUE,strata.sp=skj.alt2019.shp),silent=TRUE)
	gc(verbose=FALSE)
	B = proc.time()
	idx_vast.NoEnviro.Q = try(fit.vast(Data_Geostat.Q,RunDir=paste0(getwd()),SaveDir=paste0(getwd()),save.output=FALSE,Q_ik = Q_ik,vf.re = TRUE,FieldConfig=c(Omega1 = 1, Epsilon1 = 1, Omega2 = 1, Epsilon2 = 1),RhoConfig=c(Beta1 = 0, Beta2 = 0, Epsilon1 = 0, Epsilon2 = 0),ObsModel_ez = c(1,3),fine_scale=FALSE,input.grid.res=1,crop.extrap.by.data=TRUE,knot_method = "grid",n_x=150,Version="VAST_v8_3_0",Method="Mesh",ADREPORT=TRUE,normalize_idx=FALSE,Xconfig_zcp=NULL,slim.output=TRUE,strata.sp=skj.alt2019.shp),silent=TRUE)
	C = proc.time()

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

		# calc RMSE
			rmse.vec = rep(NA,ncol(true))
			for(j in 1:ncol(true))
			{
				rmse.vec[j] = sqrt(mean(  ( est[,j] - true[,j] )^2 ))
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

		return(data.frame(bias=bias.vec,mae=mae.vec,rmse = rmse.vec,cover=cover.vec))		
	}

	NoEnviro.NoQ.metrics = try(estimate.idx.metrics(true = simple.true.index,est = idx_vast.NoEnviro.NoQ$idx[,2:10],est.se = idx_vast.NoEnviro.NoQ$idx.se[,2:10]),silent=TRUE)
	NoEnviro.Q.metrics = try(estimate.idx.metrics(true = simple.true.index,est = idx_vast.NoEnviro.Q$idx[,2:10],est.se = idx_vast.NoEnviro.Q$idx.se[,2:10]),silent=TRUE)


# bundle output
	vast_output = list(NoEnviro.NoQ = idx_vast.NoEnviro.NoQ,NoEnviro.Q = idx_vast.NoEnviro.Q)
	vast_metrics = list(NoEnviro.NoQ = NoEnviro.NoQ.metrics,NoEnviro.Q = NoEnviro.Q.metrics)
	vast_diagnostics = list(NoEnviro.NoQ = try(idx_vast.NoEnviro.NoQ[c(9,11,5,12,8,15,16,3)],silent=TRUE),NoEnviro.Q = try(idx_vast.NoEnviro.Q[c(9,11,5,12,8,15,16,3)],silent=TRUE))
	vast_list = list(vast_output = vast_output,vast_metrics = vast_metrics,vast_diagnostics = vast_diagnostics)
# save
	# vast_list = Data_Geostat
	save(vast_list,file="vast_list.RData")
	

# clear workspace
	rm(list=ls())

