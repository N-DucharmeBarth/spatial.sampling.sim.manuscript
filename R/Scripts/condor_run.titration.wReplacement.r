# clear workspace
	rm(list=ls())

# load packages
	library(ndd.vast.utils)
	library(rgdal)
	library(splines)
	data(skj.alt2019.shp)

#_______________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________
# helper function
# Modification of ofp-sam-vast-utils v1.0.0
# fit.vast fn
# https://github.com/PacificCommunity/ofp-sam-vast-utils/blob/77da1ff3f2084ac107d386556c3ec1ccc451f4e0/R/fit.vast.r
# Allows user to specify input.grid

fit.vast.mod = function(Data_Geostat,RunDir,SaveDir,save.output=FALSE,Q_ik = NULL,vf.re = FALSE,FieldConfig=c(Omega1 = 1, Epsilon1 = 1, Omega2 = 1, Epsilon2 = 1),RhoConfig=c(Beta1 = 0, Beta2 = 0, Epsilon1 = 0, Epsilon2 = 0),ObsModel_ez = c(1,3),fine_scale=TRUE,input.grid.res=1,crop.extrap.by.data=TRUE,knot_method = "grid",n_x=100,Version="VAST_v8_3_0",Method="Mesh",ADREPORT=TRUE,normalize_idx=FALSE,Xconfig_zcp=NULL,slim.output=FALSE,strata.sp,enviro,grid_bounds)
{
	A = proc.time()

	if(!dir.exists(RunDir)){dir.create(RunDir, recursive = TRUE)}
		origwd = getwd()
		setwd(RunDir)
	if(save.output)
	{
		if(!dir.exists(SaveDir)){dir.create(SaveDir, recursive = TRUE)}
	}


	# Decide which post-hoc calculations to include in VAST output
		Options = c('SD_site_density'=FALSE, 'SD_site_logdensity'=FALSE, 'Calculate_Range'=FALSE, 'SD_observation_density'=FALSE, 'Calculate_effective_area'=FALSE,
    'Calculate_Cov_SE'=FALSE, 'Calculate_Synchrony'=FALSE, 'Calculate_Coherence'=FALSE, 'Calculate_proportion'=FALSE, 'normalize_GMRF_in_CPP'=TRUE,
    'Calculate_Fratio'=FALSE, 'Estimate_B0'=FALSE, 'Project_factors'=FALSE, 'treat_nonencounter_as_zero'=FALSE, 'simulate_random_effects'=TRUE )
		Kmeans_Config = list( "randomseed"=1, "nstart"=100, "iter.max"=1e3 )
		if(vf.re)
		{
			OverdispersionConfig = c("Eta1"=1, "Eta2"=1)
			v_i = Data_Geostat$Vessel
		} else {
			OverdispersionConfig = c("Eta1"=0, "Eta2"=0)
			v_i = NULL
		}

	# Determine strata within the study region 
			strata.limits = data.frame('STRATA'="All_areas")
			Region = "user"

	# Set spatial extent and extrapolation grid
			grid_size_km = (110 * input.grid.res)^2 # the distance between grid cells for the 2D AR1 grid if Method == "Grid"
			if(missing(grid_bounds))
			{
				grid_bounds = c(floor(min(Data_Geostat[,c("Lat")])/5)*5,ceiling(max(Data_Geostat[,c("Lat")])/5)*5,floor(min(Data_Geostat[,c("Lon")])/5)*5,ceiling(max(Data_Geostat[,c("Lon")])/5)*5)
			}
	# Jim Thorson uses 110 as the approximation for the number of km in a degree of lat/lon


			input.grid = as.matrix(expand.grid(Lat = seq(from=grid_bounds[1],to=grid_bounds[2],by=input.grid.res),
					Lon = seq(from=grid_bounds[3],to=grid_bounds[4],by=input.grid.res), Area_km2 = grid_size_km))
			
			crs.en = paste0("+proj=tpeqd +lat_1=",mean(grid_bounds[1:2])," +lon_1=",round(grid_bounds[3] + (1/3)*abs(diff(grid_bounds[3:4])))," +lat_2=",mean(grid_bounds[1:2])," +lon_2=",round(grid_bounds[3] + (2/3)*abs(diff(grid_bounds[3:4])))," +datum=WGS84 +ellps=WGS84 +units=km +no_defs")
			crs.ll = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
			if(missing(strata.sp))
			{
				Extrapolation_List = make_extrapolation_info.ndd(Region = Region,strata.limits = strata.limits,input.grid = input.grid, crs.en = crs.en,crs.ll = crs.ll)
			} else {
				Extrapolation_List = make_extrapolation_info.ndd(Region = Region,strata.limits = strata.limits,input.grid = input.grid, crs.en = crs.en,crs.ll = crs.ll,strata.sp)
			}		 

		# punch-out all extrapolation grid cells that are on land and outside of "data region"
			extrap.df = as.data.frame(Extrapolation_List$Data_Extrap)
			sp::coordinates(extrap.df) = ~E_km + N_km
			sp::proj4string(extrap.df) = crs.en

		# find overlap with land
		# load land shape file
			data("pacific.coast")
			pacific.coast = sp::spTransform(pacific.coast,crs.en)
			over.index.land = which(is.na(sp::over(extrap.df,pacific.coast)))
		if(crop.extrap.by.data)
		{
			# find overlap with data hull
			smooth.hull = smooth.hull.sp(Data_Geostat[,c("Lon","Lat")],crs.ll=crs.ll,buffer.ll=2.5,d.scalar = 0.15)
			# plot(smooth.hull)
			smooth.hull.trans = sp::spTransform(smooth.hull,crs.en)
			# identify number of extrapolation cells (within data hull and not on land) corresponding to each knot
				over.index.region = which(!is.na(sp::over(extrap.df,smooth.hull.trans)))
				over.index = intersect(over.index.land,over.index.region)
		} else {
			over.index = over.index.land
		}
		
		# modify Extrapolation_List
			Extrapolation_List$a_el = data.frame(All_areas = Extrapolation_List$a_el[over.index,])
			Extrapolation_List$Data_Extrap = Extrapolation_List$Data_Extrap[over.index,]
			Extrapolation_List$Area_km2_x = Extrapolation_List$Area_km2_x[over.index]

		# with the modified spatial list function be sure to pass to Lon_i and Lat_i the lat and lon transformed to N_km and E_km using Convert_LL_to_EastNorth_Fn.ndd()
			seed = 123 
			# ll_to_EN = Convert_LL_to_EastNorth_Fn.ndd( Lon=Data_Geostat[,'Lon'], Lat=Data_Geostat[,'Lat'],crs.en = crs.en,crs.ll = crs.ll)
			Spatial_List = make_spatial_info.ndd(n_x=n_x, Lon_i=Data_Geostat[,'Lon'], Lat_i=Data_Geostat[,'Lat'], Extrapolation_List = Extrapolation_List, knot_method="grid", Method="Mesh",
												  grid_size_km=grid_size_km, grid_size_LL=input.grid.res, fine_scale=fine_scale, Network_sz_LL=NULL,
												  iter.max=1000, randomseed=seed, nstart=100, DirPath=SaveDir, Save_Results=save.output,crs.en = crs.en,crs.ll = crs.ll)

			Data_Geostat = cbind(Data_Geostat, knot_i = Spatial_List$knot_i)

		# prep info for plotting spatial domain of the model
			MapDetails_List = FishStatsUtils::make_map_info( "Region"=Region, "NN_Extrap"=Spatial_List$PolygonList$NN_Extrap, 
			"spatial_list"=Spatial_List, "Extrapolation_List"=Extrapolation_List )

		# allow for environmental covariates
			if(missing(enviro))
	  		{
	  			X_gtp = NULL
	  			X_itp = NULL
	  			Xconfig_zcp = NULL
	  		} else {
	  			enviro.format = make_covariates.ndd(formula = as.formula(enviro[["formula"]]),covariate_data=enviro[["covariate_data"]], Year_i=Data_Geostat[,'Year'], spatial_list=Spatial_List, extrapolation_list=Extrapolation_List)
	  			X_gtp = enviro.format$X_gtp
	  			X_itp = enviro.format$X_itp
	  		}

		# run the model
			TmbData = VAST::make_data("Version"=Version, "FieldConfig"=FieldConfig, "OverdispersionConfig"=OverdispersionConfig, 
							"RhoConfig"=RhoConfig, "ObsModel_ez"=ObsModel_ez, 
							"c_iz"=as.numeric(Data_Geostat[,'Spp'])-1, "b_i"=Data_Geostat[,'Response_variable'], 
							"a_i"=Data_Geostat[,'AreaSwept_km2'], "v_i"=v_i, "Q_ik" = Q_ik, "X_gtp" = X_gtp, "X_itp" = X_itp, "Xconfig_zcp"=Xconfig_zcp,
							"t_iz"=Data_Geostat[,'Year'], 
							"Options"=Options, "spatial_list" = Spatial_List )

			TmbList = VAST::make_model("TmbData"=TmbData, "RunDir"=RunDir, "Version"=Version, "RhoConfig"=RhoConfig, "loc_x"=Spatial_List$loc_x, "Method"=Spatial_List$Method)

		# Check parameters
			Obj = TmbList[["Obj"]]
			Obj$fn( Obj$par )
			Obj$gr( Obj$par )

		# Estimate fixed effects and predict random effects
			if(ADREPORT)
			{
				Opt = TMBhelper::fit_tmb( obj = Obj, lower = TmbList[["Lower"]], upper = TmbList[["Upper"]],
				getsd = TRUE, savedir = NULL, bias.correct = FALSE, newtonsteps = 3 )
				Report = Obj$report()
				Sdreport = TMB::sdreport(Obj)
				B = proc.time()
				fit.time = (B-A)[3]
				idx = Report$Index_cyl[1,,]
				idx.se = array( TMB::summary.sdreport(Sdreport)[which(rownames(TMB::summary.sdreport(Sdreport))=="Index_cyl"),2], dim=c(dim(Report$Index_cyl)), dimnames=list(NULL,NULL,NULL) )[1,,]

				if(normalize_idx)
				{
					for(j in 1:ncol(idx))
					{
						j.mean = mean(idx[,j])
						idx[,j] = idx[,j]/j.mean
						idx.se[,j] = idx.se[,j]/j.mean
					}
				}
				vast_output = list("idx"=idx,"idx.se"=idx.se, "Opt"=Opt, "Report"=Report,"Sdreport" = Sdreport, "TmbData"=TmbData, "Extrapolation_List"=Extrapolation_List, "fit.time"=fit.time,"MapDetails_List"=MapDetails_List )
		
			} else {
				Opt = TMBhelper::fit_tmb( obj = Obj, lower = TmbList[["Lower"]], upper = TmbList[["Upper"]],
				getsd = FALSE, savedir = NULL, bias.correct = FALSE, newtonsteps = 3 )
				Report = Obj$report()
				B = proc.time()
				fit.time = (B-A)[3]
				idx = Report$Index_cyl[1,,]
				if(normalize_idx)
				{
					idx = apply(idx,2,function(x)x/mean(x))
				}
				vast_output = list("idx"=idx, "Opt"=Opt, "Report"=Report, "TmbData"=TmbData, "Extrapolation_List"=Extrapolation_List, "fit.time"=fit.time,"MapDetails_List"=MapDetails_List )
		
			}

			if(slim.output)
			{
				vast_output = vast_output[names(vast_output) %in% c("idx","idx.se","fit.time")]
				vast_output = c(vast_output,Opt)
				vast_output$X_gtp = X_gtp
			}

			if(save.output)
			{
				save(vast_output,file=paste0(SaveDir,"vast_output.RData"))
			}
			setwd(origwd)

		return(vast_output)
}
#_______________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________

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

# define grid_bounds
	grid_bounds = c(-21,51,101,211)

# fit indices
	A = proc.time()
	idx_vast.NoEnviro.NoQ = try(fit.vast.mod(Data_Geostat.noQ,RunDir=paste0(getwd()),SaveDir=paste0(getwd()),save.output=FALSE,Q_ik = NULL,vf.re = FALSE,FieldConfig=c(Omega1 = 1, Epsilon1 = 1, Omega2 = 1, Epsilon2 = 1),RhoConfig=c(Beta1 = 0, Beta2 = 0, Epsilon1 = 0, Epsilon2 = 0),ObsModel_ez = c(1,3),fine_scale=FALSE,input.grid.res=1,crop.extrap.by.data=FALSE,knot_method = "grid",n_x=150,Version="VAST_v8_3_0",Method="Mesh",ADREPORT=TRUE,normalize_idx=FALSE,Xconfig_zcp=NULL,slim.output=TRUE,strata.sp=skj.alt2019.shp,grid_bounds=grid_bounds),silent=TRUE)
	gc(verbose=FALSE)
	B = proc.time()
	idx_vast.NoEnviro.Q = try(fit.vast.mod(Data_Geostat.Q,RunDir=paste0(getwd()),SaveDir=paste0(getwd()),save.output=FALSE,Q_ik = Q_ik,vf.re = TRUE,FieldConfig=c(Omega1 = 1, Epsilon1 = 1, Omega2 = 1, Epsilon2 = 1),RhoConfig=c(Beta1 = 0, Beta2 = 0, Epsilon1 = 0, Epsilon2 = 0),ObsModel_ez = c(1,3),fine_scale=FALSE,input.grid.res=1,crop.extrap.by.data=FALSE,knot_method = "grid",n_x=150,Version="VAST_v8_3_0",Method="Mesh",ADREPORT=TRUE,normalize_idx=FALSE,Xconfig_zcp=NULL,slim.output=TRUE,strata.sp=skj.alt2019.shp,grid_bounds=grid_bounds),silent=TRUE)
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

