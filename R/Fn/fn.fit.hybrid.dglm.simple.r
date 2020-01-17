#' Fit a hybrid delta-lognormal GLM with spatiotemporal interaction terms. If there are missing spatial temporal strata these will be imputed using a simple delta GLM with no interactions
#' @param data Data.frame expecting the following columns: Response_variable,Year,Lon,Lat
#' @param n_x number of spatial strata to partition observations across for spatial effect
#' @param formula.stem Base formula to use for both the binomial and lognormal components of the model
#' @param data.weighting TRUE or FALSE, if TRUE then reweight observations so that all spatiotemporal strata are equal
#' @param n.yr.rng Number of total timesteps from first year through last year in data
#' @param seed Set seed for random number generator
#' @param scale TRUE or FALSE, if true then return index as having mean of 1 and SE as the CV
#' @param strata.sp [Optional] If present, a shapefile containing the strata boundaries to calculate the indicies for
#' @return returns a data.frame or list of data.frames (one for each spatial strata) with 3 columns: "Year","Index","SE"
#' @importFrom data.table as.data.table
#' @importFrom speedglm speedglm
#' @importFrom boot inv.logit
#' @importFrom rgeos gUnion
#' @importFrom sp coordinates
#' @importFrom sp proj4string
#' @importFrom sp over
#' @importFrom sp spTransform
#' @importFrom imputeTS na_ma
#' @importFrom ndd.vast.utils make_extrapolation_info.ndd
#' @importFrom ndd.vast.utils make_spatial_info.ndd 
#' @export
#' 

	fit.hybrid.dglm.simple = function(data, n_x = 10, formula.stem = " ~ Year * knot", data.weighting = TRUE,n.yr.rng = length(range(data$Year)[1]:range(data$Year)[2]), seed = 123, scale=TRUE, strata.sp, target.strata)
	{
		# add recursive call if strata.sp is provided
			if(missing(strata.sp))
			{
				# add 0 to Year if less than 10 (need to do this to properly match up the indices at the end)
					data$Year = sapply(data$Year,function(x)if(as.numeric(x)<10){paste0("0",x)}else{x})
				# add 0 to Year if less than 100 (need to do this to properly match up the indices at the end)
					data$Year = sapply(data$Year,function(x)if(as.numeric(x)<100){paste0("0",x)}else{x})


				# partition observations to each knot (this is now the new spatial strata)
				# Set spatial extent and extrapolation grid
					grid_size_km = (110 * 1)^2 # the distance between grid cells for the 2D AR1 grid if Method == "Grid"
					grid_bounds = c(floor(min(data[,c("Lat")])/5)*5,ceiling(max(data[,c("Lat")])/5)*5,floor(min(data[,c("Lon")])/5)*5,ceiling(max(data[,c("Lon")])/5)*5)	
					input.grid = as.matrix(expand.grid(Lat = seq(from=grid_bounds[1],to=grid_bounds[2],by=1),Lon = seq(from=grid_bounds[3],to=grid_bounds[4],by=1), Area_km2 = grid_size_km))
					crs.en = paste0("+proj=tpeqd +lat_1=",mean(grid_bounds[1:2])," +lon_1=",round(grid_bounds[3] + (1/3)*abs(diff(grid_bounds[3:4])))," +lat_2=",mean(grid_bounds[1:2])," +lon_2=",round(grid_bounds[3] + (2/3)*abs(diff(grid_bounds[3:4])))," +datum=WGS84 +ellps=WGS84 +units=km +no_defs")
					crs.ll = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"


					Extrapolation_List = ndd.vast.utils::make_extrapolation_info.ndd(Region = "User", strata.limits = data.frame('STRATA'="All_areas"),input.grid=input.grid,crs.en=crs.en,crs.ll=crs.ll)
					# modify extrapolation list to only include points within the region and not on land 
						extrap.df = as.data.frame(Extrapolation_List$Data_Extrap)
						sp::coordinates(extrap.df) = ~E_km + N_km
						sp::proj4string(extrap.df) = crs.en
						library(ndd.vast.utils)
						data("pacific.coast")
						pacific.coast = sp::spTransform(pacific.coast,crs.en)
						over.index.land = which(is.na(sp::over(extrap.df,pacific.coast)))

						if(missing(target.strata))
						{
							over.index = over.index.land
						} else {
							strata.trans = sp::spTransform(target.strata,crs.en)
							over.index.region = which(!is.na(sp::over(extrap.df,strata.trans)))
							over.index = intersect(over.index.land,over.index.region)
						}

						Extrapolation_List$a_el = data.frame(All_areas = Extrapolation_List$a_el[over.index,])
						Extrapolation_List$Data_Extrap = Extrapolation_List$Data_Extrap[over.index,]
						Extrapolation_List$Area_km2_x = Extrapolation_List$Area_km2_x[over.index]

						Spatial_List = ndd.vast.utils::make_spatial_info.ndd(n_x=10, Lon_i=data[,'Lon'], Lat_i=data[,'Lat'], Extrapolation_List = Extrapolation_List, knot_method="grid", Method="Mesh",
												  grid_size_km=grid_size_km, grid_size_LL=1, fine_scale=FALSE, Network_sz_LL=NULL,
												  iter.max=1000, randomseed=seed, nstart=100, DirPath="/dummy", Save_Results=FALSE,crs.en = crs.en,crs.ll = crs.ll)

						data$knot = Spatial_List$knot_i

				# add bin column
					data$bin = ifelse(data$Response_variable>0,1,0)

				# identify if there are any missing strata & if so, fill them in with results from a simple glm
					t.strat = 1:n.yr.rng
					t.strat = sapply(t.strat,function(x)if(as.numeric(x)<10){paste0("0",x)}else{x})
					t.strat = sapply(t.strat,function(x)if(as.numeric(x)<100){paste0("0",x)}else{x})

					s.strat = 1:n_x

					data$strat.id = paste0(data$Year,".",data$knot)
					u.strat = expand.grid(Year=t.strat,knot=s.strat)
					u.strat$strat.id = paste0(u.strat$Year,".",u.strat$knot)

					if(length(unique(data$strat.id))<nrow(u.strat))
					{
						# use simple glms to predict across all strata
							data.tmp = data
							data.tmp$knot = as.factor(as.character(data.tmp$knot))
							# fit glms
								bin.tmp = speedglm::speedglm(bin ~ Year + knot, family = binomial(link = "logit"), data = data.tmp, se = FALSE)
								pos.tmp = speedglm::speedglm(log(Response_variable) ~ Year + knot, family = gaussian(link = "identity"), data = data.tmp[which(data.tmp$bin==1),], se = FALSE)
							# predict
								u.strat$bin.pred = NA
								u.strat$pos.pred = NA
								for(i in 1:nrow(u.strat))
								{
									if(u.strat$Year[i] == "001")
									{
										bin.year.offset = 0
										pos.year.offset = 0
									} else {
										bin.year.offset = bin.tmp$coefficients[grep(paste0("Year",u.strat$Year[i]),names(bin.tmp$coefficients))]
										pos.year.offset = pos.tmp$coefficients[grep(paste0("Year",u.strat$Year[i]),names(pos.tmp$coefficients))]
									}
									if(u.strat$knot[i] == "1")
									{
										bin.knot.offset = 0
										pos.knot.offset = 0
									} else {
										bin.knot.offset = bin.tmp$coefficients[grep(paste0("knot",u.strat$knot[i]),names(bin.tmp$coefficients))]
										pos.knot.offset = pos.tmp$coefficients[grep(paste0("knot",u.strat$knot[i]),names(pos.tmp$coefficients))]
									}

									tmp.bin.pred = unname(boot::inv.logit(bin.tmp$coefficients[1] + bin.year.offset + bin.knot.offset))
									if(length(tmp.bin.pred)==0){tmp.bin.pred=NA}
									tmp.pos.pred = exp(pos.tmp$coefficients[1] + pos.year.offset + pos.knot.offset)
									if(length(tmp.pos.pred)==0){tmp.pos.pred=NA}


									u.strat$bin.pred[i] = tmp.bin.pred
									u.strat$pos.pred[i] = tmp.pos.pred
									# clean-up
										rm(list=c("bin.year.offset","pos.year.offset","bin.knot.offset","pos.knot.offset","tmp.bin.pred","tmp.pos.pred"))
								}

							u.strat$pred = u.strat$bin.pred * u.strat$pos.pred

							# create Walter's folly fantasy fill matrix
								walters.m = subset(u.strat,knot==1)$pred
								for(i in 2:max(u.strat$knot))
								{
									walters.m = cbind(walters.m,subset(u.strat,knot==i)$pred)
								}
								rownames(walters.m) = t.strat
								colnames(walters.m) = s.strat

							# impute across time within a knot for any strata that were unable to be predicted by the basic dglm
								walters.m = apply(walters.m,2,imputeTS::na_ma,k=3)

							# identify which strata are missing 
								missing.strata = u.strat[-which(u.strat$strat.id %in% unique(data$strat.id)),]

							# look-up values in walters.m & append missing to data
								for(i in 1:nrow(missing.strata))
								{
									data = rbind(data,data.frame(Response_variable = rep(walters.m[which(rownames(walters.m)==missing.strata$Year[i]),which(colnames(walters.m)==missing.strata$knot[i])],1),
										Year = rep(missing.strata$Year[i],1),
										Lon = rep(NA,1),
										Lat = rep(NA,1),
										knot = rep(missing.strata$knot[i],1),
										bin = rep(ifelse(walters.m[which(rownames(walters.m)==missing.strata$Year[i]),which(colnames(walters.m)==missing.strata$knot[i])]>0,1,0),1),
										strat.id = rep(missing.strata$strat.id[i],1)))
								}
					}

				if(data.weighting)
				{
					# calculate a weight vector for each observation given their spatiotemporal strata
					# need the number of observations per unique strata
						data$ts.cell = apply(data[,c("Year","knot")],1,paste0,collapse=".")
						strata.obs = as.data.frame(cbind(names(table(data$ts.cell)),table(data$ts.cell)))
						colnames(strata.obs) = c("strata","N")
						strata.obs$N = as.numeric(as.character(strata.obs$N))

						nrow.data = nrow(data)
						nrow.strata.obs = nrow(strata.obs)
						data$weight = nrow.data/nrow.strata.obs * (1/strata.obs[match(data$ts.cell,strata.obs$strata),]$N) # Equation 13a from Campbell 2015
						data$weight = data$weight/sum(data$weight)

					# split the data into binary catch and positive catch
						bin.data = data
							
						pos.index = which(bin.data$bin == 1)
						pos.data = data[pos.index,]

					# redefine strata weights for positve catch data
						pos.strata.obs = as.data.frame(cbind(names(table(pos.data$ts.cell)),table(pos.data$ts.cell)))
						colnames(pos.strata.obs) = c("strata","N")
						pos.strata.obs$N = as.numeric(as.character(pos.strata.obs$N))

						pos.nrow.data = nrow(pos.data)
						pos.nrow.strata.obs = nrow(pos.strata.obs)
						pos.data$weight = pos.nrow.data/pos.nrow.strata.obs * (1/pos.strata.obs[match(pos.data$ts.cell,pos.strata.obs$strata),]$N) # Equation 13a from Campbell 2015
						pos.data$weight = pos.data$weight/sum(pos.data$weight)
				} else {
					data$weight = 1/nrow(data)
					# split the data into binary catch and positive catch
						bin.data = data
						pos.index = which(bin.data$bin == 1)
						pos.data = data[pos.index,]
				}

				# convert appropriate columns to factor
					include.columns = strsplit(formula.stem, "\\s+")[[1]][which(nchar(strsplit(formula.stem, "\\s+")[[1]])>1)]
					for(i in 1:length(include.columns))
					{
						bin.data[,include.columns[i]] = as.factor(as.character(bin.data[,include.columns[i]]))
						pos.data[,include.columns[i]] = as.factor(as.character(pos.data[,include.columns[i]]))
					}

				# identify variables that have less than 2 factors and strip them out of the data going into the model
					# bin.factors = names(Filter(is.factor, bin.data))
					# bin.factors.nlevels = sapply(bin.factors,function(x)nlevels(bin.data[,x]))
					# bin.factors = bin.factors[which(bin.factors.nlevels>1)]

					# pos.factors = names(Filter(is.factor, pos.data))
					# pos.factors.nlevels = sapply(pos.factors,function(x)nlevels(pos.data[,x]))
					# pos.factors = pos.factors[which(pos.factors.nlevels>1)]

				# trim unneeded columns
					bin.data = bin.data[,c("bin",include.columns,"weight")]
					pos.data = pos.data[,c("Response_variable",include.columns,"weight")]

				# create model formulas and ammend if needed
					bin.formula =  as.formula(paste0("bin",formula.stem))
					pos.formula =  as.formula(paste0("log(Response_variable)",formula.stem))

				# fit both the positive and binomial components of the model
					bin.model = glm(bin.formula, family = binomial(link = "logit"), weights = bin.data$weight, data = bin.data)
					pos.model = glm(pos.formula, family = gaussian(link = "identity"), weights = pos.data$weight, data = pos.data)

				# make predictions across all strata
					new.data = u.strat[,c("Year","knot")]
					new.data$knot = as.factor(as.character(new.data$knot))
					bin.pred = predict(bin.model, newdata = new.data,type = c("response"),se.fit=TRUE)
					pos.pred = predict(pos.model, newdata = new.data,type = c("response"),se.fit=TRUE)
					pos.pred$fit = exp(pos.pred$fit)
					pos.pred$se.fit = sqrt(pos.pred$fit*pos.pred$se.fit^2)
					new.data$pred = bin.pred$fit * pos.pred$fit
					# sd(XY) = (var(X)var(Y)+var(X)E(Y)^2 + var(Y)E(X)^2)^0.5
					new.data$se = ((bin.pred$se.fit)^2*(pos.pred$se.fit)^2+(bin.pred$se.fit)^2*(pos.pred$fit)^2+(pos.pred$se.fit)^2*(bin.pred$fit)^2)^0.5
					new.data$cv = new.data$se/new.data$pred

				# get area of each knot & take spatial average to create index
					knot.area = as.vector(Spatial_List$a_xl)
					knot.area = knot.area/sum(knot.area)
					new.data$area = knot.area[as.numeric(as.character(new.data$knot))]
					new.data$pred_x_area = new.data$area * new.data$pred
					new.data$var_x_area = new.data$area^2 * new.data$se^2
				 	walters.dt = data.table::as.data.table(new.data)
				 	walters.dt = walters.dt[,.(Index = sum(pred_x_area),Var = sum(var_x_area)),by=Year][order(Year)]
				 	walters.dt$SE = sqrt(walters.dt$Var) 
				 	walters.dt$CV = walters.dt$SE/walters.dt$Index

					index.df = as.data.frame(walters.dt[,.(Year,Index,SE)])
					index.df$Year = 1:n.yr.rng
					
				if(scale)
					{
						mean.index = mean(index.df$Index)
						index.df$Index = index.df$Index/mean.index
						index.df$SE = index.df$SE/mean.index 
					}

				# return
					return(index.df)
			} else {
				n.substrata = length(strata.sp)
				full.reg = strata.sp[1]
				for(i in 2:length(strata.sp)){full.reg = rgeos::gUnion(full.reg,strata.sp[i])}
				strata.sp = rbind(full.reg,strata.sp)
				new_IDs = c("all.strata",paste0("sub.strata.",1:n.substrata))
				for (i in 1:length(slot(strata.sp, "polygons")))
				{
				  slot(slot(strata.sp, "polygons")[[i]], "ID") = new_IDs[i]
				}

				# create storage.structure
					output.list = as.list(rep(NA,length(strata.sp)))

				strata.points = data[,c("Lon", "Lat")]
				sp::coordinates(strata.points) = c("Lon", "Lat")
	  			sp::proj4string(strata.points) = sp::proj4string(strata.sp)
				for(i in 1:length(strata.sp))
				{
					strata.data = data[which(!is.na(sp::over(strata.points,strata.sp[i]))),]
					output.list[[i]] = fit.hybrid.dglm.simple(strata.data, n_x = 10, formula.stem = formula.stem, data.weighting = data.weighting, n.yr.rng=n.yr.rng,seed=seed,scale=scale,target.strata=strata.sp[i])
					rm(list=c("strata.data"))
				}
				return(output.list)
			}		 
	}


# # test
# 	setwd("C:/Users/nicholasd/HOME/SPC/SPC_SAM/Geostats/spatial.sampling.sim.manuscript/")	

# # load packages
# 	library(ndd.vast.utils)
# 	data(skj.alt2019.shp)
# 	i = 1
# 		save.id = i
# 				if(i<100){save.id = paste0("0",save.id)}
# 				if(i<10){save.id = paste0("0",save.id)}
# 			# bring in data
# 				load(paste0("SimData/Simple/Fixed/samp.dt.",save.id,".RData"))
# 			# format data
# 				data = as.data.frame(samp.dt[,c("response","ts","lon","lat")])
# 				colnames(data) = c("Response_variable","Year","Lon","Lat")
# ### Function defaults for testing 
# 	data = data
# 	n_x = 10
# 	formula.stem = " ~ Year * knot"
# 	data.weighting = TRUE
# 	n.yr.rng = length(range(data$Year)[1]:range(data$Year)[2])
# 	seed = 123
# 	scale = TRUE
# 	strata.sp = skj.alt2019.shp


# 	test.out = fit.hybrid.dglm.simple(data, n_x = 10, formula.stem = " ~ Year * knot", data.weighting = FALSE,n.yr.rng = length(range(data$Year)[1]:range(data$Year)[2]), seed = 123,scale=TRUE, strata.sp)

# 	load("Background_Data/data.dt.RData")
# 	data.dt = data.table::as.data.table(data.dt)
# 	n.substrata = length(strata.sp)
# 	full.reg = strata.sp[1]
# 	for(i in 2:length(strata.sp)){full.reg = rgeos::gUnion(full.reg,strata.sp[i])}
# 	strata.sp = rbind(full.reg,strata.sp)

# 	par(mfrow=c(3,3))
# 	for(i in 1:length(test.out))
# 	{
# 		extrap.points = as.data.frame(data.dt)
# 		sp::coordinates(extrap.points) = c("lon", "lat")
# 		sp::proj4string(extrap.points) = sp::proj4string(strata.sp)
# 		extrap.points$valid = sp::over(extrap.points,strata.sp[i])
# 		true.dt = data.table::as.data.table(extrap.points@data)
# 		true.index = true.dt[ts %in% 1:40 & !is.na(valid),.(Index=mean(skj.noise.patchy)),by=ts]

# 		plot(true.index$Index/mean(true.index$Index),pch=16,ylim=c(0.75,1.25))
# 		lines(test.out[[i]]$Index,col="red")
# 	}

