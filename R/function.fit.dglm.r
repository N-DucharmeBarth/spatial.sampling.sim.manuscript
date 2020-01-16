#' Fit a delta-lognormal GLM with no spatiotemporal interaction terms and all fixed effects
#' @param data Data.frame expecting the following columns: Response_variable,Year,Lon,Lat
#' @param agg.cell cell size in degrees for the spatial effect
#' @param formula.stem Base formula to use for both the binomial and lognormal components of the model
#' @param data.weighting TRUE or FALSE, if TRUE then reweight observations so that all spatiotemporal strata are equal
#' @param n.yr.rng Number of total timesteps from first year through last year in data
#' @param scale TRUE or FALSE, if true then return index as having mean of 1 and SE as the CV
#' @param strata.sp [Optional] If present, a shapefile containing the strata boundaries to calculate the indicies for
#' @return returns a data.frame or list of data.frames (one for each spatial strata) with 8 columns: "Year","nominal","bin.index","bin.se","pos.index","pos.se","index","se"
#' @importFrom data.table as.data.table
#' @importFrom speedglm speedglm
#' @importFrom boot inv.logit
#' @importFrom rgeos gUnion
#' @importFrom sp coordinates
#' @importFrom sp proj4string
#' @importFrom sp over
#' @export
#' 

	fit.dglm = function(data, agg.cell = 5, formula.stem = " ~ Year + agg.cell", data.weighting = TRUE,n.yr.rng = length(range(data$Year)[1]:range(data$Year)[2]),scale=TRUE, strata.sp)
	{
		# add recursive call if strata.sp is provided
			if(missing(strata.sp))
			{
				# add 0 to Year if less than 10 (need to do this to properly match up the indices at the end)
					data$Year = sapply(data$Year,function(x)if(as.numeric(x)<10){paste0("0",x)}else{x})
				# add 0 to Year if less than 100 (need to do this to properly match up the indices at the end)
					data$Year = sapply(data$Year,function(x)if(as.numeric(x)<100){paste0("0",x)}else{x})

				# add column for agg.cell
					if(agg.cell > 1)
					{
						cell.latlon.agg = apply(data[,c("Lon","Lat")],2,function(x)floor(x/agg.cell)*agg.cell)
						data$agg.cell = as.numeric(as.factor(apply(cell.latlon.agg,1,paste0,collapse=".")))
						rm(list=c("cell.latlon.agg"))
					} else {
						cell.latlon.agg = apply(data[,c("Lon","Lat")],2,function(x)floor(x))
						data$agg.cell = as.numeric(as.factor(apply(cell.latlon.agg,1,paste0,collapse=".")))
						rm(list=c("cell.latlon.agg"))
					}

				# add bin column
					data$bin = ifelse(data$Response_variable>0,1,0)

				if(data.weighting)
				{
					# calculate a weight vector for each observation given their spatiotemporal strata
					# need the number of observations per unique strata
						data$ts.cell = apply(data[,c("Year","agg.cell")],1,paste0,collapse=".")
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
					bin.model = speedglm::speedglm(bin.formula, family = binomial(link = "logit"), weights = bin.data$weight, data = bin.data, se = TRUE)
					pos.model = speedglm::speedglm(pos.formula, family = gaussian(link = "identity"), weights = pos.data$weight, data = pos.data, se = TRUE)

				# create index storage structure
					index.df = as.data.frame(matrix(NA,nrow=n.yr.rng,ncol=8))
					colnames(index.df) = c("Year","nominal","bin.index","bin.se","pos.index","pos.se","index","se")
					index.df$Year = 1:n.yr.rng
				# calculate nominal
					dt = data.table::as.data.table(data)
					dt = dt[,.(Nominal=mean(Response_variable,na.rm=TRUE)),by=Year]
					dt$Year = as.numeric(as.character(dt$Year))
					index.df$nominal[match(dt$Year,index.df$Year)] = dt$Nominal

				# create index
				# since we are rescaling the indices there is no need to bias correct
				# also can just use year effect since there are no interactions and we are using the rescaled indices
				# use 1st order taylor approximation for the variance back transformation
				# positive component 
					int.pos = pos.model$coefficients[grep("(Intercept)",names(pos.model$coefficients))]
					pred.pos = c(int.pos,int.pos+pos.model$coefficients[grep("Year",names(pos.model$coefficients))])
					pos.coef.index = match(names(pred.pos),names(pos.model$coefficients))

					names(pred.pos)[1] = paste0("Year",sort(unique(pos.data$Year))[1])
					index.df$pos.index[sapply(names(pred.pos),function(x)as.numeric(strsplit(x,"Year")[[1]][2]))] = exp(pred.pos)
					index.df$pos.se[sapply(names(pred.pos),function(x)as.numeric(strsplit(x,"Year")[[1]][2]))] = sqrt(exp(pred.pos)*summary(pos.model)$coefficients[pos.coef.index,2]^2)

				# binomial component
					int.bin = bin.model$coefficients[grep("(Intercept)",names(bin.model$coefficients))]
					pred.bin = c(int.bin,int.bin+bin.model$coefficients[grep("Year",names(bin.model$coefficients))])
					bin.coef.index = match(names(pred.bin),names(bin.model$coefficients))

					names(pred.bin)[1] = paste0("Year",sort(unique(bin.data$Year))[1])
					index.df$bin.index[sapply(names(pred.bin),function(x)as.numeric(strsplit(x,"Year")[[1]][2]))] = boot::inv.logit(pred.bin)
					index.df$bin.se[sapply(names(pred.bin),function(x)as.numeric(strsplit(x,"Year")[[1]][2]))] = sqrt(((exp(pred.bin))/((1+exp(pred.bin))^2))^2*(summary(bin.model)$coefficients[bin.coef.index,2])^2)

				# combine
					index.df$index = index.df$pos.index * index.df$bin.index
					# sd(XY) = (var(X)var(Y)+var(X)E(Y)^2 + var(Y)E(X)^2)^0.5
					index.df$se = ((index.df$bin.se)^2*(index.df$pos.se)^2+(index.df$bin.se)^2*(index.df$pos.index)^2+(index.df$pos.se)^2*(index.df$bin.index)^2)^0.5

					index.df = index.df[,c(Year,index,se)]
					colnames(index.df) = c("Year","Index","SE")
					
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
					output.list[[i]] = fit.dglm(strata.data, agg.cell = agg.cell, formula.stem = formula.stem, data.weighting = data.weighting, n.yr.rng=n.yr.rng,scale=scale)
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
# 	agg.cell = 5
# 	formula.stem = " ~ Year + agg.cell"
# 	data.weighting = TRUE
# 	n.yr.rng = length(range(data$Year)[1]:range(data$Year)[2])
# 	strata.sp = skj.alt2019.shp

# 	test.out = fit.dglm(data, agg.cell = 5, formula.stem = " ~ Year + agg.cell", data.weighting = FALSE,n.yr.rng = length(range(data$Year)[1]:range(data$Year)[2]), strata.sp)

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

# 		plot(1,1,type="n",xlim=c(0,n.yr.rng+1),ylim=c(-3,3))
# 		points(scale(true.index$Index),pch=16)
# 		lines(scale(test.out[[i]]$nominal),lty=2)
# 		lines(scale(test.out[[i]]$index),col="red")
# 	}

