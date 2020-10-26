#' Calculate the nominal CPUE given a data set
#' @param data Data.frame expecting the following columns: Response_variable,Year,Lon,Lat
#' @param n.yr.rng Number of total timesteps from first year through last year in data
#' @param scale TRUE or FALSE, if true then return index as having mean of 1 and SE as the CV
#' @param strata.sp [Optional] If present, a shapefile containing the strata boundaries to calculate the indicies for
#' @return returns a data.frame or list of data.frames (one for each spatial strata) with 3 columns: "Year","Index","SE"
#' @importFrom rgeos gUnion
#' @importFrom sp coordinates
#' @importFrom sp proj4string
#' @importFrom sp over
#' @importFrom sp spTransform
#' @importFrom data.table as.data.table
#' @export
#' 

nominal.simple = function(data,n.yr.rng = length(range(data$Year)[1]:range(data$Year)[2]),scale=TRUE, strata.sp)
{
	# add recursive call if strata.sp is provided
	if(missing(strata.sp))
	{
		# create index storage structure
			index.df = as.data.frame(matrix(NA,nrow=n.yr.rng,ncol=3))
			colnames(index.df) = c("Year","Index","SE")
			index.df$Year = 1:n.yr.rng
		# calculate nominal
			dt = data.table::as.data.table(data)
			dt = dt[,.(Index=mean(Response_variable,na.rm=TRUE),SE = sd(Response_variable,na.rm=TRUE)),by=Year]
			dt$Year = as.numeric(as.character(dt$Year))
			index.df$Index[match(dt$Year,index.df$Year)] = dt$Index
			index.df$SE[match(dt$Year,index.df$Year)] = dt$SE

			if(scale)
			{
				mean.index = mean(index.df$Index,na.rm=TRUE)
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
			output.list[[i]] = nominal.simple(strata.data,n.yr.rng=n.yr.rng,scale=scale)
			rm(list=c("strata.data"))
		}
		return(output.list)
	}	
}