

# Nicholas Ducharme-Barth
# 29/10/2020
# parse results from titration experiments
# iterate across data sets & calculate RMSE, bias, and coverage metrics
# create data.table of results
# data.set.id, replicate, max.dist, q, region, rmse, mrb, coverage # avg.sampling.rate
# plot results

#____________________________________________________________________________________________________________________________________________________
# set project directory
	project.dir = "C:/Users/nicholasd/HOME/SPC/SPC_SAM/Geostats/spatial.sampling.sim.manuscript/"
	setwd(project.dir)

#____________________________________________________________________________________________________________________________________________________
# load packages
	library(data.table)
	library(magrittr)

#____________________________________________________________________________________________________________________________________________________
# define helper function
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

		# calc MRB
			mrb.vec = rep(NA,ncol(true))
			for(j in 1:ncol(true))
			{
				mrb.vec[j] = mean(est[,j] - true[,j])
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
				UCI = est[,j] + 0.67449*est.se[,j] 
				LCI = est[,j] - 0.67449*est.se[,j] 
				for ( i in 1:length( COVER ) ) {
					if ( ( LCI[i]<=true[i,j] ) & ( true[i,j]<=UCI[i]) ){ COVER[i] = 1} 	
				}
				cover.vec[j] = ( sum( COVER )/ length( COVER ) ) * 100
				rm(list=c("COVER","UCI","LCI"))
			}

		return(data.frame(bias=bias.vec,mae=mae.vec,mrb=mrb.vec,rmse = rmse.vec,cover=cover.vec))		
	}

#____________________________________________________________________________________________________________________________________________________
# build out data.structure
	n.dists = 8
	min.dist.vec = seq(from=1000,to=9999,length.out=n.dists)
	n.reps = 10

	rep.dist.mat = matrix(NA,nrow=n.dists*n.reps,ncol=3)
	colnames(rep.dist.mat) = c("data.set.id","replicate","max.dist")
	iter = 1
	for(i in 1:10)
	{
		for(j in 1:length(min.dist.vec))
		{
			rep.dist.mat[iter,1] = iter
			rep.dist.mat[iter,2] = i
			rep.dist.mat[iter,3] = min.dist.vec[j]

			iter=iter+1
		}
	}

	titration.dt = as.data.table(expand.grid(data.set.id=1:(n.dists*n.reps),region=1:9,q=c(TRUE,FALSE))) %>%
					.[,replicate:=rep.dist.mat[match(data.set.id,rep.dist.mat[,1]),2]] %>%
					.[,max.dist:=rep.dist.mat[match(data.set.id,rep.dist.mat[,1]),3]] %>%
					.[,rmse:=as.numeric(NA)] %>% .[,bias:=as.numeric(NA)] %>% .[,coverage:=as.numeric(NA)]

#____________________________________________________________________________________________________________________________________________________
# load true.data and calc metrics
	load("SimData/simple.true.index.RData")
	simple.true.index = apply(simple.true.index,2,function(x)x/mean(x))

	# iterate across simulated data.sets
	for(i in 1:(n.dists*n.reps))
	{
		# load
			load(paste0("Index/Titration/",i,".vast_list.RData"))

		# calc
			NoQ.met = estimate.idx.metrics(true = simple.true.index,est = vast_list$vast_output$NoEnviro.NoQ$idx[,2:10],est.se = vast_list$vast_output$NoEnviro.NoQ$idx.se[,2:10])
			Q.met = estimate.idx.metrics(true = simple.true.index,est = vast_list$vast_output$NoEnviro.Q$idx[,2:10],est.se = vast_list$vast_output$NoEnviro.Q$idx.se[,2:10])

		# put in storage structure
			for(j in 1:9)
			{
				titration.dt[data.set.id==i&region==j&q==TRUE, rmse:=Q.met[j,4]]
				titration.dt[data.set.id==i&region==j&q==TRUE, bias:=Q.met[j,1]]
				titration.dt[data.set.id==i&region==j&q==TRUE, coverage:=Q.met[j,5]]
				titration.dt[data.set.id==i&region==j&q==FALSE, rmse:=NoQ.met[j,4]]
				titration.dt[data.set.id==i&region==j&q==FALSE, bias:=NoQ.met[j,1]]
				titration.dt[data.set.id==i&region==j&q==FALSE, coverage:=NoQ.met[j,5]]
			}

		# clean-up
			rm(list=c("vast_list","NoQ.met","Q.met"))
	}

#____________________________________________________________________________________________________________________________________________________
# calculate the proportion of each region sampled under each radius
	load("Background_Data/data.dt.RData")
	library(ndd.vast.utils)
	data(skj.alt2019.shp)
	
			# define union of model region
				strata.sp = skj.alt2019.shp
				n.substrata = length(strata.sp)
				full.reg = strata.sp[1]
				for(i in 2:length(strata.sp)){full.reg = rgeos::gUnion(full.reg,strata.sp[i])}
				strata.sp = rbind(full.reg,strata.sp)
				new_IDs = c("all.strata",paste0("sub.strata.",1:n.substrata))
				for (i in 1:length(slot(strata.sp, "polygons")))
				{
				  slot(slot(strata.sp, "polygons")[[i]], "ID") = new_IDs[i]
				}

	for(i in 1:8)
	{
		for(j in 1:9)
		{
			tmp = data.dt[ts==128,.(lon,lat,dist.jp.km)]
			sp::coordinates(tmp) = ~ lon + lat
			sp::proj4string(tmp) = sp::CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
			tmp$region = sp::over(tmp,strata.sp[j])
			tmp = subset(tmp,!is.na(region))
			titration.dt[region==j&max.dist==min.dist.vec[i], prop:=(nrow(subset(tmp,dist.jp.km<=min.dist.vec[i]))/nrow(tmp))*100]
			rm(list=c("tmp"))
		}
	}


#____________________________________________________________________________________________________________________________________________________
# make plots

	plot.titration.dt = titration.dt %>% .[region %in% c(1,2,9)] %>% .[,Distance:=factor(ceiling(max.dist/100)*100,levels=sort(unique(ceiling(max.dist/100)*100)))] %>%
										 .[,Catchability:=factor(ifelse(q,"With catchability (Q)","Without catchability (noQ)"))] %>%
										 .[,Region := factor(as.character(region),levels=c("1","2","9"),labels=c("WCPO","Region 1","Region 8"))] %>%
										 .[,prop.block:=floor(prop/10)*10] %>% .[prop.block==100,prop.block:=90] %>% .[,prop.block:=factor(as.character(prop.block),levels=rev(c("0","10","20","30","40","50","60","70","80","90")),labels=rev(c("0-9","10-19","20-29","30-39","40-49","50-59","60-69","70-79","80-89","90-100")))]


	p = plot.titration.dt %>% ggplot2::ggplot() + 
     		ggplot2::geom_boxplot(ggplot2::aes(x=Distance, y=rmse,fill=prop.block),outlier.color="gray60") + 
     		ggplot2::geom_hline(yintercept = 0,size=1,linetype="longdash",color="black") + ggplot2::ylab("RMSE") + ggplot2::xlab("Maximum distance from JP sampled") +
     		ggplot2::facet_grid( Catchability ~ Region) + ggthemes::theme_few() + ggplot2::discrete_scale(aesthetics="fill",scale_name="Area sampled (%)",name="Area sampled (%)",palette=colorRampPalette(c("#af4448","#e57373","#f06292","#ba68c8","#9575cd","#7986cb","#64b5f6")))
    ggplot2::ggsave(filename=paste0("rmse.titration.png"), plot = p, device = "png", path = "Plots/",
  			scale = 1.25, width = 12, height = 6.75, units = c("in"),
  			dpi = 300, limitsize = TRUE)

    p = plot.titration.dt %>% ggplot2::ggplot() + 
     		ggplot2::geom_boxplot(ggplot2::aes(x=Distance, y=coverage,fill=prop.block),outlier.color="gray60") + ggplot2::ylim(0,100) +
     		ggplot2::geom_hline(yintercept = 50,size=1,linetype="longdash",color="black") + ggplot2::ylab("Coverage") + ggplot2::xlab("Maximum distance from JP sampled") +
     		ggplot2::facet_grid( Catchability ~ Region) + ggthemes::theme_few() + ggplot2::discrete_scale(aesthetics="fill",scale_name="Area sampled (%)",name="Area sampled (%)",palette=colorRampPalette(c("#af4448","#e57373","#f06292","#ba68c8","#9575cd","#7986cb","#64b5f6")))
     ggplot2::ggsave(filename=paste0("coverage.titration.png"), plot = p, device = "png", path = "Plots/",
  			scale = 1.25, width = 12, height = 6.75, units = c("in"),
  			dpi = 300, limitsize = TRUE)

    p = plot.titration.dt %>% ggplot2::ggplot() + 
     		ggplot2::geom_boxplot(ggplot2::aes(x=Distance, y=bias,fill=prop.block),outlier.color="gray60") + ggplot2::ylim(-0.5,1.5) +
     		ggplot2::geom_hline(yintercept = 1,size=1,linetype="longdash",color="black") + ggplot2::ylab("Bias") + ggplot2::xlab("Maximum distance from JP sampled") +
     		ggplot2::facet_grid( Catchability ~ Region) + ggthemes::theme_few() + ggplot2::discrete_scale(aesthetics="fill",scale_name="Area sampled (%)",name="Area sampled (%)",palette=colorRampPalette(c("#af4448","#e57373","#f06292","#ba68c8","#9575cd","#7986cb","#64b5f6")))
     ggplot2::ggsave(filename=paste0("bias.titration.png"), plot = p, device = "png", path = "Plots/",
  			scale = 1.25, width = 12, height = 6.75, units = c("in"),
  			dpi = 300, limitsize = TRUE)

