

# Nicholas Ducharme-Barth
# 11/01/2021
# parse results from titration experiments (preferential sampling)
# iterate across data sets & calculate RMSE, bias, and coverage metrics
# create data.table of results
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
				t = try(summary( lm( est[,j] ~ true[,j] ) )$coef['true[, j]', 1],silent=TRUE)
				if(class(t)=="numeric"){bias.vec[j] = t}
				rm(list=c("t"))
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
				if(!(sum(is.na(UCI))==length(COVER)|sum(is.na(LCI))==length(COVER)))
				{
					for ( i in 1:length( COVER ) ) {
						if ( ( LCI[i]<=true[i,j] ) & ( true[i,j]<=UCI[i]) ){ COVER[i] = 1} 	
					}
					cover.vec[j] = ( sum( COVER )/ length( COVER ) ) * 100
				} 

				rm(list=c("COVER","UCI","LCI"))
			}

		return(data.frame(bias=bias.vec,mae=mae.vec,mrb=mrb.vec,rmse = rmse.vec,cover=cover.vec))		
	}


	pref.samp.titration = function(data.dt, ts.sampled = c(1,40), n.samp = 20000, cv = 0.15, random.seed = 123, prob.weight.exponent=0.5)
	# sample preferentially with respect to percieved abundance aka fisheries dependent
	{
		# set seed
			set.seed(random.seed)

		# subset to desired ts
			valid.ts = ts.sampled[1]:ts.sampled[2]	
			tmp.dt = data.dt[ts %in% valid.ts]

		# define number of samples taken per ts (total across ts must sum to n.samp)
			samp.per.ts = table(sample(valid.ts,n.samp,replace=TRUE))

		# define list to store samples from each ts
			ts.samp.list = as.list(rep(NA,length(valid.ts)))

		# iterate across ts and sample tmp.dt preferentially with respect to non-patchy abundance field (want samples to include zeros)
		# store the results in ts.samp.list
		# preferential sampling weight is sqrt(abundance) aka proportional to a power (prob.weight.exponent) of abundance
			for(i in 1:length(valid.ts))
			{
				# define sampling probability
					ts.dt = tmp.dt[ts == valid.ts[i]]
					prob.vector = (ts.dt$skj.noise)^prob.weight.exponent

				# sample rows
					ts.dt = ts.dt[sample(1:nrow(ts.dt),size=samp.per.ts[i],replace=TRUE,prob=prob.vector),][order(ts,lon,lat)]

				# add to ts.samp.list
					ts.samp.list[[i]] = ts.dt
			}

		# combine ts samples
			tmp.dt = do.call(rbind,ts.samp.list)

		# add observation error
			tmp.dt$response = tmp.dt$skj.noise.patchy + rnorm(length(tmp.dt$skj.noise.patchy),0,cv*tmp.dt$skj.noise.patchy)
			tmp.dt$response = ifelse(tmp.dt$response<0,0,tmp.dt$response)

		# trim columns
			tmp.dt = tmp.dt[,.(ts,yy,qq,yyqq,lon,lat,response,id.data,id.5x5,id.raster)]

		return(tmp.dt)

	}

#____________________________________________________________________________________________________________________________________________________
# build out data.structure
	n.weight = 8
	weight.vec = c(0,0.25,0.5,1,1.5,2,4,8)
	n.reps = 10

	rep.weight.mat = matrix(NA,nrow=n.weight*n.reps,ncol=3)
	colnames(rep.weight.mat) = c("data.set.id","replicate","pref.weight")
	iter = 1
	for(i in 1:10)
	{
		for(j in 1:length(weight.vec))
		{
			rep.weight.mat[iter,1] = iter
			rep.weight.mat[iter,2] = i
			rep.weight.mat[iter,3] = weight.vec[j]

			iter=iter+1
		}
	}

	titration.dt = as.data.table(expand.grid(data.set.id=1:(n.weight*n.reps),region=1:9,q=c(TRUE,FALSE))) %>%
					.[,replicate:=rep.weight.mat[match(data.set.id,rep.weight.mat[,1]),2]] %>%
					.[,pref.weight:=rep.weight.mat[match(data.set.id,rep.weight.mat[,1]),3]] %>%
					.[,rmse:=as.numeric(NA)] %>% .[,bias:=as.numeric(NA)] %>% .[,coverage:=as.numeric(NA)]

#____________________________________________________________________________________________________________________________________________________
# load true.data and calc metrics
	load("SimData/simple.true.index.RData")
	simple.true.index = apply(simple.true.index,2,function(x)x/mean(x))

	# iterate across simulated data.sets
	for(i in 1:(n.weight*n.reps))
	{
		# load
			load(paste0("Index/Titration_Preferential_wReplacement_v2/",i,".vast_list.RData"))

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
# calculate the proportion of each region sampled under each preferential sampling weight
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
		samp.dt = pref.samp.titration(data.dt, ts.sampled = c(1,120), n.samp = 60000, cv = 0.15, random.seed = 1, prob.weight.exponent=weight.vec[i])
		for(j in 1:9)
		{
			tmp = data.dt[ts==128,.(lon,lat)]
			sp::coordinates(tmp) = ~ lon + lat
			sp::proj4string(tmp) = sp::CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
			tmp$region = sp::over(tmp,strata.sp[j])
			tmp = subset(tmp,!is.na(region))

			tmp.samp = samp.dt[,.(lon,lat,ts)]
			sp::coordinates(tmp.samp) = ~ lon + lat
			sp::proj4string(tmp.samp) = sp::CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
			tmp.samp$region = sp::over(tmp.samp,strata.sp[j])
			tmp.samp = subset(tmp.samp,!is.na(region))
			tmp.samp = data.table::as.data.table(as.data.frame(tmp.samp))
			tmp.samp =  unique(tmp.samp[,.(lon,lat)])


			titration.dt[region==j&pref.weight==weight.vec[i], prop:=(nrow(tmp.samp)/nrow(tmp))*100]
			rm(list=c("tmp","tmp.samp"))
		}
		rm(list=c("samp.dt"))
	}


#____________________________________________________________________________________________________________________________________________________
# make plots

	plot.titration.dt = titration.dt %>% .[region %in% c(1,2,9)] %>% 
										 .[,Catchability:=factor(ifelse(q,"With catchability (Q)","Without catchability (noQ)"))] %>%
										 .[,Region := factor(as.character(region),levels=c("1","2","9"),labels=c("WCPO","Region 1","Region 8"))] %>%
										 .[,pref.weight.factor	:= factor(as.character(pref.weight),levels=as.character(sort(unique(c(pref.weight,seq(from=0,to=8,by=0.25))))))] %>%
										 .[,prop.block:=floor(prop/10)*10] %>% .[prop.block==100,prop.block:=90] %>% .[,prop.block:=factor(as.character(prop.block),levels=rev(c("0","10","20","30","40","50","60","70","80","90")),labels=rev(c("0-9","10-19","20-29","30-39","40-49","50-59","60-69","70-79","80-89","90-100")))]



	p = plot.titration.dt %>% ggplot2::ggplot() + 
     		ggplot2::geom_boxplot(ggplot2::aes(x=pref.weight.factor, y=rmse,fill=prop.block),outlier.color="gray60") + ggplot2::scale_x_discrete(drop=FALSE,breaks=0:8) +
     		ggplot2::geom_hline(yintercept = 0,size=1,linetype="longdash",color="black") + ggplot2::ylab("RMSE") + ggplot2::xlab("Preferential sampling weight") +
     		ggplot2::facet_grid( Catchability ~ Region) + ggthemes::theme_few(base_size = 20) + ggplot2::discrete_scale(drop=FALSE,aesthetics="fill",scale_name="Area sampled (%)",name="Area sampled (%)",palette=colorRampPalette(c("#af4448","#e57373","#f06292","#ba68c8","#9575cd","#7986cb","#64b5f6")))
    ggplot2::ggsave(filename=paste0("rmse.titration.prefv2.png"), plot = p, device = "png", path = "Plots/",
  			scale = 1.25, width = 16, height = 9, units = c("in"),
  			dpi = 300, limitsize = TRUE)

    p = plot.titration.dt %>% ggplot2::ggplot() + 
     		ggplot2::geom_boxplot(ggplot2::aes(x=pref.weight.factor, y=coverage,fill=prop.block),outlier.color="gray60") + ggplot2::ylim(0,100) + ggplot2::scale_x_discrete(drop=FALSE,breaks=0:8) +
     		ggplot2::geom_hline(yintercept = 50,size=1,linetype="longdash",color="black") + ggplot2::ylab("Coverage") + ggplot2::xlab("Preferential sampling weight") +
     		ggplot2::facet_grid( Catchability ~ Region) + ggthemes::theme_few(base_size = 20) + ggplot2::discrete_scale(aesthetics="fill",scale_name="Area sampled (%)",name="Area sampled (%)",palette=colorRampPalette(c("#af4448","#e57373","#f06292","#ba68c8","#9575cd","#7986cb","#64b5f6")))
     ggplot2::ggsave(filename=paste0("coverage.titration.prefv2.png"), plot = p, device = "png", path = "Plots/",
  			scale = 1.25, width = 12, height = 6.75, units = c("in"),
  			dpi = 300, limitsize = TRUE)

    p = plot.titration.dt %>% ggplot2::ggplot() + 
     		ggplot2::geom_boxplot(ggplot2::aes(x=pref.weight.factor, y=bias,fill=prop.block),outlier.color="gray60") + ggplot2::ylim(-0.5,1.5) + ggplot2::scale_x_discrete(drop=FALSE,breaks=0:8) +
     		ggplot2::geom_hline(yintercept = 1,size=1,linetype="longdash",color="black") + ggplot2::ylab("Bias") + ggplot2::xlab("Preferential sampling weight") +
     		ggplot2::facet_grid( Catchability ~ Region) + ggthemes::theme_few(base_size = 20) + ggplot2::discrete_scale(aesthetics="fill",scale_name="Area sampled (%)",name="Area sampled (%)",palette=colorRampPalette(c("#af4448","#e57373","#f06292","#ba68c8","#9575cd","#7986cb","#64b5f6")))
     ggplot2::ggsave(filename=paste0("bias.titration.prefv2.png"), plot = p, device = "png", path = "Plots/",
  			scale = 1.25, width = 12, height = 6.75, units = c("in"),
  			dpi = 300, limitsize = TRUE)

