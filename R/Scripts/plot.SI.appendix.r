

# Nicholas Ducharme-Barth
# 07/01/2021
# Plot SEAPODYM output for Supplementary Information section


#______________________________________________________________________________________________________________________________________________________________________________
# set working directory
	setwd("C:/Users/nicholasd/HOME/SPC/SPC_SAM/Geostats/spatial.sampling.sim.manuscript/")

#______________________________________________________________________________________________________________________________________________________________________________
# load packages
	library(magrittr)	

#______________________________________________________________________________________________________________________________________________________________________________
# load data
	load("Background_Data/data.dt.RData")
	data.dt = data.table::as.data.table(data.dt)
	load("Background_Data/sst.storage.df.RData")
	load("Background_Data/nino.df.RData")

#______________________________________________________________________________________________________________________________________________________________________________
# plot average skipjack abundance field from seapodym
	tmp.dt = data.dt %>% .[,.(skj=mean(skj)),by=.(lon,lat)] %>% .[,skj:=skj/mean(skj)]

	p = ggplot2::ggplot(tmp.dt, ggplot2::aes(x=lon,y=lat,fill=skj)) + ggplot2::coord_fixed() +  
		 ggplot2::xlab("Longitude") + ggplot2::ylab("Latitude") + ggplot2::ggtitle("Avg. SEAPODYM abundance field") +
	     ggplot2::geom_tile() + ggplot2::scale_fill_viridis_c("Relative Skipjack biomass") +
	     ggthemes::theme_few(base_size = 20)
    ggplot2::ggsave(filename=paste0("avg.skj.seapodym.png"), plot = p, device = "png", path = "Plots/",
  			scale = 1.25, width = 16, height = 9, units = c("in"),
  			dpi = 300, limitsize = TRUE)


#______________________________________________________________________________________________________________________________________________________________________________
# plot center of seapodym gravity versus enso
	nino.dt = data.table::as.data.table(nino.df)

	png(filename = "Plots/center.of.gravity.seapodym.png", width = 9, height = 9, units = "in", res = 300)
	par(mar = c(5, 5, 1, 5),mfrow=c(2,1))
	tmp.dt = data.dt[,.(skj=sum(skj),lon_skj=sum(skj*lon)),by=yyqq] %>% .[,mean_lon:=lon_skj/skj]
	plot(nino.dt[,.(yrqtr,index.scale)],type="l",lwd=2,xlim=range(nino.dt$yrqtr),ylab="ENSO",xlab="")
	legend("topleft",c("ENSO","SEAPODYM center of gravity"),col=c("black","red"),lwd=3,bty="n")
	par(new=TRUE)
	plot(tmp.dt[,.(yyqq,mean_lon)],col="red",lwd=2, type = "l", xaxt = "n", yaxt = "n",ylab = "", xlab = "",xlim=range(nino.dt$yrqtr))
	axis(side = 4)
	mtext("Center of gravity (Longitude)", side = 4, line = 3)

	tmp.dt = data.dt[,.(skj=sum(skj),lat_skj=sum(skj*lat)),by=yyqq] %>% .[,mean_lat:=lat_skj/skj]
	plot(nino.dt[,.(yrqtr,index.scale)],type="l",lwd=2,xlim=range(nino.dt$yrqtr),ylab="ENSO",xlab="Year")
	par(new=TRUE)
	plot(tmp.dt[,.(yyqq,mean_lat)],col="red",lwd=2, type = "l", xaxt = "n", yaxt = "n",ylab = "", xlab = "",xlim=range(nino.dt$yrqtr))
	axis(side = 4)
	mtext("Center of gravity (Latitude)", side = 4, line = 3)
	dev.off()


#______________________________________________________________________________________________________________________________________________________________________________
# plot sst vs seapodym skipjack
	sst.dt = data.table::as.data.table(sst.storage.df) %>% .[x >= min(data.dt$lon) & x <= max(data.dt$lon) & y >= min(data.dt$lat) & y <= max(data.dt$lat)] %>%
	         .[,yy:=as.numeric(substr(time,1,4))] %>% .[,qq:=as.numeric(substr(time,5,6))] %>% .[,yyqq:=yy+(qq-1)/4] %>% .[yyqq %in% unique(data.dt$yyqq)]

	ts.vec = sort(unique(data.dt$yyqq))
	dt.list = as.list(rep(NA,length(ts.vec)))

	for(i in 1:length(ts.vec))
	{
		tmp.dt = data.dt[yyqq==ts.vec[i]]
		tmp.sst = sst.dt[yyqq==ts.vec[i]]
		NN = RANN::nn2( data=tmp.sst[,c("x","y")], query=tmp.dt[,c("lon","lat")], k=1 )
		tmp.dt$sst = tmp.sst$sst[NN$nn.idx[,1]]
		dt.list[[i]] = tmp.dt
		rm(list=c("tmp.dt","tmp.sst","NN"))
	}

	tmp.dt=data.table::rbindlist(dt.list) %>% .[,skj:=skj/mean(skj)]

	p = ggplot2::ggplot(tmp.dt, ggplot2::aes(x=sst,y=skj)) +  
		 ggplot2::xlab("SST") + ggplot2::ylab("Relative skipjack abundance") + 
	     # ggplot2::geom_point(alpha=1/100) + 
	     ggplot2::geom_hex(bins = 50) + ggplot2::scale_fill_viridis_c("Count") +
	     ggplot2::geom_smooth(color="hotpink",size=2) +
	     ggthemes::theme_few(base_size = 20)
    ggplot2::ggsave(filename=paste0("sst.skj.seapodym.png"), plot = p, device = "png", path = "Plots/",
  			scale = 1.25, width = 9, height = 9, units = c("in"),
  			dpi = 300, limitsize = TRUE)

