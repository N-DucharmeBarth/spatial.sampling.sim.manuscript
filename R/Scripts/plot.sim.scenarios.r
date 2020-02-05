

# Nicholas Ducharme-Barth
# 04/02/2020
# Plot the effort scenarios, spatial structure & abundance

# set working directory
	project.dir = "C:/Users/nicholasd/HOME/SPC/SPC_SAM/Geostats/spatial.sampling.sim.manuscript/"
	setwd(project.dir)

# load packages
	library(ndd.vast.utils)
	library(data.table)
	library(sp)
	data(skj.alt2019.shp)
	data(pacific.coast)

# load data
	load("Background_Data/data.dt.RData")

# define plotting space
	png(filename = "Plots/sim.scenarios.png", width = 7, height = 10, units = "in", res = 300)
	layout(matrix(c(rep(c( 1, 1, 2, 2),1),
					rep(c( 3, 3, 4, 4),1),
					rep(c( 5, 6, 9, 10),1),
					rep(c( 8, 7,12,11),1),
					rep(c(13,13,13,13),1),
					rep(c(14,14,15,15),1)),nrow=6,ncol=4,byrow=TRUE),heights=c(0.825,0.1,0.60,0.60,0.05,0.825),widths=c(1,1,1,1))
	# layout.show(14)

# 1) Plot abundance
	# convert to rasters & plot
	data.abundance = data.dt[ts==1,.(Value=mean(skj.noise.patchy,na.rm=TRUE)),by=.(lon,lat)]
	data.abundance$index = ceiling(sqrt(data.abundance$Value)) - min(ceiling(sqrt(data.abundance$Value))) + 1
	abundance.cols = c("white",viridis::viridis(max(data.abundance$index)-1))

	par(mar=c(0.5,0.5,0.5,0))
	plot(data.abundance[,.(lon,lat)],type="n",axes=FALSE,xlab="",ylab="",frame.plot=FALSE,asp=1)
	points(data.abundance[,.(lon,lat)],pch=15,cex=1,col=abundance.cols[data.abundance$index])
	plot(skj.alt2019.shp,lwd=2.5,add=TRUE)
	plot(pacific.coast,border="gray50",col="gray90",add=TRUE)
	text(100,57,"Skipjack",adj=c(0,1),cex=1.5)
	text(100,52,"abundance",adj=c(0,1),cex=1.5)
	# make background circles so region numbers are readable
	for(i in 1:8)
	{
		radius = 3.5
		center_x = rgeos::gCentroid(skj.alt2019.shp,byid=TRUE)@coords[i,1]
		center_y = rgeos::gCentroid(skj.alt2019.shp,byid=TRUE)@coords[i,2]
		theta = seq(0, 2 * pi, length = 360)
		x = radius * cos(theta) + center_x
		y = radius * sin(theta) + center_y
		x = c(x,x[1])
		y = c(y,y[1])
		polygon(x,y,border=NA,col=scales::alpha("white",0.95))
	}
	text(rgeos::gCentroid(skj.alt2019.shp,byid=TRUE)@coords[,1],rgeos::gCentroid(skj.alt2019.shp,byid=TRUE)@coords[,2],1:8,cex=1.5)

# 2) Plot Preferential Effort distribution
	load("SimData/Simple120/Preferential/samp.dt.001.RData")
	pref.pts = samp.dt[ts %in% 1:12,.(lon,lat)]
	# calc number other samples within 500km, this will be color index

	neighbor.vec = rep(NA,nrow(pref.pts))
	neighbor.radius = 111*5.64
	for(i in 1:length(neighbor.vec))
	{
		neighbor.vec[i] = length(which(geosphere::distHaversine(as.vector(as.matrix(pref.pts)[i,]),as.matrix(pref.pts))/1000 <= neighbor.radius))-1
	}
	pref.pts$neighbors = neighbor.vec
	quant.probs = seq(from=0,to=0.99,length.out=15)
	neighbor.quants = quantile(pref.pts$neighbors,probs=quant.probs)
	smooth.hull.list = as.list(rep(NA,length(quant.probs)))
	for(i in 1:length(quant.probs))
	{
		smooth.hull.list[[i]] = smooth.hull.sp(as.data.frame(pref.pts[neighbors>neighbor.quants[i],.(lon,lat)]),crs.ll="+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0",buffer.ll=0.05,d.scalar = 0.15)
		smooth.hull.list[[i]]@polygons[[1]]@ID = paste0("smooth.hull.",i)
	}

	smooth.hull = do.call(rbind,smooth.hull.list)
	pref.cols = colorRampPalette(c("#64b5f6","#1976d2","#283593"))(length(quant.probs))
	

	par(mar=c(0.5,0,0.5,0.5))
	plot(data.abundance[,.(lon,lat)],type="n",axes=FALSE,xlab="",ylab="",frame.plot=FALSE,asp=1)
	plot(smooth.hull,border=NA,col=scales::alpha(pref.cols,0.15),add=TRUE)
	plot(pacific.coast,border="gray50",col="gray90",add=TRUE)
	text(210,45,"Preferential",adj=c(1,1),cex=1.5)
	text(210,40,"sampling scenario",adj=c(1,1),cex=1.5)
	text(210,35,"effort distribution",adj=c(1,1),cex=1.5)

# 3 & 4) Titles for Fixed and Rotational closures
	par(mar=c(0,0.5,0.5,0))
	plot(1,1,type="n",axes=FALSE,xlab="",ylab="",frame.plot=FALSE)
	text(1,1,"Fixed closure scenario",cex=1.5)
	abline(h=1.40,lwd=2)

	par(mar=c(0,0,0.5,0.5))
	plot(1,1,type="n",axes=FALSE,xlab="",ylab="",frame.plot=FALSE)
	text(1,1,"Rotating closure scenario",cex=1.5)  
	abline(h=1.40,lwd=2)

# 5,6,7,8) Fixed closure
	load("SimData/Simple120/Fixed/samp.dt.001.RData")

	# par(mfrow=c(2,2))
	for(i in 1:4)
	{
		if(i ==1){
			par(mar=c(0.25,1.5,0.25,0.25))
		} else if(i ==2){
			par(mar=c(0.25,0.25,0.25,1.5))
		} else if(i ==4){
			par(mar=c(0.25,1.5,0.25,0.25))
		} else {
			par(mar=c(0.25,0.25,0.25,1.5))
		}
		plot(samp.dt[ts==1,.(lon,lat)],type="n",axes=FALSE,xlab="",ylab="",frame.plot=FALSE,asp=1)
		if(i == 3)
		{
			polygon(x=c(100,210,210,100,100),y=c(-20,-20,20,20,-20),border="#c62828",density=20,angle=45,col="#c62828")
		}
		plot(pacific.coast,border="gray50",col="gray75",add=TRUE)
		points(samp.dt[ts==i,.(lon,lat)][sample(1:nrow(samp.dt[ts==i]),200)],pch=16,col="#1976d2",cex=0.5)
		text(100,60,paste0("Q",i),adj=c(0,1),cex=1.5)	
	}

# 9,10,11,12) Fixed closure
	load("SimData/Simple120/Rotating/samp.dt.001.RData")

	# par(mfrow=c(2,2))
	# Q1 the NE quadrant is closed (>155 & >15)
	# Q2 the SE quadrant is closed (>155 & <15)
	# Q3 the SW quadrant is closed (<155 & <15)
	# Q4 the NW quadrant is closed (<155 & >15)
	for(i in 1:4)
	{
		if(i ==1){
			par(mar=c(0.25,1.5,0.25,0.25))
		} else if(i ==2){
			par(mar=c(0.25,0.25,0.25,1.5))
		} else if(i ==4){
			par(mar=c(0.25,1.5,0.25,0.25))
		} else {
			par(mar=c(0.25,0.25,0.25,1.5))
		}
		plot(samp.dt[ts%in%1:4,.(lon,lat)],type="n",axes=FALSE,xlab="",ylab="",frame.plot=TRUE,asp=1)
		if(i == 1){
			polygon(x=c(155,210,210,155,155),y=c(15,15,50,50,15),border="#c62828",density=20,angle=45,col="#c62828")
		} else if(i == 2){
			polygon(x=c(155,210,210,155,155),y=c(-20,-20,15,15,-20),border="#c62828",density=20,angle=45,col="#c62828")
		} else if(i == 3){
			polygon(x=c(100,155,155,100,100),y=c(-20,-20,15,15,-20),border="#c62828",density=20,angle=45,col="#c62828")
		} else {
			polygon(x=c(100,155,155,100,100),y=c(15,15,50,50,15),border="#c62828",density=20,angle=45,col="#c62828")
		}
		plot(pacific.coast,border="gray50",col="gray90",add=TRUE)
		points(samp.dt[ts==i,.(lon,lat)][sample(1:nrow(samp.dt[ts==i]),200)],pch=16,col="#1976d2",cex=0.5)
		text(100,60,paste0("Q",i),adj=c(0,1),cex=1.5)	
	}

# 13) Divider
	par(mar=c(0.5,0.5,0.5,0.5))
	plot(1,1,type="n",axes=FALSE,xlab="",ylab="",frame.plot=FALSE)
	abline(h=1,lwd=2)

# 14) Distance from Japan
	distance.vec = c(2500,5000,7500)
	smooth.hull.list = as.list(rep(NA,length(distance.vec)))
	for(i in 1:length(distance.vec))
	{
		smooth.hull.list[[i]] = smooth.hull.sp(as.data.frame(data.dt[dist.jp.km<distance.vec[i],.(lon,lat)]),crs.ll="+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0",buffer.ll=0.01,d.scalar = 0.05)
		smooth.hull.list[[i]]@polygons[[1]]@ID = paste0("smooth.hull.",i)
	}

	smooth.hull = do.call(rbind,smooth.hull.list)

	par(mar=c(0.5,0.5,0.5,0.5))
	plot(data.dt[,.(lon,lat)],type="n",axes=FALSE,xlab="",ylab="",frame.plot=FALSE,asp=1)
	# points(data.dt[,.(lon,lat)],pch=15,cex=1,col=abundance.cols[data.abundance$index])
	plot(smooth.hull,border=NA,col=scales::alpha(rev(c("#64b5f6","#1976d2","#283593")),0.25),add=TRUE)
	plot(pacific.coast,border="gray50",col="gray90",add=TRUE)
	text(100,57,"Distance",adj=c(0,1),cex=1.5)
	text(100,52,"from Japan",adj=c(0,1),cex=1.5)
	text(150,25,paste0(2500," km"),cex=1.5)
	text(165,10,paste0(5000," km"),cex=1.5)
	text(180,-5,paste0(7500," km"),cex=1.5)

# 15) Plot expansion and contraction
	r.mat = matrix(NA,nrow=100,ncol=120)
	for(i in 1:100)
	{
		set.seed(i)
		min.dist = 1000
		max.dist = 10000
		valid.ts = 1:120
		r.vec = rep(NA,length(valid.ts))
		ts.quantiles =  floor(quantile(valid.ts,probs = seq(0, 1, 0.125)))
		r.vec[which(valid.ts == ts.quantiles[1]):which(valid.ts == ts.quantiles[2])] = max.dist
		r.vec[which(valid.ts == ts.quantiles[8]):which(valid.ts == ts.quantiles[9])] = min.dist
		bridge.index = which(is.na(r.vec))

		bridge.n = length(bridge.index)
		bridge.times = seq(0, 1, length.out=bridge.n)
		bridge.target = max.dist - min.dist # Constraint at time=1
		bridge.dW = rnorm(bridge.n,0,0.5*min.dist)
		bridge.W = cumsum(bridge.dW)
		bridge.B = bridge.W + bridge.times * (bridge.target - bridge.W[bridge.n])   # The Brownian bridge from (0,0) to (1,target)
		bridge.range = 2*(bridge.target - 0)
		bridge.B = (bridge.B - 0) %% bridge.range
		bridge.B = pmin(bridge.B, bridge.range-bridge.B) + 0
		r.vec[bridge.index] = rev(min.dist + bridge.B)
		r.mat[i,] = r.vec
	}

	quant.r = apply(r.mat/1000,2,quantile,probs=c(0.01,0.5,0.9))
	par(mar=c(3,5,0.5,0.5))
	plot(1,1,type="n",xlab="",ylab="Max distance km (1000s)", xlim=range(pretty(seq(from=1979,length.out=120,by=0.25))),ylim=range(pretty(c(0,13))),cex.lab=1.5,cex.axis=1.5,las=1,yaxt="n")
	axis(2,c(0,2.5,5,7.5,10),as.character(c(0,2.5,5,7.5,10)),las=1,cex.axis=1.5,cex=1.5)
	polygon(c(seq(from=1979,length.out=120,by=0.25),rev(seq(from=1979,length.out=120,by=0.25))),c(rev(quant.r[1,]),quant.r[3,]),border=NA,col=scales::alpha("#e57373",0.5))
	lines(seq(from=1979,length.out=120,by=0.25),rev(quant.r[2,]),lwd=3,col="#c62828")
	polygon(c(seq(from=1979,length.out=120,by=0.25),rev(seq(from=1979,length.out=120,by=0.25))),c(quant.r[1,],rev(quant.r[3,])),border=NA,col=scales::alpha("#64b5f6",0.5))
	lines(seq(from=1979,length.out=120,by=0.25),quant.r[2,],lwd=3,col="#283593")
	abline(h=0)
	abline(h=c(2.5,5,7.5),col="black",lwd=2,lty=3)
	legend("top",c("Contraction","Expansion"),lwd=3,col=c("#283593","#c62828"),bty="n",ncol=2,cex=1.5,title="Sampling scenario")

	dev.off()

