

# Nicholas Ducharme-Barth
# 05/02/2020
# plot vessel catchability


# set working directory
	project.dir = "C:/Users/nicholasd/HOME/SPC/SPC_SAM/Geostats/spatial.sampling.sim.manuscript/"
	setwd(project.dir)

# load packages
	library(data.table)
	library(magrittr)
	library(ndd.vast.utils)
	data(skj.alt2019.shp)

# source functions & load background data
	source("R/Fn/fn.add.catchability.r")
	load("Background_Data/data.dt.RData")

# prep inputs
	s = "Random"
	r = 4
	save.id = r
	if(r<100){save.id = paste0("0",save.id)}
	if(r<10){save.id = paste0("0",save.id)}
	load(paste0("SimData/Simple120/",s,"/samp.dt.",save.id,".RData"))
	samp.dt$True_Abundance = data.dt$skj.noise.patchy[samp.dt$id.data]
	samp.dt = as.data.frame(samp.dt[,c("True_Abundance","ts","lon","lat")])
	colnames(samp.dt) = c("True_Abundance","Year","Lon","Lat")

# plot
	png(filename = "Plots/vessel.catchability.png", width = 10, height = 10, units = "in", res = 300)
	layout(matrix(c(1,2,4,1,3,4),nrow=3,ncol=2),heights=c(0.4,0.3,0.3))
	layout.show(4)
	tmp = add.catchability(samp.dt,seed = r,n.vessel.target = 90,new.entry.target=30,cv = 0.15,plot=TRUE)

	# plot dist poles per class
		class.vec = c("OS","DW")[rbinom(60000,1,0.5)+1]
		Poles = rpois(length(class.vec),ifelse(class.vec=="OS",15,25))
		hist(Poles, breaks=seq(from=0,to=50,by=2), col="#e53935",main="",xlab="Number of poles fished",ylab=c("Frequency (1000s)"),cex.lab=1.5,cex.axis=1.5,yaxt="n")
		hist(Poles[class.vec!="OS"], breaks=seq(from=0,to=50,by=2), col="#2196f3", add=TRUE)
		legend("topright",c("OS", "DW (1 & 2)"),fill=c("#e53935","#2196f3"),bty="n",cex=1.5)
		axis(2, c(0,1,2,3,4,5,6)*1000,c(0,1,2,3,4,5,6),las=1,cex.axis=1.5)

	# plot Q effect by npoles per class
		q.poles = ifelse(class.vec == "OS",((-(0.25*Poles-3.75)^2)/(1+(0.25*Poles-3.75)^2) + 0.5)*0.05,((-(0.25*Poles-6.25)^2)/(1+(0.25*Poles-6.25)^2) + 0.75)*0.1)
		q.dt = as.data.table(unique(data.frame(Class=class.vec,Poles=Poles,Q=q.poles*100)))[order(Class,Poles)]
		par(mar=c(5,1,1,5))
		plot(1,1,type="n",frame.plot=FALSE,xlim=c(0,50),ylim=c(-2.5,8),xlab="Number of poles fished",ylab="",cex.axis=1.5,cex.lab=1.5,las=1,yaxt="n")
		axis(4,at=c(-2,0,2,4,6,8),las=1,cex.axis=1.5)
		axis(4,at=mean(c(-2.5,8)),"Catchability effect (x 1e-02)",tick=FALSE,cex.axis=1.5,line=2)
		abline(h=0,col="gray70",lty=3,lwd=2)
		lines(q.dt[Class=="OS",.(Poles,Q)],col="#e53935",lwd=2)
		lines(q.dt[Class=="DW",.(Poles,Q)],col="#2196f3",lwd=2)

		par(mar=c(5,5,1,5))
		q.dt = as.data.table(tmp)[,.(q=mean(q)),by=Year] %>% .[,q:=scale(q)]
		plot(1, 1, type = "n", axes = TRUE, xlab = "Year", ylab = "",xlim=c(-0.1*max(q.dt$Year),1.1*max(q.dt$Year)),ylim=range(pretty(q.dt$q)),cex.lab=1.5,cex.axis=1.5,las=1)
		lines(q.dt,lwd=2)
		text(0,2,"Average set-specific catchability (relative)",adj=c(0,1),cex=1.5)

	dev.off()
