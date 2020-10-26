

#' Nicholas Ducharme-Barth
#' 23/01/2020
#' Add catchability column to simulated data 
#' @param data Data.frame expecting the following columns: True_Abundance,Year,Lon,Lat
#' @param seed random seed
#' @param n.vessel.target desired number of random vessels to simulate (will not equal this exactly)
#' @param new.entry.target desired number of vessels that enter the fishery at atime (will not equal this exactly)
#' @param cv The CV for the observation error
#' @param plot True or False, plot the vessel characteristics?
#' @return returns a data.frame  columns: Response_variable,Year,Lon,Lat,Vessel,Gear_Config,Class,Poles
#' @export

add.catchability = function(data,seed = 123,n.vessel.target = 20,new.entry.target=3,cv = 0.15,plot=FALSE)
{
	# Response_Variable = Q*E*A
	# q = a + b*vessel + b*gear_config + b*class * b*bs(Poles)
	# Q = q + e
	# e = N(0,cv*q)
	set.seed(seed)

	# generate vessels and activity profiles
		# initialize
		yr.rng = diff(range(data$Year)) + 1
		vessel.life = (yr.rng/n.vessel.target)*new.entry.target

		new.entries = rpois(1,new.entry.target)
		new.entries = ifelse(new.entries==0,1,new.entries)
		term.yr = rpois(new.entries,vessel.life)
		term.yr = ifelse(term.yr==0,1,term.yr)
		start.yr = rep(1,length(term.yr))

		# define storage structure
		vessel.df = as.matrix(data.frame(ID = NA, start.yr = start.yr, term.yr = term.yr, effect = NA, class = NA, gear_config = NA))

		# iteratively generate more vessels
		while(max(vessel.df[,"term.yr"])<max(data$Year))
		{
			new.entries = rpois(1,new.entry.target)
			new.entries = ifelse(new.entries==0,1,new.entries)
			# add some variability to start year
				if(length(term.yr)==1)
				{
					start.yr = term.yr
				} else {
					start.yr = floor(mean(term.yr))
					start.year.adj = sample(-ceiling(vessel.life*0.25):ceiling(vessel.life*0.25),new.entries,replace=TRUE)
					start.yr = floor(mean(term.yr)) + start.year.adj
					start.yr = ifelse(start.yr<=0,1,start.yr)
					if(min(start.yr)>max(term.yr))
					{
						start.yr[sample(1:length(start.yr),1)] = term.yr
					}
				}

			term.yr = rpois(new.entries,vessel.life)
			term.yr = ifelse(term.yr==0,1,term.yr)
			term.yr = start.yr + term.yr - 1
			term.yr = ifelse(term.yr>max(data$Year),max(data$Year),term.yr)
			# start.yr = rep(start.yr,length(term.yr))
			term.yr = term.yr[order(start.yr)]
			start.yr = start.yr[order(start.yr)]

			vessel.df = rbind(vessel.df,cbind(rep(NA,length(term.yr)),start.yr,term.yr,rep(NA,length(term.yr)),rep(NA,length(term.yr)),rep(NA,length(term.yr))))
		}
		vessel.df = as.data.frame(vessel.df)
		vessel.df$ID = 1:nrow(vessel.df)
		vessel.df$effect = sort(rnorm(nrow(vessel.df),0,0.05))
		vessel.df$class = c("OS","DW")[rbinom(nrow(vessel.df),size=1,prob=0.5)+1] 
		prob.gear_config = 1/(1+exp(-0.001*yr.rng*(vessel.df$start.yr-0.5*yr.rng)))
		vessel.df$gear_config = rbinom(length(prob.gear_config),size=1,prob=prob.gear_config)+1
		vessel.df$gear_config = ifelse(vessel.df$class == "OS",0,vessel.df$gear_config)

		if(plot)
		{
			class.cols = c("#e53935","#64b5f6","#1976d2")
			par(mar=c(5,5,1,5))
			plot(1, 1, type = "n", axes = FALSE, xlab = "", ylab = "",xlim=c(-0.1*max(data$Year),1.1*max(data$Year)),ylim=c(0,nrow(vessel.df)+1),cex.lab=1.5,cex.axis=1.5,las=1)
			axis(side=4, at = pretty(c(0,nrow(vessel.df)+1)),cex.axis=1.5,las=1)
			mtext("Number of unique vessels", side=4, line=3,cex=1.5)
			
			tot.activ.v = sapply(1:max(data$Year),function(x)sum(cbind( vessel.df[,2] <= x & vessel.df[,3] >= x )))
			polygon(c(1:max(data$Year),rev(1:max(data$Year))),c(tot.activ.v,rep(0,length(tot.activ.v))),border="gray65",col="gray90")

			par(new = TRUE)
			plot(1, 1, type = "n", axes = TRUE,frame.plot=TRUE, xlab = "Year", ylab = "Relative vessel power",xlim=c(-0.1*max(data$Year),1.1*max(data$Year)),ylim=c(-2.5,2.5),cex.lab=1.5,cex.axis=1.5,las=1)
			legend("topleft",c("OS","DW 1","DW 2","# unique vessels"),title=paste0("N = ",nrow(vessel.df)),pch=c(16,16,16,NA),lwd=c(2,2,2,NA),lty=c(2,1,1,NA),col=c(class.cols,NA),fill=c(NA,NA,NA,"gray90"),border=c(NA,NA,NA,"gray65"),cex=1.5,bty="n")
			# text(0,nrow(vessel.df)+1,paste0("N = ",nrow(vessel.df)),adj=c(0,1),cex=1.5)
			vessel.df$power = vessel.df$effect + ifelse(vessel.df$class == "OS",0.02,0.08) + c(-0.05,0,0.1)[vessel.df$gear_config+1]
			vessel.df$power.scale = scale(vessel.df$power)
			for(v in 1:nrow(vessel.df))
			{
				power = vessel.df$power.scale[v]
				cbind( vessel.df[,2] <= 1 & vessel.df[,3] >= 1 )
				# points(vessel.df[v,2:3],rep(power,2),pch=16,col=class.cols[vessel.df$gear_config[v]+1],cex=1.5*exp(5*vessel.df$effect[v]))
				points(vessel.df[v,2:3],rep(power,2),pch=16,col=class.cols[vessel.df$gear_config[v]+1],cex=1.5)
				# lines(vessel.df[v,2:3],rep(power,2),lty=ifelse(vessel.df$gear_config[v]==0,2,1),lwd=3*exp(5*vessel.df$effect[v]),col=class.cols[vessel.df$gear_config[v]+1])
				lines(vessel.df[v,2:3],rep(power,2),lty=ifelse(vessel.df$gear_config[v]==0,2,1),lwd=2,col=class.cols[vessel.df$gear_config[v]+1])
				rm(list=c("power"))
			}
		}
		
	# assign each sample to a vessel based on if it was active during that year
		data$Vessel = NA
		for(t in 1:max(data$Year))
		{
			activ.v = which(cbind( vessel.df[,2] <= t & vessel.df[,3] >= t )==TRUE)
			data.target = which(data$Year == t)
			data$Vessel[data.target] = sample(activ.v,length(data.target),replace=TRUE)
			rm(list=c("activ.v","data.target"))
		}
		data$Gear_Config = vessel.df$gear_config[data$Vessel]
		data$Class = vessel.df$class[data$Vessel]
		data$Poles = rpois(nrow(data),ifelse(data$Class=="OS",15,25))

	# define q data.frame
		q.df = data.frame(int=rep(1,nrow(data)),vessel = vessel.df$effect[data$Vessel], gear_config = c(-0.05,0,0.1)[data$Gear_Config+1],class=rep(0,nrow(data)))
		q.df$poles = ifelse(data$Class == "OS",((-(0.25*data$Poles-3.75)^2)/(1+(0.25*data$Poles-3.75)^2) + 0.5)*0.05,((-(0.25*data$Poles-6.25)^2)/(1+(0.25*data$Poles-6.25)^2) + 0.75)*0.1)
 		q.df$q = rowSums(q.df)
 		q.df$Q = q.df$q + rnorm(nrow(q.df),0,cv*q.df$q)

 	# sample
 		data$Response_variable = data$True_Abundance*q.df$Q

 	# return
 		return(data[,c("Response_variable","Year","Lon","Lat","Vessel","Gear_Config","Class","Poles")])
}
