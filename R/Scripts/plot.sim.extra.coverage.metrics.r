

# Nicholas Ducharme-Barth
# 03/02/2020
# plot the results

#____________________________________________________________________________________________________________________________________________________________________________
# load packages
	library(magrittr)

#____________________________________________________________________________________________________________________________________________________________________________
# set working directory
	project.dir = "C:/Users/nicholasd/HOME/SPC/SPC_SAM/Geostats/spatial.sampling.sim.manuscript/"
	setwd(project.dir)

#____________________________________________________________________________________________________________________________________________________________________________
# source
	source("R/Fn/fn.plot.sim.results.summary.r")

#____________________________________________________________________________________________________________________________________________________________________________
# load
	load("SimData/simple.true.index.RData")
	load("Index/ResultsDF/ts.df.MRwR.RData")
	load("Index/ResultsDF/nominal.df.wR.RData")

#____________________________________________________________________________________________________________________________________________________________________________
# plot extra coverage metric boxplots
	load("Index/ResultsDF/metric.MRwR.df.RData")

	metric.dt = data.table::as.data.table(metric.df) %>% .[Region=="all"&mgc<=1e-04] %>% .[Model=="Enviro",Model:="E"] %>%
			 .[Model=="NoEnviro",Model:="noE"] %>% .[Model=="EnviroSVC",Model:="ESVC"] %>%
			 .[Model=="NoEnviroSVC",Model:="noESVC"] %>% .[Catchability=="Q",Model := paste0(Catchability,"-",Model)] %>% .[,Model:=factor(Model,levels=c("noE","noESVC","E","ESVC","Q-noE","Q-noESVC","Q-E","Q-ESVC"))] %>%
			 .[Metric %in% c("cover.50","cover.70","cover.95")] %>% .[,Metric:=factor(Metric,levels=c("cover.50","cover.70","cover.95"),labels=c("Coverage: 50% CI","Coverage: 70% CI","Coverage: 95% CI"))]

			u.scenario = sort(unique(metric.dt$Scenario))
			u.metric = rev(sort(unique(metric.dt$Metric)))
			u.model = sort(unique(metric.dt$Model))
			# ylim.list = list(RMSD=range(pretty(metric.dt[Metric=="RMSD"]$Value,na.rm=TRUE)),Bias=range(pretty(metric.dt[Metric=="Bias"]$Value,na.rm=TRUE)),Coverage=range(pretty(metric.dt[Metric=="Coverage"]$Value,na.rm=TRUE)))
			yaxis.list = list(Coverage.50=c(0,50,100),Coverage.70=c(0,70,100),Coverage.50=c(0,95,100))
			x.offset = 0.2
			fill.vec = c("#ffebee","#ffcdd2","#ef9a9a","#e57373","#e3f2fd","#bbdefb","#90caf9","#64b5f6") #100,400,300,600
			line.vec = c(rep("#7f0000",4),rep("#002f6c",4))

			png(filename = "Plots/sim.metrics.extra.coverage.wR.png", width = 16, height = 9, units = "in", res=300,pointsize = 12, bg = "white")
			layout(matrix(c(34,28,29,30,31,32,33,19,1,2,3,4,5,6,20,7,8,9,10,11,12,21,13,14,15,16,17,18,35,22,23,24,25,26,27),nrow=5,ncol=7,byrow = TRUE),widths=c(0.25,rep(1,6)),heights=c(0.1,rep(1,3),0.4))
			layout.show(35)
			for(j in 1:length(u.metric))
			{
				for(i in 1:length(u.scenario))
				{
					xlim = c(0,10)
					ylim = range(yaxis.list[[u.metric[j]]])
					if(j==3)
					{
						ylim.text = ylim[1] - 0.075*abs(diff(ylim))
						ylim[1] = ylim[1] - 0.1*abs(diff(ylim))
					}

					par(mar=c(0.5,0.25,0.5,0.25))
					if(j==3)
					{
						plot(1,1,type="n",xlab="",xlim=xlim,ylab="",ylim=ylim,xaxt="n",yaxt="n",xaxs="i",yaxs="i",frame.plot=FALSE)
						polygon(x=c(0,10,10,0),y=c(range(yaxis.list[[u.metric[j]]])[1],range(yaxis.list[[u.metric[j]]])[1],range(yaxis.list[[u.metric[j]]])[2],range(yaxis.list[[u.metric[j]]])[2]))
						axis(1,at=c(1,2,3,4,6,7,8,9),u.model,las=2,cex.axis=1.25, col = NA, col.ticks = 1)
						abline(h=yaxis.list[[u.metric[j]]][2],lty=3)
					} else {
						plot(1,1,type="n",xlab="",xlim=xlim,ylab="",ylim=ylim,xaxt="n",yaxt="n",xaxs="i",yaxs="i")
						abline(h=yaxis.list[[u.metric[j]]][2],lty=3)
					}
					if(i==1)
					{
						axis(2,yaxis.list[[u.metric[j]]],las=1,cex.axis=1.25)
						axis(2,mean(yaxis.list[[u.metric[j]]]),u.metric[j],las=0,cex.axis=1.5,tick=FALSE,line=2)
					}
					if(j==1){axis(3,at=mean(xlim),u.scenario[i],tick=FALSE,cex.axis=1.5,line=-0.5)}
					if(j==2){abline(h=1,lty=3)}

					for(m in 1:length(u.model))
					{
						tmp.dt = na.omit(metric.dt[Scenario == u.scenario[i] & Metric == u.metric[j] & Model == u.model[m]])
						tmp.val = tmp.dt$Value
						tmp.25 = quantile(tmp.val,probs=c(0.25))
						tmp.50 = quantile(tmp.val,probs=c(0.5))
						tmp.75 = quantile(tmp.val,probs=c(0.75))
						tmp.iqr = c(tmp.25-1.5*(tmp.75-tmp.25),tmp.75+1.5*(tmp.75-tmp.25))
						if(min(tmp.val)>tmp.iqr[1]){tmp.iqr[1]=min(tmp.val)}
						if(max(tmp.val)<tmp.iqr[2]){tmp.iqr[2]=max(tmp.val)}
						tmp.out = tmp.val[which(tmp.val<tmp.iqr[1]|tmp.val>tmp.iqr[2])]
						if(m>4){m.idx=m;m=m+1} else {m.idx=m}
						points(rep(m,length(tmp.out)),tmp.out,pch=16,cex=0.75,col="gray70")
						lines(rep(m,2),tmp.iqr,col=line.vec[m.idx])
						polygon(x = c(m-x.offset,m+x.offset,m+x.offset,m-x.offset), y = c(tmp.25,tmp.25,tmp.75,tmp.75), density = NULL, angle = 45, border = line.vec[m.idx], col = fill.vec[m.idx])
						lines(c(m-x.offset,m+x.offset),rep(tmp.50,2),col=line.vec[m.idx])
						if(j==3){text(m,ylim.text,length(tmp.val),adj=c(0.5,0),cex=1.25)}
					}
				}
			}
			dev.off()	
