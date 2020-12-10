

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
	load("Index/ResultsDF/ts.df.MR.RData")
	load("Index/ResultsDF/nominal.df.RData")

#____________________________________________________________________________________________________________________________________________________________________________
# plot example time series 
	tmp.df = na.omit(ts.df)
	tmp.dt = data.table::as.data.table(tmp.df) %>% .[Region=="all"&mgc<=1e-04] %>% .[Model=="Enviro",Model:="E"] %>%
			 .[Model=="NoEnviro",Model:="noE"] %>% .[Model=="EnviroSVC",Model:="ESVC"] %>%
			 .[Model=="NoEnviroSVC",Model:="noESVC"] %>% .[Catchability=="Q",Model := paste0(Catchability,"-",Model)] %>% .[,Model:=as.factor(Model)] %>%
			 .[,Type := "Estimated"]

	# tmp.dt = tmp.dt[,.(Index = median(Index),Index.10 = quantile(Index,probs=0.1),Index.90 = quantile(Index,probs=0.9)),by=.(Scenario,Catchability,Model,Region,Year)]
	# tmp.dt$Model = "Estimated"
	true.df = data.frame(Year=seq(from=1979,length.out=120,by=0.25),True=simple.true.index[,1])
	true.df$Type = "True"

	nom.df = na.omit(subset(nominal.df,Region=="all"))
	nom.dt = data.table::as.data.table(nom.df) %>% .[,Type := "Nominal"]
	tmp.nom.1 = data.table::copy(nom.dt); tmp.nom.1 = tmp.nom.1 %>% .[,Model:="E"] %>% .[Catchability=="Q",Model := paste0(Catchability,"-","E")]
	tmp.nom.2 = data.table::copy(nom.dt); tmp.nom.2 = tmp.nom.2 %>% .[,Model:="noE"] %>% .[Catchability=="Q",Model := paste0(Catchability,"-","noE")]
	tmp.nom.3 = data.table::copy(nom.dt); tmp.nom.3 = tmp.nom.3 %>% .[,Model:="ESVC"] %>% .[Catchability=="Q",Model := paste0(Catchability,"-","ESVC")]
	tmp.nom.4 = data.table::copy(nom.dt); tmp.nom.4 = tmp.nom.4 %>% .[,Model:="noESVC"] %>% .[Catchability=="Q",Model := paste0(Catchability,"-","noESVC")]
	nom.dt = rbind(tmp.nom.1,tmp.nom.2,tmp.nom.3,tmp.nom.4)

	tmp.dt = rbind(tmp.dt[,.(Scenario,Catchability,Replicate,Model,Region,TS,Year,Index,SE,Type)],nom.dt[,.(Scenario,Catchability,Replicate,Model,Region,TS,Year,Index,SE,Type)])
	

	# extract replicates for given scenarios
	# replicates were randomly selected from set of replicates where all 4 models were fit
		list.dt = as.list(rep(NA,12))
		list.dt[[1]] = tmp.dt[Scenario == "Contraction" & Catchability == "noQ" & Replicate == 97]
		list.dt[[2]] = tmp.dt[Scenario == "Contraction" & Catchability == "Q" & Replicate == 56]
		list.dt[[3]] = tmp.dt[Scenario == "Expansion" & Catchability == "noQ" & Replicate == 8]
		list.dt[[4]] = tmp.dt[Scenario == "Expansion" & Catchability == "Q" & Replicate == 46]
		list.dt[[5]] = tmp.dt[Scenario == "Fixed" & Catchability == "noQ" & Replicate == 62]
		list.dt[[6]] = tmp.dt[Scenario == "Fixed" & Catchability == "Q" & Replicate == 54]
		list.dt[[7]] = tmp.dt[Scenario == "Preferential" & Catchability == "noQ" & Replicate == 32]
		list.dt[[8]] = tmp.dt[Scenario == "Preferential" & Catchability == "Q" & Replicate == 76]
		list.dt[[9]] = tmp.dt[Scenario == "Random" & Catchability == "noQ" & Replicate == 12]
		list.dt[[10]] = tmp.dt[Scenario == "Random" & Catchability == "Q" & Replicate == 15]
		list.dt[[11]] = tmp.dt[Scenario == "Rotating" & Catchability == "noQ" & Replicate == 100]
		list.dt[[12]] = tmp.dt[Scenario == "Rotating" & Catchability == "Q" & Replicate == 19]
		plot.dt = data.table::rbindlist(list.dt) %>% .[Type=="Estimated",l95 := Index-SE*1.96] %>% .[Type=="Estimated",u95 := Index+SE*1.96]

	p = plot.dt %>% ggplot2::ggplot() +  
	    ggplot2::facet_grid(Model~Scenario) +  ggplot2::ggtitle("Region: WCPO") +
	   	ggplot2::geom_line(data = plot.dt[Type=="Nominal"], ggplot2::aes(x=Year,y=Index),color="gray70") +
	    ggplot2::geom_ribbon(data = plot.dt[Type=="Estimated"],  ggplot2::aes(x=Year, ymin=l95,ymax=u95),color=NA,fill="#90caf9") +
	    ggplot2::geom_line(data = plot.dt[Type=="Estimated"], ggplot2::aes(x=Year,y=Index),color="#005cb2") +
		ggplot2::geom_line(data = true.df, ggplot2::aes(x=Year,y=True)) +
	    ggthemes::theme_few() +  ggplot2::scale_fill_manual(values=scales::alpha(c("#03a9f4","gray70","black"),0.5),name = "Index Type") +  ggplot2::scale_color_manual(values = c("#0277bd","gray30","black"),name = "Index Type")
	ggsave(filename=paste0("sim.ts.example.png"), plot = p, device = "png", path = "Plots/",
  			scale = 1, width = 16, height = 9, units = c("in"),
  			dpi = 300, limitsize = TRUE)

#____________________________________________________________________________________________________________________________________________________________________________
# plot metric boxplots
	load("Index/ResultsDF/metric.MR.df.RData")

	metric.dt = data.table::as.data.table(metric.df) %>% .[Region=="all"&mgc<=1e-04] %>% .[Model=="Enviro",Model:="E"] %>%
			 .[Model=="NoEnviro",Model:="noE"] %>% .[Model=="EnviroSVC",Model:="ESVC"] %>%
			 .[Model=="NoEnviroSVC",Model:="noESVC"] %>% .[Catchability=="Q",Model := paste0(Catchability,"-",Model)] %>% .[,Model:=factor(Model,levels=c("noE","noESVC","E","ESVC","Q-noE","Q-noESVC","Q-E","Q-ESVC"))] %>%
			 .[Metric %in% c("bias","rmsd","cover.50")] %>% .[,Metric:=factor(Metric,levels=c("rmsd","bias","cover.50"),labels=c("RMSD","Bias","Coverage"))]

			u.scenario = sort(unique(metric.dt$Scenario))
			u.metric = sort(unique(metric.dt$Metric))
			u.model = sort(unique(metric.dt$Model))
			# ylim.list = list(RMSD=range(pretty(metric.dt[Metric=="RMSD"]$Value,na.rm=TRUE)),Bias=range(pretty(metric.dt[Metric=="Bias"]$Value,na.rm=TRUE)),Coverage=range(pretty(metric.dt[Metric=="Coverage"]$Value,na.rm=TRUE)))
			yaxis.list = list(RMSD=c(0,0.05,0.1,0.15,0.2),Bias=c(-0.5,0,0.5,1,1.5),Coverage=c(0,25,50,75,100))
			x.offset = 0.2
			fill.vec = c("#ffebee","#ffcdd2","#ef9a9a","#e57373","#e3f2fd","#bbdefb","#90caf9","#64b5f6") #100,400,300,600
			line.vec = c(rep("#7f0000",4),rep("#002f6c",4))

			png(filename = "Plots/sim.metrics.png", width = 16, height = 9, units = "in", res=300,pointsize = 12, bg = "white")
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
						abline(h=50,lty=3)
					} else {
						plot(1,1,type="n",xlab="",xlim=xlim,ylab="",ylim=ylim,xaxt="n",yaxt="n",xaxs="i",yaxs="i")
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
