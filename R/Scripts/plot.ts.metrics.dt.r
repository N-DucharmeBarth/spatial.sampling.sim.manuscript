

# load packages
	library(ggplot2)
	library(ggthemes)
	library(data.table)

# set project directory
	project.dir = "C:/Users/nicholasd/HOME/SPC/SPC_SAM/Geostats/spatial.sampling.sim.manuscript/"
	setwd(project.dir)
	Save.Dir = "Plots/"
# load metrics
	load("Index/ResultsDF/ts.metrics.dt.RData")

# plot - type 1
	tmp.dt = ts.metrics.dt[mgc<1e-04&Scenario %in% c("Contraction") & Catchability == "noQ" & Model == "NoEnviro", .(rmse.025=quantile(rmse,probs=c(0.025), na.rm = TRUE),rmse.25=quantile(rmse,probs=c(0.25), na.rm = TRUE),rmse.50=quantile(rmse,probs=0.5, na.rm = TRUE),rmse.75=quantile(rmse,probs=c(0.75), na.rm = TRUE),rmse.975=quantile(rmse,probs=c(0.975), na.rm = TRUE),sample.prop=median(sample.prop,na.rm=TRUE)),by=.(Scenario,Catchability,Model,Region,Year)]
	p = ggplot(data = tmp.dt, aes(x=Year, y=rmse.50)) + geom_hline(yintercept = 0,size=1,linetype="solid",colour="black") + geom_vline(xintercept = 1982.5,size=0.5,linetype="dashed",colour="gray70") + 
     	geom_ribbon(aes(ymin=rmse.25,ymax=rmse.75,x=Year),alpha=0.25,fill="blue") + geom_ribbon(aes(ymin=rmse.025,ymax=rmse.975,x=Year),alpha=0.25,fill="blue") +
     	geom_line(aes(y=rmse.50,x=Year),alpha=0.5,colour="blue",size=1.25) +
     	geom_smooth(aes(y=sample.prop/4), method="loess", colour="red",span=0.25,se=FALSE) +
		coord_cartesian(ylim=c(0, 0.25)) + 
     	scale_y_continuous(name="RMSE", sec.axis=sec_axis(~.*4, name="Proportion knots sampled")) +     	
     	facet_wrap( ~ Region) + theme_few() + theme(axis.title.y = element_text(colour = "blue"), axis.title.y.right = element_text(colour = "red"))
   
    ggsave(filename=paste0("Contraction_RMSE_Sampling_v1.png"), plot = p, device = "png", path = Save.Dir,
  			scale = 1.25, width = 9, height = 6, units = c("in"),
  			dpi = 300, limitsize = TRUE)

# plot - type 2
	tmp.dt = ts.metrics.dt[mgc<1e-04&Scenario %in% c("Contraction") & Catchability == "noQ" & Model == "NoEnviro"]
	p = ggplot(data = tmp.dt,aes(y=rmse,x=Year)) + geom_hline(yintercept = 0,size=1,colour="black") + geom_vline(xintercept = 1982.5,size=0.5,linetype="dashed",colour="gray70") +
		geom_point(size=0.15,alpha=0.1) + 
		geom_smooth( method="loess", col="blue",span=0.05,se=FALSE) +
		geom_smooth(aes(y=sample.prop/4), method="loess", col="red",span=0.25,se=FALSE) +
		coord_cartesian(ylim=c(0, 0.25)) + 
     	scale_y_continuous(name="RMSE", sec.axis=sec_axis(~.*4, name="Proportion knots sampled")) +
     	facet_wrap( ~ Region) + theme_few() + theme(axis.title.y = element_text(colour = "blue"), axis.title.y.right = element_text(colour = "red"))
     ggsave(filename=paste0("Contraction_RMSE_Sampling_v2.png"), plot = p, device = "png", path = Save.Dir,
  			scale = 1.25, width = 9, height = 6, units = c("in"),
  			dpi = 300, limitsize = TRUE)

