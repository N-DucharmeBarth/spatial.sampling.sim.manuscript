

#' Nicholas Ducharme-Barth
#' 31/01/2020
#' Define plotting functions
#' Plot metrics (bias,rmsd,mae,cover) by Model & Scenario. Data subsetted by Region & Catchability
#' Plot time series by Region & Catchability. Data subsetted by Model & Scenario
#' @param metric.df metric.df created by parse.output.r
#' @param ts.df ts.df created by parse.output.r
#' @param nominal.df bring in the nominal index
#' @param simple.true.index Bring in the true index
#' @param plot.type For time series plots: 1 = Model.Region & 2 = Model.Scenario
#' @param r Subset data to only this Region
#' @param q Subset data to only these Catchability scenarios
#' @param m Subset data to only this Model
#' @param s Subset data to only this Scenario
#' @param Save True or False
#' @param Save.Dir Path to directory where output is saved if Save == TRUE
#' @import ggplot2
#' @import ggthemes
#' @importFrom scales alpha
#' @import data.table
#' @export

# load packages
	library(ggplot2)
	library(ggthemes)

# define for testing
	# r = "all"
	# q = "Q"
	# Save = TRUE
	# Save.Dir = "Plots/Metric/"

# define functions

	# plot.bias(metric.df,r="all",q="noQ",Save=TRUE,Save.Dir="Plots/Metric/")
	# plot.rmse(metric.df,r="all",q="noQ",Save=TRUE,Save.Dir="Plots/Metric/")
	# plot.mae(metric.df,r="all",q="noQ",Save=TRUE,Save.Dir="Plots/Metric/")
	# plot.cover(metric.df,r="all",q="noQ",Save=TRUE,Save.Dir="Plots/Metric/")

	plot.bias.trend = function(metric.df,r,q,Save,Save.Dir)
	{
		# prep data
			tmp.df = na.omit(subset(metric.df,Metric=="bias"&Region==r&Catchability==q&mgc<=1e-04))

		# define helper function for calculating sample size
			n_bias <- function(y, upper_limit = 1.5) {
			  return( 
			    data.frame(
			      y = 0.9 * upper_limit,
			      label = paste(length(y), '\n')
			    )
			  )
			}

		# define plot
			p = ggplot(data = tmp.df, aes(x=Scenario, y=Value,fill=Model)) + geom_hline(yintercept = c(0,0.5,1.5),size=0.5,linetype="longdash",colour="gray70") +
     		geom_boxplot(outlier.color="gray60") + ylim(0,1.5) + stat_summary(fun.data = n_bias,geom="text",position = position_dodge(width = 0.75), aes(group=Model)) +
     		geom_hline(yintercept = 1,size=1,linetype="longdash",colour="black") + xlab(paste0("Region: ",r," Simulation: ",q)) + ylab("Bias")+
     		facet_wrap( ~ Scenario, scales="free_x") + theme_few() + scale_fill_manual(values=economist_pal()(4))

		if(Save)
		{
			ggsave(filename=paste0("bias.trend.",r,".",q,".png"), plot = p, device = "png", path = Save.Dir,
  			scale = 1.25, width = 6, height = 6, units = c("in"),
  			dpi = 300, limitsize = TRUE)
		} else {
			p
		}
	}

	plot.bias.coefficient = function(metric.df,r,q,Save,Save.Dir)
	{
		# prep data
			tmp.df = na.omit(subset(metric.df,Metric=="bias.coefficient"&Region==r&Catchability==q&mgc<=1e-04))

		# define helper function for calculating sample size
			n_bias <- function(y, upper_limit = 1) {
			  return( 
			    data.frame(
			      y = 0.9 * upper_limit,
			      label = paste(length(y), '\n')
			    )
			  )
			}

		# define plot
			p = ggplot(data = tmp.df, aes(x=Scenario, y=Value,fill=Model)) + geom_hline(yintercept = c(-0.5,0.5),size=0.5,linetype="longdash",colour="gray70") +
     		geom_boxplot(outlier.color="gray60") + ylim(-1,1) + stat_summary(fun.data = n_bias,geom="text",position = position_dodge(width = 0.75), aes(group=Model)) +
     		geom_hline(yintercept = 0,size=1,linetype="longdash",colour="black") + xlab(paste0("Region: ",r," Simulation: ",q)) + ylab(expression(paste("Bias coefficient (",kappa,")")))+
     		facet_wrap( ~ Scenario, scales="free_x") + theme_few() + scale_fill_manual(values=economist_pal()(4))

		if(Save)
		{
			ggsave(filename=paste0("bias.coefficient.",r,".",q,".png"), plot = p, device = "png", path = Save.Dir,
  			scale = 1.25, width = 6, height = 6, units = c("in"),
  			dpi = 300, limitsize = TRUE)
		} else {
			p
		}
	}

	plot.rmse = function(metric.df,r,q,Save,Save.Dir)
	{
		# prep data
			tmp.df = na.omit(subset(metric.df,Metric=="rmsd"&Region==r&Catchability==q&mgc<=1e-04))

		# define helper function for calculating sample size
			n_bias <- function(y, upper_limit = 0.3) {
			  return( 
			    data.frame(
			      y = 0.9 * upper_limit,
			      label = paste(length(y), '\n')
			    )
			  )
			}

		# define plot
			p = ggplot(data = tmp.df, aes(x=Scenario, y=Value,fill=Model)) + geom_hline(yintercept = c(0.1,0.2,0.3),size=0.5,linetype="longdash",colour="gray70") +
     		geom_boxplot(outlier.color="gray60") + ylim(0,0.3) + stat_summary(fun.data = n_bias,geom="text",position = position_dodge(width = 0.75), aes(group=Model)) +
     		geom_hline(yintercept = 0,size=1,linetype="solid",colour="black") + xlab(paste0("Region: ",r," Simulation: ",q)) + ylab("RMSE")+
     		facet_wrap( ~ Scenario, scales="free_x") + theme_few() + scale_fill_manual(values=economist_pal()(4))

		if(Save)
		{
			ggsave(filename=paste0("rmse.",r,".",q,".png"), plot = p, device = "png", path = Save.Dir,
  			scale = 1.25, width = 6, height = 6, units = c("in"),
  			dpi = 300, limitsize = TRUE)
		} else {
			p
		}
	}

	plot.mae = function(metric.df,r,q,Save,Save.Dir)
	{
		# prep data
			tmp.df = na.omit(subset(metric.df,Metric=="mae"&Region==r&Catchability==q&mgc<=1e-04))

		# define helper function for calculating sample size
			n_bias <- function(y, upper_limit = 0.3) {
			  return( 
			    data.frame(
			      y = 0.9 * upper_limit,
			      label = paste(length(y), '\n')
			    )
			  )
			}

		# define plot
			p = ggplot(data = tmp.df, aes(x=Scenario, y=Value,fill=Model)) + geom_hline(yintercept = c(0.1,0.2,0.3),size=0.5,linetype="longdash",colour="gray70") +
     		geom_boxplot(outlier.color="gray60") + ylim(0,0.3) + stat_summary(fun.data = n_bias,geom="text",position = position_dodge(width = 0.75), aes(group=Model)) +
     		geom_hline(yintercept = 0,size=1,linetype="solid",colour="black") + xlab(paste0("Region: ",r," Simulation: ",q)) + ylab("MAE")+
     		facet_wrap( ~ Scenario, scales="free_x") + theme_few() + scale_fill_manual(values=economist_pal()(4))

		if(Save)
		{
			ggsave(filename=paste0("mae.",r,".",q,".png"), plot = p, device = "png", path = Save.Dir,
  			scale = 1.25, width = 6, height = 6, units = c("in"),
  			dpi = 300, limitsize = TRUE)
		} else {
			p
		}
	}

	plot.cover = function(metric.df,r,q,Save,Save.Dir)
	{
		# prep data
			tmp.df = na.omit(subset(metric.df,Metric=="cover.50"&Region==r&Catchability==q&mgc<=1e-04))

		# define helper function for calculating sample size
			n_bias <- function(y, upper_limit = 100) {
			  return( 
			    data.frame(
			      y = 0.05 * upper_limit,
			      label = paste(length(y), '\n')
			    )
			  )
			}

		# define plot
			p = ggplot(data = tmp.df, aes(x=Scenario, y=Value,fill=Model)) + geom_hline(yintercept = c(0,25,75,100),size=0.5,linetype="longdash",colour="gray70") +
     		geom_boxplot(outlier.color="gray60") + ylim(0,100) + stat_summary(fun.data = n_bias,geom="text",position = position_dodge(width = 0.75), aes(group=Model)) +
     		geom_hline(yintercept = 50,size=1,linetype="longdash",colour="black") + xlab(paste0("Region: ",r," Simulation: ",q)) + ylab("Cover")+
     		facet_wrap( ~ Scenario, scales="free_x") + theme_few() + scale_fill_manual(values=economist_pal()(4))

		if(Save)
		{
			ggsave(filename=paste0("cover.",r,".",q,".png"), plot = p, device = "png", path = Save.Dir,
  			scale = 1.25, width = 6, height = 6, units = c("in"),
  			dpi = 300, limitsize = TRUE)
		} else {
			p
		}
	}

	plot.ts = function(ts.df,nominal.df,simple.true.index,plot.type,r,s,q,m,Save,Save.Dir)
	{
		if(plot.type == 1)
		{
			tmp.df = na.omit(subset(ts.df,Region==r&Catchability==q&Model==m&mgc<=1e-04))
			tmp.dt = data.table::as.data.table(tmp.df)
			tmp.dt = tmp.dt[,.(Index = median(Index),Index.10 = quantile(Index,probs=0.1),Index.90 = quantile(Index,probs=0.9)),by=.(Scenario,Catchability,Model,Region,Year)]
			tmp.dt$Model = "Estimated"
			true.df = data.frame(Year=seq(from=1979,length.out=120,by=0.25),True=simple.true.index[,r])
			true.df$Model = "True"

			nom.df = na.omit(subset(nominal.df,Region==r&Catchability==q&Model=="Nominal"))
			nom.dt = data.table::as.data.table(nom.df)
			nom.dt = nom.dt[,.(Index = median(Index),Index.10 = quantile(Index,probs=0.1),Index.90 = quantile(Index,probs=0.9)),by=.(Scenario,Catchability,Model,Region,Year)]

			tmp.dt = rbind(tmp.dt,nom.dt)
			
			p = ggplot(data = tmp.dt, aes(x=Year, y=Index,color=Model,fill=Model)) + facet_wrap(~Scenario) + ggtitle(paste0("Catchability: ",q," - Region: ",r," - Model: ",m)) +
				geom_ribbon( aes(x=Year, ymin=Index.10,ymax=Index.90,color=NA)) + geom_line(aes(x=Year,y=Index)) + geom_point(data = true.df,aes(x=Year,y=True)) +
			    theme_few() + scale_fill_manual(values=scales::alpha(c("#03a9f4","gray70","black"),0.5),name = "Index Type") + scale_color_manual(values = c("#0277bd","gray30","black"),name = "Index Type")
			if(Save)
			{
				ggsave(filename=paste0("ts.",q,".",m,".",r,".png"), plot = p, device = "png", path = Save.Dir,
	  			scale = 1.25, width = 9, height = 6, units = c("in"),
	  			dpi = 300, limitsize = TRUE)
			} else {
				p
			}
		} else {
			tmp.df = na.omit(subset(ts.df,Scenario==s&Catchability==q&Model==m&mgc<=1e-04))
			tmp.dt = data.table::as.data.table(tmp.df)
			tmp.dt = tmp.dt[,.(Index = median(Index),Index.10 = quantile(Index,probs=0.1),Index.90 = quantile(Index,probs=0.9)),by=.(Scenario,Catchability,Model,Region,Year)]
			tmp.dt$Model = "Estimated"
			
			true.df = data.frame(Year=rep(seq(from=1979,length.out=120,by=0.25),9),True=do.call(c,lapply(1:9,function(x)simple.true.index[,x])),Region=do.call(c,lapply(c("all",1:8),function(x)rep(x,120))))
			true.df$Model = "True"

			nom.df = na.omit(subset(nominal.df,Scenario==s&Catchability==q&Model=="Nominal"))
			nom.dt = data.table::as.data.table(nom.df)
			nom.dt = nom.dt[,.(Index = median(Index),Index.10 = quantile(Index,probs=0.1),Index.90 = quantile(Index,probs=0.9)),by=.(Scenario,Catchability,Model,Region,Year)]

			tmp.dt = rbind(tmp.dt,nom.dt)
			
			p = ggplot(data = tmp.dt, aes(x=Year, y=Index,color=Model,fill=Model)) + facet_wrap(~Region) + ggtitle(paste0("Catchability: ",q," - Effort scenario: ",s," - Model: ",m)) +
				geom_ribbon( aes(x=Year, ymin=Index.10,ymax=Index.90,color=NA)) + geom_line(aes(x=Year,y=Index)) + geom_point(data = true.df,aes(x=Year,y=True)) +
			    theme_few() + scale_fill_manual(values=scales::alpha(c("#03a9f4","gray70","black"),0.5),name = "Index Type") + scale_color_manual(values = c("#0277bd","gray30","black"),name = "Index Type")
			if(Save)
			{
				ggsave(filename=paste0("ts.",q,".",m,".",s,".png"), plot = p, device = "png", path = Save.Dir,
	  			scale = 1.25, width = 9, height = 6, units = c("in"),
	  			dpi = 300, limitsize = TRUE)
			} else {
				p
			}
		}
	}

	plot.mpe = function(metric.df,r,q,Save,Save.Dir)
	{
		# prep data
			tmp.df = na.omit(subset(metric.df,Metric=="MPE"&Region==r&Catchability==q&mgc<=1e-04))

		# define helper function for calculating sample size
			n_bias <- function(y, upper_limit = 10) {
			  return( 
			    data.frame(
			      y = 0.75 * upper_limit,
			      label = paste(length(y), '\n')
			    )
			  )
			}

		# define plot
			p = ggplot(data = tmp.df, aes(x=Scenario, y=Value,fill=Model)) + geom_hline(yintercept = c(-10,-5,5,10),size=0.5,linetype="longdash",colour="gray70") +
     		geom_boxplot(outlier.color="gray60") + ylim(-10,10) + stat_summary(fun.data = n_bias,geom="text",position = position_dodge(width = 0.75), aes(group=Model)) +
     		geom_hline(yintercept = c(0),size=0.5,linetype="longdash",colour="black") + xlab(paste0("Region: ",r," Simulation: ",q)) + ylab("MPE (%)")+
     		facet_wrap( ~ Scenario, scales="free_x") + theme_few() + scale_fill_manual(values=economist_pal()(4))

		if(Save)
		{
			ggsave(filename=paste0("mpe.",r,".",q,".png"), plot = p, device = "png", path = Save.Dir,
  			scale = 1.25, width = 6, height = 6, units = c("in"),
  			dpi = 300, limitsize = TRUE)
		} else {
			p
		}
	}