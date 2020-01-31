

#' Nicholas Ducharme-Barth
#' 31/01/2020
#' Define plotting functions
#' Plot metrics (bias,rmsd,mae,cover) by Model & Scenario. Data subsetted by Region & Catchability
#' Plot time series by Region & Catchability. Data subsetted by Model & Scenario
#' @param metric.df metric.df created by parse.output.r
#' @param ts.df ts.df created by parse.output.r
#' @param r Subset data to only this Region
#' @param q Subset data to only these Catchability scenarios
#' @param m Subset data to only this Model
#' @param s Subset data to only this Scenario
#' @param Save True or False
#' @param Save.Dir Path to directory where output is saved if Save == TRUE
#' @import ggplot2
#' @import ggthemes
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

	plot.bias = function(metric.df,r,q,Save,Save.Dir)
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
			ggsave(filename=paste0("bias.",r,".",q,".png"), plot = p, device = "png", path = Save.Dir,
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
			tmp.df = na.omit(subset(metric.df,Metric=="cover"&Region==r&Catchability==q&mgc<=1e-04))

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
