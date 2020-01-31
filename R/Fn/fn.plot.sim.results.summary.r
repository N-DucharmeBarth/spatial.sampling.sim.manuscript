

# Nicholas Ducharme-Barth
# 31/01/2020
# Define plotting functions
# Plot metrics (bias,rmsd,mae,cover) by Model & Scenario. Data subsetted by Region & Catchability
# Plot time series by Region & Catchability. Data subsetted by Model & Scenario

library(ggplot2)
library(ggthemes)

n_bias <- function(y, upper_limit = 1.5) {
  return( 
    data.frame(
      y = 0.9 * upper_limit,
      label = paste(length(y), '\n')
    )
  )
}


tmp.df = na.omit(subset(metric.df,Metric=="bias"&Region=="all"&Catchability=="noQ"&mgc<=1e-04))
colnames(tmp.df)[which(colnames(tmp.df)=="Value")] = "Bias"
colnames(tmp.df)[which(colnames(tmp.df)=="Scenario")] = "Region.All"
tmp.df$Combo = as.factor(paste0(tmp.df$Scenario,".",tmp.df$Model))
p <- ggplot(data = tmp.df, aes(x=Region.All, y=Bias,fill=Model)) + geom_hline(yintercept = c(0,0.5,1.5),size=0.5,linetype="longdash",colour="gray70") +
     geom_boxplot() + ylim(0,1.5) + stat_summary(fun.data = n_bias,geom="text",position = position_dodge(width = 0.75), aes(group=Model)) +
     geom_hline(yintercept = 1,size=1,linetype="longdash",colour="black") +
     facet_wrap( ~ Region.All, scales="free_x") + theme_few() + scale_fill_manual(values=economist_pal()(4))
p 



test.mat = cbind(as.vector(X_gtp[,1,1]),as.vector(X_gtp[,1,2]),as.vector(X_gtp[,1,3]))
test.response = test.mat %*% vast_list$vast_output$Enviro$par[249:251]
plot(Extrapolation_List$Data_Extrap[,1],test.response)
