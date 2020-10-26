

# define quarterly nino4 index 
# https://esrl.noaa.gov/psd/enso/dashboard.html

setwd("C:/Users/nicholasd/HOME/SPC/SPC_SAM/Geostats/spatial.sampling.sim.manuscript/")	
nino.mat = as.matrix(read.csv("Background_Data/nino4monthly.csv",sep=";"))
nino.df = data.frame(yrqtr = rep(NA,nrow(nino.mat)*4),yr = rep(NA,nrow(nino.mat)*4),qtr = rep(NA,nrow(nino.mat)*4),index = rep(NA,nrow(nino.mat)*4))

iter = 1
index.list = list(c(2,3,4),c(5,6,7),c(8,9,10),c(11,12,13))
for(i in 1:nrow(nino.mat))
{
	for(j in 1:4)
	{
		nino.df$yr[iter] = nino.mat[i,1]
		nino.df$qtr[iter] = j
		nino.df$index[iter] = mean(nino.mat[i,index.list[[j]]])

		iter=iter+1
	}
}

nino.df$yrqtr = nino.df$yr + c(0,0.25,0.5,0.75)[nino.df$qtr]
nino.df$index.scale = scale(nino.df$index)

save(nino.df,file="Background_Data/nino.df.RData")