

# Nicholas Ducharme-Barth
# 03/02/2020
# Calculate the nominal CPUE for all replicates within each effort scenario and each sampling scenario

# set working directory
	project.dir = "C:/Users/nicholasd/HOME/SPC/SPC_SAM/Geostats/spatial.sampling.sim.manuscript/"
	setwd(project.dir)

# load packages
	library(ndd.vast.utils)
	data(skj.alt2019.shp)

# source functions
	source("R/Fn/fn.add.catchability.r")
	source("R/Fn/fn.nominal.simple.r")
	load("Background_Data/data.dt.RData")

# create storage structure
	nominal.df = expand.grid(Scenario = c("Contraction", "Expansion", "Fixed", "Preferential", "Random", "Rotating"),
						Catchability = c("noQ","Q"),
						Replicate = 1:100,
						Model = c("Nominal"),
						TS = 1:120,
						Region = c("all",1:8))

	nominal.df$Year = seq(from=1979,length.out=120,by=0.25)[nominal.df$TS]
	nominal.df$Index = NA
	nominal.df$SE = NA

# iterate
	for(s in c("Contraction", "Expansion", "Fixed", "Preferential", "Random", "Rotating"))
	{

		for(r in 1:100)
		{
			save.id = r
			if(r<100){save.id = paste0("0",save.id)}
			if(r<10){save.id = paste0("0",save.id)}
			# load replicate
				load(paste0("SimData/Simple120_wReplacement/",s,"/samp.dt.",save.id,".RData"))


				# no Catchability
				tmp.data = samp.dt
				colnames(tmp.data)[c(1,5:7)] = c("Year","Lon","Lat","Response_variable")
				tmp.noQ = nominal.simple(as.data.frame(tmp.data),n.yr.rng = 120,scale=TRUE, strata.sp=skj.alt2019.shp)
					
				# with Catchability
				samp.dt$True_Abundance = data.dt$skj.noise.patchy[samp.dt$id.data]
				samp.dt = as.data.frame(samp.dt[,c("True_Abundance","ts","lon","lat")])
				colnames(samp.dt) = c("True_Abundance","Year","Lon","Lat")
				Data_Geostat = add.catchability(samp.dt,seed = r,n.vessel.target = 90,new.entry.target=30,cv = 0.15,plot=FALSE)
				tmp.Q = nominal.simple(Data_Geostat,n.yr.rng = 120,scale=TRUE, strata.sp=skj.alt2019.shp)

			# assign back to data.frame
				pnt.noQ = which(nominal.df$Scenario == s & nominal.df$Catchability == "noQ" & nominal.df$Replicate == r)
				tmp.noQ = do.call(rbind,tmp.noQ)
				nominal.df$Index[pnt.noQ] = tmp.noQ$Index
				nominal.df$SE[pnt.noQ] = tmp.noQ$SE

				pnt.Q = which(nominal.df$Scenario == s & nominal.df$Catchability == "Q" & nominal.df$Replicate == r)
				tmp.Q = do.call(rbind,tmp.Q)
				nominal.df$Index[pnt.Q] = tmp.Q$Index
				nominal.df$SE[pnt.Q] = tmp.Q$SE

			# clean-up
				rm(list=c("save.id","samp.dt","tmp.data","tmp.noQ","Data_Geostat","tmp.Q","pnt.noQ","pnt.Q"))
		}
	}

# change order of columns TS & Region
	nominal.df = nominal.df[,c(1,2,3,4,6,5,7,8,9)]
	save(nominal.df,file="Index/ResultsDF/nominal.df.wR.RData")

