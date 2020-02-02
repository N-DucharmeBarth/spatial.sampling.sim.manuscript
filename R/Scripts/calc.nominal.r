

# Nicholas Ducharme-Barth
# 03/02/2020
# Calculate the nominal CPUE for all replicates within each effort scenario and each sampling scenario

# set working directory
	project.dir = "C:/Users/nicholasd/HOME/SPC/SPC_SAM/Geostats/spatial.sampling.sim.manuscript/"
	setwd(project.dir)

# source functions
	source("R/Fn/fn.add.catchability.r")
	source("R/Fn/fn.nominal.simple.r")

# create storage structure
	nominal.df = expand.grid(Scenario = c("Contraction", "Expansion", "Fixed", "Preferential", "Random", "Rotating"),
						Catchability = c("noQ","Q"),
						Replicate = 1:100,
						Model = c("Nominal"),
						Region = c("all",1:8),
						TS = 1:120)

	nominal.df$Year = seq(from=1979,length.out=120,by=0.25)[nominal.df$TS]
	nominal.df$Index = NA
	nominal.df$SE = NA

# iterate
	for(s in c("Contraction", "Expansion", "Fixed", "Preferential", "Random", "Rotating"))
	{
		for(q in c("noQ","Q"))
		{
			for(r in 1:100)
			{
				# load replicate
					load(paste0("Index/Simple120/",s,"_",q,"/vast/",r,".vast_list.RData"))
				if(q == "noQ")
				{

				} else {

				}
			}
		}
	}