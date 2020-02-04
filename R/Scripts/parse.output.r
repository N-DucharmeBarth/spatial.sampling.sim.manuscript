

# Nicholas Ducharme-Barth
# 30/01/2020
# Parse diagnostics & metrics from condor runs
# store results in a data.frame
# diag.df: Scenario, Catchability, Replicate, Model, Index, Machine, Condor.time, Fit.time, mgc, error
# metric.df: Scenario, Catchability, Replicate, Model, Region, Metric, Value
# ts.df: Scenario, Catchability, Replicate, Model, Region, TS, Year, Index, SE

# 1) iterate over Q and noQ
# 2) iterate over Contraction, Expansion, Fixed, Preferential, Random & Rotating
# 3) Get diagnostics from vast_list & log

project.dir = "C:/Users/nicholasd/HOME/SPC/SPC_SAM/Geostats/spatial.sampling.sim.manuscript/"
setwd(project.dir)

# define reps to read
	reps = 61:80

# define storage structures
	diag.df = expand.grid(Scenario = c("Contraction", "Expansion", "Fixed", "Preferential", "Random", "Rotating"),
						  Catchability = c("noQ","Q"),
						  Replicate = reps,
						  Model = c("Enviro","NoEnviro","EnviroSVC","NoEnviroSVC"))
	diag.df$Index = NA
	diag.df$Machine = NA
	diag.df$condor.time = NA
	diag.df$fit.time = NA
	diag.df$mgc = NA
	diag.df$error.type = NA

	metric.df = expand.grid(Scenario = c("Contraction", "Expansion", "Fixed", "Preferential", "Random", "Rotating"),
						  Catchability = c("noQ","Q"),
						  Replicate = reps,
						  Model = c("Enviro","NoEnviro","EnviroSVC","NoEnviroSVC"),
						  Region = c("all",1:8),
						  Metric = c("bias","mae","rmsd","cover"))
	metric.df$Value = NA
	metric.df$mgc = NA

	ts.df = expand.grid(Scenario = c("Contraction", "Expansion", "Fixed", "Preferential", "Random", "Rotating"),
						  Catchability = c("noQ","Q"),
						  Replicate = reps,
						  Model = c("Enviro","NoEnviro","EnviroSVC","NoEnviroSVC"),
						  Region = c("all",1:8),
						  TS = 1:120)

	ts.df$Year = seq(from=1979,length.out=120,by=0.25)[ts.df$TS]
	ts.df$Index = NA
	ts.df$SE = NA
	ts.df$mgc = NA

# open connection
session = ssh::ssh_connect("nicholasd@noumultifancl02")

for(q in c("noQ","Q"))
{
	for(s in c("Contraction", "Expansion", "Fixed", "Preferential", "Random", "Rotating"))
	{
		for(r in reps)
		{
			pnt = which(diag.df$Scenario == s & diag.df$Catchability == q & diag.df$Replicate == r)
			# read vast_output & read from log
				try(load(paste0("Index/Simple120/",s,"_",q,"/vast/",r,".vast_list.RData")),silent=TRUE)
				log = try(readLines(paste0("Index/Simple120/",s,"_",q,"/log/",list.files(paste0("Index/Simple120/",s,"_",q,"/log"))[grep(paste0(".",r,".condor_R.log"),list.files(paste0("Index/Simple120/","/",s,"_",q,"/log")),fixed=TRUE)])),silent=TRUE)
				if("vast_list" %in% ls() == FALSE){vast_list = NULL}
			# parse log
				if(length(log)>1)
				{
					if(length(log[grep("Job terminated",log)])>0)
					{
						tmp.hex = ssh::ssh_exec_internal(session,command=paste0("nslookup ",strsplit(strsplit(log[grep("Job executing on host",log)],"<")[[1]][2],":")[[1]][1]))$stdout
						# 0a is the unix/linux line ending
						tmp.var = rawToChar(tmp.hex[c(1,(which(tmp.hex == "0a") + 1)[1:(length(which(tmp.hex == "0a"))-1)])[1]:((which(tmp.hex == "0a") - 1)[1])])
						Machine = trimws(strsplit(strsplit(tmp.var,"=")[[1]][2],"[.]")[[1]][1])
						condor.time = as.numeric(difftime(strptime(trimws(strsplit(strsplit(log[grep("Job terminated",log)],")")[[1]][2],"Job")[[1]][1]),format="%m/%d %H:%M:%S"),strptime(trimws(strsplit(strsplit(log[grep("Job submitted from host",log)],")")[[1]][2],"Job")[[1]][1]),format="%m/%d %H:%M:%S"),units="hours"))
						
						# add to df
						diag.df$Machine[pnt] = Machine
						diag.df$condor.time[pnt] = condor.time
						diag.df$error.type[pnt] = "none"

						# clean-up
						rm(list=c("tmp.hex","tmp.var","Machine","condor.time"))
					} else {
						diag.df$error.type[pnt] = "still_running"
					}
				} else {
					# add to df
					diag.df$error.type[pnt] = "missing_log"
				}
			# get model diagnostics and metrics
			if(length(vast_list)==3)
			{
				# iterate across models
				for(m in c("Enviro","NoEnviro","EnviroSVC","NoEnviroSVC"))
				{
					pntm.diag = which(diag.df$Scenario == s & diag.df$Catchability == q & diag.df$Replicate == r & diag.df$Model == m)
					pntm.met = which(metric.df$Scenario == s & metric.df$Catchability == q & metric.df$Replicate == r & metric.df$Model == m)

					# get diagnostics and metrics
					if(length(vast_list$vast_output[[m]])>1)
					{
						fit.time = unname(vast_list$vast_output[[m]]$fit.time/(60*60))
						if("opt" %in% names(vast_list$vast_output[[m]]))
						{
							mgc = vast_list$vast_output[[m]]$opt$max_gradient
						} else {
							mgc = vast_list$vast_output[[m]]$max_gradient
						}
						diag.df$fit.time[pntm.diag] = fit.time
						diag.df$mgc[pntm.diag] = mgc
						metric.df$mgc[pntm.met] = mgc
						if(length(as.vector(as.matrix(vast_list$vast_metric[[m]])))==36)
						{
							metric.df$Value[pntm.met] =  as.vector(as.matrix(vast_list$vast_metric[[m]]))
						}
						if(try(length(as.vector(as.matrix(vast_list$vast_output[[m]]$idx[,-1]))),silent=TRUE)==1080)
						{
							idx.raw = as.matrix(vast_list$vast_output[[m]]$idx[,-1])
							idx.se = try(as.matrix(vast_list$vast_output[[m]]$idx.se[,-1]),silent=TRUE)
							idx.raw.means = colMeans(idx.raw)
							for(a in 1:length(idx.raw.means))
							{
								pntm.ts = which(ts.df$Scenario == s & ts.df$Catchability == q & ts.df$Replicate == r & ts.df$Model == m & ts.df$Region == c("all",1:8)[a])
								ts.df$mgc[pntm.ts] =  mgc
								ts.df$Index[pntm.ts] =  idx.raw[,a]/idx.raw.means[a]
								if(length(idx.se)==1080)
								{
									ts.df$SE[pntm.ts] =  idx.se[,a]/idx.raw.means[a]
								}
							}
							# clean-up
							rm(list=c("idx.raw","idx.se","idx.raw.means","pntm.ts"))
						}
						rm(list=c("mgc","fit.time"))
					} else {
						# error type
						if(grepl("singular",vast_list$vast_output[[m]])){
							error = "singular"
						} else if (grepl("manifold",vast_list$vast_output[[m]])) {
							error = "manifold"
						} else if (grepl("compile",vast_list$vast_output[[m]])){
							error = "compile"
						} else {
							error = "other"
						}
						diag.df$error.type[pntm.diag] = error
						rm(list = c("error"))
					}

					# clean-up
					rm(list=c("pntm.diag","pntm.met"))
				}
			} else {
				diag.df$Index[pnt] = FALSE
				diag.df$error.type[pnt] = "missing_results"
			}

			# clean-up
				rm(list=c("pnt","log","vast_list"))
		}
	}
}

ssh::ssh_disconnect(session)

dir.create(paste0("Index/ResultsDF/"),recursive=TRUE,showWarnings=FALSE)

# create temporary df
	tmp.diag.df = diag.df
	tmp.metric.df = metric.df
	tmp.ts.df = ts.df

# load existing df
	load("Index/ResultsDF/diag.df.RData")
	load("Index/ResultsDF/metric.df.RData")
	load("Index/ResultsDF/ts.df.RData")

# append to existing
	diag.df = rbind(diag.df[-which(diag.df$Replicate %in% reps),],tmp.diag.df)
	metric.df = rbind(metric.df[-which(metric.df$Replicate %in% reps),],tmp.metric.df)
	ts.df = rbind(ts.df[-which(ts.df$Replicate %in% reps),],tmp.ts.df)

# save
	save(diag.df,file="Index/ResultsDF/diag.df.RData")
	save(metric.df,file="Index/ResultsDF/metric.df.RData")
	save(ts.df,file="Index/ResultsDF/ts.df.RData")

	