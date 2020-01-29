

# Nicholas Ducharme-Barth
# 30/01/2020
# Download jobs from condor

# 1) iterate over Q and noQ
# 2) iterate over Contraction, Expansion, Fixed, Preferential, Random & Rotating
# 3) Download the jobs

project.dir = "C:/Users/nicholasd/HOME/SPC/SPC_SAM/Geostats/spatial.sampling.sim.manuscript/"
setwd(project.dir)	
library(ssh)

session = ssh_connect("nicholasd@noumultifancl02")
launch_machine.stem = "/home/nicholasd/spatial.sampling.sim.manuscript/Index/Simple120/"

condor.files.dir = paste0(project.dir,"condor/condor.files/")
condor.loading.dock.dir = paste0(project.dir,"condor/condor_loading_dock/")

for(q in c("noQ","Q"))
{
	for(s in c("Contraction", "Expansion", "Fixed", "Preferential", "Random", "Rotating"))
	{
		# build sub directory for log files & indices
			dir.create(paste0("Index/Simple120/","/",s,"_",q,"/log"),recursive=TRUE,showWarnings=FALSE)
			dir.create(paste0("Index/Simple120/","/",s,"_",q,"/vast"),recursive=TRUE,showWarnings=FALSE)

		# identify files on condor
			tmp = ssh_exec_internal(session,command=paste0('cd ',launch_machine.stem,'/',s,'_',q,'; ls'))
			tmp.hex = tmp$stdout
			var.vec = rep(NA,length(which(tmp.hex == "0a"))) # 0a is the unix/linux line ending and indicates the start of a new name
			end.read = which(tmp.hex == "0a") - 1
			start.read = c(1,(which(tmp.hex == "0a") + 1)[1:(length(which(tmp.hex == "0a"))-1)])
			for(i in 1:length(which(tmp.hex == "0a")))
			{
				var.vec[i]=rawToChar(tmp.hex[start.read[i]:end.read[i]])
			}  

		# download *.vast_list.RData
			download.var = grep("vast_list",var.vec,value=TRUE)
			for(i in 1:length(download.var)){scp_download(session,files=paste0(paste0(launch_machine.stem,'/',s,'_',q),'/',download.var[i]),to= paste0("Index/Simple120/","/",s,"_",q,"/vast"))}

		# download *.log files
			download.var = grep("condor_R.log",var.vec,value=TRUE)
			for(i in 1:length(download.var)){scp_download(session,files=paste0(paste0(launch_machine.stem,'/',s,'_',q),'/',download.var[i]),to= paste0("Index/Simple120/","/",s,"_",q,"/log"))}

		# clean-up workspace
			rm(list=c("tmp","tmp.hex","var.vec","end.read","start.read","download.var"))
	}
}

ssh_disconnect(session)


