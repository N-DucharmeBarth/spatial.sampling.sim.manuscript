

# Nicholas Ducharme-Barth
# 29/01/2020
# Launch jobs to condor

# 1) iterate over Q and noQ
# 2) iterate over Contraction, Expansion, Fixed, Preferential, Random & Rotating
# 3) Estimate the indices for the preceeding simulations (10 reps each, though this can be modified)

project.dir = "C:/Users/nicholasd/HOME/SPC/SPC_SAM/Geostats/spatial.sampling.sim.manuscript/"
setwd(project.dir)	
library(ssh)

session = ssh_connect("nicholasd@suvofpsubmit")
launch_machine.stem = "/home/nicholasd/spatial.sampling.sim.manuscript/Index/Simple120_wReplacement/"
ssh_exec_wait(session, command = paste0("mkdir -p ",launch_machine.stem))

condor.files.dir = paste0(project.dir,"condor/condor.files/")
condor.loading.dock.dir = paste0(project.dir,"condor/condor_loading_dock/")
reps = 1:100 # which reps to use to calculate the indices

# upload R portable one time
    scp_upload(session,files=paste0(project.dir,"condor/condor_files/r361port.tar.gz"),to="/home/nicholasd/")
	
# tar & upload all other files one time
	file.copy(from=paste0(project.dir,"Background_Data/data.dt.RData"),to=condor.loading.dock.dir, overwrite = TRUE, recursive = FALSE, copy.mode=TRUE)
	file.copy(from=paste0(project.dir,"Background_Data/sst.storage.df.RData"),to=condor.loading.dock.dir, overwrite = TRUE, recursive = FALSE, copy.mode=TRUE)
	file.copy(from=paste0(project.dir,"Background_Data/nino.df.RData"),to=condor.loading.dock.dir, overwrite = TRUE, recursive = FALSE, copy.mode=TRUE)
	file.copy(from=paste0(project.dir,"R/Fn/fn.add.catchability.r"),to=condor.loading.dock.dir, overwrite = TRUE, recursive = FALSE, copy.mode=TRUE)
	file.copy(from=paste0(project.dir,"SimData/simple.true.index.RData"),to=condor.loading.dock.dir, overwrite = TRUE, recursive = FALSE, copy.mode=TRUE)
	file.copy(from=paste0(project.dir,"condor/condor_files/rm_except"),to=condor.loading.dock.dir, overwrite = TRUE, recursive = FALSE, copy.mode=TRUE)
	file.copy(from=paste0(project.dir,"condor/condor_files/VAST_v8_3_0.cpp"),to=condor.loading.dock.dir, overwrite = TRUE, recursive = FALSE, copy.mode=TRUE)
	file.copy(from=paste0(project.dir,"condor/condor_files/VAST_v8_3_0.o"),to=condor.loading.dock.dir, overwrite = TRUE, recursive = FALSE, copy.mode=TRUE)
	file.copy(from=paste0(project.dir,"condor/condor_files/VAST_v8_3_0.so"),to=condor.loading.dock.dir, overwrite = TRUE, recursive = FALSE, copy.mode=TRUE)

	TarList=c("data.dt.RData", "fn.add.catchability.r", "nino.df.RData", "rm_except", "simple.true.index.RData", "sst.storage.df.RData", "VAST_v8_3_0.*" )
	shell(paste0("cd ",condor.loading.dock.dir,"& C:/cygwin64/bin/tar.exe -czvf Background.tar.gz ",paste(TarList,collapse=' ')),translate=TRUE)
    scp_upload(session,files=paste0(condor.loading.dock.dir,c("Background.tar.gz")),to="/home/nicholasd/")

# clean-up
   	file.remove(paste0(condor.loading.dock.dir,list.files(condor.loading.dock.dir)))
	rm(list=c("TarList"))


for(q in c("noQ","Q"))
{
	for(s in c("Contraction", "Expansion", "Fixed", "Preferential", "Random", "Rotating"))
	{
		# grab sim data and move to loading dock
			file.copy(from=paste0(project.dir,"SimData/Simple120_wReplacement/",s,"/",list.files(paste0(project.dir,"SimData/Simple120_wReplacement/",s,"/"))[reps]), to=condor.loading.dock.dir, overwrite = TRUE, recursive = FALSE, copy.mode=TRUE)

		# grab files and move to loading dock: condor_*.sub, samp.dt.names.txt, runR_*.bat, condor_run.*.r, data.dt.RData, fn.add.catchability.r, nino.df.RData, r361port.tar.gz, rm_except, simple.true.index.RData, sst.storage.df.RData, VAST_v8_3_0.*
			if(q == "Q")
			{
				file.copy(from=paste0(project.dir,"R/Condor/",c("condor_Q.sub","runR_Q.bat","samp.dt.names.txt")), to=condor.loading.dock.dir, overwrite = TRUE, recursive = FALSE, copy.mode=TRUE)
				file.copy(from=paste0(project.dir,"R/Scripts/",c("condor_run.Q.r")), to=condor.loading.dock.dir, overwrite = TRUE, recursive = FALSE, copy.mode=TRUE)
			} else {
				file.copy(from=paste0(project.dir,"R/Condor/",c("condor_noQ.sub","runR_noQ.bat","samp.dt.names.txt")), to=condor.loading.dock.dir, overwrite = TRUE, recursive = FALSE, copy.mode=TRUE)
				file.copy(from=paste0(project.dir,"R/Scripts/",c("condor_run.noQ.r")), to=condor.loading.dock.dir, overwrite = TRUE, recursive = FALSE, copy.mode=TRUE)
			}

		# dos2unix: runR_*.bat
			if(q == "Q")
			{
				shell(paste0("cd ",condor.loading.dock.dir,"& C:/cygwin64/bin/dos2unix.exe condor_Q.sub"),translate=TRUE)
				shell(paste0("cd ",condor.loading.dock.dir,"& C:/cygwin64/bin/dos2unix.exe runR_Q.bat"),translate=TRUE)

			} else {
				shell(paste0("cd ",condor.loading.dock.dir,"& C:/cygwin64/bin/dos2unix.exe condor_noQ.sub"),translate=TRUE)
				shell(paste0("cd ",condor.loading.dock.dir,"& C:/cygwin64/bin/dos2unix.exe runR_noQ.bat"),translate=TRUE)
			}
			shell(paste0("cd ",condor.loading.dock.dir,"& C:/cygwin64/bin/dos2unix.exe samp.dt.names.txt"),translate=TRUE)

		# create Start.tar.gz: data.dt.RData, fn.add.catchability.r, nino.df.RData, r361port.tar.gz, rm_except, simple.true.index.RData, sst.storage.df.RData, VAST_v8_3_0.* 
			if(q == "Q")
			{
				TarList=c("condor_run.Q.r")
				shell(paste0("cd ",condor.loading.dock.dir,"& C:/cygwin64/bin/tar.exe -czvf Start.tar.gz ",paste(TarList,collapse=' ')),translate=TRUE)
			} else {
				TarList=c("condor_run.noQ.r")
				shell(paste0("cd ",condor.loading.dock.dir,"& C:/cygwin64/bin/tar.exe -czvf Start.tar.gz ",paste(TarList,collapse=' ')),translate=TRUE)
			}

		# create new directories for results: local & remote
			dir.create(paste0("Index/Simple120_wReplacement/","/",s,"_",q),recursive=TRUE,showWarnings=FALSE)
			ssh_exec_wait(session, command = paste0("mkdir ",launch_machine.stem,"/",s,"_",q))

		# send files to remote: condor_*.sub, samp.dt.names.txt, samp.dt.*.RData, runR_*.bat, Start.tar.gz
			if(q == "Q")
			{
       			scp_upload(session,files=paste0(condor.loading.dock.dir,c("condor_Q.sub","samp.dt.names.txt","runR_Q.bat","Start.tar.gz",list.files(paste0(project.dir,"SimData/Simple120_wReplacement/",s,"/"))[reps])),to=paste0(launch_machine.stem,"/",s,"_",q))
			} else {
       			scp_upload(session,files=paste0(condor.loading.dock.dir,c("condor_noQ.sub","samp.dt.names.txt","runR_noQ.bat","Start.tar.gz",list.files(paste0(project.dir,"SimData/Simple120_wReplacement/",s,"/"))[reps])),to=paste0(launch_machine.stem,"/",s,"_",q))
			}

		# submit job
       		if(q == "Q")
       		{
       			ssh_exec_wait(session,command=paste0('cd ',paste0(launch_machine.stem,"/",s,"_",q),'; chmod 777 runR_Q.bat'))
       			ssh_exec_wait(session,command=paste0('cd ',paste0(launch_machine.stem,"/",s,"_",q),'; condor_submit condor_Q.sub'))
       		} else {
       			ssh_exec_wait(session,command=paste0('cd ',paste0(launch_machine.stem,"/",s,"_",q),'; chmod 777 runR_noQ.bat'))
       			ssh_exec_wait(session,command=paste0('cd ',paste0(launch_machine.stem,"/",s,"_",q),'; condor_submit condor_noQ.sub'))
       		}

		# remove files from loading dock
        	file.remove(paste0(condor.loading.dock.dir,list.files(condor.loading.dock.dir)))

		# clean-up workspace
        	rm(list=c("TarList"))
	}
}

ssh_disconnect(session)


