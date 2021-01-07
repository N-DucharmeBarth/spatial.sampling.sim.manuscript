

# Nicholas Ducharme-Barth
# 07/01/2021
# Launch jobs to condor

project.dir = "C:/Users/nicholasd/HOME/SPC/SPC_SAM/Geostats/spatial.sampling.sim.manuscript/"
setwd(project.dir)	
library(ssh)

session = ssh_connect("nicholasd@suvofpsubmit")
launch_machine.stem = "/home/nicholasd/spatial.sampling.sim.manuscript/Titration_Preferential/"
ssh_exec_wait(session, command = paste0("mkdir -p ",launch_machine.stem))

condor.files.dir = paste0(project.dir,"condor/condor.files/")
condor.loading.dock.dir = paste0(project.dir,"condor/condor_loading_dock/")
reps = 1:80 # which reps to use to calculate the indices

		# grab sim data and move to loading dock
			file.copy(from=paste0(project.dir,"SimData/Titration_Preferential/",list.files(paste0(project.dir,"SimData/Titration_Preferential/"))[reps]), to=condor.loading.dock.dir, overwrite = TRUE, recursive = FALSE, copy.mode=TRUE)

		# grab files and move to loading dock: condor_*.sub, samp.dt.names.txt, runR_*.bat, condor_run.*.r, data.dt.RData, fn.add.catchability.r, nino.df.RData, r361port.tar.gz, rm_except, simple.true.index.RData, sst.storage.df.RData, VAST_v8_3_0.*
			file.copy(from=paste0(project.dir,"R/Condor/",c("condor_titration.sub","runR_titration.bat","titration.samp.dt.names.txt")), to=condor.loading.dock.dir, overwrite = TRUE, recursive = FALSE, copy.mode=TRUE)
			file.copy(from=paste0(project.dir,"R/Scripts/",c("condor_run.titration.r")), to=condor.loading.dock.dir, overwrite = TRUE, recursive = FALSE, copy.mode=TRUE)
			file.copy(from=paste0(project.dir,"Background_Data/",c("data.dt.RData")), to=condor.loading.dock.dir, overwrite = TRUE, recursive = FALSE, copy.mode=TRUE)
			file.copy(from=paste0(project.dir,"R/Fn/",c("fn.add.catchability.r")), to=condor.loading.dock.dir, overwrite = TRUE, recursive = FALSE, copy.mode=TRUE)
			file.copy(from=paste0(project.dir,"SimData/",c("simple.true.index.RData")), to=condor.loading.dock.dir, overwrite = TRUE, recursive = FALSE, copy.mode=TRUE)
			file.copy(from=paste0(project.dir,"condor/condor_files/",list.files(paste0(project.dir,"condor/condor_files/"))), to=condor.loading.dock.dir, overwrite = TRUE, recursive = FALSE, copy.mode=TRUE)

		# dos2unix: runR_*.bat
			shell(paste0("cd ",condor.loading.dock.dir,"& C:/cygwin64/bin/dos2unix.exe condor_titration.sub"),translate=TRUE)
			shell(paste0("cd ",condor.loading.dock.dir,"& C:/cygwin64/bin/dos2unix.exe runR_titration.bat"),translate=TRUE)
			shell(paste0("cd ",condor.loading.dock.dir,"& C:/cygwin64/bin/dos2unix.exe titration.samp.dt.names.txt"),translate=TRUE)

		# create Start.tar.gz: data.dt.RData, fn.add.catchability.r, nino.df.RData, r361port.tar.gz, rm_except, simple.true.index.RData, sst.storage.df.RData, VAST_v8_3_0.* 
				TarList=c("condor_run.titration.r","data.dt.RData", "fn.add.catchability.r", "r361port.tar.gz", "rm_except", "simple.true.index.RData", "VAST_v8_3_0.*")
				shell(paste0("cd ",condor.loading.dock.dir,"& C:/cygwin64/bin/tar.exe -czvf Start.tar.gz ",paste(TarList,collapse=' ')),translate=TRUE)

		# create new directories for results: local & remote
			dir.create(paste0("Index/Titration_Preferential/"),recursive=TRUE,showWarnings=FALSE)

		# send files to remote: condor_*.sub, samp.dt.names.txt, samp.dt.*.RData, runR_*.bat, Start.tar.gz
       		scp_upload(session,files=paste0(condor.loading.dock.dir,c("condor_titration.sub","titration.samp.dt.names.txt","runR_titration.bat","Start.tar.gz",list.files(paste0(project.dir,"SimData/Titration_Preferential/"))[reps])),to=paste0(launch_machine.stem,"/"))


		# submit job
       			ssh_exec_wait(session,command=paste0('cd ',paste0(launch_machine.stem,"/"),'; chmod 777 runR_titration.bat'))
       			ssh_exec_wait(session,command=paste0('cd ',paste0(launch_machine.stem,"/"),'; condor_submit condor_titration.sub'))

		# remove files from loading dock
        	file.remove(paste0(condor.loading.dock.dir,list.files(condor.loading.dock.dir)))

		# clean-up workspace
        	rm(list=c("TarList"))

ssh_disconnect(session)


