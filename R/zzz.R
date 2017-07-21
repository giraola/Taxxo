.onLoad<-function(...){
		
	# Dependencies
	
	packages.needed<-c('ape','seqinr','doMC','foreach','plyr')
	#packages.needed<-c('msa','ape','seqinr','doMC','foreach','plyr')
	packages.installed<-installed.packages()[,'Package']
	
	dif<-setdiff(packages.needed,packages.installed)
	
	if (length(dif)>0){
		
		le<-length(dif)
				
		for (d in dif){
			
			install.packages(d)
		}
	}
	
	# External software
	
	os<-.getOS()

	if (os=='linux'){

		clustalo_linux<-paste(system.file('clustalo',package='taxxo'),'/linux/',sep='')
		hmmer_linux<-paste(system.file('hmmer3',package='taxxo'),'/linux/',sep='')
		blast_linux<-paste(system.file('blast',package='taxxo'),'/linux/',sep='')
		prodigal_linux<-paste(system.file('prodigal',package='taxxo'),'/linux/',sep='')
		barrnap_all<-paste(system.file('barrnap',package='taxxo'),'/',sep='')
		fasttree_linux<-paste(system.file('fasttree',package='taxxo'),'/linux/',sep='')
		
		list_blast_linux<-list.files(blast_linux,pattern='.zip',full.names=T)
	
		if (length(list_blast_linux)>0){
		
			lbl<-gsub('//','/',list_blast_linux)
			
			for (l in lbl){
			
				system(paste('unzip ',l,' -d ',blast_linux,sep=''),ignore.stdout=T)
			}

			system(paste('rm -rf ',blast_linux,'*.zip',sep=''),ignore.stdout=T)
		}
		
		list_barrnap_linux<-list.files(barrnap_all,pattern='.zip',full.names=T)
		
		if (length(list_barrnap_linux)>0){
			
			lpl<-gsub('//','/',list_barrnap_linux)
			system(paste('unzip ',lpl,' -d ',barrnap_all,sep=''),ignore.stdout=T)
			system(paste('rm -rf ',lpl,sep=''),ignore.stdout=T)
		}
		
		system(paste('chmod +x ',clustalo_linux,'*',sep=''),ignore.stdout=T,ignore.stderr=T)
		system(paste('chmod +x ',hmmer_linux,'*',sep=''),ignore.stdout=T,ignore.stderr=T)
		system(paste('chmod +x ',blast_linux,'*',sep=''),ignore.stdout=T,ignore.stderr=T)
		system(paste('chmod +x ',prodigal_linux,'prodigal',sep=''),ignore.stdout=T,ignore.stderr=T)
		system(paste('chmod +x ',barrnap_all,'common/bin/barrnap',sep=''),ignore.stdout=T,ignore.stderr=T)
		system(paste('chmod +x ',barrnap_all,'common/binaries/linux/nhmmer',sep=''),ignore.stdout=T,ignore.stderr=T)
		system(paste('chmod +x ',fasttree_linux,'FastTree',sep=''),ignore.stdout=T,ignore.stderr=T)
		
	} else if (os=='darwin'){
	
		clustalo_darwin<-paste(system.file('clustalo',package='taxxo'),'/darwin',sep='')
		hmmer_darwin<-paste(system.file('hmmer3',package='taxxo'),'/darwin/',sep='')
		blast_darwin<-paste(system.file('blast',package='taxxo'),'/darwin/',sep='')
		prodigal_darwin<-paste(system.file('prodigal',package='taxxo'),'/darwin/',sep='')
		barrnap_all<-paste(system.file('barrnap',package='taxxo'),'/',sep='')
		fasttree_darwin<-paste(system.file('fasttree',package='taxxo'),'/darwin/',sep='')

		list_blast_darwin<-list.files(blast_darwin,pattern='.zip',full.names=T)
	
		if (length(list_blast_darwin)>0){
		
			lbd<-gsub('//','/',list_blast_darwin)
			
			for (l in lbd){
			
				system(paste('unzip ',l,' -d ',blast_darwin,sep=''),ignore.stdout=T)
			}
			
			system(paste('rm -rf ',blast_darwin,'*.zip',sep=''),ignore.stdout=T)
		}
		
		list_barrnap_darwin<-list.files(barrnap_all,pattern='.zip',full.names=T)
		
		if (length(list_barrnap_darwin)>0){
			
			lpd<-gsub('//','/',list_barrnap_darwin)
			system(paste('unzip ',lpd,' -d ',barrnap_all,sep=''),ignore.stdout=T)
			system(paste('rm -rf ',lpd,sep=''),ignore.stdout=T)
		}
		
		system(paste('chmod +x ',clustalo_darwin,'*',sep=''),ignore.stdout=T,ignore.stderr=T)
		system(paste('chmod +x ',hmmer_darwin,'*',sep=''),ignore.stdout=T,ignore.stderr=T)
		system(paste('chmod +x ',blast_darwin,'*',sep=''),ignore.stdout=T,ignore.stderr=T)
		system(paste('chmod +x ',prodigal_darwin,'prodigal',sep=''),ignore.stdout=T,ignore.stderr=T)
		system(paste('chmod +x ',barrnap_all,'common/bin/barrnap',sep=''),ignore.stdout=T,ignore.stderr=T)
		system(paste('chmod +x ',barrnap_all,'common/binaries/darwin/nhmmer',sep=''),ignore.stdout=T,ignore.stderr=T)
		system(paste('chmod +x ',fasttree_darwin,'FastTree',sep=''),ignore.stdout=T,ignore.stderr=T)

	}	
}
