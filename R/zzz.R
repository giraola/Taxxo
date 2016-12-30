.onLoad<-function(){
		
	# Dependencies
	
	print('[1] Checking for dependencies...')
	
	packages.needed<-c('msa','phangorn','seqinr','doMC','foreach','plyr')
	packages.installed<-installed.packages()[,'Package']
	
	dif<-setdiff(packages.needed,packages.installed)
	
	if (length(dif)==0){
		
		print('   All dependencies available')
		
	} else {
		
		le<-length(dif)
		
		print(paste('   Installing',le,'packages'))
		
		for (d in dif){
			
			install.packages(d)
		}
	}
	
	# External software

	print('[2] Setting up external dependencies...')
	
	os<-.getOS()

	if (os=='linux'){

		hmmer_linux<-paste(system.file('hmmer3',package='taxxo'),'/linux/',sep='')
		blast_linux<-paste(system.file('blast',package='taxxo'),'/linux/',sep='')
		prodigal_linux<-paste(system.file('prodigal',package='taxxo'),'/linux/',sep='')
		barrnap_all<-paste(system.file('barrnap',package='taxxo'),'/',sep='')
		
		list_hmmer_linux<-list.files(hmmer_linux,pattern='.zip',full.names=T)
	
		if (length(list_hmmer_linux)>0){
		
			lhl<-gsub('//','/',list_hmmer_linux)
			system(paste('unzip ',lhl,' -d ',hmmer_linux,sep=''))
			system(paste('rm -rf ',hmmer_linux,'*.zip',sep=''))
		}
		
		list_blast_linux<-list.files(blast_linux,pattern='.zip',full.names=T)
	
		if (length(list_blast_linux)>0){
		
			lbl<-gsub('//','/',list_blast_linux)
			system(paste('unzip ',lbl,' -d ',blast_linux,sep=''))
			system(paste('rm -rf ',blast_linux,'*.zip',sep=''))
		}
		
		list_barrnap_linux<-list.files(barrnap_all,pattern='.zip',full.names=T)
		
		if (length(list_barrnap_linux)>0){
			
			lpl<-gsub('//','/',list_barranp_linux)
			system(paste('unzip ',lpl,' -d ',barrnap_all,sep=''))
			system(paste('rm -rf ',barrnap_all,'*.zip',sep=''))
		}
		
		system(paste('chmod +x ',hmmer_linux,'*',sep=''))
		system(paste('chmod +x ',blast_linux,'*',sep=''))
		system(paste('chmod +x ',prodigal_linux,'prodigal'))
		system(paste('chmod +x ',barrnap_all,'common/bin/barrnap',sep=''))
		system(paste('chmod +x ',barrnap_all,'common/binaries/linux/nhmmer',sep=''))
		
	} else if (os=='darwin'){
	
		hmmer_darwin<-paste(system.file('hmmer3',package='taxxo'),'/darwin/',sep='')
		blast_darwin<-paste(system.file('blast',package='taxxo'),'/darwin/',sep='')
		prodigal_darwin<-paste(system.file('prodigal',package='taxxo'),'/darwin/',sep='')
		barrnap_all<-paste(system.file('barrnap',package='taxxo'),'/',sep='')
		
		list_hmmer_darwin<-list.files(hmmer_darwin,pattern='.zip',full.names=T)
	
		if (length(list_hmmer_darwin)>0){
		
			lhd<-gsub('//','/',list_hmmer_darwin)
			system(paste('unzip ',lhd,' -d ',hmmer_darwin,sep=''))
			system(paste('rm -rf ',hmmer_darwin,'*.zip',sep=''))
		}
		
		list_blast_darwin<-list.files(blast_darwin,pattern='.zip',full.names=T)
	
		if (length(list_blast_darwin)>0){
		
			lbd<-gsub('//','/',list_blast_darwin)
			system(paste('unzip ',lbd,' -d ',blast_darwin,sep=''))
			system(paste('rm -rf ',blast_darwin,'*.zip',sep=''))
		}
		
		list_barrnap_darwin<-list.files(barrnap_all,pattern='.zip',full.names=T)
		
		if (length(list_barrnap_darwin)>0){
			
			lpd<-gsub('//','/',list_barranp_darwin)
			system(paste('unzip ',lpd,' -d ',barrnap_all,sep=''))
			system(paste('rm -rf ',barrnap_all,'*.zip',sep=''))
		}
		
		system(paste('chmod +x ',hmmer_darwin,'*',sep=''))
		system(paste('chmod +x ',blast_darwin,'*',sep=''))
		system(paste('chmod +x ',prodigal_darwin,'prodigal'))
		system(paste('chmod +x ',barrnap_all,'common/bin/barrnap',sep=''))
		system(paste('chmod +x ',barrnap_all,'common/binaries/darwin/nhmmer',sep=''))
	}
	
	print('   External software checked')
}	