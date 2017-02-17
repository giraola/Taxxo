#' prodigal
#'
#' Predicts protein coding genes from the genome using Prodigal.
#' @param path is the full path to where genomes in fasta format are placed.
#' @param pattern a pattern (like '.fasta') for recognizing the genome files.
#' @param outdir a name for the output directory that will be created for results.
#' @keywords gene prediction Prodigal
#' @export
#' @examples
#' rrna(path='.',pattern='.fasta',outdir='prodigal_output')


prodigal<-function(path='.',
				   pattern='.fasta',
				   outdir='prodigal_output',
				   proc=2
				   )
				   
				   {
	# Options #
	
	options(getClass.msg=FALSE)
	os<-.getOS()
	
	# Select OS #
	
	if (os=='linux'){
	
		prodigal<-paste(system.file('prodigal',package='taxxo'),'/linux/prodigal',sep='')
		
	} else if (os=='darwin'){
		
		prodigal<-paste(system.file('prodigal',package='taxxo'),'/darwin/prodigal',sep='')

	} else {
		
		stop('ERROR: Unknown OS, unable to proceed.')
	}
		
	# Dependencies #
	
	suppressPackageStartupMessages(library(seqinr,quietly=T))
	suppressPackageStartupMessages(library(foreach,quietly=T))
	suppressPackageStartupMessages(library(doMC,quietly=T))

	# Internal functions #
	
	chunker<-function(m,n){

		s<-split(m,cut(seq_along(m),n,labels=F))
		return(s)
	}
	
	system('touch prodigal.err')
	
	prodigal_parallel<-function(z){
		
		chun<-chunks[[z]]
		
		for (h in chun){
			
			ingenome<-fullpath[h]
			outgenome<-paste(fullpath[h],'.out',sep='')
			logfile<-paste(fullpath[h],'.log',sep='')
			faa<-paste(fullpath[h],'.faa',sep='')
			ffn<-paste(fullpath[h],'.ffn',sep='')
			
			cmd<-paste(prodigal,' -i ',ingenome,' -a ',faa,' -d ',ffn,' -o ',outgenome,' &> ',logfile,sep='')
		
			system(cmd)
		}
	}
	
	parse_annot<-function(y){
		
		chun<-chunks[[y]]
		
		for (h in chun){
			
			infaa<-paste(fullpath[h],'.faa',sep='')
			inffn<-paste(fullpath[h],'.ffn',sep='')
		
			if (file.info(infaa)$size>0){
			
				faa<-read.fasta(infaa)
				ffn<-read.fasta(inffn)
		
				sfa<-lapply(getSequence(faa),toupper)
				sfn<-lapply(getSequence(ffn),toupper)
			
				len<-length(sfa)
			
				nam<-paste(flist[h],seq(1,len),sep='_')
			
				write.fasta(sfa,names=nam,file=paste(fullpath[h],'.faa',sep=''))
				write.fasta(sfn,names=nam,file=paste(fullpath[h],'.ffn',sep=''))
			
			} else {
				
				mssg<-paste('Files',infaa,'and',inffn,'are empty!')
				warning(mssg)
				
				cat(mssg,sep='\n',append=T,file='prodigal.err')
			}
		}
	}
				
	# Start analysis #

	flist<-list.files(path=path,pattern=pattern)
	fullpath<-gsub('//','/',list.files(path=path,pattern=pattern,full.names=T))
		
	chunks<-chunker(1:length(flist),proc)
	
	# Run Prodigal in parallel #

	registerDoMC(proc)
	
	foreach(z=1:proc) %dopar% {
		
		prodigal_parallel(z)
	}
			
	registerDoMC(proc)
	
	foreach(z=1:proc) %dopar% {
		
		parse_annot(z)
	}
	
	system(paste('mkdir',outdir))	
	
	system(paste('mv ',path,'/*.faa ',outdir,sep=''))
	system(paste('mv ',path,'/*.ffn ',outdir,sep=''))
	system(paste('mv ',path,'/*.out ',outdir,sep=''))
	system(paste('mv ',path,'/*.log ',outdir,sep=''))
	system(paste('mv ',path,'/prodigal.err',outdir,sep=''))
}
