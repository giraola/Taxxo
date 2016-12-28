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
	
	require(seqinr,quietly=T)
	require(foreach,quietly=T)
	require(doMC,quietly=T)

	# Internal functions #
	
	chunker<-function(m,n){

		s<-split(m,cut(seq_along(m),n,labels=F))
		return(s)
	}
	
	prodigal_parallel<-function(z){
		
		chun<-chunks[[z]]
		
		for (h in chun){
			
			ingenome<-fullpath[h]
			outgenome<-gsub(pattern,'.out',fullpath[h])
			logfile<-gsub(pattern,'.log',fullpath[h])
			faa<-gsub(pattern,'.faa',fullpath[h])
			ffn<-gsub(pattern,'.ffn',fullpath[h])
			
			cmd<-paste(prodigal,' -i ',ingenome,' -a ',faa,' -d ',ffn,' -o ',outgenome,' &> ',logfile,sep='')
		
			system(cmd)
		}
	}
	
	parse_annot<-function(y){
		
		chun<-chunks[[y]]
		
		for (h in chun){
			
			infaa<-gsub(pattern,'.faa',fullpath[h])
			inffn<-gsub(pattern,'.ffn',fullpath[h])
		
			faa<-read.fasta(infaa)
			ffn<-read.fasta(inffn)
		
			sfa<-lapply(getSequence(faa),toupper)
			sfn<-lapply(getSequence(ffn),toupper)
			
			len<-length(sfa)
			
			nam<-paste(gsub(pattern,'',flist[h]),seq(1,len),sep='_')
			
			write.fasta(sfa,names=nam,file=gsub(pattern,'.faa',fullpath[h]))
			write.fasta(sfn,names=nam,file=gsub(pattern,'.ffn',fullpath[h]))
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
	
	system(paste('mv ',path,'*.faa ',outdir,sep=''))
	system(paste('mv ',path,'*.ffn ',outdir,sep=''))
	system(paste('mv ',path,'*.out ',outdir,sep=''))
	system(paste('mv ',path,'*.log ',outdir,sep=''))	
}