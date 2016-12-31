#' uprot
#'
#' Identifies and extracts a set of 40 prokaryotic universal proteins used as phylogenetic markers.
#' @param path is the full path to where genome annotations in fasta format are placed.
#' @param pattern a pattern (like '.faa') for recognizing the annotation files.
#' @param outdir a name for the output directory that will be created for results.
#' @param proc is the number of threads used for protein searches.
#' @keywords ribosomal proteins universal markers
#' @export
#' @examples
#' uprot(path='./prodigal_output',pattern='.faa',outdir='uprot_output',proc=2)

uprot<-function(pattern='.faa',
			    path='.',
			    outdir='uprot_output',
			    proc=2)
			   
			    {

	# Options #
	
	options(getClass.msg=FALSE)
	gw<-getwd()
	os<-.getOS()
	
	# Dependencies #
	
	require(seqinr,quietly=T)
	require(foreach,quietly=T)
	require(doMC,quietly=T)
	require(msa,quietly=T)
	require(plyr,quietly=T)

	# Internal functions #
	
	chunker<-function(m,n){

		s<-split(m,cut(seq_along(m),n,labels=F))
		return(s)
	}
				    	
	# Select OS #
	
	if (os=='linux'){
	
		hmmsearch<-paste(system.file('hmmer3',package='taxxo'),'/linux/hmmsearch',sep='')
		bldbcd<-paste(system.file('blast',package='taxxo'),'/linux/blastdbcmd',sep='')
		mkbldb<-paste(system.file('blast',package='taxxo'),'/linux/makeblastdb',sep='')
		
	} else if (os=='darwin'){
		
		hmmsearch<-paste(system.file('hmmer3',package='taxxo'),'/darwin/hmmsearch',sep='')
		bldbcd<-paste(system.file('blast',package='taxxo'),'/darwin/blastdbcmd',sep='')
		mkbldb<-paste(system.file('blast',package='taxxo'),'/darwin/makeblastdb',sep='')

	} else {
		
		stop('ERROR: Unknown OS, unable to proceed.')
	}

	# List HMMs and genome annotations
	
	phmm<-paste(system.file('exdata',package='taxxo'),'/uprot',sep='')
	thrs<-paste(system.file('exdata',package='taxxo'),'/uprot/tresholds.csv',sep='')
	
	flist<-gsub('//','/',list.files(path=path,pattern=pattern,full.names=T))
	fnams<-list.files(path=path,pattern=pattern)		
	hmms<-list.files(phmm,pattern='.hmm',full.names=T)
	hmm2<-list.files(phmm,pattern='.hmm')
	
	# Find proteins
	
	result<-matrix(ncol=length(hmms),nrow=length(flist),NA)
	colnames(result)<-gsub('.hmm','',hmm2)
	rownames(result)<-gsub(pattern,'',fnams)
	
	tresholds<-read.table(thrs,sep='\t',header=F,skip=1)
	
	for (f in 1:length(flist)){
		
			genome<-gsub(pattern,'',fnams[f])
					
		for (h in 1:length(hmms)){
			
			model<-gsub('.hmm','',hmm2[h])
			
			if (model%in%tresholds$V1){
			
				thr<-tresholds[which(tresholds[,1]==model),2]
			
				cmd<-paste(hmmsearch,
				           ' --noali -T ',thr,
				           ' --cpu ',proc,
				           ' -o search_tmp',
			    	       ' --domtblout out_tmp',
			        	   ' ',hmms[h],' ',flist[f],sep='')

			} else {
				
				cmd<-paste(hmmsearch,
				           ' --noali',
				           ' --cpu ',proc,
				           ' -o search_tmp',
			    	       ' --domtblout out_tmp',
			        	   ' ',hmms[h],' ',flist[f],sep='')
			}
			
			system(cmd)
			
			rl<-readLines('out_tmp')
			rl<-rl[which(!grepl("\\#",rl))]
			rl<-gsub('[ ]+',' ',rl)
			
			if (length(rl)>0){
			
				lst<-strsplit(rl,' ')
					
				hit<-sapply(lst,function(x){x[1]})
  				pfmID<-sapply(lst,function(x){x[2]})
  				query<-sapply(lst,function(x){x[4]})
  				evalu<-as.numeric(sapply(lst,function(x){x[13]}))
  				score<-as.numeric(sapply(lst,function(x){x[14]}))
  			
				htab<-data.frame(Query=query,
								 Hit=hit,
								 PfamID=pfmID,
								 Evalue=evalu,
								 Score=score,
        	                  	 stringsAsFactors=F)
                          	 
        	    dimi<-dim(htab)[1]
  			
  				if (dimi>1){
  				
  					maxi<-max(htab$Score)
  					gene<-as.vector(htab[which(htab$Score==maxi),2])
  				
  				} else if (dimi==1){
  				
  					gene<-as.vector(htab[1,2])
  				}
  				
  				result[f,h]<-gene
  				
  				system('rm -rf search_tmp out_tmp')
			}
		}
	}
  			
  	system(paste('mkdir',outdir))
  	setwd(outdir)
  			
  	presence_absence_uprot<-result
  			
  	save(presence_absence_uprot,file='presence_absence_uprot.Rdata')				
 		
 	# Make databases
 		 		
 	system('mkdir databases')		
 	setwd('databases')
 		
 	for (f in 1:length(flist)){
 		
 		dn<-gsub(pattern,'',fnams[f])
 			
 		cmd<-paste(mkbldb,
 				   ' -in ../.',flist[f],
 				   ' -dbtype prot -hash_index -parse_seqids -title ',
 				   dn,' -out ',dn,sep='')
 					   
		system(cmd,ignore.stderr=T,ignore.stdout=T)
 	}
 		
 	setwd('../')
 		 		
 	# Extract sequences
 		 		
 	uprots<-colnames(presence_absence_uprot)
 	genoms<-rownames(presence_absence_uprot)
 		
 	for (u in 1:length(uprots)){
 			
 		outfile<-paste(uprots[u],'.uprot.faa',sep='')
 			
 		prots<-as.vector(presence_absence_uprot[,u])
			
		for (p in 1:length(prots)){
				
			prot<-prots[p]
				
			if (is.na(prot)==F){
					
				cmd<-paste(bldbcd,' -entry ',prot,' -db ./databases/',genoms[p],' >> ',outfile,sep='')
					
				system(cmd)
					
			} else {
					
				cat(paste('>',genoms[p],'_absent',sep=''),
					file=outfile,
					sep='\n',
					append=T)
			}
		}
	}
 		
	# Align sequences
 		
 	multi<-list.files(pattern='.uprot.faa')
 		
 	for (m in multi){
 			
 		out<-gsub('.faa','.ali',m)
 		
 		alignment<-msa(inputSeqs=m,method='Muscle',type='protein')	
		aliconver<-msaConvert(alignment,type='seqinr::alignment')
			
		seqs<-lapply(aliconver$seq,s2c)
		nams<-gsub(' ','',gsub('>','',system(paste("grep '>' ",m,sep=''),intern=T)))
				
		write.fasta(seqs,names=nams,file=out)
	}
 		
 	# Concatenate alignment
 		
 	alis<-list.files(pattern='.uprot.ali')
 		
 	catmat<-NULL
 		
 	for (a in alis){
 			
		fasta<-read.fasta(a)
 		sequs<-lapply(getSequence(fasta),toupper)
 		smatx<-do.call(rbind,sequs)
 			
 		catmat<-cbind(catmat,smatx)
 		catlis<-alply(catmat,1)
 		
 		nams<-getName(fasta)
 		nams<-unlist(lapply(nams,function(x){strsplit(x,'_')[[1]][1]}))
 			
 		write.fasta(catlis,names=nams,file='universal_proteins.ali')
 	}
 	
 	setwd(gw)
}
