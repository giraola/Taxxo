#' anib
#'
#' Calculates ANIb (Average Nucleotide Identity using Blast) values for all genomes present in a directory.
#' @param path is the full path to where genomes assemblies in fasta format are placed.
#' @param pattern a pattern (like '.fna') for recognizing the genome files.
#' @param outdir a name for the output directory that will be created for results.
#' @param proc is the number of threads for parallel computing.
#' @param reference is a genome name used to compare one vs. the rest. If set to 'all' (default), all vs. all is computed.
#' @keywords ANI species blastn
#' @export
#' @examples
#' anib(path='.',pattern='.fna',reference='all',outdir='anib_result',proc=2)

anib<-function(pattern='.fna',
			   reference='all',
			   path='.',
			   outdir='anib_output',
			   proc=2)
			   
			   {

	# Options #
	
	options(getClass.msg=FALSE)
	os<-.getOS()
		
	# Select OS #
	
	if (os=='linux'){
	
		blastn<-paste(system.file('blast',package='taxxo'),'/linux/blastn',sep='')
		
	} else if (os=='darwin'){
		
		blastn<-paste(system.file('blast',package='taxxo'),'/darwin/blastn',sep='')

	} else {
		
		stop('ERROR: Unknown OS, unable to proceed.')
	}
	
	# Dependencies #

	library(foreach,quietly=T)
	library(doMC,quietly=T)
	library(seqinr,quietly=T)
		
	# Define internal function #

	chunker<-function(m,n){

		s<-split(m,cut(seq_along(m),n,labels=F))
		return(s)
	}

	# List genomes #
	
	genomes<-gsub('//','/',list.files(path=path,pattern=pattern,full.names=T))
	g.names<-gsub(pattern,'',list.files(path=path,pattern=pattern))

	# ANIb parallel #
	
	aniparallel<-function(y){

		chun<-chunks[[y]]
		lchun<-length(chun)
		comb<-combina[,chun]
		
		if (lchun>1){
			
			lcomb<-dim(comb)[2]
			query1<-comb[1,i]
			subjt1<-comb[2,i]
			query2<-comb[2,i]
			subjt2<-comb[1,i]
			
		} else if (lchun==1){
			
			lcomb<-lchun
			query1<-comb[1]
			subjt1<-comb[2]
			query2<-comb[2]
			subjt2<-comb[1]
		}

		result<-NULL

		for (i in 1:lcomb){

			# BLAST 1

			#query1<-comb[1,i]
			#subjt1<-comb[2,i]
			
			blout1<-paste('blout1',y,i,sep='.')
			query<-paste('query1',y,i,'fasta',sep='.')

			qfasta<-read.fasta(query1)
			qseq<-unlist(lapply(getSequence(qfasta),toupper))
				
			rnd<-round(length(qseq)/1020)
			fragments<-chunker(qseq,rnd)
				
			headers<-paste('fragment',seq(1:length(fragments)),sep='_')
				
			write.fasta(fragments,names=headers,file=query)

			blcmd1<-paste(blastn," -query ",query," -subject ",subjt1,
						  " -xdrop_gap_final 150 ",
						  " -outfmt '6 qseqid pident length qlen gaps' -out ",
						  blout1,sep='')
						 
			system(blcmd1,ignore.stderr=T,ignore.stdout=T)
			
			if (file.info(blout1)$size>0){
				
				tab1<-read.csv(blout1,sep='\t',header=F)
				tab1[,6]<-(tab1[,3]-tab1[,5])/tab1[,4]
				
				ani1<-mean(tab1[which(tab1[,6]>=.70 & tab1[,2]>=30),2])
				
			} else {
				
				ani1<-0
			}
			
			system(paste('rm -rf',blout1,query))
				
			# BLAST 2
							
			#query2<-comb[2,i]
			#subjt2<-comb[1,i]

			qfasta<-read.fasta(query2)
			qseq<-unlist(lapply(getSequence(qfasta),toupper))
	
			blout2<-paste('blout2',y,i,sep='.')
			query<-paste('query2',y,i,'fasta',sep='.')

			qfasta<-read.fasta(query2)
			qseq<-unlist(lapply(getSequence(qfasta),toupper))
				
			rnd<-round(length(qseq)/1020)
			fragments<-chunker(qseq,rnd)
				
			headers<-paste('fragment',seq(1:length(fragments)),sep='_')
				
			write.fasta(fragments,names=headers,file=query)
	
			blcmd2<-paste(blastn," -query ",query," -subject ",subjt2,
						  " -xdrop_gap_final 150 ",
						  " -outfmt '6 qseqid pident length qlen gaps' -out ",
						  blout2,sep='')
						 
			system(blcmd2,ignore.stderr=T,ignore.stdout=T)
						
			if (file.info(blout2)$size>0){
				
				tab2<-read.csv(blout2,sep='\t',header=F)
				tab2[,6]<-(tab2[,3]-tab2[,5])/tab2[,4]
				
				ani2<-mean(tab2[which(tab2[,6]>=.70 & tab2[,2]>=30),2])
				
			} else {
				
				ani2<-0
			}
						
			system(paste('rm -rf',blout2,query))
			
			# RESULT

			anib<-(ani1+ani2)/2
			
			u1<-unlist(strsplit(comb[1,i],'/'))
			u2<-unlist(strsplit(comb[2,i],'/'))
			
			l1<-length(u1)
			l2<-length(u2)
			
			n1<-gsub(pattern,'',u1[l1])
			n2<-gsub(pattern,'',u2[l2])

			result<-rbind(result,c(n1,n2,round(anib,digits=2)))
		}

		return(result)
	}
	
	if (reference=='all'){
	
		# Combinations #
		
		combina<-combn(genomes,2)
		combdim<-dim(combina)[2]

		chunks<-chunker(c(1:combdim),proc)

		registerDoMC(proc)

		final<-foreach(y=1:proc) %dopar% {
		
			aniparallel(y)	
		}
	
		final<-do.call('rbind',final)
		anib_result<-as.data.frame(final)
		colnames(anib_result)<-c('Genome.A','Genome.B','ANIb')
		
		save(anib_result,file='anib_result.Rdata')
		
	} else if (reference!='all'){
		
		if (reference%in%genomes){
			
			genomes2<-genomes[-which(genomes==reference)]
			refset<-rep(reference,length(genomes2))
			
			combina<-rbind(genomes2,refset)
			rownames(combina)<-NULL
			combdim<-dim(combina)[2]
			
			chunks<-chunker(c(1:combdim),proc)
			
			registerDoMC(proc)

			final<-foreach(y=1:proc) %dopar% {
		
				aniparallel(y)	
			}
	
			final<-do.call('rbind',final)
			anib_result<-as.data.frame(final)
			colnames(anib_result)<-c('Genome.A','Genome.B','ANIb')
	
			save(anib_result,file='anib_result.Rdata')
					
		} else {
			
			stop('ERROR: The reference genome is not in the list')
		}
	}
	
	system(paste('mkdir',outdir))
	system(paste('mv anib_result.Rdata',outdir))
}




