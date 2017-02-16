#' aai
#'
#' Calculates AAI (Average Amino acid Identity) values for all genomes present in the folder.
#' @param path is the full path to where genome annotations in fasta format are placed.
#' @param pattern a pattern (like '.faa') for recognizing the genome files.
#' @param outdir a name for the output directory that will be created for results.
#' @param proc is the number of threads for parallel computing.
#' @param reference is a genome name used to compare one vs. the rest. If set to 'all' (default), all vs. all is computed.
#' @keywords AAI species
#' @export
#' @examples
#' aai(path='./prodigal_output',pattern='.faa',reference='all',outdir='aai_result',proc=2)

aai<-function(

	pattern='.faa',
	reference='all',
	path='.',
	outdir='aai_output',
	proc=2

	)
			   
	{

	# Options #
	
	options(getClass.msg=FALSE)
	os<-.getOS()

	# Dependencies #

	suppressPackageStartupMessages(library(seqinr,quietly=T))

	# Define internal function #

	chunker<-function(m,n){

		s<-split(m,cut(seq_along(m),n,labels=F))
		return(s)
	}
	
	# Select OS #
	
	if (os=='linux'){
	
		blastp<-paste(system.file('blast',package='taxxo'),'/linux/blastp',sep='')
		mkbldb<-paste(system.file('blast',package='taxxo'),'/linux/makeblastdb',sep='')
		
	} else if (os=='darwin'){
		
		blastp<-paste(system.file('blast',package='taxxo'),'/darwin/blastp',sep='')
		mkbldb<-paste(system.file('blast',package='taxxo'),'/darwin/makeblastdb',sep='')
		
	} else {
		
		stop('ERROR: Unknown OS, unable to proceed.')
	}

	# List genomes #
	
	genomes<-gsub('//','/',list.files(path=path,pattern=pattern,full.names=T))
	g.names<-gsub(pattern,'',list.files(path=path,pattern=pattern))
	
	# Make databases #
	
	system('mkdir tmpaaidbs')
	
	setwd('./tmpaaidbs')
	
	for (g in 1:length(genomes)){
		
		cmd<-paste(mkbldb,' -in ../',genomes[g],' -dbtype prot -title ',g.names[g],' -out ',g.names[g],sep='')
		system(cmd,ignore.stderr=T,ignore.stdout=T)
	}
	
	setwd('../')
	
	# Define internal function #

	fastasplit<-function(f,n){

		fasta<-read.fasta(f)
		seqs<-lapply(getSequence(fasta),toupper)
		vecu<-seq(1:length(seqs))
		
		s<-split(vecu,cut(seq_along(vecu),n,labels=F))

		for (x in 1:length(s)){
		
			write.fasta(seqs[s[[x]]],names=paste('name',s[[x]],sep=''),file.out=paste('splitted_',x,'.faa',sep=''))
		}
	}

	# Run reciprocal AAI #

	if (reference=='all'){

		result<-NULL
	
		# Combinations #
		
		combina<-combn(genomes,2)
		combdim<-dim(combina)[2]

		for (i in 1:combdim){

			# BLAST 1
			
			query1<-combina[1,i]
				
			aux1<-strsplit(combina[2,i],'/')
			len1<-length(aux1[[1]])
			subjt1<-gsub(pattern,'',aux1[[1]][len1])

			fastasplit(f=query1,n=proc)

			cmd1.1<-paste("seq 1 ",proc,
				    " | xargs -P ",proc," -I {} ",blastp," -query splitted_{}.faa -db ./tmpaaidbs/",subjt1,
				    " -evalue 0.001 -outfmt '6 qseqid sseqid pident gaps length qlen' -out partial_{}",sep='')

			cmd1.2<-'cat partial_* > blout1'
	
			system(cmd1.1,wait=T)
			system(cmd1.2,wait=T)

			blout1<-read.csv('blout1',sep='\t',header=F)
			blout1[,7]<-(blout1[,5]-blout1[,4])/blout1[,6]
			blout1<-blout1[which(blout1[,3]>=50 & blout1[,7]>=.7),]

			system('rm -rf partial_* splitted_* blout1')

			# BLAST 2

			query2<-combina[2,i]
				
			aux2<-strsplit(combina[1,i],'/')
			len2<-length(aux2[[1]])
			subjt2<-gsub(pattern,'',aux2[[1]][len2])	

			fastasplit(f=query2,n=proc)

			cmd2.1<-paste("seq 1 ",proc,
				    " | xargs -P ",proc," -I {} ",blastp," -query splitted_{}.faa -db ./tmpaaidbs/",subjt2,
				    " -evalue 0.001 -outfmt '6 qseqid sseqid pident gaps length qlen' -out partial_{}",sep='')
			
			cmd2.2<-'cat partial_* > blout2'
	
			system(cmd2.1,wait=T)
			system(cmd2.2,wait=T)

			blout2<-read.csv('blout2',sep='\t',header=F)
			blout2[,7]<-(blout2[,5]-blout2[,4])/blout2[,6]
			blout2<-blout2[which(blout2[,3]>=50 & blout2[,7]>=.7),]

			system('rm -rf partial_* splitted_* blout2')
			
			# RESULT
			
			pairs1<-paste(blout1[,1],blout1[,2],sep='_')
			pairs2<-paste(blout2[,2],blout2[,1],sep='_')
			
			allpairs<-table(c(pairs1,pairs2))
			brh<-names(allpairs[which(allpairs>1)])
			
			aai1<-mean(c(blout1[which(pairs1%in%brh),3],blout2[which(pairs2%in%brh),3]))
			
			n1<-subjt2
			n2<-subjt1
			
			result<-rbind(result,c(n1,n2,round(aai1,digits=2)))
		}

		aai_result<-as.data.frame(result)
		colnames(aai_result)<-c('Genome A','Genome B','AAI')
	
		save(aai_result,file='aai_result.Rdata')

		
	} else if (reference!='all'){
		
		genaux<-do.call(rbind,strsplit(genomes,'/'))
		genomes2<-as.vector(genaux[,dim(genaux)[2]])
		
		if (reference%in%genomes2){
		
			result<-NULL
			
			wh<-which(genomes2==reference)
			genomes3<-genomes[-wh]
			refset<-rep(genomes[wh],length(genomes3))
			
			combina<-rbind(genomes3,refset)
			rownames(combina)<-NULL
			combdim<-dim(combina)[2]

			for (i in 1:combdim){

				# BLAST 1
			
				query1<-combina[1,i]
				
				aux1<-strsplit(combina[2,i],'/')
				len1<-length(aux1[[1]])
				subjt1<-gsub(pattern,'',aux1[[1]][len1])

				fastasplit(f=query1,n=proc)

				cmd1.1<-paste("seq 1 ",proc,
		                              " | xargs -P ",proc," -I {} ",blastp," -query splitted_{}.faa -db ./tmpaaidbs/",subjt1,
					      " -evalue 0.001 -outfmt '6 qseqid sseqid pident gaps length qlen' -out partial_{}",sep='')

				cmd1.2<-'cat partial_* > blout1'
	
				system(cmd1.1,wait=T)
				system(cmd1.2,wait=T)

				blout1<-read.csv('blout1',sep='\t',header=F)
				blout1[,7]<-(blout1[,5]-blout1[,4])/blout1[,6]
				blout1<-blout1[which(blout1[,3]>=50 & blout1[,7]>=.7),]

				system('rm -rf partial_* splitted_* blout1')

				# BLAST 2

				query2<-combina[2,i]
				
				aux2<-strsplit(combina[1,i],'/')
				len2<-length(aux2[[1]])
				subjt2<-gsub(pattern,'',aux2[[1]][len2])	

				fastasplit(f=query2,n=proc)

				cmd2.1<-paste("seq 1 ",proc,
				              " | xargs -P ",proc," -I {} ",blastp," -query splitted_{}.faa -db ./tmpaaidbs/",subjt2,
				              " -evalue 0.001 -outfmt '6 qseqid sseqid pident gaps length qlen' -out partial_{}",sep='')
			
				cmd2.2<-'cat partial_* > blout2'
	
				system(cmd2.1,wait=T)
				system(cmd2.2,wait=T)

				blout2<-read.csv('blout2',sep='\t',header=F)
				blout2[,7]<-(blout2[,5]-blout2[,4])/blout2[,6]
				blout2<-blout2[which(blout2[,3]>=50 & blout2[,7]>=.7),]

				system('rm -rf partial_* splitted_* blout2')
			
				# RESULT
			
				pairs1<-paste(blout1[,1],blout1[,2],sep='_')
				pairs2<-paste(blout2[,2],blout2[,1],sep='_')
			
				allpairs<-table(c(pairs1,pairs2))
				brh<-names(allpairs[which(allpairs>1)])
			
				aai1<-mean(c(blout1[which(pairs1%in%brh),3],blout2[which(pairs2%in%brh),3]))
			
				n1<-subjt2
				n2<-subjt1
			
				result<-rbind(result,c(n1,n2,round(aai1,digits=2)))
			}

			aai_result<-as.data.frame(result)
			colnames(aai_result)<-c('Genome A','Genome B','AAI')
	
			save(aai_result,file='aai_result.Rdata')
							
		} else {
			
			stop('ERROR: The reference genome is not in the list')
		}
	}
	
	system(paste('mkdir',outdir))
	system(paste('mv aai_result.Rdata',outdir))
	system('rm -rf tmpaaidbs blout1* blout2*')
}
