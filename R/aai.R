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
		nams<-getName(fasta)
		vecu<-seq(1:length(seqs))
		
		s<-split(vecu,cut(seq_along(vecu),n,labels=F))

		for (x in 1:length(s)){
		
			write.fasta(seqs[s[[x]]],names=nams[s[[x]]],file.out=paste('splitted_',x,'.faa',sep=''))
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

			# BLAST 2

			query2<-combina[2,i]
				
			aux2<-strsplit(combina[1,i],'/')
			len2<-length(aux2[[1]])
			subjt2<-gsub(pattern,'',aux2[[1]][len2])	

			fastasplit(f=query2,n=proc)

			cmd2.1<-paste("seq 1 ",proc,
				    " | xargs -P ",proc," -I {} ",blastp," -query splitted_{}.faa -db ./tmpaaidbs/",subjt2,
				    " -evalue 0.001 -outfmt '6 sseqid qseqid pident gaps length qlen' -out partial_{}",sep='')
			
			cmd2.2<-'cat partial_* > blout2'
	
			system(cmd2.1,wait=T)
			system(cmd2.2,wait=T)
			
			# RESULT

			cmd3<-'cat blout1 blout2 > blout3'

			system(cmd3)
			system('rm -rf partial_* splitted_* blout*')

			blout<-read.csv('blout3',sep='\t',header=F)
			blout[,7]<-(blout[,5]-blout[,4])/blout[,6]
			blout<-blout[which(blout[,3]>=50 & blout[,7]>=.7),]
			
			pairs<-table(paste(blout[,1],blout[,2],sep='_'))
			brh<-names(pairs[which(pairs>1)])
			
			aai<-mean(blout[which(pairs%in%brh),3])
			
			n1<-subjt2
			n2<-subjt1
			
			result<-rbind(result,c(n1,n2,round(aai,digits=2)))
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
			
				# RESULT
			
				cmd3<-'cat blout1 blout2 > blout3'

				system(cmd3)
				system('rm -rf partial_* splitted_* blout*')

				blout<-read.csv('blout3',sep='\t',header=F)
				blout[,7]<-(blout[,5]-blout[,4])/blout[,6]
				blout<-blout[which(blout[,3]>=50 & blout[,7]>=.7),]
			
				pairs<-table(paste(blout[,1],blout[,2],sep='_'))
				brh<-names(pairs[which(pairs>1)])
			
				aai<-mean(blout[which(pairs%in%brh),3])
			
				n1<-subjt2
				n2<-subjt1
			
				result<-rbind(result,c(n1,n2,round(aai,digits=2)))
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
	system('rm -rf tmpaaidbs')
}
