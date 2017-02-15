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

aai<-function(pattern='.faa',
			   reference='all',
			   path='.',
			   outdir='aai_output',
			   proc=2)
			   
			   {

	# Options #
	
	options(getClass.msg=FALSE)
	os<-.getOS()

	# Dependencies #

	suppressPackageStartupMessages(library(foreach,quietly=T))
	suppressPackageStartupMessages(library(doMC,quietly=T))
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
	
	# Check multi cores #
	
	if(proc>length(genomes)){
		
		proc<-length(genomes)
	}

	# AAI parallel #
	
	aaiparallel<-function(y){

		chun<-chunks[[y]]	
		
		if (length(chun)>1){
		
			comb<-combina[,chun]

			result<-NULL

			for (i in 1:dim(comb)[2]){

				# BLAST 1

				query1<-comb[1,i]
				
				aux1<-strsplit(comb[2,i],'/')
				len1<-length(aux1[[1]])
				subjt1<-gsub(pattern,'',aux1[[1]][len1])

				blcmd1<-paste(blastp," -query ",query1," -db ./tmpaaidbs/",subjt1,
							 " -outfmt '6 qseqid sseqid pident length qlen gaps' -out blout1.",
							 y,".",i,sep='')
                        
				system(blcmd1,ignore.stderr=T,ignore.stdout=T)
			
				blout1<-read.csv(paste('blout1.',y,'.',i,sep=''),sep='\t',header=F)
				blout1[,7]<-(blout1[,4]-blout1[,6])/blout1[,5]
				blout1<-blout1[which(blout1[,3]>=50 & blout1[,7]>=.7),]
								
				# BLAST 2
							
				query2<-comb[2,i]
				
				aux2<-strsplit(comb[1,i],'/')
				len2<-length(aux2[[1]])
				subjt2<-gsub(pattern,'',aux2[[1]][len2])

				blcmd2<-paste(blastp," -query ",query2," -db ./tmpaaidbs/",subjt2,
							 " -outfmt '6 qseqid sseqid pident length qlen gaps' -out blout2.",
							 y,".",i,sep='')
                        
				system(blcmd2,ignore.stderr=T,ignore.stdout=T)
			
				blout2<-read.csv(paste('blout2.',y,'.',i,sep=''),sep='\t',header=F)
				blout2[,7]<-(blout2[,4]-blout2[,6])/blout2[,5]
				blout2<-blout2[which(blout2[,3]>=50 & blout2[,7]>=.7),]
			
				# RESULT
			
				pairs1<-paste(blout1[,1],blout1[,2],sep='_')
				pairs2<-paste(blout2[,2],blout2[,1],sep='_')
			
				allpairs<-table(c(pairs1,pairs2))
				brh<-names(allpairs[which(allpairs>1)])
			
				aai1<-mean(c(blout1[which(pairs1%in%brh),3],blout2[which(pairs2%in%brh),3]))
			
				u1<-unlist(strsplit(comb[1,i],'/'))
				u2<-unlist(strsplit(comb[2,i],'/'))
			
				l1<-length(u1)
				l2<-length(u2)
			
				n1<-gsub(pattern,'',u1[l1])
				n2<-gsub(pattern,'',u2[l2])
			
				result<-rbind(result,c(n1,n2,round(aai1,digits=2)))
			}
			
		} else if (length(chun)==1){

			comb<-combina[,chun]
			
			result<-NULL

			# BLAST 1

			query1<-comb[1]
			
			aux1<-strsplit(comb[2],'/')
			len1<-length(aux1[[1]])
			subjt1<-gsub(pattern,'',aux1[[1]][len1])

			blcmd1<-paste(blastp," -query ",query1," -db ./tmpaaidbs/",subjt1,
						 " -outfmt '6 qseqid sseqid pident length qlen gaps' -out blout1.",
						 y,".",i,sep='')
                        
			system(blcmd1,ignore.stderr=T,ignore.stdout=T)
			
			blout1<-read.csv(paste('blout1.',y,'.',i,sep=''),sep='\t',header=F)
			blout1[,7]<-(blout1[,4]-blout1[,6])/blout1[,5]
			blout1<-blout1[which(blout1[,3]>=50 & blout1[,7]>=.7),]

			# BLAST 2
							
			query2<-comb[2]
			
			aux2<-strsplit(comb[1],'/')
			len2<-length(aux2[[1]])
			subjt2<-gsub(pattern,'',aux2[[1]][len2])

			blcmd2<-paste(blastp," -query ",query2," -db ./tmpaaidbs/",subjt2,
						 " -outfmt '6 qseqid sseqid pident length qlen gaps' -out blout2.",
						 y,".",i,sep='')
                        
			system(blcmd2,ignore.stderr=T,ignore.stdout=T)
			
			blout2<-read.csv(paste('blout2.',y,'.',i,sep=''),sep='\t',header=F)
			blout2[,7]<-(blout2[,4]-blout2[,6])/blout2[,5]
			blout2<-blout2[which(blout2[,3]>=50 & blout2[,7]>=.7),]

			# RESULT
			
			pairs1<-paste(blout1[,1],blout1[,2],sep='_')
			pairs2<-paste(blout2[,2],blout2[,1],sep='_')
			
			allpairs<-table(c(pairs1,pairs2))
			brh<-names(allpairs[which(allpairs>1)])
			
			aai1<-mean(c(blout1[which(pairs1%in%brh),3],blout2[which(pairs2%in%brh),3]))
			
			u1<-unlist(strsplit(comb[1],'/'))
			u2<-unlist(strsplit(comb[2],'/'))
			
			l1<-length(u1)
			l2<-length(u2)
			
			n1<-gsub(pattern,'',u1[l1])
			n2<-gsub(pattern,'',u2[l2])
			
			result<-rbind(result,c(n1,n2,round(aai1,digits=2)))
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
		
			aaiparallel(y)	
		}
	
		final<-do.call('rbind',final)
		aai_result<-as.data.frame(final)
		colnames(aai_result)<-c('Genome A','Genome B','AAI')
	
		save(aai_result,file='aai_result.Rdata')
		
	} else if (reference!='all'){
		
		genaux<-do.call(rbind,strsplit(genomes,'/'))
		genomes2<-as.vector(genaux[,dim(genaux)[2]])
		
		if (reference%in%genomes2){
			
			wh<-which(genomes2==reference)
			genomes3<-genomes[-wh]
			refset<-rep(genomes[wh],length(genomes3))
			
			combina<-rbind(genomes3,refset)
			rownames(combina)<-NULL
			combdim<-dim(combina)[2]
			
			chunks<-chunker(c(1:combdim),proc)
			
			registerDoMC(proc)

			final<-foreach(y=1:proc) %dopar% {
		
				aaiparallel(y)	
			}
	
			final<-do.call('rbind',final)
			aai_result<-as.data.frame(final)
			colnames(aai_result)<-c('Genome.A','Genome.B','AAI')
	
			save(aai_result,file='aai_result.Rdata')
					
		} else {
			
			stop('ERROR: The reference genome is not in the list')
		}
	}
	
	system(paste('mkdir',outdir))
	system(paste('mv aai_result.Rdata',outdir))
	system('rm -rf tmpaaidbs blout1* blout2*')
}
