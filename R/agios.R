#' agios
#'
#' Calculates AGIOS (Average Genome Identity of Orthologous Sequences) values for all genomes present in the folder.
#' @param path is the full path to where genome annotations in fasta format are placed.
#' @param pattern a pattern (like '.ffn') for recognizing the genome annotation files in nucleotides.
#' @param outdir a name for the output directory that will be created for results.
#' @param proc is the number of threads for parallel computing.
#' @param reference is a genome name used to compare one vs. the rest. If set to 'all' (default), all vs. all is computed.
#' @keywords AGIOS species
#' @export
#' @examples
#' agios(path='./prodigal_output',pattern='.ffn',reference='all',outdir='agios_result',proc=2)

agios<-function(

	pattern='.ffn',
	reference='all',
	path='.',
	outdir='agios_output',
	proc=2

	)
			   
	{

	# Options #
	
	options(getClass.msg=FALSE)
	os<-.getOS()

	# Dependencies #

	suppressPackageStartupMessages(library(seqinr,quietly=T))
	suppressPackageStartupMessages(library(msa,quietly=T))

	# Select OS #
	
	if (os=='linux'){
	
		blastp<-paste(system.file('blast',package='taxxo'),'/linux/blastp',sep='')
		mkbldb<-paste(system.file('blast',package='taxxo'),'/linux/makeblastdb',sep='')
		bldbcm<-paste(system.file('blast',package='taxxo'),'/linux/blastdbcmd',sep='')		

	} else if (os=='darwin'){
		
		blastp<-paste(system.file('blast',package='taxxo'),'/darwin/blastp',sep='')
		mkbldb<-paste(system.file('blast',package='taxxo'),'/darwin/makeblastdb',sep='')
		bldbcm<-paste(system.file('blast',package='taxxo'),'/darwin/blastdbcmd',sep='')
		
	} else {
		
		stop('ERROR: Unknown OS, unable to proceed.')
	}

	# List genomes #
	
	genomes<-gsub('//','/',list.files(path=path,pattern=pattern,full.names=T))
	g.names<-gsub(pattern,'',list.files(path=path,pattern=pattern))

	# Translate ffns #

	system('mkdir agiosfaa')

	setwd('./agiosfaa')

	for (g in 1:length(genomes)){

		f<-read.fasta(paste('../',genomes[g],sep=''))
		s<-lapply(getSequence(f),toupper)
		d<-lapply(s,seqinr::translate)
		n<-getName(f)

		write.fasta(d,names=n,file=paste(g.names[g],'.faa',sep=''))
	}

	path2<-getwd()

	genomes2<-gsub('//','/',list.files(path=path2,pattern='.faa',full.names=T))

	setwd('../')
	
	# Make databases #
	
	system('mkdir tmpagiosdbs')
	
	setwd('./tmpagiosdbs')
	
	for (g in 1:length(genomes2)){
		
		cmd<-paste(mkbldb,' -in ',genomes2[g],' -dbtype prot -title ',g.names[g],' -out ',g.names[g],sep='')
		system(cmd,ignore.stderr=T,ignore.stdout=T)
	}
	
	system(paste('cat ',paste(paste('../',genomes,sep=''),collapse=' '),' > agiosall.ffn',sep=''))
	system(paste(mkbldb,' -in agiosall.ffn -dbtype nucl -hash_index -parse_seqids -title ffns -out ffns',sep=''),ignore.stderr=T,ignore.stdout=T)

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

	#

	chunker<-function(m,n){

		s<-split(m,cut(seq_along(m),n,labels=F))
		return(s)
	}

	# Run reciprocal AGIOS #

	if (reference=='all'){

		result<-NULL
	
		# Combinations #
		
		combina<-combn(genomes2,2)
		combdim<-dim(combina)[2]

		for (i in 1:combdim){

			# BLAST 1
			
			query1<-combina[1,i]
				
			aux1<-strsplit(combina[2,i],'/')
			len1<-length(aux1[[1]])
			subjt1<-gsub('.faa','',aux1[[1]][len1])

			fastasplit(f=query1,n=proc)

			cmd1.1<-paste("seq 1 ",proc,
				    " | xargs -P ",proc," -I {} ",blastp," -query splitted_{}.faa -db ./tmpagiosdbs/",subjt1,
				    " -evalue 0.001 -outfmt '6 qseqid sseqid pident gaps length qlen' -out partial_{}",sep='')

			cmd1.2<-'cat partial_* > blout1'
	
			system(cmd1.1,wait=T)
			system(cmd1.2,wait=T)

			# BLAST 2

			query2<-combina[2,i]
				
			aux2<-strsplit(combina[1,i],'/')
			len2<-length(aux2[[1]])
			subjt2<-gsub('.faa','',aux2[[1]][len2])	

			fastasplit(f=query2,n=proc)

			cmd2.1<-paste("seq 1 ",proc,
				    " | xargs -P ",proc," -I {} ",blastp," -query splitted_{}.faa -db ./tmpagiosdbs/",subjt2,
				    " -evalue 0.001 -outfmt '6 sseqid qseqid pident gaps length qlen' -out partial_{}",sep='')
			
			cmd2.2<-'cat partial_* > blout2'
	
			system(cmd2.1,wait=T)
			system(cmd2.2,wait=T)
			
			# RESULT

			cmd3<-'cat blout1 blout2 > blout3'

			system(cmd3)

			blout<-read.csv('blout3',sep='\t',header=F)

			system('rm -rf partial_* splitted_* blout*')

			blout[,7]<-(blout[,5]-blout[,4])/blout[,6]
			blout<-blout[which(blout[,3]>=30 & blout[,7]>=.50),]
			
			pairs<-table(paste(blout[,1],blout[,2],sep='@'))
			brh<-names(pairs[which(pairs>1)])
			
			system('mkdir agiosbrhs')

			setwd('./agiosbrhs')

			for (b in 1:length(brh)){

				ub<-unlist(strsplit(brh[b],'@'))
				cat(ub,sep='\n',file=paste('brh_',b,'.txt',sep=''))
			}

			texts<-list.files(pattern='.txt')

			for (e in texts){

				cmd5<-paste(bldbcm,' -entry_batch ',e,' -db ../tmpagiosdbs/ffns > ',gsub('.txt','.fasta',e),sep='')
				system(cmd5,ignore.stderr=T)
			}

			fastas<-list.files(pattern='.fasta')

			for (f in fastas){

				partial<-NULL
				
				aux<-capture.output(
				alignment<-msa(inputSeqs=f,method='ClustalOmega',type='dna'))
				aliconver<-msaConvert(alignment,type='seqinr::alignment')
				
				dista<-dist.alignment(aliconver)^2
				ident<-round((1-dista)*100,digits=2)

				partial<-c(partial,ident)
			}
			
			agios<-mean(partial)
			
			n1<-subjt2
			n2<-subjt1
			
			result<-rbind(result,c(n1,n2,round(agios,digits=2)))

			setwd('../')
		}

		agios_result<-as.data.frame(result)
		colnames(agios_result)<-c('Genome A','Genome B','AGIOS')
	
		save(agios_result,file='agios_result.Rdata')

	} else if (reference!='all'){
		
		genaux<-do.call(rbind,strsplit(genomes2,'/'))
		genomes3<-as.vector(genaux[,dim(genaux)[2]])
		
		if (reference%in%genomes3){
		
			result<-NULL
			
			wh<-which(genomes3==reference)
			genomes4<-genomes2[-wh]
			refset<-rep(genomes2[wh],length(genomes4))
			
			combina<-rbind(genomes4,refset)
			rownames(combina)<-NULL
			combdim<-dim(combina)[2]

			for (i in 1:combdim){

					query1<-combina[1,i]
				
			aux1<-strsplit(combina[2,i],'/')
			len1<-length(aux1[[1]])
			subjt1<-gsub('.faa','',aux1[[1]][len1])

			fastasplit(f=query1,n=proc)

			cmd1.1<-paste("seq 1 ",proc,
				    " | xargs -P ",proc," -I {} ",blastp," -query splitted_{}.faa -db ./tmpagiosdbs/",subjt1,
				    " -evalue 0.001 -outfmt '6 qseqid sseqid pident gaps length qlen' -out partial_{}",sep='')

			cmd1.2<-'cat partial_* > blout1'
	
			system(cmd1.1,wait=T)
			system(cmd1.2,wait=T)

			# BLAST 2

			query2<-combina[2,i]
				
			aux2<-strsplit(combina[1,i],'/')
			len2<-length(aux2[[1]])
			subjt2<-gsub('.faa','',aux2[[1]][len2])	

			fastasplit(f=query2,n=proc)

			cmd2.1<-paste("seq 1 ",proc,
				    " | xargs -P ",proc," -I {} ",blastp," -query splitted_{}.faa -db ./tmpagiosdbs/",subjt2,
				    " -evalue 0.001 -outfmt '6 sseqid qseqid pident gaps length qlen' -out partial_{}",sep='')
			
			cmd2.2<-'cat partial_* > blout2'
	
			system(cmd2.1,wait=T)
			system(cmd2.2,wait=T)
			
			# RESULT

			cmd3<-'cat blout1 blout2 > blout3'

			system(cmd3)

			blout<-read.csv('blout3',sep='\t',header=F)

			system('rm -rf partial_* splitted_* blout*')

			blout[,7]<-(blout[,5]-blout[,4])/blout[,6]
			blout<-blout[which(blout[,3]>=30 & blout[,7]>=.50),]
			
			pairs<-table(paste(blout[,1],blout[,2],sep='@'))
			brh<-names(pairs[which(pairs>1)])
			
			system('mkdir agiosbrhs')

			setwd('./agiosbrhs')

			for (b in 1:length(brh)){

				ub<-unlist(strsplit(brh[b],'@'))
				cat(ub,sep='\n',file=paste('brh_',b,'.txt',sep=''))
			}

			texts<-list.files(pattern='.txt')

			for (e in texts){

				cmd5<-paste(bldbcm,' -entry_batch ',e,' -db ../tmpagiosdbs/ffns > ',gsub('.txt','.fasta',e),sep='')
				system(cmd5,ignore.stderr=T)
			}

			fastas<-list.files(pattern='.fasta')

			for (f in fastas){

				partial<-NULL
			
				aux<-capture.output(
				alignment<-msa(inputSeqs=f,method='ClustalOmega',type='dna'))
				aliconver<-msaConvert(alignment,type='seqinr::alignment')
				
				dista<-dist.alignment(aliconver)^2
				ident<-round((1-dista)*100,digits=2)

				partial<-c(partial,ident)
			}
			
			agios<-mean(partial)
			
			n1<-subjt2
			n2<-subjt1
			
			result<-rbind(result,c(n1,n2,round(agios,digits=2)))

			setwd('../')
		}

		agios_result<-as.data.frame(result)
		colnames(agios_result)<-c('Genome A','Genome B','AGIOS')
	
		save(agios_result,file='agios_result.Rdata')
							
		} else {
			
			stop('ERROR: The reference genome is not in the list')
		}
	}
	
	system(paste('mkdir',outdir))
	system(paste('mv agios_result.Rdata',outdir))
	system('rm -rf tmpagiosdbs agiosbrhs agiosfaa')
}
