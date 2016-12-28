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

agios<-function(pattern='.ffn',
			    reference='all',
			    path='.',
			    outdir='agios_output',
			    proc=2)
			   
			   {
			   	
	# Options #
	
	options(getClass.msg=FALSE)
	os<-.getOS()
	
	# Dependencies #

	library(foreach,quietly=T)
	library(doMC,quietly=T)
	library(seqinr,quietly=T)
	library(msa,quietly=T)

	# Define internal function #

	chunker<-function(m,n){

		s<-split(m,cut(seq_along(m),n,labels=F))
		return(s)
	}
	
	# Select OS #
	
	if (os=='linux'){
	
		tblastx<-paste(system.file('blast',package='taxxo'),'/linux/tblastx',sep='')
		bldbcmd<-paste(system.file('blast',package='taxxo'),'/linux/blastdbcmd',sep='')
		
	} else if (os=='darwin'){
		
		tblastx<-paste(system.file('blast',package='taxxo'),'/darwin/tblastx',sep='')
		bldbcmd<-paste(system.file('blast',package='taxxo'),'/darwin/blastdbcmd',sep='')

	} else {
		
		stop('ERROR: Unknown OS, unable to proceed.')
	}

	# List genomes #
	
	genomes<-gsub('//','/',list.files(path=path,pattern=pattern,full.names=T))
	g.names<-gsub(pattern,'',list.files(path=path,pattern=pattern))
	
	# Make databases #
	
	system('mkdir databases')
	setwd('./databases')

	system(paste('cat ',paste(genomes,collapse=' '),' > all.ffn',sep=''))
	system('makeblastdb -in all.ffn -dbtype nucl -hash_index -parse_seqids -title ffns -out ffns',
		   ignore.stderr=T,ignore.stdout=T)

	setwd('../')

	# AGIOS parallel #
	
	agiosparallel<-function(y){

		chun<-chunks[[y]]	
		comb<-combina[,chun]

		result<-NULL

		for (i in 1:dim(comb)[2]){

			# BLAST 1

			query1<-comb[1,i]
			subjt1<-comb[2,i]

			blcmd1<-paste(tblastx," -query ",query1," -subject ",subjt1,
						 " -outfmt '6 qseqid sseqid pident qcovs sframe qframe' -out blout1.",
						 y,".",i,sep=''
                        )
                        
			system(blcmd1,ignore.stderr=T,ignore.stdout=T)
			
			blout1<-read.csv(paste('blout1.',y,'.',i,sep=''),sep='\t',header=F)
			blout1<-blout1[which(blout1[,5]==1 & blout1[,6]==1),]
			blout1<-blout1[which(blout1[,3]>=50 & blout1[,4]>=70),]
								
			# BLAST 2
							
			query2<-comb[2,i]
			subjt2<-comb[1,i]

			blcmd2<-paste(tblastx," -query ",query2," -subject ",subjt2,
						 " -outfmt '6 qseqid sseqid pident qcovs sframe qframe' -out blout2.",
						 y,".",i,sep=''
                        )
                        
			system(blcmd2,ignore.stderr=T,ignore.stdout=T)
			
			blout2<-read.csv(paste('blout2.',y,'.',i,sep=''),sep='\t',header=F)
			blout2<-blout2[which(blout2[,5]==1 & blout2[,6]==1),]
			blout2<-blout2[which(blout2[,3]>=50 & blout2[,4]>=70),]
			
			# RESULT
			
			pairs1<-paste(blout1[,1],blout1[,2],sep='@')
			pairs2<-paste(blout2[,2],blout2[,1],sep='@')
			
			allpairs<-table(c(pairs1,pairs2))
			brh<-names(allpairs[which(allpairs>1)])
			
			partial<-NULL
			
			for (b in 1:length(brh)){
				
				brhout<-paste('brh',y,'.',b,'.tmp',sep='')
				brhfas<-paste('brh',y,'.',b,'.tmp2.fasta',sep='')
				
				cat(unlist(strsplit(brh[b],'@',fixed=T)),sep='\n',file=brhout)
				
				cmd<-paste(bldbcmd,' -entry_batch ',brhout,' -db ./databases/ffns > ',brhfas,sep='')
				
				system(cmd)
				
				alignment<-msa(inputSeqs=brhfas,method='Muscle',type='dna')	
				aliconver<-msaConvert(alignment,type='seqinr::alignment')
				
				dista<-dist.alignment(aliconver)^2
				ident<-round((1-dista)*100,digits=2)

				partial<-c(partial,ident)
			}
			
			u1<-unlist(strsplit(comb[1,i],'/'))
			u2<-unlist(strsplit(comb[2,i],'/'))
			
			l1<-length(u1)
			l2<-length(u2)
			
			n1<-gsub(pattern,'',u1[l1])
			n2<-gsub(pattern,'',u2[l2])
			
			result<-rbind(result,c(n1,n2,round(mean(partial),digits=2)))
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
		
			agiosparallel(y)	
		}
	
		final<-do.call('rbind',final)
		agios_result<-as.data.frame(final)
		colnames(agios_result)<-c('Genome.A','Genome.B','AGIOS')
	
		save(agios_result,file='agios_result.Rdata')
		
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
		
				agiosparallel(y)	
			}
	
			final<-do.call('rbind',final)
			agios_result<-as.data.frame(final)
			colnames(agios_result)<-c('Genome.A','Genome.B','AGIOS')
	
			save(agios_result,file='agios_result.Rdata')
					
		} else {
			
			stop('ERROR: The reference genome is not in the list')
		}
	}
	
	system(paste('mkdir',outdir))
	system(paste('mv agios_result.Rdata',outdir))
}