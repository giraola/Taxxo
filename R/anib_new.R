#' anib_new
#'
#' Calculates ANIb (Average Nucleotide Identity using Blast) values for all genomes present in a directory.
#' @param path is the full path to where genomes assemblies in fasta format are placed.
#' @param pattern a pattern (like '.fna') for recognizing the genome files.
#' @param outdir a name for the output directory that will be created for results.
#' @param proc is the number of threads for parallel computing.
#' @param reference is a genome name used to compare one vs. the rest. If set to 'all' (default), all vs. all is computed.
#' @param win is the size of the fragments to compare (default = 1020).
#' @param winst is the step size to generate fragments (default = 0, no overlap).
#' @param imin is the identity cut-off to report a hit (defualt = 30).
#' @param cmin is the alignment query coverage cut-off to report a hit (defaut = 0.7).

#' @keywords ANI species blastn
#' @export
#' @examples
#' anib_new(path='.',pattern='.fna',reference='all',outdir='anib_result',proc=2)

anib_new <- function(

			pattern='.fna',
			reference='all',
			path='.',
	
			win=1020,
			winst=0,
			imin=30,
			cmin=0.7,
			
			outdir='anib_result',
			proc=2
			
		) {

	
	# Options #
	
	options(getClass.msg=FALSE)
	os <- .getOS()
		
	# Select OS #
	
	if (os=='linux'){
	
		blastn   <- paste(system.file('blast',package='taxxo'),'/linux/blastn',sep='')
		mkbldb   <- paste(system.file('blast',package='taxxo'),'/linux/makeblastdb',sep='')
		
	} else if (os=='darwin'){
		
		blastn   <- paste(system.file('blast',package='taxxo'),'/darwin/blastn',sep='')
		mkbldb   <- paste(system.file('blast',package='taxxo'),'/darwin/makeblastdb',sep='')

	} else {
		
		stop('ERROR: Unknown OS, unable to proceed.')
	}
	
	# Dependencies #

	suppressMessages(library(foreach,quietly=T))
	suppressMessages(library(doMC,quietly=T))
	suppressMessages(library(seqinr,quietly=T))
	suppressMessages(library(Biostrings,quietly=T))
	suppressMessages(library(parallel,quietly=T))

	
	# Define internal function #

	chunker <- function(m,n){

		s <- split(m,cut(seq_along(m),n,labels=F))
	
		return(s)
	}
	
	# spliter
	
	spliter <- function(
	
				w=win,
				s=winst,
				v
				
	) {
	
		n  <- length(v)
		x  <- 1
		x2 <- w
		o  <- 1
		l  <- list()

		if ( n>w & x2<=(n-s) ) {

			while ( x2<=(n-s) ) {
	
				ini    <- x
				fin    <- x2
				l[[o]] <- v[ini:fin]
				x      <- x+s
				x2     <- x2+s
				o      <- o+1

			}

			l[[o]] <- v[(fin+1):n]
			
		} else {
						
			l[[1]] <- v
		}
	
		return(l)
	}
	
	# List genomes #
	
	genomes <- gsub('//','/',list.files(path=path,pattern=pattern,full.names=T))
	g.names <- gsub(pattern,'',list.files(path=path,pattern=pattern))
	
	if( proc>length(genomes) ){
		
		proc <- length(genomes)
	}
	
	# Make databases #
	
	system('mkdir blast_databases')	
	setwd('blast_databases')
		
	fdatabase <- function(f) {
			
		cmd1 <- paste(mkbldb,
			      ' -in .',genomes[f],
			      ' -dbtype nucl ',
			      ' -title ',g.names[f],
			      ' -out ',g.names[f],
			      sep='')
			
		system(cmd1,ignore.stderr=T,ignore.stdout=T)
	}
	
		mclapply(1:length(genomes),fdatabase,mc.cores=proc) -> aux1
		
		rm(aux1)
		
		setwd('../')
			
	# ANIb parallel #
	
	aniparallel <- function(y) {
		
		chun <- chunks[[y]]

			comb <- combina[,chun]
		
			result <- NULL

			for (i in 1:dim(comb)[2]) {

				# BLAST 1

				query1 <- comb[1,i]
				
				blout1 <- paste('blout1',y,i,sep='.')
				query  <- paste('query1',y,i,'fasta',sep='.')

				qfasta <- read.fasta(query1)
				qseq   <- lapply(getSequence(qfasta),toupper)
	
				lseq   <- lapply(qseq,function(s){spliter(w=win,s=winst,v=s)})
				lseq2  <- unlist(lseq,recursive=F)
				llen   <- length(lseq2)
				
				fragms <- paste('fragment',seq(1:llen),sep='_')
				
				s1     <- lapply(lseq2,c2s)
				s2     <- lapply(s1,DNAString)
				s3     <- DNAStringSet(s2)
				
				names(s3) <- fragms
				
				writeXStringSet(s3,file=query,format='fasta')
				
							
				dbase1 <- gsub('//','/',
					      paste(path,'/blast_databases/',gsub('.fna','',comb[2,i]),sep=''))

				blcmd1 <- paste(blastn," -query ",query," -db ",dbase1,
							  " -max_target_seqs 1 -penalty -1 -xdrop_gap_final 150",
							  " -outfmt '6 qseqid pident length qlen gaps' -out ",
							  blout1,sep='')
						 
				system(blcmd1,ignore.stderr=T,ignore.stdout=T)
			
				if (file.info(blout1)$size>0) {
				
					tab1     <- read.csv(blout1,sep='\t',header=F)
					tab1[,6] <- (tab1[,3]-tab1[,5])/tab1[,4]
					ani1     <- mean(tab1[which(tab1[,6]>=cmin & tab1[,2]>=imin),2])
				
				} else {
				
					ani1 <- 0
				}
			
				system(paste('rm -rf',blout1,query))
			

				# BLAST 2
							
				query2<-comb[2,i]
				
				blout2 <- paste('blout2',y,i,sep='.')
				query  <- paste('query2',y,i,'fasta',sep='.')
				
				qfasta <- read.fasta(query2)
				qseq   <- lapply(getSequence(qfasta),toupper)

				lseq   <- lapply(qseq,function(s){spliter(w=win,s=winst,v=s)})
				lseq2  <- unlist(lseq,recursive=F)
				llen   <- length(lseq2)
				
				fragms <- paste('fragment',seq(1:llen),sep='_')
				
				s1     <- lapply(lseq2,c2s)
				s2     <- lapply(s1,DNAString)
				s3     <- DNAStringSet(s2)
				
				names(s3) <- fragms
				
				writeXStringSet(s3,file=query,format='fasta')
				
					
				dbase2 <- gsub('//','/',
					      paste(path,'/blast_databases/',gsub('.fna','',comb[2,i]),sep=''))

				blcmd2 <- paste(blastn," -query ",query," -db ",dbase2,
						 	  " -max_target_seqs 1 -penalty -1 -xdrop_gap_final 150",
						 	  " -outfmt '6 qseqid pident length qlen gaps' -out ",
							  blout2,sep='')

				system(blcmd2,ignore.stderr=T,ignore.stdout=T)
		
				if (file.info(blout2)$size>0) {
			
					tab2     <- read.csv(blout2,sep='\t',header=F)
					tab2[,6] <- (tab2[,3]-tab2[,5])/tab2[,4]	
					ani2     <- mean(tab2[which(tab2[,6]>=cmin & tab2[,2]>=imin),2])
				
				} else {
				
					ani2 <- 0
				}
			
				system(paste('rm -rf',blout2,query))
	
				# RESULT

				anib <- (ani1+ani2)/2	
				u1   <- unlist(strsplit(comb[1,i],'/'))
				u2   <- unlist(strsplit(comb[2,i],'/'))
			
				l1   <- length(u1)
				l2   <- length(u2)
			
				n1   <- gsub(pattern,'',u1[l1])
				n2   <- gsub(pattern,'',u2[l2])

				result <- rbind(result,c(n1,n2,round(anib,digits=2)))
			}
		
		return(result)
	}
	
	
	### RUN ANIb ###
	
	
	if (reference=='all') {
	
		# Combinations #
		
		combina <- combn(genomes,2)
		combdim <- dim(combina)[2]
		chunks  <- chunker(c(1:combdim),proc)
		
		registerDoMC(proc)

		final <- foreach(y=1:proc) %dopar% {
		
			aniparallel(y)	
		}
	
		final                 <- do.call('rbind',final)
		anib_result           <- as.data.frame(final)
		colnames(anib_result) <- c('Genome.A','Genome.B','ANIb')
		
		save(anib_result,file='anib_result.Rdata')
		
	} else if (reference!='all') {
		
		genaux   <- do.call(rbind,strsplit(genomes,'/'))
		genomes2 <- as.vector(genaux[,dim(genaux)[2]])
		
		if (reference%in%genomes2) {
			
			wh       <- which(genomes2==reference)
			genomes3 <- genomes[-wh]
			refset   <- rep(genomes[wh],length(genomes3))
			
			combina           <- rbind(genomes3,refset)
			rownames(combina) <- NULL
			combdim           <- dim(combina)[2]
			
			chunks <- chunker(c(1:combdim),proc)
			
			registerDoMC(proc)

			final <- foreach(y=1:proc) %dopar% {
		
				aniparallel(y)	
			}
	
			final                 <- do.call('rbind',final)
			anib_result           <- as.data.frame(final)
			colnames(anib_result) <- c('Genome.A','Genome.B','ANIb')
	
			save(anib_result,file='anib_result.Rdata')
					
		} else {
			
			stop('ERROR: The reference genome is not in the list')
		}
	}
	
	system('rm -rf blast_databases')
	system(paste('mkdir',outdir))
	system(paste('mv anib_result.Rdata',outdir))
}
