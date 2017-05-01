#' uprot
#'
#' Identifies and extracts a set of 40 prokaryotic universal proteins used as phylogenetic markers.
#' @param path is the full path to where genome annotations in fasta format are placed.
#' @param pattern a pattern (like '.faa') for recognizing the annotation files.
#' @param outdir a name for the output directory that will be created for results.
#' @param align takes a logical (default, TRUE) indicating if multiple sequence alignment is performed.
#' @param phylogeny takes 'NJ' (default) for a Neighbor-Joining tree, 'ML' for a Maximum-Likelihood tree or 'NO' for avoiding tree reconstruction.
#' @param proc is the number of threads used for protein searches.
#' @keywords ribosomal proteins universal markers
#' @export
#' @examples
#' uprot(path='./prodigal_output',pattern='.faa',outdir='uprot_output',proc=2)

uprot<-function(pattern='.faa',
			    path='.',
			    outdir='uprot_output',
			    align=TRUE,
			    phylogeny='NJ',
			    proc=2)
			   
			    {

	# Options #
	
	options(getClass.msg=FALSE)
	gw<-getwd()
	os<-.getOS()
	
	# Dependencies #
	
	suppressMessages(require(seqinr,quietly=T))
	suppressMessages(require(foreach,quietly=T))
	suppressMessages(require(doMC,quietly=T))
	suppressMessages(require(msa,quietly=T))
	suppressMessages(require(plyr,quietly=T))
	suppressMessages(require(ape,quietly=T))


	
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
		fasttree<-paste(system.file('fasttree',package='taxxo'),'/linux/FastTree',sep='')
		
	} else if (os=='darwin'){
		
		hmmsearch<-paste(system.file('hmmer3',package='taxxo'),'/darwin/hmmsearch',sep='')
		bldbcd<-paste(system.file('blast',package='taxxo'),'/darwin/blastdbcmd',sep='')
		mkbldb<-paste(system.file('blast',package='taxxo'),'/darwin/makeblastdb',sep='')
		fasttree<-paste(system.file('fasttree',package='taxxo'),'/darwin/FastTree',sep='')
		
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
			
			system('rm -rf out_tmp search_tmp')
			
			lr<-length(rl)
			
			if (lr>0){
			
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
  				
  					sc<-htab$Score
  				
  					if ((sc[1]!=sc[2]) & dimi==2){
  				
  						maxi<-max(htab$Score)
  						gene<-as.vector(htab[which(htab$Score==maxi),2])
  						
  					} else {
  						
  						gene<-as.vector(htab[1,2])
  					}
  				
  				} else if (dimi==1){
  				
  					gene<-as.vector(htab[1,2])
  				}
  				
  				result[f,h]<-gene
  			}
		}
	}
  			
  	system(paste('mkdir',outdir))
  			
  	presence_absence_uprot<-result		
  	save(presence_absence_uprot,file='presence_absence_uprot.Rdata')				
 	system(paste('mv presence_absence_uprot.Rdata',outdir),ignore.stdout=T)
	
 	# Make databases
 	
	cmd2<-paste('cat ',paste(flist,collapse=' '),' > all_genomes_uprot.faa',sep='')
	system(cmd2)
	
	cmd3<-paste(mkbldb,
		    ' -in all_genomes_uprot.faa',
		    ' -dbtype prot -hash_index -parse_seqids -title',
		    ' uprotdb -out uprotdb',sep='')
	system(cmd3,ignore.stdout=T)
	system('rm -rf all_genomes_uprot.faa')
	
 	# Extract sequences
 		 		
 	uprots<-colnames(presence_absence_uprot)
 	genoms<-rownames(presence_absence_uprot)
 		
 	for (u in 1:length(uprots)){
 			
 		outfile<-paste(uprots[u],'.uprot.faa',sep='')
		absents<-paste(uprots[u],'.absent',sep='')
 			
 		prots<-as.vector(presence_absence_uprot[,u])
			
		for (p in 1:length(prots)){
				
			prot<-prots[p]
				
			if (is.na(prot)==F){
					
				cmd<-paste(bldbcd,' -entry ',prot,' -db uprotdb >> ',outfile,sep='')
					
				system(cmd)
					
			} else {
				
				system(paste('touch',outfile))
				
				cat(paste(genoms[p],'_absent','\t',p,sep=''),sep='\n',file=absents,append=T)
			}
		}
	}
 		
	# Align sequences
	
	if (align==TRUE){
 		
 		multi<-list.files(pattern='.uprot.faa')
 		
 		for (m in multi){
			
			numse<-system(paste('grep -c ">" ',m,sep=''),intern=T)
			fsize<-file.info(m)$size
			
			if (numse > 1 & fsize > 0){
 			
 				out<-gsub('.faa','.ali',m)
 		
 				aux<-capture.output(
				alignment<-msa(inputSeqs=m,method='ClustalOmega',type='protein',order='input'))
				aliconver<-msaConvert(alignment,type='seqinr::alignment')
			
				seqs<-lapply(aliconver$seq,s2c)
				nams<-gsub(' ','',gsub('>','',system(paste("grep '>' ",m,sep=''),intern=T)))
				
				write.fasta(seqs,names=nams,file=out)
			
			} else if (numse == 1 & fsize > 0){
				
				cat(readLines(m),sep='\n',file=gsub('.faa','.ali',m))
				
			} else {
				
				out<-gsub('.faa','.ali',m)
				system(paste('touch',out))
			}
		}
		
		# Add absents
		
		abfiles<-list.files(pattern='.absent')
		
		if (length(abfiles)>0){
			
			for (a in abfiles){
				
				# Open corresponding alignment
				
				lfiles<-length(flist)
				lwc<-strsplit(system(paste('wc -l',a),intern=T),' ')[[1]][1]
				
				if (lfiles!=lwc){
				
					badfile<-gsub('.absent','.uprot.ali',a)
					fasta<-read.fasta(badfile)
					sequs<-lapply(getSequence(fasta),toupper)
					namos<-getName(fasta)
				
					gaps<-list()
					gaps[[1]]<-rep('-',length(sequs[[1]]))
				
					# Read file with absents
				
					ab<-read.table(a,sep='\t',header=F)
				
					counter<-0
				
					for (b in 1:dim(ab)[1]){
							
						nam<-as.vector(ab[b,1])
						pos<-as.vector(ab[b,2])+counter				
		
						counter<-counter+1
					
						sequs1<-sequs[1:pos]
						sequs2<-sequs[(pos+1):length(sequs)]

						namos1<-namos[1:pos]
						namos2<-namos[(pos+1):length(namos)]
					
						sequs<-c(sequs1,gaps,sequs2)
						namos<-c(namos1,nam,namos2)
					}
				
					write.fasta(sequs,names=namos,file=badfile)
				}
			}
		}

 		# Concatenate alignment
 		
 		alis<-list.files(pattern='.uprot.ali')
 		
 		catmat<-NULL
 		
 		for (a in alis){
			
			if (file.info(a)$size>0){ 
 			
				fasta<-read.fasta(a)
 				sequs<-lapply(getSequence(fasta),toupper)
 				smatx<-do.call(rbind,sequs)
 			
 				catmat<-cbind(catmat,smatx)
 				catlis<-alply(catmat,1)
 		
				nams<-gsub('>','',system(paste('grep ">" ,a',sep=''),intern=T))
 				#nams<-getName(fasta)
 				#nams<-unlist(lapply(nams,function(x){strsplit(x,'_')[[1]][1]}))
 			
 				write.fasta(catlis,names=nams,file='universal_proteins.ali')
			}
 		}
 	
		system(paste('mv *uprot.faa',outdir))
		system(paste('mv *uprot.ali',outdir))
		system(paste('mv universal_proteins.ali',outdir))
		system('rm -rf uprotdb*')
		system('rm -rf *.absent')
	
	} else {
		
		system(paste('mv uprot.faa',outdir))
		system('rm -rf uprotdb*')
	}
			
	if (align==TRUE & phylogeny!='NO'){
		
		setwd(outdir)
		
		if (phylogeny=='NJ'){
		
			upali<-read.fasta('universal_proteins.ali')
			asali<-as.alignment(upali)
			dista<-dist.alignment(asali)^2			
			tre<-nj(dista)

			write.tree(tre,file='NJ.uprot.tree.nwk')
					
		} else if (phylogeny=='ML'){
			
			cmd<-paste(fasttree,'universal_proteins.ali > ML.uprot.tree.nwk')
			
			system(cmd,ignore.stderr=T)
		}
	
	} else if (align==FALSE & phylogeny==TRUE){
		
		stop('For tree building both "align" and "phylogeny" must be set!')
	}
	
 	setwd(gw)
}
