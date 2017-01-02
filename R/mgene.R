#' mgene
#'
#' Identifies and extracts a desired gene used as taxonomic/phylogenetic marker.
#' @param path is the full path to where genome annotations in fasta format are placed.
#' @param pattern is a pattern (like '.faa') for recognizing the annotation files.
#' @param outdir is a name for the output directory that will be created for results.
#' @param marker is the multifasta file in amino acids for the desired marker gene.
#' @param align takes a logical (default, TRUE) indiciating if producing a multiple sequence alignment.
#' @param phylogeny takes a logical (default, TRUE) indicating if building NJ phylogeny.
#' @param proc is the number of threads used for protein searches.
#' @keywords marker gene
#' @export
#' @examples
#' mgene(path='./prodigal_output',pattern='.faa',marker='gyrB.faa',outdir='mgene_output',proc=2)

mgene<-function(pattern='.faa',
			    path='.',
			    outdir='mgene_output',
			    marker,
			    align=TRUE,
			    phylogeny=TRUE,
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
		hmmbuild<-paste(system.file('hmmer3',package='taxxo'),'/linux/hmmbuild',sep='')
		bldbcd<-paste(system.file('blast',package='taxxo'),'/linux/blastdbcmd',sep='')
		mkbldb<-paste(system.file('blast',package='taxxo'),'/linux/makeblastdb',sep='')
		
	} else if (os=='darwin'){
		
		hmmsearch<-paste(system.file('hmmer3',package='taxxo'),'/darwin/hmmsearch',sep='')
		hmmbuild<-paste(system.file('hmmer3',package='taxxo'),'/darwin/hmmbuild',sep='')
		bldbcd<-paste(system.file('blast',package='taxxo'),'/darwin/blastdbcmd',sep='')
		mkbldb<-paste(system.file('blast',package='taxxo'),'/darwin/makeblastdb',sep='')

	} else {
		
		stop('ERROR: Unknown OS, unable to proceed.')
	}
	
	umarker<-unlist(strsplit(marker,'.'))
	ulength<-length(umarker)
	nmarker<-paste(umarker[-ulength],collapse='_')
	aliname<-paste(nmarker,'mgene.ali',sep='.')
	hmmname<-paste(nmarker,'mgene.hmm',sep='.')

	# Align marker
		
	aux<-capture.output(
	alignment<-msa(inputSeqs=marker,method='ClustalOmega',type='protein'))
	aliconver<-msaConvert(alignment,type='seqinr::alignment')
			
	seqs<-lapply(aliconver$seq,s2c)
	nams<-gsub(' ','',gsub('>','',system(paste("grep '>' ",marker,sep=''),intern=T)))
		
	write.fasta(seqs,names=nams,file=aliname)
	
	# Build HMM
	
	cmd<-paste(hmmbuild,hmmname,aliname)
	
	system(cmd)
	
	flist<-list.files(path=path,pattern=pattern,full.names=T)
	fnams<-gsub(pattern,'',list.files(path=path,pattern=pattern))
	
	# Find proteins
		
	result<-matrix(ncol=1,nrow=length(flist),NA)
	colnames(result)<-nmarker
	rownames(result)<-fnams
	
	for (f in 1:length(flist)){
								
		cmd<-paste(hmmsearch,
		           ' --noali --cpu ',proc,
		           ' -o search_tmp',
		           ' --domtblout out_tmp',
		           ' ',hmmname,' ',flist[f],sep='')

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
  				
  			result[f,1]<-gene
  				
  			system('rm -rf search_tmp out_tmp')
		}
	}
		
  	system(paste('mkdir',outdir))
  	setwd(outdir)
  	
  	paname<-paste('presence_absence',nmarker,sep='_')

	assign(pname,result)
	  			
  	save(get(pname),file=paste(pname,'Rdata',sep='.'))				
 		
 	# Make databases
 		 		
	cmd2<-paste('cat ',paste(flist,collapse=' '),' > all_genomes_mgene.faa',sep='')
	system(cmd2)
	
	cmd3<-paste(mkbldb,
		    ' -in all_genomes_mgene.faa',
		    ' -dbtype prot -hash_index -parse_seqids -title',
		    ' mgenedb -out mgenedb',sep='')
	
	system(cmd3,ignore.stdout=T)
	
	system('rm -rf all_genomes_mgene.faa') 	
		 		
 	# Extract sequences
 		 		
 	prots<-as.vector(result[,1])
 	genoms<-rownames(result)
 			
	outfile<-paste(nmarker,'_all.faa',sep='')
 				
	for (p in 1:length(prots)){
				
		prot<-prots[p]
				
		if (is.na(prot)==F){
					
			cmd<-paste(bldbcd,' -entry ',prot,' -db mgenedb >> ',outfile,sep='')
					
			system(cmd)
					
		} else {
					
			cat(paste('>',genoms[p],'_absent',sep=''),
				file=outfile,
				sep='\n',
				append=T)
		}
	}
 		
	# Align sequences
 		
 	if (align==TRUE){	
 	
 		m<-list.files(pattern='_all.faa')
 				
 		out<-gsub('.faa','.ali',m)
 		
 		aux<-capture.output(
 		alignment<-msa(inputSeqs=m,method='ClustalOmega',type='protein'))
		aliconver<-msaConvert(alignment,type='seqinr::alignment')
			
		seqs<-lapply(aliconver$seq,s2c)
		nams<-gsub(' ','',gsub('>','',system(paste("grep '>' ",m,sep=''),intern=T)))
	
		write.fasta(seqs,names=nams,file=out)
		
		system(paste('mv *_all.faa',outdir))
		system(paste('mv *_all.ali',outdir))
		system(paste('mv *mgene.ali',outdir))
		system(paste('mv *mgene.hmm',outdir))
		system('rm -rf mgenedb*')

				
 	} else {
 		
 		system(paste('mv *_all.faa',outdir))
 		system(paste('mv *mgene.ali',outdir))
 		system(paste('mv *mgene.hmm',outdir))
 		system('rm -rf mgenedb*')
 	}
 	
 	# Phylogeny
 	
 	if (align==TRUE & phylogeny==TRUE){
		
		phydat<-msaConvert(alignment,type='phangorn::phyDat')
		distan<-dist.ml(phydat,model='JTT')
		tre<-NJ(distan)
		
		write.tree(tre,file=paste('NJ.',nmarker,'.tree.nwk',sep=''))
		
		system(paste('mv NJ.*.tree.nwk',outdir))
	
	} else if (align==FALSE & phylogeny==TRUE){
		
		stop('For tree building both "align" and "phylogeny" must be set TRUE')
	}
	
	setwd(gw)
 }
