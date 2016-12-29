#' rrna
#'
#' Outputs a individual rRNA gene sequences and a multiple sequence alignment.
#' @param path is the full path to where genomes in fasta format are placed.
#' @param pattern a pattern (like '.fasta') for recognizing the genome files.
#' @param outdir a name for the output directory that will be created for results.
#' @param subunit takes '16S' (default) or '23S'.
#' @param multiple takes a logical (default, FALSE) for managing multiple copies of rRNA genes. If more than one copy is found, FALSE will return one random copy while TRUE will return all of them.
#' @param distance takes a logical (default, TRUE) indicating if a distance matrix is calculated and outputed or not.
#' @param phylogeny takes a logical (default, TRUE) indicating if a Neighbor-Joining phylogeny is constructed or not.
#' @param kingdom takes 'bacteria' (default) or 'archaea'.
#' @param align takes a logical (default, TRUE) indicating if sequence alignment is performed.
#' @keywords rRNA 16S 23S
#' @export
#' @examples
#' rrna(align=T,kingdom='bacteria',multiple=F,outdir='./test_genomes/rrna_result',path='./test_genomes/',pattern='.fna',subunit='16S',phylogeny=TRUE)


rrna<-function(path,
			   pattern,
			   outdir='rrna_out',
			   subunit='16S',
			   kingdom='bacteria',
			   align=TRUE,
			   multiple=FALSE,
			   distance=TRUE,
			   phylogeny=TRUE,
			   ...)
			   
			   {
	
	# Options #
	
	options(getClass.msg=FALSE)
	gw<-getwd()
	barrnap<-paste(system.file('barrnap',package='taxxo'),'/common/bin/barrnap',sep='')

	# Dependencies #
	
	require(seqinr,quietly=T)
	require(msa,quietly=T)
	require(phangorn,quietly=T)
	
	# Check input #
	
	if (kingdom=='bacteria'){
		
		king<-'bac'
	
	} else if (kingdom=='archaea'){
		
		king<-'arc'
		
	} else {
		
		stop("ERROR: Parameter kingdom must be 'archaea' or 'bacteria'")
	}
	
	if (!subunit%in%c('16S','23S','both')){
		
		stop("ERROR: Parameter subunit must be '16S', '23S' or 'both'")
	}
	
	# Get genomes #
	
	flist<-list.files(path=path,pattern=pattern,full.names=T)
	fnams<-list.files(path=path,pattern=pattern)
	
	outfiles<-paste(gsub(pattern,'',fnams),'.',subunit,'.fasta',sep='')
		
	# Create output directory #
	
	system(paste('mkdir',outdir))
	system('touch rrna.err')
	
	# Perform rRNA search with barrnap #
	
	for (f in 1:length(flist)){
		
		outgff<-paste(gsub('.fasta','',fnams[f]),'.tmp.gff',sep='')
		
		cmd<-paste(barrnap,' --kingdom ',king,' ',flist[f],' > ',outgff,sep='')
		system(cmd,ignore.stderr=T)
		
		genome<-read.fasta(flist[f])
		sequen<-lapply(getSequence(genome),toupper)
		snames2<-gsub('>','',system(paste("grep '>' ",flist[f],sep=''),intern=T))
		snames<-unlist(lapply(snames2,function(x){strsplit(x,' ')[[1]][1]}))
		
		if (file.info(outgff)$size==0){
			
			mssg<-paste('No ',subunit,' gene found in genome ',fnams[f],sep='')
			cat(mssg,file='rrna.err',append=T,sep='\n')
			
			stop(mssg)
			
		} else {

			gff<-read.table(outgff,sep='\t',skip=1)
	
			sun<-as.vector(gff[,9])
		
			if (subunit=='16S'){

				grp<-grep(subunit,sun)
					
				if (length(grp)==0){
					
					mssg<-paste('No ',subunit,' gene found in genome ',fnams[f],sep='')
					cat(mssg,file='rrna.err',append=T,sep='\n')
			
					stop(mssg)
					
				} else if (length(grp)==1){
						
					ini<-gff[grp,4]
					fin<-gff[grp,5]
					std<-gff[grp,7]
					nam<-as.vector(gff[grp,1])
			
					contig<-which(snames==nam)
				
					if (std=='+'){
				
						gene<-sequen[[contig]][ini:fin]
				
					} else if (std=='-'){
					
						gene<-toupper(rev(comp(sequen[[contig]][ini:fin])))
					}
				
					write.fasta(gene,names=gsub('.fasta','',outfiles[f]),file=outfiles[f])

					system(paste('mv *.16S.fasta',outdir))

				} else if (length(grp)>1){
				
					if (multiple==F){
						
						ini<-gff[grp[1],4]
						fin<-gff[grp[1],5]
						std<-gff[grp[1],7]
						nam<-as.vector(gff[grp[1],1])
				
						contig<-which(snames==nam)
				
						if (std=='+'){
				
							gene<-sequen[[contig]][ini:fin]
				
						} else if (std=='-'){
					
							gene<-toupper(rev(comp(sequen[[contig]][ini:fin])))
						}
				
						write.fasta(gene,names=gsub('.fasta','',outfiles[f]),file=outfiles[f])

						system(paste('mv *.16S.fasta',outdir))
						
					} else if (multiple==T){
						
						gff2<-gff[grp,]
						
						ini<-gff2[,4]
						fin<-gff2[,5]
						std<-as.vector(gff2[,7])
						nam<-as.vector(gff2[,1])
						
						for (g in 1:length(grp)){
							
							contig<-which(snames==nam[g])
							
							if (std[g]=='+'){
				
								gene<-sequen[[contig]][ini[g]:fin[g]]
				
							} else {
					
								gene<-toupper(rev(comp(sequen[[contig]][ini[g]:fin[g]])))
							}
							
							write.fasta(gene,names=gsub('fasta',g,outfiles[f]),file=outfiles[f],open='a')
						}

						system(paste('mv *.16S.fasta',outdir))
					}
				}
				
			} else if (subunit=='23S'){
			
				if (length(grp)==0){
					
					mssg<-paste('No ',subunit,' gene found in genome ',fnams[f],sep='')
					cat(mssg,file='rrna.err',append=T,sep='\n')
			
					stop(mssg)
					
				} else if (length(grp)==1){
						
					ini<-gff[grp,4]
					fin<-gff[grp,5]
					std<-gff[grp,7]
					nam<-as.vector(gff[grp,1])
			
					contig<-which(snames==nam)

					if (std=='+'){
				
						gene<-sequen[[contig]][ini:fin]
				
					} else if (std=='-'){
					
						gene<-toupper(rev(comp(sequen[[contig]][ini:fin])))
					}
				
					write.fasta(gene,names=gsub('.fasta','',outfiles[f]),file=outfiles[f])

					system(paste('mv *.23S.fasta',outdir))

				} else if (length(grp)>1){
				
					if (multiple==F){
						
						ini<-gff[grp[1],4]
						fin<-gff[grp[1],5]
						std<-gff[grp[1],7]
						nam<-as.vector(gff[grp[1],1])
							
						contig<-which(snames==nam)
								
						if (std=='+'){
							
							gene<-sequen[[contig]][ini:fin]
				
						} else if (std=='-'){
					
							gene<-toupper(rev(comp(sequen[[contig]][ini:fin])))
						}
				
						write.fasta(gene,names=gsub('.fasta','',outfiles[f]),file=outfiles[f])

						system(paste('mv *.23S.fasta',outdir))
						
					} else if (multiple==T){
						
						gff2<-gff[grp,]

						ini<-gff2[grp,4]
						fin<-gff2[grp,5]
						std<-as.vector(gff2[grp,7])
						nam<-as.vector(gff2[grp,1])
						
						for (g in 1:length(grp)){
							
							contig<-which(snames==nam[g])

							if (std[g]=='+'){
				
								gene<-sequen[[contig]][ini[g]:fin[g]]
				
							} else {
					
								gene<-toupper(rev(comp(sequen[[contig]][ini[g]:fin[g]])))
							}
							
							write.fasta(gene,names=gsub('fasta',g,outfiles[f]),file=outfiles[f],open='a')
						}

						system(paste('mv *.23S.fasta',outdir))
					}
				}
				
			} else {
				
				stop('Parameter subunit must be "16S" or "23S"')
			}
		}
	}	
					
	system('rm -rf *.tmp.gff')
	system(paste('mv rrna.err',outdir))
	
	# Align sequences #
	
	if (align==TRUE){
		
		setwd(outdir)
		
		if (subunit=='16S'){
			
			system('cat *.16S.fasta > all.16S.fasta')
			namali<-gsub('>','',system("grep '>' all.16S.fasta",intern=T))
			
			aux<-capture.output(
			alignment<-msa(inputSeqs='all.16S.fasta',method='ClustalOmega',type='dna'))
			aliconver<-msaConvert(alignment,type='seqinr::alignment')
			
			seqs<-lapply(aliconver$seq,s2c)
			
			write.fasta(seqs,names=namali,file='alignment.16S.fasta')
			
		} else if (subunit=='23S'){
			
			system('cat *.23S.fasta > all.23S.fasta')
			namali<-gsub('>','',system("grep '>' all.23S.fasta",intern=T))
			
			aux<-capture.output(
			alignment<-msa(inputSeqs='all.23S.fasta',method='ClustalO',type='dna'))
			aliconver<-msaConvert(alignment,type='seqinr::alignment')
			
			seqs<-lapply(aliconver$seq,s2c)
			
			write.fasta(seqs,names=namali,file='alignment.23S.fasta')
		}			
	}
	
	# Distance matrix #
	
	if (distance==TRUE){
		
		aliconver$nam<-namali
		
		dista1<-dist.alignment(aliconver)
		dista2<-round(dista1^2,digits=3)
		
		dnam<-paste('distance_matrix_',subunit,sep='')
		
		assign(dnam,dista2)
		
		save(list=dnam,file=paste(dnam,'Rdata',sep='.'))
	}
	
	# Phylogeny #
	
	if (phylogeny==TRUE){
		
		phydat<-msaConvert(alignment,'phangorn::phyDat')
		distma<-dist.ml(phydat,model='F81')
		tre<-NJ(distma)
		write.tree(tre,file=paste('NJ.',subunit,'.tree.nwk',sep=''))		
	}
	
	setwd(gw)
}