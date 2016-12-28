#' rrna
#'
#' Outputs a individual rRNA gene sequences and a multiple sequence alignment.
#' @param path is the full path to where genomes in fasta format are placed.
#' @param pattern a pattern (like '.fasta') for recognizing the genome files.
#' @param outdir a name for the output directory that will be created for results.
#' @param subunit takes '16S' (default), '23S' or 'both'.
#' @param kingdom takes 'bacteria' (default) or 'archaea'.
#' @param align takes a logical (default, TRUE) indicating if sequence alignment is performed.
#' @keywords rRNA 16S 23S
#' @export
#' @examples
#' rrna(path='.',pattern='.fasta',outdir='rrna_output',subunit='16S',kingdom='bacteria')

rrna<-function(path,
			   pattern,
			   outdir='rrna_out',
			   subunit='16S',
			   kingdom='bacteria',
			   align=TRUE)
			   
			   {
	
	# Options #
	
	options(getClass.msg=FALSE)
					    	
	barrnap<-paste(system.file('barrnap',package='taxxo'),'/common/bin/barrnap',sep='')

	# Dependencies #
	
	require(seqinr,quietly=T)
	require(msa,quietly=T)
	
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
	system('touch rrna.log')
	
	# Perform rRNA search with barrnap #
	
	for (f in 1:length(flist)){
		
		outgff<-paste(gsub('.fasta','',fnams[f]),'.tmp.gff',sep='')
		
		cmd<-paste(barrnap,' --kingdom ',king,' ',flist[f],' > ',outgff,sep='')
		system(cmd,ignore.stderr=T)
		
		genome<-read.fasta(flist[f])
		sequen<-lapply(getSequence(genome),toupper)
		snames<-gsub('>','',system(paste("grep '>' ",flist[f],sep=''),intern=T))
		
		gff<-read.table(outgff,sep='\t',skip=1)
	
		sun<-as.vector(gff[,9])
		
		if (subunit=='16S'){

			grp<-grep(subunit,sun)
					
			if (length(grp)>0){
			
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

			} else {
				
				mssg<-paste('No ',subunit,' gene found in genome ',fnams[f],sep='')
				cat(mssg,file='rrna.log',append=T,sep='\n')
			}
			
		} else if (subunit=='23S'){
			
			grp<-grep(subunit,sun)
				
			if (length(grp)>0){
			
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
				
			} else {
				
				mssg<-paste('No ',subunit,' gene found in genome ',fnams[f],sep='')
				cat(mssg,file='rrna.log',append=T,sep='\n')
			}
			
		} else if (subunit=='both'){
			
			grps<-grep('16S',sun)
			grpl<-grep('23S',sun)
			
			if (length(grps)>0){
				
				inis<-gff[grps,4]
				fins<-gff[grps,5]
				stds<-gff[grps,7]
				nams<-as.vector(gff[grps,1])
				
				contigs<-which(snames==nams)
				
				if (stds=='+'){
				
					genes<-sequen[[contigs]][inis:fins]
				
				} else if (stds=='-'){
					
					genes<-toupper(rev(comp(sequen[[contigs]][inis:fins])))
				}
				
				write.fasta(genes,names=gsub('both.fasta','16S',outfiles[f]),file=gsub('both','16S',outfiles[f]))
				
				system(paste('mv *.16S.fasta',outdir))
				
			} else {
				
				mssg<-paste('No 16S gene found in genome ',fnams[f],sep='')
				cat(mssg,file='rrna.log',append=T,sep='\n')
			}
				
			if (length(grpl)>0){
				
				inil<-gff[grpl,4]
				finl<-gff[grpl,5]
				stdl<-gff[grpl,7]
				naml<-as.vector(gff[grpl,1])
				
				contigl<-which(snames==naml)
				
				if (stdl=='+'){
				
					genel<-sequen[[contigl]][inil:finl]
				
				} else if (stdl=='-'){
					
					genel<-toupper(rev(comp(sequen[[contigl]][inil:finl])))
				}
				
				write.fasta(genel,names=gsub('both.fasta','23S',outfiles[f]),file=gsub('both','23S',outfiles[f]))
				
				system(paste('mv *.23S.fasta',outdir))
				
			} else {
				
				mssg<-paste('No 23S gene found in genome ',fnams[f],sep='')
				cat(mssg,file='rrna.log',append=T,sep='\n')
			}
		}		
	}
	
	system('rm -rf *.tmp.gff')
	system(paste('mv rrna.log',outdir))
	
	# Align sequences #
	
	if (align==TRUE){
		
		setwd(outdir)
		
		if (subunit=='16S'){
			
			system('cat *.16S.fasta > all.16S.fasta')
			namali<-gsub('>','',system("grep '>' all.16S.fasta",intern=T))
			
			alignment<-msa(inputSeqs='all.16S.fasta',method='Muscle',type='dna')	
			aliconver<-msaConvert(alignment,type='seqinr::alignment')
			
			seqs<-lapply(aliconver$seq,s2c)
			
			write.fasta(seqs,names=namali,file='alignment.16S.fasta')
			
		} else if (subunit=='23S'){
			
			system('cat *.23S.fasta > all.23S.fasta')
			namali<-gsub('>','',system("grep '>' all.23S.fasta",intern=T))
			
			alignment<-msa(inputSeqs='all.23S.fasta',method='Muscle',type='dna')	
			aliconver<-msaConvert(alignment,type='seqinr::alignment')
			
			seqs<-lapply(aliconver$seq,s2c)
			
			write.fasta(seqs,names=namali,file='alignment.23S.fasta')
			
		} else if (subunit=='both'){
			
			system('cat *.16S.fasta > all.16S.fasta')
			namalis<-gsub('>','',system("grep '>' all.16S.fasta",intern=T))
			
			alignments<-msa(inputSeqs='all.16S.fasta',method='Muscle',type='dna')	
			aliconvers<-msaConvert(alignments,type='seqinr::alignment')
			
			seqs<-lapply(aliconvers$seq,s2c)
			
			write.fasta(seqs,names=namalis,file='alignment.16S.fasta')

			system('cat *.23S.fasta > all.23S.fasta')
			namalil<-gsub('>','',system("grep '>' all.23S.fasta",intern=T))
			
			alignmentl<-msa(inputSeqs='all.23S.fasta',method='Muscle',type='dna')	
			aliconverl<-msaConvert(alignmentl,type='seqinr::alignment')
			
			seql<-lapply(aliconverl$seq,s2c)
			
			write.fasta(seql,names=namalil,file='alignment.23S.fasta')
		}
	}
	
	setwd('../')
}