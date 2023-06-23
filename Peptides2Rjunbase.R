# Libraries
library(GenomicAlignments)
library(dplyr)
library(chromoMap)

# Data
linear<-read.delim("C:/Users/bea_f/OneDrive/Escritorio/HematoLaFe/Neoantigens-RNA-splicing/RJunBase/detail_LS_annotation.txt")
non.can.bam<-readGAlignments("Peptides/non-canonical.sorted.bam",use.names = TRUE)
junctions<-read.delim("RJunBase/RJunbase_annotated.txt",sep="\t")

######################## BED file: do not execute again ########################
chr.order<-c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY","chrM")

bed.file <- file("./RJunBase/RJunbase.txt", "wb")
df<-data.frame(chrom=sub(":.*","",linear$Junction.location),
               start=sub(".*: *(.*?) *\\|.*", "\\1", linear$Junction.location),
               end=sub(".*\\| *(.*?) *:.*", "\\1", linear$Junction.location),
               strand=sub(".*:","",linear$Junction.location))
df<-df[order(as.numeric(df$start)),]
df<-df[order(match(df$chrom,chr.order)),]
write.table(x=df,file=bed.file,row.names = FALSE,quote=FALSE,col.names=TRUE)
close(bed.file)


bed.annotated.file <- file("./RJunBase/RJunbase_annotated.txt", "wb")
df<-data.frame(chrom=sub(":.*","",linear$Junction.location),
               start=sub(".*: *(.*?) *\\|.*", "\\1", linear$Junction.location),
               end=sub(".*\\| *(.*?) *:.*", "\\1", linear$Junction.location),
               strand=sub(".*:","",linear$Junction.location))
df<-cbind(df,linear)
df<-df[order(as.numeric(df$start)),]
df<-df[order(match(df$chrom,chr.order)),]
write.table(x=df,file=bed.annotated.file,row.names = FALSE,quote=FALSE,col.names=TRUE,sep="\t")
close(bed.annotated.file)


########################### Overlapping regions ################################
junctions$chrom[junctions$chrom=="chrM"]<-"chrMT"

##### Junction frequency in tumor and normal tissues #####
pattern <- "\\(\\s*(.*?)\\s*\\)"
nt.freq<- regmatches(junctions$NT.frequency.sample.number., regexec(pattern, junctions$NT.frequency.sample.number.))
nt.freq.indiv<-unlist(lapply(nt.freq,function(x){x[[2]]}))

normal.freq<- regmatches(junctions$Normal.frequency.sample.number., regexec(pattern, junctions$Normal.frequency.sample.number.))
normal.freq.indiv<-unlist(lapply(normal.freq,function(x){x[[2]]}))

tumor.freq<- regmatches(junctions$Tumor.frequency.sample.number., regexec(pattern, junctions$Tumor.frequency.sample.number.))
tumor.freq.indiv<-unlist(lapply(tumor.freq,function(x){x[[2]]}))

##### Convert data to genomic ranges #####
gr<-as.data.frame(non.can.bam) 
gr<-makeGRangesFromDataFrame(gr)

gr.UE<-GRanges(seqnames=junctions$chrom,
               ranges = IRanges(start=junctions$start-36,end=junctions$start),
               strand=junctions$strand)
gr.DE<-GRanges(seqnames=junctions$chrom,
               ranges = IRanges(start=junctions$end,end=junctions$end+36),
               strand=junctions$strand)

##### Find overlaps #####
result.UE<-findOverlaps(gr,gr.UE,type="within")
result.DE<-findOverlaps(gr,gr.DE,type="within")

## Write output: do not execute again
gr.df<-as.data.frame(gr)
colnames(gr.df)<-paste0("peptide_",colnames(gr.df))

write.table(x=cbind(gr.df[queryHits(result.UE),],junctions[subjectHits(result.UE),]),
            file="peptidesUE.txt",quote=FALSE)

write.table(x=cbind(gr.df[queryHits(result.DE),],junctions[subjectHits(result.DE),]),
            file="peptidesDE.txt",quote=FALSE)

##### Analysis #####
queries<-c(queryHits(result.UE),queryHits(result.DE)) # All queries with hits
unique.queries<-queries[!duplicated(queries) & !duplicated(queries,fromLast=TRUE)] # 2186 queries with only 1 hit
UE.DE.queries<-unique(queryHits(result.UE)[queryHits(result.UE) %in% queryHits(result.DE)]) # 103 queries with hit in UE and DE

UE.unique<-result.UE[queryHits(result.UE) %in% unique.queries] # 931
DE.unique<-result.DE[queryHits(result.DE) %in% unique.queries] # 1255

UE.duplicated<-result.UE[!queryHits(result.UE) %in% unique.queries & 
                           !queryHits(result.UE) %in% UE.DE.queries] # 1014 unique
DE.duplicated<-result.DE[!queryHits(result.DE) %in% unique.queries & 
                           !queryHits(result.DE) %in% UE.DE.queries] # 1351
UE.common<-result.UE[queryHits(result.UE) %in% UE.DE.queries]
DE.common<-result.DE[queryHits(result.DE) %in% UE.DE.queries]

View(cbind(as.data.frame(UE.unique),junctions[subjectHits(UE.unique),]))
View(cbind(as.data.frame(DE.unique),junctions[subjectHits(DE.unique),]))

View(cbind(as.data.frame(UE.duplicated),junctions[subjectHits(UE.duplicated),]))
View(cbind(as.data.frame(DE.duplicated),junctions[subjectHits(DE.duplicated),]))

View(cbind(as.data.frame(UE.common),junctions[subjectHits(UE.common),]))
View(cbind(as.data.frame(DE.common),junctions[subjectHits(DE.common),]))


## Search for more tumor-specific peptides (normal frequency=0 independent of tumor frequency)
# Queries qith unique hits
sum(nt.freq.indiv[subjectHits(UE.unique)]==0 & normal.freq.indiv[subjectHits(UE.unique)]==0)
sum(nt.freq.indiv[subjectHits(DE.unique)]==0 & normal.freq.indiv[subjectHits(DE.unique)]==0)

# Queries with common hits in UE and DE from different junctions
UE.DE.unique<-vector()
for (i in 1:length(UE.DE.queries)){
  ue<-grep(UE.DE.queries[i],queryHits(UE.common))
  de<-grep(UE.DE.queries[i],queryHits(DE.common))
  if(all(c(normal.freq.indiv[subjectHits(UE.common)[ue]]==0 &
           nt.freq.indiv[subjectHits(UE.common)[ue]]==0 &
           normal.freq.indiv[subjectHits(DE.common)[de]]==0 &
           nt.freq.indiv[subjectHits(DE.common)[de]]==0))){
    UE.DE.unique[i]<-UE.DE.queries[i]
  }
  print(i)
}

UE.DE.unique[!is.na(UE.DE.unique)] # 4 (all non tumor-specific)

## Duplicated queries with all hits in UE having freq 0 in normal tissue
UE.duplicated.unique<-vector()
for (i in 1:(length(unique(queryHits(UE.duplicated))))){
  index<-subjectHits(UE.duplicated)[queryHits(UE.duplicated)==unique(queryHits(UE.duplicated))[i]]
  if (all(c(normal.freq.indiv[index]==0 &
            nt.freq.indiv[index]==0))){
    UE.duplicated.unique[i]<-unique(queryHits(UE.duplicated))[i]
  }
}

UE.duplicated.unique[!is.na(UE.duplicated.unique)] # 1 (non tumor-specific)
View(cbind(as.data.frame(UE.duplicated[queryHits(UE.duplicated) %in% UE.duplicated.unique[!is.na(UE.duplicated.unique)]]),
           junctions[subjectHits(UE.duplicated[queryHits(UE.duplicated) %in% UE.duplicated.unique[!is.na(UE.duplicated.unique)]]),]))

## Duplicated queries with all hits in DE having freq 0 in normal tissue
DE.duplicated.unique<-vector()
for (i in 1:(length(unique(queryHits(DE.duplicated))))){
  index<-subjectHits(DE.duplicated)[queryHits(DE.duplicated)==unique(queryHits(DE.duplicated))[i]]
  if (all(c(normal.freq.indiv[index]==0 &
            nt.freq.indiv[index]==0))){
    DE.duplicated.unique[i]<-unique(queryHits(DE.duplicated))[i]
  }
}

DE.duplicated.unique[!is.na(DE.duplicated.unique)] # 4 (all non tumor-specific)
View(cbind(as.data.frame(DE.duplicated[queryHits(DE.duplicated) %in% DE.duplicated.unique[!is.na(DE.duplicated.unique)]]),
           junctions[subjectHits(DE.duplicated[queryHits(DE.duplicated) %in% DE.duplicated.unique[!is.na(DE.duplicated.unique)]]),]))


##### Table following Rjunbase criteria #####
# For duplicated peptides all junctions must be tumor specific to classify peptides as tumor specific
UE.duplicated.unique<-vector()
for (i in 1:(length(unique(queryHits(UE.duplicated))))){
  index<-subjectHits(UE.duplicated)[queryHits(UE.duplicated)==unique(queryHits(UE.duplicated))[i]]
  if (all(junctions$Tumor.specific.types[index]!="non tumor-specific")){
    UE.duplicated.unique[i]<-index
  }
}


DE.duplicated.unique<-vector()
for (i in 1:(length(unique(queryHits(DE.duplicated))))){
  index<-subjectHits(DE.duplicated)[queryHits(DE.duplicated)==unique(queryHits(DE.duplicated))[i]]
  if (all(junctions$Tumor.specific.types[index]!="non tumor-specific")){
    DE.duplicated.unique[i]<-index
  }
}

UE.DE.unique<-vector()
DE.UE.unique<-vector()
for (i in 1:length(UE.DE.queries)){
  ue<-grep(UE.DE.queries[i],queryHits(UE.common))
  de<-grep(UE.DE.queries[i],queryHits(DE.common))
  if(all(c(junctions$Tumor.specific.types[subjectHits(UE.common)[ue]]!="non tumor-specific",
           junctions$Tumor.specific.types[subjectHits(DE.common)[de]]!="non tumor-specific"))){
    UE.DE.unique[i]<-subjectHits(UE.common)[ue]
    DE.UE.unique[i]<-subjectHits(DE.common)[de]
  }
  print(i)
}


df<-rbind(data.frame(peptide=names(gr[queryHits(UE.unique)]),
               n_hits="1",
               exon_hit="UE",
               tumor_specific=ifelse(junctions$Tumor.specific.types[subjectHits(UE.unique)]=="non tumor-specific","no","yes"),
               chromosome=junctions$chrom[subjectHits(UE.unique)],
               start=junctions$start[subjectHits(UE.unique)],
               end=junctions$end[subjectHits(UE.unique)],
               strand=junctions$strand[subjectHits(UE.unique)],
               origin=junctions$Alternative.splicing.types[subjectHits(UE.unique)],
               tumor_type=junctions$Tumor.specific.types[subjectHits(UE.unique)]),
      data.frame(peptide=names(gr[queryHits(DE.unique)]),
                 n_hits="1",
                 exon_hit="DE",
                 tumor_specific=ifelse(junctions$Tumor.specific.types[subjectHits(DE.unique)]=="non tumor-specific","no","yes"),
                 chromosome=junctions$chrom[subjectHits(DE.unique)],
                 start=junctions$start[subjectHits(DE.unique)],
                 end=junctions$end[subjectHits(DE.unique)],
                 strand=junctions$strand[subjectHits(DE.unique)],
                 origin=junctions$Alternative.splicing.types[subjectHits(DE.unique)],
                 tumor_type=junctions$Tumor.specific.types[subjectHits(DE.unique)]),
      data.frame(peptide=names(gr[queryHits(UE.duplicated)]),
                 n_hits=">1",
                 exon_hit="UE",
                 tumor_specific=ifelse(queryHits(UE.duplicated) %in% UE.duplicated.unique,"yes","no"),
                 chromosome=junctions$chrom[subjectHits(UE.duplicated)],
                 start=junctions$start[subjectHits(UE.duplicated)],
                 end=junctions$end[subjectHits(UE.duplicated)],
                 strand=junctions$strand[subjectHits(UE.duplicated)],
                 origin=junctions$Alternative.splicing.types[subjectHits(UE.duplicated)],
                 tumor_type=junctions$Tumor.specific.types[subjectHits(UE.duplicated)]),
      data.frame(peptide=names(gr[queryHits(DE.duplicated)]),
                 n_hits=">1",
                 exon_hit="DE",
                 tumor_specific=ifelse(queryHits(DE.duplicated) %in% DE.duplicated.unique,"yes","no"),
                 chromosome=junctions$chrom[subjectHits(DE.duplicated)],
                 start=junctions$start[subjectHits(DE.duplicated)],
                 end=junctions$end[subjectHits(DE.duplicated)],
                 strand=junctions$strand[subjectHits(DE.duplicated)],
                 origin=junctions$Alternative.splicing.types[subjectHits(DE.duplicated)],
                 tumor_type=junctions$Tumor.specific.types[subjectHits(DE.duplicated)]), 
      data.frame(peptide=names(gr[queryHits(UE.common)]),
                 n_hits=">1",
                 exon_hit="UE;DE",
                 tumor_specific=ifelse(queryHits(UE.common) %in% UE.DE.unique,"yes","no"),
                 chromosome=junctions$chrom[subjectHits(UE.common)],
                 start=junctions$start[subjectHits(UE.common)],
                 end=junctions$end[subjectHits(UE.common)],
                 strand=junctions$strand[subjectHits(UE.common)],
                 origin=junctions$Alternative.splicing.types[subjectHits(UE.common)],
                 tumor_type=junctions$Tumor.specific.types[subjectHits(UE.common)]),
      data.frame(peptide=names(gr[queryHits(DE.common)]),
                 n_hits=">1",
                 exon_hit="UE;DE",
                 tumor_specific=ifelse(queryHits(DE.common) %in% DE.UE.unique,"yes","no"),
                 chromosome=junctions$chrom[subjectHits(DE.common)],
                 start=junctions$start[subjectHits(DE.common)],
                 end=junctions$end[subjectHits(DE.common)],
                 strand=junctions$strand[subjectHits(DE.common)],
                 origin=junctions$Alternative.splicing.types[subjectHits(DE.common)],
                 tumor_type=junctions$Tumor.specific.types[subjectHits(DE.common)])
)

length(unique(df$peptide[df$tumor_specific=="yes"])) # 15

pep.file <- file("./RJunBase/Peptides2RjunbaseC1.txt", "wb")
write.table(x=df,file=pep.file,row.names = FALSE,quote=FALSE,col.names=TRUE)
close(pep.file)

##### Table following my criteria: normal frequency =0 #####
UE.DE.unique<-vector()
for (i in 1:length(UE.DE.queries)){
  ue<-grep(UE.DE.queries[i],queryHits(UE.common))
  de<-grep(UE.DE.queries[i],queryHits(DE.common))
  if(all(c(normal.freq.indiv[subjectHits(UE.common)[ue]]==0 &
           nt.freq.indiv[subjectHits(UE.common)[ue]]==0 &
           normal.freq.indiv[subjectHits(DE.common)[de]]==0 &
           nt.freq.indiv[subjectHits(DE.common)[de]]==0))){
    UE.DE.unique[i]<-UE.DE.queries[i]
  }
  print(i)
}

DE.duplicated.unique<-vector()
for (i in 1:(length(unique(queryHits(DE.duplicated))))){
  index<-subjectHits(DE.duplicated)[queryHits(DE.duplicated)==unique(queryHits(DE.duplicated))[i]]
  if (all(junctions$Tumor.specific.types[index]!="non tumor-specific")){
    DE.duplicated.unique[i]<-index
  }
}

UE.DE.unique<-vector()
DE.UE.unique<-vector()
for (i in 1:length(UE.DE.queries)){
  ue<-grep(UE.DE.queries[i],queryHits(UE.common))
  de<-grep(UE.DE.queries[i],queryHits(DE.common))
  if(all(c(junctions$Tumor.specific.types[subjectHits(UE.common)[ue]]!="non tumor-specific",
           junctions$Tumor.specific.types[subjectHits(DE.common)[de]]!="non tumor-specific"))){
    UE.DE.unique[i]<-subjectHits(UE.common)[ue]
    DE.UE.unique[i]<-subjectHits(DE.common)[de]
  }
  print(i)
}

df$tumor_specific[which(c(nt.freq.indiv[subjectHits(UE.unique)]==0 & normal.freq.indiv[subjectHits(UE.unique)]==0,
                    nt.freq.indiv[subjectHits(DE.unique)]==0 & normal.freq.indiv[subjectHits(DE.unique)]==0))]<-"yes"

df$tumor_specific[df$peptide %in% names(gr)[UE.DE.unique[!is.na(UE.DE.unique)]]]<-"yes"
df$tumor_specific[df$peptide %in% names(gr)[UE.duplicated.unique[!is.na(UE.duplicated.unique)]]]<-"yes"
df$tumor_specific[df$peptide %in% names(gr)[DE.duplicated.unique[!is.na(DE.duplicated.unique)]]]<-"yes"

length(unique(df$peptide[df$tumor_specific=="yes"]))

pep.file <- file("./RJunBase/Peptides2RjunbaseC2.txt", "wb")
write.table(x=df,file=pep.file,row.names = FALSE,quote=FALSE,col.names=TRUE)
close(pep.file)

##### Barplots #####
# https://cran.r-project.org/web/packages/chromoMap/vignettes/chromoMap.html
# All overlaps by chromosome and by start coordinate
# Frequency as percentage. Tumor: 10,283 samples. Normal (GTEx): 17,382. NT:¿?
# Also we have median expression
# chromosome coverage¿?
# position must relate to UE or DE or both

index<-unique(c(subjectHits(result.UE),subjectHits(result.DE))) # 2006 unique junctions
pos<-vector()
for(i in 1:length(index)){
  if(index[i] %in% subjectHits(UE.unique) | index[i] %in% subjectHits(UE.duplicated) | index[i] %in% subjectHits(UE.common)){
    pos[i]<-as.numeric(sub("\\(.*","",junctions$start[index[i]]))
  }
  if(index[i] %in% subjectHits(DE.unique) | index[i] %in% subjectHits(DE.duplicated) | index[i] %in% subjectHits(DE.common)){
    pos[i]<-as.numeric(sub("\\(.*","",junctions$end[index[i]]))
  }
}


y.normal<-as.numeric(sub("\\(.*","",junctions$Normal.frequency.sample.number.[which(junctions$chrom=="chr1")[which(junctions$chrom=="chr1") %in% index]]))
y.tumor<-as.numeric(sub("\\(.*","",junctions$Tumor.frequency.sample.number.[which(junctions$chrom=="chr1")[which(junctions$chrom=="chr1") %in% index]]))

pos<-as.numeric(sub("\\(.*","",junctions$start[which(junctions$chrom=="chr1")[which(junctions$chrom=="chr1") %in% index]]))
chr1.annot<-cbind(paste0("junc",1:length(y.normal)),rep("chr1",length(y.normal)),pos,pos+1,y.normal)
chr1.annot.tumor<-cbind(paste0("junc",1:length(y.tumor)),rep("chr1",length(y.tumor)),pos,pos+1,y.tumor)

chr_file_1<-cbind("chr1",1,248956422)
write.table(x=chr1.annot,file="chr1_annot.txt",quote=FALSE,row.names = FALSE,col.names = FALSE,sep="\t")
write.table(x=chr1.annot.tumor,file="chr1_annot_tumor.txt",quote=FALSE,row.names = FALSE,col.names = FALSE,sep="\t")
write.table(x=chr_file_1,file="chr_file_1.txt",quote=FALSE,row.names = FALSE,col.names = FALSE,sep="\t")
par(mfrow=c(2,2))
chromoMap("chr_file_1.txt",
          "chr1_annot.txt",
          data_based_color_map = TRUE,
          data_type="numeric",
          plots="bar",
          plot_color="blue")
chromoMap("chr_file_1.txt",
          "chr1_annot_tumor.txt",
          data_based_color_map = TRUE,
          data_type="numeric",
          plots="bar",
          plot_color="red")

plot.junctions<-function(junctions,index,chr){
  y.normal<-as.numeric(sub("\\(.*","",junctions$Normal.frequency.sample.number.[which(junctions$chrom==chr)[which(junctions$chrom==chr) %in% index]]))
  y.nt<-as.numeric(sub("\\(.*","",junctions$NT.frequency.sample.number.[which(junctions$chrom==chr)[which(junctions$chrom==chr) %in% index]]))
  y.tumor<-as.numeric(sub("\\(.*","",junctions$Tumor.frequency.sample.number.[which(junctions$chrom==chr)[which(junctions$chrom==chr) %in% index]]))
  barplot(matrix(c(y.normal,y.nt,y.tumor),nrow=3),col=c("blue","green","red"),beside=FALSE,border=c("blue","green","red"))
  print(length(y.normal))
}

plot.junctions(junctions,index,"chr2")


plot.junctions<-function(junctions,index,chr){
  y.normal<-as.numeric(sub("\\(.*","",junctions$Normal.frequency.sample.number.[which(junctions$chrom==chr)[which(junctions$chrom==chr) %in% index]]))
  y.nt<-as.numeric(sub("\\(.*","",junctions$NT.frequency.sample.number.[which(junctions$chrom==chr)[which(junctions$chrom==chr) %in% index]]))
  y.tumor<-as.numeric(sub("\\(.*","",junctions$Tumor.frequency.sample.number.[which(junctions$chrom==chr)[which(junctions$chrom==chr) %in% index]]))
  df<-data.frame(pos,y.normal,y.nt,y.tumor)
  ggplot(data=df,aes(x=pos,y=y.normal))+geom_area(aes(fill="blue"))
  print(df)
  print(length(y.normal))
}

plot.junctions<-function(junctions,index,chr){
  pos<-as.numeric(sub("\\(.*","",junctions$start[which(junctions$chrom==chr)[which(junctions$chrom==chr) %in% index]]))
  y.normal<-as.numeric(sub("\\(.*","",junctions$Normal.frequency.sample.number.[which(junctions$chrom==chr)[which(junctions$chrom==chr) %in% index]]))
  y.nt<-as.numeric(sub("\\(.*","",junctions$NT.frequency.sample.number.[which(junctions$chrom==chr)[which(junctions$chrom==chr) %in% index]]))
  y.tumor<-as.numeric(sub("\\(.*","",junctions$Tumor.frequency.sample.number.[which(junctions$chrom==chr)[which(junctions$chrom==chr) %in% index]]))
  gr.normal<-GRanges(seqnames=chr,ranges=IRanges())
  ggplot(data=df,aes(x=pos,y=y.normal))+geom_area(aes(fill="blue"))
  print(df)
  print(length(y.normal))
}


gr.normal<-GRanges(seqnames="chr1",ranges=IRanges(start=pos,end=pos))
kp <- plotKaryotype(plot.type=2, chromosomes = "chr1",genome="hg38",
                    zoom=GRanges(seqnames="chr1",ranges=IRanges(start=0,end=7000000)))
kpPlotCoverage(kp, data=gr.normal[1:5])


pos.chr<-c(1:248956422)
pos.chr<-pos.chr[!pos.chr %in% pos]
pos.final<-c(pos,pos.chr)
length(pos.chr)+length(unique(pos))
y<-c(y.normal,rep(0,length(pos.chr)))
plot(pos.final,y)
barplot(y[order(pos.final,decreasing = FALSE)])
