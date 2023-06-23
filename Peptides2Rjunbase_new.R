# Libraries
library(GenomicAlignments)
library(dplyr)
library(chromoMap)
library(ggplot2)
library(forcats)

# Data
linear<-read.delim("C:/Users/bea_f/OneDrive/Escritorio/HematoLaFe/Neoantigens-RNA-splicing/RJunBase/detail_LS_annotation.txt")
non.can.bam<-readGAlignments("Peptides/non-canonical.sorted.bam",use.names = TRUE)
junctions<-read.delim("RJunBase/RJunbase_annotated.txt",sep="\t")

# BAM structure
sum(duplicated(paste(names(non.can.bam),start(non.can.bam),end(non.can.bam),strand(non.can.bam),sep="|")))
length(non.can.bam[!duplicated(paste(names(non.can.bam),start(non.can.bam),end(non.can.bam),strand(non.can.bam),sep="|"))]) # 27,551

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
gr<-GRanges(non.can.bam[!duplicated(paste(names(non.can.bam),start(non.can.bam),end(non.can.bam),strand(non.can.bam),sep="|"))])

gr.UE<-GRanges(seqnames=junctions$chrom,
               ranges = IRanges(start=junctions$start-36,end=junctions$start),
               strand=junctions$strand)
gr.DE<-GRanges(seqnames=junctions$chrom,
               ranges = IRanges(start=junctions$end,end=junctions$end+36),
               strand=junctions$strand)

##### Find overlaps #####
result.UE<-findOverlaps(gr,gr.UE,type="within")
result.DE<-findOverlaps(gr,gr.DE,type="within")

##### Analysis #####
length(unique(c(queryHits(result.UE),queryHits(result.DE)))) # 1,094
queries<-c(queryHits(result.UE),queryHits(result.DE)) # All queries with hits
unique.queries<-queries[!duplicated(queries) & !duplicated(queries,fromLast=TRUE)] # 529 queries with only 1 hit
UE.DE.queries<-unique(queryHits(result.UE)[queryHits(result.UE) %in% queryHits(result.DE)]) # 22 queries with hit in UE and DE

UE.unique<-result.UE[queryHits(result.UE) %in% unique.queries] # 266
DE.unique<-result.DE[queryHits(result.DE) %in% unique.queries] # 263

UE.duplicated<-result.UE[!queryHits(result.UE) %in% unique.queries & 
                           !queryHits(result.UE) %in% UE.DE.queries] # 279 unique
DE.duplicated<-result.DE[!queryHits(result.DE) %in% unique.queries & 
                           !queryHits(result.DE) %in% UE.DE.queries] # 264 unique

UE.common<-result.UE[queryHits(result.UE) %in% UE.DE.queries] # 22
DE.common<-result.DE[queryHits(result.DE) %in% UE.DE.queries] # 22

View(cbind(as.data.frame(UE.unique),junctions[subjectHits(UE.unique),]))
View(cbind(as.data.frame(DE.unique),junctions[subjectHits(DE.unique),]))

View(cbind(as.data.frame(UE.duplicated),junctions[subjectHits(UE.duplicated),]))
View(cbind(as.data.frame(DE.duplicated),junctions[subjectHits(DE.duplicated),]))

View(cbind(as.data.frame(UE.common),junctions[subjectHits(UE.common),]))
View(cbind(as.data.frame(DE.common),junctions[subjectHits(DE.common),]))

## Search for more tumor-specific peptides (normal frequency=0 independent of tumor frequency)
# Queries qith unique hits
sum(nt.freq.indiv[subjectHits(UE.unique)]==0 & normal.freq.indiv[subjectHits(UE.unique)]==0) # 31
sum(nt.freq.indiv[subjectHits(DE.unique)]==0 & normal.freq.indiv[subjectHits(DE.unique)]==0) # 31

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

UE.DE.unique[!is.na(UE.DE.unique)] # 1 (all junctions non tumor-specific)
subjectHits(UE.common)[queryHits(UE.common) %in% UE.DE.unique]
subjectHits(DE.common)[queryHits(DE.common) %in% UE.DE.unique]
junctions[376866,]
junctions[376865,]

## Duplicated queries with all hits in UE having freq 0 in normal tissue
UE.duplicated.unique<-vector()
for (i in 1:(length(unique(queryHits(UE.duplicated))))){
  index<-subjectHits(UE.duplicated)[queryHits(UE.duplicated)==unique(queryHits(UE.duplicated))[i]]
  if (all(c(normal.freq.indiv[index]==0 &
            nt.freq.indiv[index]==0))){
    UE.duplicated.unique[i]<-unique(queryHits(UE.duplicated))[i]
  }
}

UE.duplicated.unique[!is.na(UE.duplicated.unique)] # 1 (all junctions non tumor-specific)
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

DE.duplicated.unique[!is.na(DE.duplicated.unique)] # 1 (all junctions non tumor-specific)
View(cbind(as.data.frame(DE.duplicated[queryHits(DE.duplicated) %in% DE.duplicated.unique[!is.na(DE.duplicated.unique)]]),
           junctions[subjectHits(DE.duplicated[queryHits(DE.duplicated) %in% DE.duplicated.unique[!is.na(DE.duplicated.unique)]]),]))


##### Table following Rjunbase criteria #####
# For duplicated peptides all junctions must be tumor specific to classify peptides as tumor specific
UE.duplicated.unique<-vector()
for (i in 1:(length(unique(queryHits(UE.duplicated))))){
  index<-subjectHits(UE.duplicated)[queryHits(UE.duplicated)==unique(queryHits(UE.duplicated))[i]]
  if (all(junctions$Tumor.specific.types[index]!="non tumor-specific")){
    UE.duplicated.unique[i]<-unique(queryHits(UE.duplicated))[i]
  }
}


DE.duplicated.unique<-vector()
for (i in 1:(length(unique(queryHits(DE.duplicated))))){
  index<-subjectHits(DE.duplicated)[queryHits(DE.duplicated)==unique(queryHits(DE.duplicated))[i]]
  if (all(junctions$Tumor.specific.types[index]!="non tumor-specific")){
    DE.duplicated.unique[i]<-unique(queryHits(UE.duplicated))[i]
  }
}

UE.DE.unique<-vector()
DE.UE.unique<-vector()
for (i in 1:length(UE.DE.queries)){
  ue<-grep(UE.DE.queries[i],queryHits(UE.common))
  de<-grep(UE.DE.queries[i],queryHits(DE.common))
  if(all(c(junctions$Tumor.specific.types[subjectHits(UE.common)[ue]]!="non tumor-specific",
           junctions$Tumor.specific.types[subjectHits(DE.common)[de]]!="non tumor-specific"))){
    UE.DE.unique[i]<-UE.DE.queries[i]
  }
  print(i)
}



df<-cbind(c(queryHits(result.UE),queryHits(result.DE)),
          c(subjectHits(result.UE),subjectHits(result.DE)),
          names(gr[c(queryHits(result.UE),queryHits(result.DE))]),
          start(gr)[c(queryHits(result.UE),queryHits(result.DE))],
          end(gr)[c(queryHits(result.UE),queryHits(result.DE))],
           junctions[c(subjectHits(result.UE),subjectHits(result.DE)),c(1,2,3,5,6,12:20)])
colnames(df)[c(1,2,3,4,5,14,15,16)]<-c("query","hit","peptide","peptide_start","peptide_end","Tumor_frequency","NT_frequency","Normal_frequency")
colnames(df)<-gsub("\\.","_",colnames(df))
df$Tumor_frequency<-sub("\\(.*","",df$Tumor_frequency)
df$NT_frequency<-sub("\\(.*","",df$NT_frequency)
df$Normal_frequency<-sub("\\(.*","",df$Normal_frequency)

df$Tumor_samples<-tumor.freq.indiv[df$hit]
df$NT_samples<-nt.freq.indiv[df$hit]
df$Normal_samples<-normal.freq.indiv[df$hit]

df$n_hits<-ifelse(df$query %in% unique.queries,"1",">1")
df$exon_hit<-ifelse(df$query %in% c(queryHits(UE.unique),queryHits(UE.duplicated)),"UE",
                    ifelse(df$query %in% c(queryHits(DE.unique),queryHits(DE.duplicated)),"DE",
                    ifelse(df$query %in% queryHits(UE.common),"UE;DE","NA")))

df$tumor_specific<-ifelse(df$query %in% c(queryHits(UE.unique),queryHits(DE.unique)) & df$Tumor_specific_types!="non tumor-specific","yes",
ifelse(df$query %in% UE.duplicated.unique,"yes",
ifelse(df$query %in% DE.duplicated.unique,"yes",
ifelse(df$query %in% UE.DE.unique,"yes","no"))))

df<-df[,c(1:5,23,24,6:16,20:22,17,25,18,19)]


# dfC1<-rbind(data.frame(query=queryHits(UE.unique),
#                      hit=subjectHits(UE.unique),
#                      peptide=names(gr[queryHits(UE.unique)]),
#                      peptide_start=start(gr[queryHits(UE.unique)]),
#                      n_hits="1",
#                      exon_hit="UE",
#                      tumor_specific=ifelse(junctions$Tumor.specific.types[subjectHits(UE.unique)]=="non tumor-specific","no","yes"),
#                      chromosome=junctions$chrom[subjectHits(UE.unique)],
#                      start=junctions$start[subjectHits(UE.unique)],
#                      end=junctions$end[subjectHits(UE.unique)],
#                      strand=junctions$strand[subjectHits(UE.unique)],
#                      origin=junctions$Alternative.splicing.types[subjectHits(UE.unique)],
#                      tumor_type=junctions$Tumor.specific.types[subjectHits(UE.unique)],
#                      tumor_freq=sub("\\(.*","",junctions$Tumor.frequency.sample.number.[subjectHits(UE.unique)]),
#                      nt_freq=sub("\\(.*","",junctions$NT.frequency.sample.number.[subjectHits(UE.unique)]),
#                      normal_freq=sub("\\(.*","",junctions$Normal.frequency.sample.number.[subjectHits(UE.unique)])
#                      ),
#           data.frame(query=queryHits(DE.unique),
#                      hit=subjectHits(DE.unique),
#                      peptide=names(gr[queryHits(DE.unique)]),
#                      peptide_start=start(gr[queryHits(DE.unique)]),
#                      n_hits="1",
#                      exon_hit="DE",
#                      tumor_specific=ifelse(junctions$Tumor.specific.types[subjectHits(DE.unique)]=="non tumor-specific","no","yes"),
#                      chromosome=junctions$chrom[subjectHits(DE.unique)],
#                      start=junctions$start[subjectHits(DE.unique)],
#                      end=junctions$end[subjectHits(DE.unique)],
#                      strand=junctions$strand[subjectHits(DE.unique)],
#                      origin=junctions$Alternative.splicing.types[subjectHits(DE.unique)],
#                      tumor_type=junctions$Tumor.specific.types[subjectHits(DE.unique)],
#                      tumor_freq=sub("\\(.*","",junctions$Tumor.frequency.sample.number.[subjectHits(DE.unique)]),
#                      nt_freq=sub("\\(.*","",junctions$NT.frequency.sample.number.[subjectHits(DE.unique)]),
#                      normal_freq=sub("\\(.*","",junctions$Normal.frequency.sample.number.[subjectHits(DE.unique)])
#                      ),
#           data.frame(query=queryHits(UE.duplicated),
#                      hit=subjectHits(UE.duplicated),
#                      peptide=names(gr[queryHits(UE.duplicated)]),
#                      peptide_start=start(gr[queryHits(UE.duplicated)]),
#                      n_hits=">1",
#                      exon_hit="UE",
#                      tumor_specific=ifelse(queryHits(UE.duplicated) %in% UE.duplicated.unique,"yes","no"),
#                      chromosome=junctions$chrom[subjectHits(UE.duplicated)],
#                      start=junctions$start[subjectHits(UE.duplicated)],
#                      end=junctions$end[subjectHits(UE.duplicated)],
#                      strand=junctions$strand[subjectHits(UE.duplicated)],
#                      origin=junctions$Alternative.splicing.types[subjectHits(UE.duplicated)],
#                      tumor_type=junctions$Tumor.specific.types[subjectHits(UE.duplicated)],
#                      tumor_freq=sub("\\(.*","",junctions$Tumor.frequency.sample.number.[subjectHits(UE.duplicated)]),
#                      nt_freq=sub("\\(.*","",junctions$NT.frequency.sample.number.[subjectHits(UE.duplicated)]),
#                      normal_freq=sub("\\(.*","",junctions$Normal.frequency.sample.number.[subjectHits(UE.duplicated)])
#                      ),
#           data.frame(query=queryHits(DE.duplicated),
#                      hit=subjectHits(DE.duplicated),
#                      peptide=names(gr[queryHits(DE.duplicated)]),
#                      peptide_start=start(gr[queryHits(DE.duplicated)]),
#                      n_hits=">1",
#                      exon_hit="DE",
#                      tumor_specific=ifelse(queryHits(DE.duplicated) %in% DE.duplicated.unique,"yes","no"),
#                      chromosome=junctions$chrom[subjectHits(DE.duplicated)],
#                      start=junctions$start[subjectHits(DE.duplicated)],
#                      end=junctions$end[subjectHits(DE.duplicated)],
#                      strand=junctions$strand[subjectHits(DE.duplicated)],
#                      origin=junctions$Alternative.splicing.types[subjectHits(DE.duplicated)],
#                      tumor_type=junctions$Tumor.specific.types[subjectHits(DE.duplicated)],
#                      tumor_freq=sub("\\(.*","",junctions$Tumor.frequency.sample.number.[subjectHits(DE.duplicated)]),
#                      nt_freq=sub("\\(.*","",junctions$NT.frequency.sample.number.[subjectHits(DE.duplicated)]),
#                      normal_freq=sub("\\(.*","",junctions$Normal.frequency.sample.number.[subjectHits(DE.duplicated)])
#                      ), 
#           data.frame(query=queryHits(UE.common),
#                      hit=subjectHits(UE.common),
#                      peptide=names(gr[queryHits(UE.common)]),
#                      peptide_start=start(gr[queryHits(UE.common)]),
#                      n_hits=">1",
#                      exon_hit="UE;DE",
#                      tumor_specific=ifelse(queryHits(UE.common) %in% UE.DE.unique,"yes","no"),
#                      chromosome=junctions$chrom[subjectHits(UE.common)],
#                      start=junctions$start[subjectHits(UE.common)],
#                      end=junctions$end[subjectHits(UE.common)],
#                      strand=junctions$strand[subjectHits(UE.common)],
#                      origin=junctions$Alternative.splicing.types[subjectHits(UE.common)],
#                      tumor_type=junctions$Tumor.specific.types[subjectHits(UE.common)],
#                      tumor_freq=sub("\\(.*","",junctions$Tumor.frequency.sample.number.[subjectHits(UE.common)]),
#                      nt_freq=sub("\\(.*","",junctions$NT.frequency.sample.number.[subjectHits(UE.common)]),
#                      normal_freq=sub("\\(.*","",junctions$Normal.frequency.sample.number.[subjectHits(UE.common)])
#                      ),
#           data.frame(query=queryHits(DE.common),
#                      hit=subjectHits(DE.common),
#                      peptide=names(gr[queryHits(DE.common)]),
#                      peptide_start=start(gr[queryHits(DE.common)]),
#                      n_hits=">1",
#                      exon_hit="UE;DE",
#                      tumor_specific=ifelse(queryHits(DE.common) %in% DE.UE.unique,"yes","no"),
#                      chromosome=junctions$chrom[subjectHits(DE.common)],
#                      start=junctions$start[subjectHits(DE.common)],
#                      end=junctions$end[subjectHits(DE.common)],
#                      strand=junctions$strand[subjectHits(DE.common)],
#                      origin=junctions$Alternative.splicing.types[subjectHits(DE.common)],
#                      tumor_type=junctions$Tumor.specific.types[subjectHits(DE.common)],
#                      tumor_freq=sub("\\(.*","",junctions$Tumor.frequency.sample.number.[subjectHits(DE.common)]),
#                      nt_freq=sub("\\(.*","",junctions$NT.frequency.sample.number.[subjectHits(DE.common)]),
#                      normal_freq=sub("\\(.*","",junctions$Normal.frequency.sample.number.[subjectHits(DE.common)])
#                      )
# )

length(unique(df$query[df$tumor_specific=="yes"])) # 14
length(unique(df$peptide[df$tumor_specific=="yes"])) # 8

pep.file <- file("./RJunBase/Peptides2RjunbaseC1.txt", "wb")
write.table(x=df,file=pep.file,row.names = FALSE,quote=FALSE,col.names=TRUE,sep="\t")
close(pep.file)

##### Table following my criteria: normal frequency = 0 #####
UE.DE.unique<-vector() # 1
for (i in 1:length(UE.DE.queries)){
  ue<-grep(UE.DE.queries[i],queryHits(UE.common))
  de<-grep(UE.DE.queries[i],queryHits(DE.common))
  if(all(c(normal.freq.indiv[subjectHits(UE.common)[ue]]==0 &
           nt.freq.indiv[subjectHits(UE.common)[ue]]==0 &
           normal.freq.indiv[subjectHits(DE.common)[de]]==0 &
           nt.freq.indiv[subjectHits(DE.common)[de]]==0))){
    UE.DE.unique<-c(UE.DE.unique,UE.DE.queries[i])
  }
  print(i)
}

UE.duplicated.unique<-vector() # 1
for (i in 1:(length(unique(queryHits(UE.duplicated))))){
  index<-subjectHits(UE.duplicated)[queryHits(UE.duplicated)==unique(queryHits(UE.duplicated))[i]]
  if (all(c(normal.freq.indiv[index]==0 &
            nt.freq.indiv[index]==0))){
    UE.duplicated.unique<-c(UE.duplicated.unique,unique(queryHits(UE.duplicated))[i])
  }
}

DE.duplicated.unique<-vector() # 1
for (i in 1:(length(unique(queryHits(DE.duplicated))))){
  index<-subjectHits(DE.duplicated)[queryHits(DE.duplicated)==unique(queryHits(DE.duplicated))[i]]
  if (all(c(normal.freq.indiv[index]==0 &
            nt.freq.indiv[index]==0))){
    DE.duplicated.unique<-c(DE.duplicated.unique,unique(queryHits(DE.duplicated))[i])
  }
}



dfC2<-df
dfC2$tumor_specific[dfC2$query %in% c(queryHits(UE.unique),queryHits(DE.unique)) & 
                            dfC2$NT_frequency==0 &
                            dfC2$Normal_frequency==0]<-"yes"
dfC2$tumor_specific[dfC2$query %in% c(UE.DE.unique,UE.duplicated.unique,DE.duplicated.unique)]<-"yes"


# dfC2$tumor_specific[df$peptide %in% names(gr)[UE.DE.unique[!is.na(UE.DE.unique)]]]<-"yes"
# dfC2$tumor_specific[df$peptide %in% names(gr)[UE.duplicated.unique[!is.na(UE.duplicated.unique)]]]<-"yes"
# dfC2$tumor_specific[df$peptide %in% names(gr)[DE.duplicated.unique[!is.na(DE.duplicated.unique)]]]<-"yes"

length(unique(dfC2$query[dfC2$tumor_specific=="yes"])) # 65
length(unique(dfC2$peptide[dfC2$tumor_specific=="yes"])) # 45

pep.file <- file("./RJunBase/Peptides2RjunbaseC2.txt", "wb")
write.table(x=dfC2,file=pep.file,row.names = FALSE,quote=FALSE,col.names=TRUE,sep="\t")
close(pep.file)

##### Barplots #####
# https://cran.r-project.org/web/packages/chromoMap/vignettes/chromoMap.html
# All overlaps by chromosome and by start coordinate
# Frequency as percentage. Tumor: 10,283 samples. Normal (GTEx): 17,382. NT:¿?
# Also we have median expression
# chromosome coverage¿?
# position must relate to UE or DE or both

index<-unique(c(subjectHits(result.UE),subjectHits(result.DE))) # 2006 unique junctions (11 shared between result.UE and result.DE)
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
          plot_color="blue",
          title="Normal")
chromoMap("chr_file_1.txt",
          "chr1_annot_tumor.txt",
          data_based_color_map = TRUE,
          data_type="numeric",
          plots="bar",
          plot_color="red",
          title="Tumor",
          n_win.factor = 1) 






gr.normal<-GRanges(seqnames="chr1",ranges=IRanges(start=pos,end=pos))
kp <- plotKaryotype(plot.type=2, chromosomes = "chr1",genome="hg38",
                    zoom=GRanges(seqnames="chr1",ranges=IRanges(start=0,end=7000000)))
kpPlotCoverage(kp, data=gr.normal[1:5])


##### Plot mean frequencies by peptide (ggplot) #####
plot.junctions<-function(df,chr,n=0,m=0){
  y.tumor<-numeric()
  y.normal<-numeric()
  y.nt<-numeric()
  times<-vector()
  peptide.names<-vector()
  peptide.start<-vector()
  df<-df[df$chromosome==chr,]
  print(df)
  index.query<-unique(df$query)
  for (i in 1:length(index.query)){
    times[i]<-sum(df$query==index.query[i])
    t.f<-as.numeric(df$tumor_freq[df$query==index.query[i]])
    n.f<-as.numeric(df$normal_freq[df$query==index.query[i]])
    nt.f<-as.numeric(df$nt_freq[df$query==index.query[i]])
    peptide.names[i]<-unique(df$peptide[df$query==index.query[i]])
    peptide.start[i]<-unique(df$peptide_start[df$query==index.query[i]])
    
    y.tumor<-c(y.tumor,mean(t.f))
    y.normal<-c(y.normal,mean(n.f))
    y.nt<-c(y.nt,mean(nt.f))
  }
  index<-order(y.tumor,decreasing = FALSE)
  df.freq<-data.frame(start=peptide.start[index],
                      peptide=paste0(peptide.names[index],"\n","(",peptide.start[index],")"),
                      freq=c(y.tumor[index],y.nt[index],y.normal[index]),
                      tissue=c(rep("tumor",length(y.tumor)),rep("NT",length(y.nt)),rep("normal",length(y.normal))))
  if(n>0 & m>0){
  df.freq<-df.freq[c(n:m,
                  which(df.freq$tissue=="NT")[n]:(which(df.freq$tissue=="NT")[n]+(m-n)),
                  which(df.freq$tissue=="normal")[n]:(which(df.freq$tissue=="normal")[n]+(m-n))),]
  }
  print(df.freq)
  ggplot(data=df.freq,aes(x=fct_inorder(peptide),y=freq,fill=factor(tissue,levels=unique(tissue))))+
    geom_bar(stat="identity",position=position_dodge(),width=0.7)+
    labs(fill="Tissue")+
    # scale_fill_discrete(labels=c("Tumor","Adjacent","GTEx"))+
    theme(axis.text.x = element_text(size=7,angle=45))+
    ggtitle(chr)+
    xlab("Peptide (start)")+
    ylab("Junction mean frequency")+
    ylim(0,1)+
    scale_fill_manual(values=c("#E69F00","#009E73","#39BEB1"),labels=c("Tumor","Adjacent","GTEx"))
}
par(mfrow=c(2,2))
a<-plot.junctions(dfC2,"chr1",n=1,m=20)
b<-plot.junctions(dfC2,"chr1",n=21,m=30)
c<-plot.junctions(dfC2,"chr1",n=31,m=50)
d<-plot.junctions(dfC2,"chr1",n=51,m=70)

ggarrange(a+rremove("legend"),b+rremove("legend"),ncol=1,nrow=2)
plot.junctions(dfC2,"chr18")



#####################
library(AnnotationHub)
ah<-AnnotationHub()
ahub.chain <- subset(AnnotationHub(), rdataclass == "ChainFile" & species == "Homo sapiens")
query(ahub.chain, c("hg19", "hg38"))
chain <- ahub.chain[ahub.chain$title == "hg19ToHg38.over.chain.gz"]
chain <- chain[[1]]

RI<-rtracklayer::import("C:/Users/bea_f/Downloads/merge_graphs_intron_retention_C2.gff3")
RI<-RI[RI$type=="exon"]
seqlevels(RI) <- c(paste0("chr",seqlevels(RI)[1:24]),"chrM")
RI.hg38 <- liftOver(RI, chain)


seqlevels(RI.hg38) <- c(seqlevels(RI.hg38)[1:24],"chrMT")
result.RI<-findOverlaps(gr,RI.hg38,type="within")

############################
flair<-rtracklayer::import("C:/Users/bea_f/OneDrive/Escritorio/HematoLaFe/GTFs/flair_filter_transcripts.gtf")
gencode<-rtracklayer::import("C:/Users/bea_f/OneDrive/Escritorio/HematoLaFe/GTFs/gencode.v36.annotation.gtf")
beat<-rtracklayer::import("C:/Users/bea_f/OneDrive/Escritorio/HematoLaFe/GTFs/BEATv36_AML.gtf")
tcga<-rtracklayer::import("C:/Users/bea_f/OneDrive/Escritorio/HematoLaFe/GTFs/merged_TCGA.gtf")
healthy<-rtracklayer::import("C:/Users/bea_f/OneDrive/Escritorio/HematoLaFe/GTFs/BEATv36_healthy.gtf")


seqlevels(flair)
seqlevels(gencode)
seqlevels(gr) <- c(seqlevels(gr)[1:24],"chrM")

flair<-flair[flair$type=="exon"]
gencode<-gencode[gencode$type=="exon" | gencode$type=="UTR"]
beat<-beat[beat$type=="exon"]
tcga<-tcga[tcga$type=="exon"]

result.flair<-findOverlaps(gr,flair,type="within")
result.gencode<-findOverlaps(gr,gencode,type="within")


length(unique(c(queryHits(result.flair),queryHits(result.gencode))))
View(as.data.frame(gr[-unique(c(queryHits(result.flair),queryHits(result.gencode)))]))
gr.target<-gr[-unique(c(queryHits(result.flair),queryHits(result.gencode)))]
View(as.data.frame(gr.target,row.names = NULL))

result.beat<-findOverlaps(gr.target,beat,type="within")
result.tcga<-findOverlaps(gr.target,tcga,type="within")
result.healthy<-findOverlaps(gr.target,healthy,type="within")
View(as.data.frame(gr.target[unique(c(queryHits(result.beat)[!queryHits(result.beat) %in% queryHits(result.healthy)],
                                      queryHits(result.tcga)[!queryHits(result.tcga) %in% queryHits(result.healthy)]))],row.names = NULL))
