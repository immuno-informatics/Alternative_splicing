# Libraries
library(seqinr)
library(Biostrings)
library(BSgenome.Hsapiens.UCSC.hg38)
library(TAPseq)
library(data.table)

# GTFs
gencode<-rtracklayer::import("gencode.v36.annotation.gtf")
beat<-rtracklayer::import("BEATv36_AML.gtf")
tcga<-rtracklayer::import("merged_TCGA.gtf")
flair<-rtracklayer::import("flair_filter_transcripts.gtf")
hl60.1<-rtracklayer::import("GSM5331246_ENCFF770VLX_transcriptome_annotations_GRCh38.gtf")
hl60.2<-rtracklayer::import("GSM5331247_ENCFF376DFZ_transcriptome_annotations_GRCh38.gtf")

# RJunBase
junctions<-read.delim("RJunBase/RJunbase_annotated.txt",sep="\t")

###################### RJunBase: translate junctions ###########################
##### Translate 18 bp upstream and dowstream the junctions #####
UE.gr<-GRanges(seqnames=junctions$chrom,
               ranges=IRanges(start=junctions$start-17,end=junctions$start),
               strand=junctions$strand)
DE.gr<-GRanges(seqnames=junctions$chrom,
               ranges=IRanges(start=junctions$end,end=junctions$end+17),
               strand=junctions$strand)

UE.rna<-getSeq(BSgenome.Hsapiens.UCSC.hg38,UE.gr)
DE.rna<-getSeq(BSgenome.Hsapiens.UCSC.hg38,DE.gr)
target.seq<-ifelse(junctions$strand=="+",paste0(UE.rna,DE.rna),paste0(DE.rna,UE.rna))
names(target.seq)<-junctions$JunctionID

# Get 3 reading frames from each sequence
# Differentiate between std and mt genomes for next step
target.seq.std<-lapply(1:3,function(i){substring(target.seq[junctions$chrom!="chrM"],
                                                 i,
                                                 width(target.seq[junctions$chrom!="chrM"]))})
target.seq.mt<-lapply(1:3,function(i){substring(target.seq[junctions$chrom=="chrM"],
                                                i,
                                                width(target.seq[junctions$chrom=="chrM"]))})


# Check available genetic codes in Biostrings
SGC0<-getGeneticCode("SGC0")
SGC1<-getGeneticCode("SGC1")

# Translate (use no.init.codon=TRUE to avoid the use of alternative initiation codons)
# std: from ATG to TAA/TAG/TGA
orf1.std<-as.character(Biostrings::translate(DNAStringSet(target.seq.std[[1]]),
                                             genetic.code = SGC0,
                                             if.fuzzy.codon = "solve",
                                             no.init.codon = TRUE))
orf2.std<-as.character(Biostrings::translate(DNAStringSet(target.seq.std[[2]]),
                                             genetic.code = SGC0,
                                             if.fuzzy.codon = "solve",
                                             no.init.codon = TRUE))
orf3.std<-as.character(Biostrings::translate(DNAStringSet(target.seq.std[[3]]),
                                             genetic.code = SGC0,
                                             if.fuzzy.codon = "solve",
                                             no.init.codon = TRUE))
# mt: from ATA/ATG to TAA/TAG/AGA/AGG)
orf1.mt<-as.character(Biostrings::translate(DNAStringSet(target.seq.mt[[1]]),
                                            genetic.code = SGC1,
                                            if.fuzzy.codon = "solve",
                                            no.init.codon = TRUE))
orf2.mt<-as.character(Biostrings::translate(DNAStringSet(target.seq.mt[[2]]),
                                            genetic.code = SGC1,
                                            if.fuzzy.codon = "solve",
                                            no.init.codon = TRUE))
orf3.mt<-as.character(Biostrings::translate(DNAStringSet(target.seq.mt[[3]]),
                                            genetic.code = SGC1,
                                            if.fuzzy.codon = "solve",
                                            no.init.codon = TRUE))

##### Keep ORFs that are >=7aa #####
# 1. Split each transcript sequence by stop codon (*)
# 2. If sequence before first stop codon is >=7 aa keep it
# 3. For the rest, we need to find an ORF (i.e. beginning with M)
# 4. Report only sequences >=7aa

# ORF1 std
orf1.std.split<-sapply(orf1.std,function(x){
  a<-unlist(strsplit(x,"\\*"))
  nchar(a[1])>=7
  b<-a[-1]
  b<-b[grep("M",b)]
  pos<-sapply(b,function(x){unlist(gregexpr("M",b))[1]})
  c<-substring(b,pos,nchar(b))
  return(c(a[1],c)[nchar(c(a[1],c))>=7])
})
orf1.std.stack<-stack(orf1.std.split)
length(names(orf1.std.split)[which(!names(orf1.std.split) %in% orf1.std.stack$ind)]) # 126,597 transcripts do not generate sequences >=7aa; drop this levels
orf1.std.stack<-stack(orf1.std.split,drop=TRUE)
sum(duplicated(orf1.std.stack$ind)) # 0
orf1.std.stack$ID<-paste0(orf1.std.stack$ind,"_ORF1")
sum(duplicated(orf1.std.stack$values)) # 5,007,531 are duplicated sequences: aggregate
orf1.std.stack<-aggregate(ID~values,orf1.std.stack,paste,collapse=";") # 530,543

# ORF2 std
orf2.std.split<-sapply(orf2.std,function(x){
  a<-unlist(strsplit(x,"\\*"))
  nchar(a[1])>=7
  b<-a[-1]
  b<-b[grep("M",b)]
  pos<-sapply(b,function(x){unlist(gregexpr("M",b))[1]})
  c<-substring(b,pos,nchar(b))
  return(c(a[1],c)[nchar(c(a[1],c))>=7])
})
orf2.std.stack<-stack(orf2.std.split)
length(names(orf2.std.split)[which(!names(orf2.std.split) %in% orf2.std.stack$ind)]) # 159,296 transcripts do not generate sequences >=7aa; drop this levels
orf2.std.stack<-stack(orf2.std.split,drop=TRUE)
sum(duplicated(orf2.std.stack$ind)) # 0
orf2.std.stack$ID<-paste0(orf2.std.stack$ind,"_ORF2")
sum(duplicated(orf2.std.stack$values)) # 21,879 are duplicated sequences: aggregate
orf2.std.stack<-aggregate(ID~values,orf2.std.stack,paste,collapse=";") # 500,098

# ORF3 std
orf3.std.split<-sapply(orf3.std,function(x){
  a<-unlist(strsplit(x,"\\*"))
  nchar(a[1])>=7
  b<-a[-1]
  b<-b[grep("M",b)]
  pos<-sapply(b,function(x){unlist(gregexpr("M",b))[1]})
  c<-substring(b,pos,nchar(b))
  return(c(a[1],c)[nchar(c(a[1],c))>=7])
})
orf3.std.stack<-stack(orf3.std.split)
length(names(orf3.std.split)[which(!names(orf3.std.split) %in% orf3.std.stack$ind)]) # 147,921 transcripts do not generate sequences >=7aa; drop this levels
orf3.std.stack<-stack(orf3.std.split,drop=TRUE)
sum(duplicated(orf3.std.stack$ind)) # 0
orf3.std.stack$ID<-paste0(orf3.std.stack$ind,"_ORF3")
sum(duplicated(orf3.std.stack$values)) # 24,242 are duplicated sequences: aggregate
orf3.std.stack<-aggregate(ID~values,orf3.std.stack,paste,collapse=";") # 509,110

# ORF1 mt
orf1.mt.split<-sapply(orf1.mt,function(x){
  a<-unlist(strsplit(x,"\\*"))
  nchar(a[1])>=7
  b<-a[-1]
  b<-b[grep("M",b)]
  pos<-sapply(b,function(x){unlist(gregexpr("M",b))[1]})
  c<-substring(b,pos,nchar(b))
  return(c(a[1],c)[nchar(c(a[1],c))>=7])
})
orf1.mt.stack<-stack(orf1.mt.split)
length(names(orf1.mt.split)[which(!names(orf1.mt.split) %in% orf1.mt.stack$ind)]) # 408 transcripts do not generate sequences >=7aa; drop this levels
orf1.mt.stack<-stack(orf1.mt.split,drop=TRUE)
sum(duplicated(orf1.mt.stack$ind)) # 0
orf1.mt.stack$ID<-paste0(orf1.mt.stack$ind,"_ORF1")
sum(duplicated(orf1.mt.stack$values)) # 0 are duplicated sequences: aggregate
orf1.mt.stack<-aggregate(ID~values,orf1.mt.stack,paste,collapse=";") # 336

# ORF2 mt
orf2.mt.split<-sapply(orf2.mt,function(x){
  a<-unlist(strsplit(x,"\\*"))
  nchar(a[1])>=7
  b<-a[-1]
  b<-b[grep("M",b)]
  pos<-sapply(b,function(x){unlist(gregexpr("M",b))[1]})
  c<-substring(b,pos,nchar(b))
  return(c(a[1],c)[nchar(c(a[1],c))>=7])
})
orf2.mt.stack<-stack(orf2.mt.split)
length(names(orf2.mt.split)[which(!names(orf2.mt.split) %in% orf2.mt.stack$ind)]) # 545 transcripts do not generate sequences >=7aa; drop this levels
orf2.mt.stack<-stack(orf2.mt.split,drop=TRUE)
sum(duplicated(orf2.mt.stack$ind)) # 0
orf2.mt.stack$ID<-paste0(orf2.mt.stack$ind,"_ORF2")
sum(duplicated(orf2.mt.stack$values)) # 1 are duplicated sequences: aggregate
orf2.mt.stack<-aggregate(ID~values,orf2.mt.stack,paste,collapse=";") # 198

# ORF3 mt
orf3.mt.split<-sapply(orf3.mt,function(x){
  a<-unlist(strsplit(x,"\\*"))
  nchar(a[1])>=7
  b<-a[-1]
  b<-b[grep("M",b)]
  pos<-sapply(b,function(x){unlist(gregexpr("M",b))[1]})
  c<-substring(b,pos,nchar(b))
  return(c(a[1],c)[nchar(c(a[1],c))>=7])
})
orf3.mt.stack<-stack(orf3.mt.split)
length(names(orf3.mt.split)[which(!names(orf3.mt.split) %in% orf3.mt.stack$ind)]) # 376 transcripts do not generate sequences >=7aa; drop this levels
orf3.mt.stack<-stack(orf3.mt.split,drop=TRUE)
sum(duplicated(orf3.mt.stack$ind)) # 0
orf3.mt.stack$ID<-paste0(orf3.mt.stack$ind,"_ORF3")
sum(duplicated(orf3.mt.stack$values)) # 0 are duplicated sequences: aggregate
orf3.mt.stack<-aggregate(ID~values,orf3.mt.stack,paste,collapse=";") # 368

# Rbind the data frames and aggregate to get unique sequences
rjunbase.fa<-rbind(orf1.std.stack,orf2.std.stack,orf3.std.stack,
                   orf1.mt.stack,orf2.mt.stack,orf3.mt.stack)
rjunbase.fa<-aggregate(ID~values,rjunbase.fa,paste,collapse=";") # 1,540,536
rjunbase.fa<-data.frame(ID=rjunbase.fa$ID,sequence=rjunbase.fa$values)

# Write only unannotated junctions
sum(grepl("^UN_",junctions$JunctionID)) # 324,651
write.table(x=rjunbase.fa[grep("^UN_",rjunbase.fa$ID),],
            file="Fragpipe/inputs/RJunBase_pep.tsv",row.names = FALSE,quote=FALSE,col.names=TRUE,sep="\t")


############################ rMATS Beat AML ####################################
################################## SE ##########################################
SE<-read.delim("SE_annotated.txt",sep="\t")

# Select events
i1<-which(SE$FDR<0.05 & SE$IncLevelDifference>0) # inc (which() returns row position)
i2<-which(SE$FDR<0.05 & SE$IncLevelDifference<0) # exc

i3<-which(SE$all.0.healthy==19 & SE$gt.0.aml>=1) # 57,005 (inc)
i4<-which(SE$all.1.healthy==19 & SE$lt.1.aml>=1) # 189,937 (exc)


length(unique(c(i1,i2,i3,i4))) # 280,210
SE<-SE[unique(c(i1,i2,i3,i4)),] # rownames are row original positions
sum(unique(c(i1,i3)) %in% unique(c(i2,i4))) # 3,244 that could be inclusion or exclusion depending on IncLevelDifference


UE.gr<-GRanges(seqnames=SE$chr,
               ranges=IRanges(start=SE$upstreamES+1,end=SE$upstreamEE),
               strand=SE$strand)
SE.gr<-GRanges(seqnames=SE$chr,
               ranges=IRanges(start=SE$exonStart_0base+1,end=SE$exonEnd),
               strand=SE$strand)
DE.gr<-GRanges(seqnames=SE$chr,
               ranges=IRanges(start=SE$downstreamES+1,end=SE$downstreamEE),
               strand=SE$strand)

# Get sequence por upstream, downstream and skipping exons
UE.rna<-getSeq(BSgenome.Hsapiens.UCSC.hg38,UE.gr)
SE.rna<-getSeq(BSgenome.Hsapiens.UCSC.hg38,SE.gr)
DE.rna<-getSeq(BSgenome.Hsapiens.UCSC.hg38,DE.gr)

# Inclusion and exclusion isoforms
rna.inc<-ifelse(SE$strand[rownames(SE) %in% unique(c(i1,i3))]=="+",
                paste0(UE.rna,SE.rna,DE.rna)[rownames(SE) %in% unique(c(i1,i3))],
                paste0(DE.rna,SE.rna,UE.rna)[rownames(SE) %in% unique(c(i1,i3))])
rna.exc<-ifelse(SE$strand[rownames(SE) %in% unique(c(i2,i4))]=="+",
                paste0(UE.rna,DE.rna)[rownames(SE) %in% unique(c(i2,i4))],
                paste0(DE.rna,UE.rna)[rownames(SE) %in% unique(c(i2,i4))])

names(rna.inc)<-paste0("AMLvsHealthy_SEinc_ID",SE$ID[rownames(SE) %in% unique(c(i1,i3))])
names(rna.exc)<-paste0("AMLvsHealthy_SEexc_ID",SE$ID[rownames(SE) %in% unique(c(i2,i4))])


length(rna.inc)+length(rna.exc) # 280,210

inc.seq<-lapply(1:3,function(i){substring(rna.inc,i,nchar(rna.inc))})
exc.seq<-lapply(1:3,function(i){substring(rna.exc,i,nchar(rna.exc))})

orf1.inc<-as.character(Biostrings::translate(DNAStringSet(inc.seq[[1]]),
                                             if.fuzzy.codon = "solve",
                                             no.init.codon = TRUE))
orf2.inc<-as.character(Biostrings::translate(DNAStringSet(inc.seq[[2]]),
                                             if.fuzzy.codon = "solve",
                                             no.init.codon = TRUE))
orf3.inc<-as.character(Biostrings::translate(DNAStringSet(inc.seq[[3]]),
                                             if.fuzzy.codon = "solve",
                                             no.init.codon = TRUE))

orf1.exc<-as.character(Biostrings::translate(DNAStringSet(exc.seq[[1]]),
                                             if.fuzzy.codon = "solve",
                                             no.init.codon = TRUE))
orf2.exc<-as.character(Biostrings::translate(DNAStringSet(exc.seq[[2]]),
                                             if.fuzzy.codon = "solve",
                                             no.init.codon = TRUE))
orf3.exc<-as.character(Biostrings::translate(DNAStringSet(exc.seq[[3]]),
                                             if.fuzzy.codon = "solve",
                                             no.init.codon = TRUE))

##### Keep ORFs that are >=7aa #####
# ORF1
orf1.inc.split<-sapply(orf1.inc,function(x){
  a<-unlist(strsplit(x,"\\*"))
  nchar(a[1])>=7
  b<-a[-1]
  b<-b[grep("M",b)]
  pos<-sapply(b,function(x){unlist(gregexpr("M",b))[1]})
  c<-substring(b,pos,nchar(b))
  return(c(a[1],c)[nchar(c(a[1],c))>=7])
})
orf1.inc.stack<-stack(orf1.inc.split)
length(names(orf1.inc.split)[which(!names(orf1.inc.split) %in% orf1.inc.stack$ind)]) # 2,476 events do not generate sequences >=7aa; drop this levels
orf1.inc.stack<-stack(orf1.inc.split,drop=TRUE)
sum(duplicated(orf1.inc.stack$ind))
orf1.inc.stack$ID<-paste0(orf1.inc.stack$ind,"_ORF1.",unlist(sapply(table(orf1.inc.stack$ind),function(x){1:x})))
sum(duplicated(orf1.inc.stack$values)) # 42,960 are duplicated sequences: aggregate
orf1.inc.stack<-aggregate(ID~values,orf1.inc.stack,paste,collapse=";")

# ORF1
orf1.exc.split<-sapply(orf1.exc,function(x){
  a<-unlist(strsplit(x,"\\*"))
  nchar(a[1])>=7
  b<-a[-1]
  b<-b[grep("M",b)]
  pos<-sapply(b,function(x){unlist(gregexpr("M",b))[1]})
  c<-substring(b,pos,nchar(b))
  return(c(a[1],c)[nchar(c(a[1],c))>=7])
})
orf1.exc.stack<-stack(orf1.exc.split)
length(names(orf1.exc.split)[which(!names(orf1.exc.split) %in% orf1.exc.stack$ind)]) # 12,781 events do not generate sequences >=7aa; drop this levels
orf1.exc.stack<-stack(orf1.exc.split,drop=TRUE)
orf1.exc.stack$ID<-paste0(orf1.exc.stack$ind,"_ORF1.",unlist(sapply(table(orf1.exc.stack$ind),function(x){1:x})))
sum(duplicated(orf1.exc.stack$values)) # 182,727 are duplicated sequences: aggregate
orf1.exc.stack<-aggregate(ID~values,orf1.exc.stack,paste,collapse=";")

# ORF2
orf2.inc.split<-sapply(orf2.inc,function(x){
  a<-unlist(strsplit(x,"\\*"))
  nchar(a[1])>=7
  b<-a[-1]
  b<-b[grep("M",b)]
  pos<-sapply(b,function(x){unlist(gregexpr("M",b))[1]})
  c<-substring(b,pos,nchar(b))
  return(c(a[1],c)[nchar(c(a[1],c))>=7])
})
orf2.inc.stack<-stack(orf2.inc.split)
length(names(orf2.inc.split)[which(!names(orf2.inc.split) %in% orf2.inc.stack$ind)]) # 3,481 events do not generate sequences >=7aa; drop this levels
orf2.inc.stack<-stack(orf2.inc.split,drop=TRUE)
orf2.inc.stack$ID<-paste0(orf2.inc.stack$ind,"_ORF2.",unlist(sapply(table(orf2.inc.stack$ind),function(x){1:x})))
sum(duplicated(orf2.inc.stack$values)) # 45,595 are duplicated sequences: aggregate
orf2.inc.stack<-aggregate(ID~values,orf2.inc.stack,paste,collapse=";")


# ORF2
orf2.exc.split<-sapply(orf2.exc,function(x){
  a<-unlist(strsplit(x,"\\*"))
  nchar(a[1])>=7
  b<-a[-1]
  b<-b[grep("M",b)]
  pos<-sapply(b,function(x){unlist(gregexpr("M",b))[1]})
  c<-substring(b,pos,nchar(b))
  return(c(a[1],c)[nchar(c(a[1],c))>=7])
})
orf2.exc.stack<-stack(orf2.exc.split)
length(names(orf2.exc.split)[which(!names(orf2.exc.split) %in% orf2.exc.stack$ind)]) # 18,610 events do not generate sequences >=7aa; drop this levels
orf2.exc.stack<-stack(orf2.exc.split,drop=TRUE)
orf2.exc.stack$ID<-paste0(orf2.exc.stack$ind,"_ORF2.",unlist(sapply(table(orf2.exc.stack$ind),function(x){1:x})))
sum(duplicated(orf2.exc.stack$values)) # 186,941 are duplicated sequences: aggregate
orf2.exc.stack<-aggregate(ID~values,orf2.exc.stack,paste,collapse=";")

# ORF3
orf3.inc.split<-sapply(orf3.inc,function(x){
  a<-unlist(strsplit(x,"\\*"))
  nchar(a[1])>=7
  b<-a[-1]
  b<-b[grep("M",b)]
  pos<-sapply(b,function(x){unlist(gregexpr("M",b))[1]})
  c<-substring(b,pos,nchar(b))
  return(c(a[1],c)[nchar(c(a[1],c))>=7])
})
orf3.inc.stack<-stack(orf3.inc.split)
length(names(orf3.inc.split)[which(!names(orf3.inc.split) %in% orf3.inc.stack$ind)]) # 3,611 events do not generate sequences >=7aa; drop this levels
orf3.inc.stack<-stack(orf3.inc.split,drop=TRUE)
orf3.inc.stack$ID<-paste0(orf3.inc.stack$ind,"_ORF3.",unlist(sapply(table(orf3.inc.stack$ind),function(x){1:x})))
sum(duplicated(orf3.inc.stack$values)) # 42,599 are duplicated sequences: aggregate
orf3.inc.stack<-aggregate(ID~values,orf3.inc.stack,paste,collapse=";")

# ORF3
orf3.exc.split<-sapply(orf3.exc,function(x){
  a<-unlist(strsplit(x,"\\*"))
  nchar(a[1])>=7
  b<-a[-1]
  b<-b[grep("M",b)]
  pos<-sapply(b,function(x){unlist(gregexpr("M",b))[1]})
  c<-substring(b,pos,nchar(b))
  return(c(a[1],c)[nchar(c(a[1],c))>=7])
})
orf3.exc.stack<-stack(orf3.exc.split)
length(names(orf3.exc.split)[which(!names(orf3.exc.split) %in% orf3.exc.stack$ind)]) # 20,225 events do not generate sequences >=7aa; drop this levels
orf3.exc.stack<-stack(orf3.exc.split,drop=TRUE)
orf3.exc.stack$ID<-paste0(orf3.exc.stack$ind,"_ORF3.",unlist(sapply(table(orf3.exc.stack$ind),function(x){1:x})))
sum(duplicated(orf3.exc.stack$values)) # 177,244 are duplicated sequences: aggregate
orf3.exc.stack<-aggregate(ID~values,orf3.exc.stack,paste,collapse=";")

# Rbind and aggregate duplicated sequences
SE.fa<-rbind(orf1.inc.stack,orf1.exc.stack,orf2.inc.stack,orf2.exc.stack,orf3.inc.stack,orf3.exc.stack)
SE.fa<-aggregate(ID~values,SE.fa,paste,collapse=";") # 812,933
SE.fa<-data.frame(ID=SE.fa$ID,sequence=SE.fa$values)

##### Write .tsv with unique sequences #####
write.table(x=SE.fa,file="Fragpipe/inputs/BeatAMLvsHealthy_SE.tsv",row.names = FALSE,quote=FALSE,col.names=TRUE,sep="\t")

################################# A5SS #########################################
A5SS<-read.delim("A5SS_annotated.txt")
# Select events
i1<-which(A5SS$FDR<0.05 & A5SS$IncLevelDifference>0) # 1,934 inc (which returns row position)
i2<-which(A5SS$FDR<0.05 & A5SS$IncLevelDifference<0) # 3,325 exc

i3<-which(A5SS$all.0.healthy==19 & A5SS$gt.0.aml>=1) # 14,490 (inc)
i4<-which(A5SS$all.1.healthy==19 & A5SS$lt.1.aml>=1) # 35,971 (exc)

length(unique(c(i1,i2,i3,i4))) # 55,700
A5SS<-A5SS[unique(c(i1,i2,i3,i4)),] # rownames are row original positions
sum(unique(c(i1,i3)) %in% unique(c(i2,i4))) # 0 that could be inclusion or exclusion depending on rMATS PSI or my PSI


long.gr<-GRanges(seqnames=A5SS$chr,
               ranges=IRanges(start=A5SS$longExonStart_0base+1,end=A5SS$longExonEnd),
               strand=A5SS$strand)
short.gr<-GRanges(seqnames=A5SS$chr,
               ranges=IRanges(start=A5SS$shortES+1,end=A5SS$shortEE),
               strand=A5SS$strand)
flanking.gr<-GRanges(seqnames=A5SS$chr,
               ranges=IRanges(start=A5SS$flankingES+1,end=A5SS$flankingEE),
               strand=A5SS$strand)
# Get sequence for long, short and flanking exons
long.rna<-getSeq(BSgenome.Hsapiens.UCSC.hg38,long.gr)
short.rna<-getSeq(BSgenome.Hsapiens.UCSC.hg38,short.gr)
flanking.rna<-getSeq(BSgenome.Hsapiens.UCSC.hg38,flanking.gr)

# Inclusion and exclusion isoforms
rna.inc<-paste0(long.rna,flanking.rna)[rownames(A5SS) %in% unique(c(i1,i3))]
rna.exc<-paste0(short.rna,flanking.rna)[rownames(A5SS) %in% unique(c(i2,i4))]

names(rna.inc)<-paste0("AMLvsHealthy_A5SSinc_ID",A5SS$ID[rownames(A5SS) %in% unique(c(i1,i3))])
names(rna.exc)<-paste0("AMLvsHealthy_A5SSexc_ID",A5SS$ID[rownames(A5SS) %in% unique(c(i2,i4))])


length(rna.inc)+length(rna.exc) # 55,700, OK

inc.seq<-lapply(1:3,function(i){substring(rna.inc,i,nchar(rna.inc))})
exc.seq<-lapply(1:3,function(i){substring(rna.exc,i,nchar(rna.exc))})

orf1.inc<-as.character(Biostrings::translate(DNAStringSet(inc.seq[[1]]),
                                             if.fuzzy.codon = "solve",
                                             no.init.codon = TRUE))
orf2.inc<-as.character(Biostrings::translate(DNAStringSet(inc.seq[[2]]),
                                             if.fuzzy.codon = "solve",
                                             no.init.codon = TRUE))
orf3.inc<-as.character(Biostrings::translate(DNAStringSet(inc.seq[[3]]),
                                             if.fuzzy.codon = "solve",
                                             no.init.codon = TRUE))

orf1.exc<-as.character(Biostrings::translate(DNAStringSet(exc.seq[[1]]),
                                             if.fuzzy.codon = "solve",
                                             no.init.codon = TRUE))
orf2.exc<-as.character(Biostrings::translate(DNAStringSet(exc.seq[[2]]),
                                             if.fuzzy.codon = "solve",
                                             no.init.codon = TRUE))
orf3.exc<-as.character(Biostrings::translate(DNAStringSet(exc.seq[[3]]),
                                             if.fuzzy.codon = "solve",
                                             no.init.codon = TRUE))

##### Keep ORFs that are >=7aa #####
# ORF1
orf1.inc.split<-sapply(orf1.inc,function(x){
  a<-unlist(strsplit(x,"\\*"))
  nchar(a[1])>=7
  b<-a[-1]
  b<-b[grep("M",b)]
  pos<-sapply(b,function(x){unlist(gregexpr("M",b))[1]})
  c<-substring(b,pos,nchar(b))
  return(c(a[1],c)[nchar(c(a[1],c))>=7])
})
orf1.inc.stack<-stack(orf1.inc.split)
length(names(orf1.inc.split)[which(!names(orf1.inc.split) %in% orf1.inc.stack$ind)]) # 693 event do not generate sequences >=7aa; drop this levels
orf1.inc.stack<-stack(orf1.inc.split,drop=TRUE)
orf1.inc.stack$ID<-paste0(orf1.inc.stack$ind,"_ORF1.",unlist(sapply(table(orf1.inc.stack$ind),function(x){1:x})))
sum(duplicated(orf1.inc.stack$values)) # 3,060 are duplicated sequences: aggregate
orf1.inc.stack<-aggregate(ID~values,orf1.inc.stack,paste,collapse=";")

# ORF1
orf1.exc.split<-sapply(orf1.exc,function(x){
  a<-unlist(strsplit(x,"\\*"))
  nchar(a[1])>=7
  b<-a[-1]
  b<-b[grep("M",b)]
  pos<-sapply(b,function(x){unlist(gregexpr("M",b))[1]})
  c<-substring(b,pos,nchar(b))
  return(c(a[1],c)[nchar(c(a[1],c))>=7])
})
orf1.exc.stack<-stack(orf1.exc.split)
length(names(orf1.exc.split)[which(!names(orf1.exc.split) %in% orf1.exc.stack$ind)]) # 2,868 events do not generate sequences >=7aa; drop this levels
orf1.exc.stack<-stack(orf1.exc.split,drop=TRUE)
orf1.exc.stack$ID<-paste0(orf1.exc.stack$ind,"_ORF1.",unlist(sapply(table(orf1.exc.stack$ind),function(x){1:x})))
sum(duplicated(orf1.exc.stack$values)) # 15,954 are duplicated sequences: aggregate
orf1.exc.stack<-aggregate(ID~values,orf1.exc.stack,paste,collapse=";")

# ORF2
orf2.inc.split<-sapply(orf2.inc,function(x){
  a<-unlist(strsplit(x,"\\*"))
  nchar(a[1])>=7
  b<-a[-1]
  b<-b[grep("M",b)]
  pos<-sapply(b,function(x){unlist(gregexpr("M",b))[1]})
  c<-substring(b,pos,nchar(b))
  return(c(a[1],c)[nchar(c(a[1],c))>=7])
})
orf2.inc.stack<-stack(orf2.inc.split)
length(names(orf2.inc.split)[which(!names(orf2.inc.split) %in% orf2.inc.stack$ind)]) # 870 events do not generate sequences >=7aa; drop this levels
orf2.inc.stack<-stack(orf2.inc.split,drop=TRUE)
orf2.inc.stack$ID<-paste0(orf2.inc.stack$ind,"_ORF2.",unlist(sapply(table(orf2.inc.stack$ind),function(x){1:x})))
sum(duplicated(orf2.inc.stack$values)) # 3,178 are duplicated sequences: aggregate
orf2.inc.stack<-aggregate(ID~values,orf2.inc.stack,paste,collapse=";")

# ORF2
orf2.exc.split<-sapply(orf2.exc,function(x){
  a<-unlist(strsplit(x,"\\*"))
  nchar(a[1])>=7
  b<-a[-1]
  b<-b[grep("M",b)]
  pos<-sapply(b,function(x){unlist(gregexpr("M",b))[1]})
  c<-substring(b,pos,nchar(b))
  return(c(a[1],c)[nchar(c(a[1],c))>=7])
})
orf2.exc.stack<-stack(orf2.exc.split)
length(names(orf2.exc.split)[which(!names(orf2.exc.split) %in% orf2.exc.stack$ind)]) # 3,908 events do not generate sequences >=7aa; drop this levels
orf2.exc.stack<-stack(orf2.exc.split,drop=TRUE)
orf2.exc.stack$ID<-paste0(orf2.exc.stack$ind,"_ORF2.",unlist(sapply(table(orf2.exc.stack$ind),function(x){1:x})))
sum(duplicated(orf2.exc.stack$values)) # 16,409 are duplicated sequences: aggregate
orf2.exc.stack<-aggregate(ID~values,orf2.exc.stack,paste,collapse=";")

# ORF3
orf3.inc.split<-sapply(orf3.inc,function(x){
  a<-unlist(strsplit(x,"\\*"))
  nchar(a[1])>=7
  b<-a[-1]
  b<-b[grep("M",b)]
  pos<-sapply(b,function(x){unlist(gregexpr("M",b))[1]})
  c<-substring(b,pos,nchar(b))
  return(c(a[1],c)[nchar(c(a[1],c))>=7])
})
orf3.inc.stack<-stack(orf3.inc.split)
length(names(orf3.inc.split)[which(!names(orf3.inc.split) %in% orf3.inc.stack$ind)]) # 899 events do not generate sequences >=7aa; drop this levels
orf3.inc.stack<-stack(orf3.inc.split,drop=TRUE)
orf3.inc.stack$ID<-paste0(orf3.inc.stack$ind,"_ORF3.",unlist(sapply(table(orf3.inc.stack$ind),function(x){1:x})))
sum(duplicated(orf3.inc.stack$values)) # 3,032 are duplicated sequences: aggregate
orf3.inc.stack<-aggregate(ID~values,orf3.inc.stack,paste,collapse=";")

# ORF3
orf3.exc.split<-sapply(orf3.exc,function(x){
  a<-unlist(strsplit(x,"\\*"))
  nchar(a[1])>=7
  b<-a[-1]
  b<-b[grep("M",b)]
  pos<-sapply(b,function(x){unlist(gregexpr("M",b))[1]})
  c<-substring(b,pos,nchar(b))
  return(c(a[1],c)[nchar(c(a[1],c))>=7])
})
orf3.exc.stack<-stack(orf3.exc.split)
length(names(orf3.exc.split)[which(!names(orf3.exc.split) %in% orf3.exc.stack$ind)]) # 3,810 events do not generate sequences >=7aa; drop this levels
orf3.exc.stack<-stack(orf3.exc.split,drop=TRUE)
orf3.exc.stack$ID<-paste0(orf3.exc.stack$ind,"_ORF3.",unlist(sapply(table(orf3.exc.stack$ind),function(x){1:x})))
sum(duplicated(orf3.exc.stack$values)) # 16,080 are duplicated sequences: aggregate
orf3.exc.stack<-aggregate(ID~values,orf3.exc.stack,paste,collapse=";")

# Rbind and aggregate duplicated sequences
A5SS.fa<-rbind(orf1.inc.stack,orf1.exc.stack,orf2.inc.stack,orf2.exc.stack,orf3.inc.stack,orf3.exc.stack)
A5SS.fa<-aggregate(ID~values,A5SS.fa,paste,collapse=";") # 219,973
A5SS.fa<-data.frame(ID=A5SS.fa$ID,sequence=A5SS.fa$values)

##### Write .tsv with unique sequences #####
write.table(x=A5SS.fa,file="Fragpipe/inputs/BeatAMLvsHealthy_A5SS.tsv",row.names = FALSE,quote=FALSE,col.names=TRUE,sep="\t")



################################# A3SS #########################################
A3SS<-read.delim("A3SS_annotated.txt")
i1<-which(A3SS$FDR<0.05 & A3SS$IncLevelDifference>0) # 3,073 inc (which returns row position)
i2<-which(A3SS$FDR<0.05 & A3SS$IncLevelDifference<0) # 5,777 exc

i3<-which(A3SS$all.0.healthy==19 & A3SS$gt.0.aml>=1) # 17,340 (inc)
i4<-which(A3SS$all.1.healthy==19 & A3SS$lt.1.aml>=1) # 94,428 (exc)

length(unique(c(i1,i2,i3,i4))) # 120,419
A3SS<-A3SS[unique(c(i1,i2,i3,i4)),] # rownames are row original positions
sum(unique(c(i1,i3)) %in% unique(c(i2,i4))) # 0 that could be inclusion or exclusion depending on rMATS PSI or my PSI


long.gr<-GRanges(seqnames=A3SS$chr,
                 ranges=IRanges(start=A3SS$longExonStart_0base+1,end=A3SS$longExonEnd),
                 strand=A3SS$strand)
short.gr<-GRanges(seqnames=A3SS$chr,
                  ranges=IRanges(start=A3SS$shortES+1,end=A3SS$shortEE),
                  strand=A3SS$strand)
flanking.gr<-GRanges(seqnames=A3SS$chr,
                     ranges=IRanges(start=A3SS$flankingES+1,end=A3SS$flankingEE),
                     strand=A3SS$strand)

# Get sequence for long, short and flanking exons
long.rna<-getSeq(BSgenome.Hsapiens.UCSC.hg38,long.gr)
short.rna<-getSeq(BSgenome.Hsapiens.UCSC.hg38,short.gr)
flanking.rna<-getSeq(BSgenome.Hsapiens.UCSC.hg38,flanking.gr)

# Inclusion and exclusion isoforms
rna.inc<-paste0(flanking.rna,long.rna)[rownames(A3SS) %in% unique(c(i1,i3))]
rna.exc<-paste0(flanking.rna,short.rna)[rownames(A3SS) %in% unique(c(i2,i4))]

names(rna.inc)<-paste0("AMLvsHealthy_A3SSinc_ID",A3SS$ID[rownames(A3SS) %in% unique(c(i1,i3))])
names(rna.exc)<-paste0("AMLvsHealthy_A3SSexc_ID",A3SS$ID[rownames(A3SS) %in% unique(c(i2,i4))])


length(rna.inc)+length(rna.exc) # 120,419, OK

inc.seq<-lapply(1:3,function(i){substring(rna.inc,i,nchar(rna.inc))})
exc.seq<-lapply(1:3,function(i){substring(rna.exc,i,nchar(rna.exc))})

orf1.inc<-as.character(Biostrings::translate(DNAStringSet(inc.seq[[1]]),
                                             if.fuzzy.codon = "solve",
                                             no.init.codon = TRUE))
orf2.inc<-as.character(Biostrings::translate(DNAStringSet(inc.seq[[2]]),
                                             if.fuzzy.codon = "solve",
                                             no.init.codon = TRUE))
orf3.inc<-as.character(Biostrings::translate(DNAStringSet(inc.seq[[3]]),
                                             if.fuzzy.codon = "solve",
                                             no.init.codon = TRUE))

orf1.exc<-as.character(Biostrings::translate(DNAStringSet(exc.seq[[1]]),
                                             if.fuzzy.codon = "solve",
                                             no.init.codon = TRUE))
orf2.exc<-as.character(Biostrings::translate(DNAStringSet(exc.seq[[2]]),
                                             if.fuzzy.codon = "solve",
                                             no.init.codon = TRUE))
orf3.exc<-as.character(Biostrings::translate(DNAStringSet(exc.seq[[3]]),
                                             if.fuzzy.codon = "solve",
                                             no.init.codon = TRUE))

##### Keep ORFs that are >=7aa #####
# ORF1
orf1.inc.split<-sapply(orf1.inc,function(x){
  a<-unlist(strsplit(x,"\\*"))
  nchar(a[1])>=7
  b<-a[-1]
  b<-b[grep("M",b)]
  pos<-sapply(b,function(x){unlist(gregexpr("M",b))[1]})
  c<-substring(b,pos,nchar(b))
  return(c(a[1],c)[nchar(c(a[1],c))>=7])
})
orf1.inc.stack<-stack(orf1.inc.split)
length(names(orf1.inc.split)[which(!names(orf1.inc.split) %in% orf1.inc.stack$ind)]) # 907 events do not generate sequences >=7aa; drop this levels
orf1.inc.stack<-stack(orf1.inc.split,drop=TRUE)
orf1.inc.stack$ID<-paste0(orf1.inc.stack$ind,"_ORF1.",unlist(sapply(table(orf1.inc.stack$ind),function(x){1:x})))
sum(duplicated(orf1.inc.stack$values)) # 3,812 are duplicated sequences: aggregate
orf1.inc.stack<-aggregate(ID~values,orf1.inc.stack,paste,collapse=";")

# ORF1
orf1.exc.split<-sapply(orf1.exc,function(x){
  a<-unlist(strsplit(x,"\\*"))
  nchar(a[1])>=7
  b<-a[-1]
  b<-b[grep("M",b)]
  pos<-sapply(b,function(x){unlist(gregexpr("M",b))[1]})
  c<-substring(b,pos,nchar(b))
  return(c(a[1],c)[nchar(c(a[1],c))>=7])
})
orf1.exc.stack<-stack(orf1.exc.split)
length(names(orf1.exc.split)[which(!names(orf1.exc.split) %in% orf1.exc.stack$ind)]) # 7,408 events do not generate sequences >=7aa; drop this levels
orf1.exc.stack<-stack(orf1.exc.split,drop=TRUE)
orf1.exc.stack$ID<-paste0(orf1.exc.stack$ind,"_ORF1.",unlist(sapply(table(orf1.exc.stack$ind),function(x){1:x})))
sum(duplicated(orf1.exc.stack$values)) # 54,728 are duplicated sequences: aggregate
orf1.exc.stack<-aggregate(ID~values,orf1.exc.stack,paste,collapse=";")



# ORF2
orf2.inc.split<-sapply(orf2.inc,function(x){
  a<-unlist(strsplit(x,"\\*"))
  nchar(a[1])>=7
  b<-a[-1]
  b<-b[grep("M",b)]
  pos<-sapply(b,function(x){unlist(gregexpr("M",b))[1]})
  c<-substring(b,pos,nchar(b))
  return(c(a[1],c)[nchar(c(a[1],c))>=7])
})
orf2.inc.stack<-stack(orf2.inc.split)
length(names(orf2.inc.split)[which(!names(orf2.inc.split) %in% orf2.inc.stack$ind)]) # 1,393 events do not generate sequences >=7aa; drop this levels
orf2.inc.stack<-stack(orf2.inc.split,drop=TRUE)
orf2.inc.stack$ID<-paste0(orf2.inc.stack$ind,"_ORF2.",unlist(sapply(table(orf2.inc.stack$ind),function(x){1:x})))
sum(duplicated(orf2.inc.stack$values)) # 4,368 are duplicated sequences: aggregate
orf2.inc.stack<-aggregate(ID~values,orf2.inc.stack,paste,collapse=";")


# ORF2
orf2.exc.split<-sapply(orf2.exc,function(x){
  a<-unlist(strsplit(x,"\\*"))
  nchar(a[1])>=7
  b<-a[-1]
  b<-b[grep("M",b)]
  pos<-sapply(b,function(x){unlist(gregexpr("M",b))[1]})
  c<-substring(b,pos,nchar(b))
  return(c(a[1],c)[nchar(c(a[1],c))>=7])
})
orf2.exc.stack<-stack(orf2.exc.split)
length(names(orf2.exc.split)[which(!names(orf2.exc.split) %in% orf2.exc.stack$ind)]) # 10,862 events do not generate sequences >=7aa; drop this levels
orf2.exc.stack<-stack(orf2.exc.split,drop=TRUE)
orf2.exc.stack$ID<-paste0(orf2.exc.stack$ind,"_ORF2.",unlist(sapply(table(orf2.exc.stack$ind),function(x){1:x})))
sum(duplicated(orf2.exc.stack$values)) # 56,679 are duplicated sequences: aggregate
orf2.exc.stack<-aggregate(ID~values,orf2.exc.stack,paste,collapse=";")

# ORF3
orf3.inc.split<-sapply(orf3.inc,function(x){
  a<-unlist(strsplit(x,"\\*"))
  nchar(a[1])>=7
  b<-a[-1]
  b<-b[grep("M",b)]
  pos<-sapply(b,function(x){unlist(gregexpr("M",b))[1]})
  c<-substring(b,pos,nchar(b))
  return(c(a[1],c)[nchar(c(a[1],c))>=7])
})
orf3.inc.stack<-stack(orf3.inc.split)
length(names(orf3.inc.split)[which(!names(orf3.inc.split) %in% orf3.inc.stack$ind)]) # 1,516 events do not generate sequences >=7aa; drop this levels
orf3.inc.stack<-stack(orf3.inc.split,drop=TRUE)
orf3.inc.stack$ID<-paste0(orf3.inc.stack$ind,"_ORF3.",unlist(sapply(table(orf3.inc.stack$ind),function(x){1:x})))
sum(duplicated(orf3.inc.stack$values)) # 4,009 are duplicated sequences: aggregate
orf3.inc.stack<-aggregate(ID~values,orf3.inc.stack,paste,collapse=";")


# ORF3
orf3.exc.split<-sapply(orf3.exc,function(x){
  a<-unlist(strsplit(x,"\\*"))
  nchar(a[1])>=7
  b<-a[-1]
  b<-b[grep("M",b)]
  pos<-sapply(b,function(x){unlist(gregexpr("M",b))[1]})
  c<-substring(b,pos,nchar(b))
  return(c(a[1],c)[nchar(c(a[1],c))>=7])
})
orf3.exc.stack<-stack(orf3.exc.split)
length(names(orf3.exc.split)[which(!names(orf3.exc.split) %in% orf3.exc.stack$ind)]) # 11,830 events do not generate sequences >=7aa; drop this levels
orf3.exc.stack<-stack(orf3.exc.split,drop=TRUE)
orf3.exc.stack$ID<-paste0(orf3.exc.stack$ind,"-ORF3.",unlist(sapply(table(orf3.exc.stack$ind),function(x){1:x})))
sum(duplicated(orf3.exc.stack$values)) # 53,238 are duplicated sequences: aggregate
orf3.exc.stack<-aggregate(ID~values,orf3.exc.stack,paste,collapse=";")

# Rbind and aggregate duplicated sequences
A3SS.fa<-rbind(orf1.inc.stack,orf1.exc.stack,orf2.inc.stack,orf2.exc.stack,orf3.inc.stack,orf3.exc.stack)
A3SS.fa<-aggregate(ID~values,A3SS.fa,paste,collapse=";") # 296,819
A3SS.fa<-data.frame(ID=A3SS.fa$ID,sequence=A3SS.fa$values)


##### Write .tsv with unique sequences #####
write.table(x=A3SS.fa,file="Fragpipe/inputs/BeatAMLvsHealthy_A3SS.tsv",row.names = FALSE,quote=FALSE,col.names=TRUE,sep="\t")


################################## RI ##########################################
RI<-read.delim("RI_annotated.txt")
# Select events
i1<-which(RI$FDR<0.05 & RI$IncLevelDifference>0) # 2,016 inc (which returns row position)
i2<-which(RI$FDR<0.05 & RI$IncLevelDifference<0) # 5,597 exc

i3<-which(RI$all.0.healthy==19 & RI$gt.0.aml>=1) # 43 (inc)
i4<-which(RI$all.1.healthy==19 & RI$lt.1.aml>=1) # 11,537 (exc)

length(unique(c(i1,i2,i3,i4))) # 19,138
RI<-RI[unique(c(i1,i2,i3,i4)),] # rownames are row original positions
sum(unique(c(i1,i3)) %in% unique(c(i2,i4))) # 0 that could be inclusion or exclusion depending on rMATS PSI or my PSI

# Get sequence for upstream exon, downstream exon and retained intron
UE.gr<-GRanges(seqnames=RI$chr,
                 ranges=IRanges(start=RI$upstreamES+1,end=RI$upstreamEE),
                 strand=RI$strand)
DE.gr<-GRanges(seqnames=RI$chr,
                  ranges=IRanges(start=RI$downstreamES+1,end=RI$downstreamEE),
                  strand=RI$strand)
RI.gr<-GRanges(seqnames=RI$chr,
                     ranges=IRanges(start=RI$riExonStart_0base+1,end=RI$riExonEnd),
                     strand=RI$strand)


UE.rna<-getSeq(BSgenome.Hsapiens.UCSC.hg38,UE.gr)
DE.rna<-getSeq(BSgenome.Hsapiens.UCSC.hg38,DE.gr)
RI.rna<-getSeq(BSgenome.Hsapiens.UCSC.hg38,RI.gr)

# Inclusion and exclusion isoforms
rna.inc<-RI.rna[rownames(RI) %in% unique(c(i1,i3))] # RI goes from UEstart to DEend
rna.exc<-ifelse(RI$strand[rownames(RI) %in% unique(c(i2,i4))]=="+",
                paste0(UE.rna,DE.rna)[rownames(RI) %in% unique(c(i2,i4))],
                paste0(DE.rna,UE.rna)[rownames(RI) %in% unique(c(i2,i4))])

names(rna.inc)<-paste0("AMLvsHealthy_RIinc_ID",RI$ID[rownames(RI) %in% unique(c(i1,i3))])
names(rna.exc)<-paste0("AMLvsHealthy_RIexc_ID",RI$ID[rownames(RI) %in% unique(c(i2,i4))])


length(rna.inc)+length(rna.exc) # 19,138, OK

inc.seq<-lapply(1:3,function(i){substring(rna.inc,i,nchar(rna.inc))})
exc.seq<-lapply(1:3,function(i){substring(rna.exc,i,nchar(rna.exc))})

orf1.inc<-as.character(Biostrings::translate(DNAStringSet(inc.seq[[1]]),
                                             if.fuzzy.codon = "solve",
                                             no.init.codon = TRUE))
orf2.inc<-as.character(Biostrings::translate(DNAStringSet(inc.seq[[2]]),
                                             if.fuzzy.codon = "solve",
                                             no.init.codon = TRUE))
orf3.inc<-as.character(Biostrings::translate(DNAStringSet(inc.seq[[3]]),
                                             if.fuzzy.codon = "solve",
                                             no.init.codon = TRUE))

orf1.exc<-as.character(Biostrings::translate(DNAStringSet(exc.seq[[1]]),
                                             if.fuzzy.codon = "solve",
                                             no.init.codon = TRUE))
orf2.exc<-as.character(Biostrings::translate(DNAStringSet(exc.seq[[2]]),
                                             if.fuzzy.codon = "solve",
                                             no.init.codon = TRUE))
orf3.exc<-as.character(Biostrings::translate(DNAStringSet(exc.seq[[3]]),
                                             if.fuzzy.codon = "solve",
                                             no.init.codon = TRUE))

##### Keep ORFs that are >=7aa #####
# ORF1
orf1.inc.split<-sapply(orf1.inc,function(x){
  a<-unlist(strsplit(x,"\\*"))
  nchar(a[1])>=7
  b<-a[-1]
  b<-b[grep("M",b)]
  pos<-sapply(b,function(x){unlist(gregexpr("M",b))[1]})
  c<-substring(b,pos,nchar(b))
  return(c(a[1],c)[nchar(c(a[1],c))>=7])
})
orf1.inc.stack<-stack(orf1.inc.split)
length(names(orf1.inc.split)[which(!names(orf1.inc.split) %in% orf1.inc.stack$ind)]) # 99 events do not generate sequences >=7aa; drop this levels
orf1.inc.stack<-stack(orf1.inc.split,drop=TRUE)
orf1.inc.stack$ID<-paste0(orf1.inc.stack$ind,"_ORF1.",unlist(sapply(table(orf1.inc.stack$ind),function(x){1:x})))
sum(duplicated(orf1.inc.stack$values)) # 211 are duplicated sequences: aggregate
orf1.inc.stack<-aggregate(ID~values,orf1.inc.stack,paste,collapse=";")


# ORF1
orf1.exc.split<-sapply(orf1.exc,function(x){
  a<-unlist(strsplit(x,"\\*"))
  nchar(a[1])>=7
  b<-a[-1]
  b<-b[grep("M",b)]
  pos<-sapply(b,function(x){unlist(gregexpr("M",b))[1]})
  c<-substring(b,pos,nchar(b))
  return(c(a[1],c)[nchar(c(a[1],c))>=7])
})
orf1.exc.stack<-stack(orf1.exc.split)
length(names(orf1.exc.split)[which(!names(orf1.exc.split) %in% orf1.exc.stack$ind)]) # 1,364 events do not generate sequences >=7aa; drop this levels
orf1.exc.stack<-stack(orf1.exc.split,drop=TRUE)
orf1.exc.stack$ID<-paste0(orf1.exc.stack$ind,"_ORF1.",unlist(sapply(table(orf1.exc.stack$ind),function(x){1:x})))
sum(duplicated(orf1.exc.stack$values)) # 2,851 are duplicated sequences: aggregate
orf1.exc.stack<-aggregate(ID~values,orf1.exc.stack,paste,collapse=";")

# ORF2
orf2.inc.split<-sapply(orf2.inc,function(x){
  a<-unlist(strsplit(x,"\\*"))
  nchar(a[1])>=7
  b<-a[-1]
  b<-b[grep("M",b)]
  pos<-sapply(b,function(x){unlist(gregexpr("M",b))[1]})
  c<-substring(b,pos,nchar(b))
  return(c(a[1],c)[nchar(c(a[1],c))>=7])
})
orf2.inc.stack<-stack(orf2.inc.split)
length(names(orf2.inc.split)[which(!names(orf2.inc.split) %in% orf2.inc.stack$ind)]) # 118 event do not generate sequences >=7aa; drop this levels
orf2.inc.stack<-stack(orf2.inc.split,drop=TRUE)
orf2.inc.stack$ID<-paste0(orf2.inc.stack$ind,"_ORF2.",unlist(sapply(table(orf2.inc.stack$ind),function(x){1:x})))
sum(duplicated(orf2.inc.stack$values)) # 229 are duplicated sequences: aggregate
orf2.inc.stack<-aggregate(ID~values,orf2.inc.stack,paste,collapse=";")

# ORF2
orf2.exc.split<-sapply(orf2.exc,function(x){
  a<-unlist(strsplit(x,"\\*"))
  nchar(a[1])>=7
  b<-a[-1]
  b<-b[grep("M",b)]
  pos<-sapply(b,function(x){unlist(gregexpr("M",b))[1]})
  c<-substring(b,pos,nchar(b))
  return(c(a[1],c)[nchar(c(a[1],c))>=7])
})
orf2.exc.stack<-stack(orf2.exc.split)
length(names(orf2.exc.split)[which(!names(orf2.exc.split) %in% orf2.exc.stack$ind)]) # 1,780 event do not generate sequences >=7aa; drop this levels
orf2.exc.stack<-stack(orf2.exc.split,drop=TRUE)
orf2.exc.stack$ID<-paste0(orf2.exc.stack$ind,"_ORF2.",unlist(sapply(table(orf2.exc.stack$ind),function(x){1:x})))
sum(duplicated(orf2.exc.stack$values)) # 3,278 are duplicated sequences: aggregate
orf2.exc.stack<-aggregate(ID~values,orf2.exc.stack,paste,collapse=";")

# ORF3
orf3.inc.split<-sapply(orf3.inc,function(x){
  a<-unlist(strsplit(x,"\\*"))
  nchar(a[1])>=7
  b<-a[-1]
  b<-b[grep("M",b)]
  pos<-sapply(b,function(x){unlist(gregexpr("M",b))[1]})
  c<-substring(b,pos,nchar(b))
  return(c(a[1],c)[nchar(c(a[1],c))>=7])
})
orf3.inc.stack<-stack(orf3.inc.split)
length(names(orf3.inc.split)[which(!names(orf3.inc.split) %in% orf3.inc.stack$ind)]) # 131 event do not generate sequences >=7aa; drop this levels
orf3.inc.stack<-stack(orf3.inc.split,drop=TRUE)
orf3.inc.stack$ID<-paste0(orf3.inc.stack$ind,"_ORF3.",unlist(sapply(table(orf3.inc.stack$ind),function(x){1:x})))
sum(duplicated(orf3.inc.stack$values)) # 201 are duplicated sequences: aggregate
orf3.inc.stack<-aggregate(ID~values,orf3.inc.stack,paste,collapse=";")

# ORF3
orf3.exc.split<-sapply(orf3.exc,function(x){
  a<-unlist(strsplit(x,"\\*"))
  nchar(a[1])>=7
  b<-a[-1]
  b<-b[grep("M",b)]
  pos<-sapply(b,function(x){unlist(gregexpr("M",b))[1]})
  c<-substring(b,pos,nchar(b))
  return(c(a[1],c)[nchar(c(a[1],c))>=7])
})
orf3.exc.stack<-stack(orf3.exc.split)
length(names(orf3.exc.split)[which(!names(orf3.exc.split) %in% orf3.exc.stack$ind)]) # 1,917 event do not generate sequences >=7aa; drop this levels
orf3.exc.stack<-stack(orf3.exc.split,drop=TRUE)
orf3.exc.stack$ID<-paste0(orf3.exc.stack$ind,"_ORF3.",unlist(sapply(table(orf3.exc.stack$ind),function(x){1:x})))
sum(duplicated(orf3.exc.stack$values)) # 2,928 are duplicated sequences: aggregate
orf3.exc.stack<-aggregate(ID~values,orf3.exc.stack,paste,collapse=";")

RI.fa<-rbind(orf1.inc.stack,orf1.exc.stack,orf2.inc.stack,orf2.exc.stack,orf3.inc.stack,orf3.exc.stack)
RI.fa<-aggregate(ID~values,RI.fa,paste,collapse=";") # 62,991
RI.fa<-data.frame(ID=RI.fa$ID,sequence=RI.fa$values)


##### Write .tsv with unique sequences #####
write.table(x=RI.fa,file="Fragpipe/inputs/BeatAMLvsHealthy_RI.tsv",row.names = FALSE,quote=FALSE,col.names=TRUE,sep="\t")



############################ BEAT-AML Arriba ###################################
# Fusions for the whole Beat AML cohort identified with Arriba
# File list (tsv files)
# Arriba pipeline was run in Afrodita
file.list<-list.files("arriba/files/",full.names=TRUE) # Directory with all .tsv files

# Read all the files  
df.list <- lapply(file.list, fread)

# Collapse all files
dt <- rbindlist(df.list,idcol="id") # Number 269 is empty
dt$fusion<-paste(dt$'#gene1',
                 dt$breakpoint1,
                 sub(".*\\/","",dt$`strand1(gene/fusion)`),
                 dt$gene2,
                 dt$breakpoint2,
                 sub(".*\\/","",dt$`strand1(gene/fusion)`),
                 sep="|")

# Unique fusions and sequence
fusions.unique<-dt[!duplicated(dt$fusion),]
length(unique(dt$fusion))==nrow(fusions.unique)
fusions.unique<-fusions.unique[!grepl("\\.",fusions.unique$`strand1(gene/fusion)`) &
                                 !grepl("\\.",fusions.unique$`strand2(gene/fusion)`),] # remove fusions with unknown translation strand
table(fusions.unique$`strand1(gene/fusion)`,fusions.unique$`strand2(gene/fusion)`)

# Translate 7 aa = 6 aa +1 aa
UE.gr<-GRanges(seqnames=sub(":.*","",fusions.unique$breakpoint1),
               ranges=IRanges(start=ifelse(sub(".*\\/","",fusions.unique$`strand1(gene/fusion)`)=="+",
                                           as.numeric(sub(".*:","",fusions.unique$breakpoint1))-17,
                                           as.numeric(sub(".*:","",fusions.unique$breakpoint1))),
                              end=ifelse(sub(".*\\/","",fusions.unique$`strand1(gene/fusion)`)=="+",
                                         as.numeric(sub(".*:","",fusions.unique$breakpoint1)),
                                         as.numeric(sub(".*:","",fusions.unique$breakpoint1))+17)),
               strand=sub(".*\\/","",fusions.unique$`strand1(gene/fusion)`))

DE.gr<-GRanges(seqnames=sub(":.*","",fusions.unique$breakpoint2),
               ranges=IRanges(start=ifelse(sub(".*\\/","",fusions.unique$`strand2(gene/fusion)`)=="+",
                                           as.numeric(sub(".*:","",fusions.unique$breakpoint2)),
                                           as.numeric(sub(".*:","",fusions.unique$breakpoint2))-17),
                              end=ifelse(sub(".*\\/","",fusions.unique$`strand2(gene/fusion)`)=="+",
                                         as.numeric(sub(".*:","",fusions.unique$breakpoint2))+17,
                                         as.numeric(sub(".*:","",fusions.unique$breakpoint2)))),
               strand=sub(".*\\/","",fusions.unique$`strand2(gene/fusion)`))

UE.rna<-getSeq(BSgenome.Hsapiens.UCSC.hg38,UE.gr)
DE.rna<-getSeq(BSgenome.Hsapiens.UCSC.hg38,DE.gr)
target.seq<-paste0(UE.rna,DE.rna)
names(target.seq)<-paste0("ARRIBA_",fusions.unique$fusion)

# There are no fusions involving chrM
# Get 3 reading frames

fusion.seq<-lapply(1:3,function(i){substring(target.seq,i,nchar(target.seq))})
orf1<-as.character(Biostrings::translate(DNAStringSet(fusion.seq[[1]]),
                                         if.fuzzy.codon = "solve",
                                         no.init.codon = TRUE))
orf2<-as.character(Biostrings::translate(DNAStringSet(fusion.seq[[2]]),
                                         if.fuzzy.codon = "solve",
                                         no.init.codon = TRUE))
orf3<-as.character(Biostrings::translate(DNAStringSet(fusion.seq[[3]]),
                                         if.fuzzy.codon = "solve",
                                         no.init.codon = TRUE))

##### Keep ORFs that are >=7aa #####
# ORF1
orf1.split<-sapply(orf1,function(x){
  a<-unlist(strsplit(x,"\\*"))
  nchar(a[1])>=7
  b<-a[-1]
  b<-b[grep("M",b)]
  pos<-sapply(b,function(x){unlist(gregexpr("M",b))[1]})
  c<-substring(b,pos,nchar(b))
  return(c(a[1],c)[nchar(c(a[1],c))>=7])
})
orf1.stack<-stack(orf1.split)
length(names(orf1.split)[which(!names(orf1.split) %in% orf1.stack$ind)]) # 1,116 fusions do not generate sequences >=7aa; drop this levels
orf1.stack<-stack(orf1.split,drop=TRUE) # drop names that are not used
sum(duplicated(orf1.stack$ind)) # 0
orf1.stack$ID<-paste0(orf1.stack$ind,"_ORF1")
sum(duplicated(orf1.stack$values)) # 426 are duplicated sequences: aggregate
orf1.stack<-aggregate(ID~values,orf1.stack,paste,collapse=";")

# ORF2
orf2.split<-sapply(orf2,function(x){
  a<-unlist(strsplit(x,"\\*"))
  nchar(a[1])>=7
  b<-a[-1]
  b<-b[grep("M",b)]
  pos<-sapply(b,function(x){unlist(gregexpr("M",b))[1]})
  c<-substring(b,pos,nchar(b))
  return(c(a[1],c)[nchar(c(a[1],c))>=7])
})
orf2.stack<-stack(orf2.split)
length(names(orf2.split)[which(!names(orf2.split) %in% orf2.stack$ind)]) # 1,573 fusions do not generate sequences >=7aa; drop this levels
orf2.stack<-stack(orf2.split,drop=TRUE) # drop names that are not used
sum(duplicated(orf2.stack$ind)) # 0
orf2.stack$ID<-paste0(orf2.stack$ind,"_ORF2")
sum(duplicated(orf2.stack$values)) # 351 are duplicated sequences: aggregate
orf2.stack<-aggregate(ID~values,orf2.stack,paste,collapse=";")

# ORF3
orf3.split<-sapply(orf3,function(x){
  a<-unlist(strsplit(x,"\\*"))
  nchar(a[1])>=7
  b<-a[-1]
  b<-b[grep("M",b)]
  pos<-sapply(b,function(x){unlist(gregexpr("M",b))[1]})
  c<-substring(b,pos,nchar(b))
  return(c(a[1],c)[nchar(c(a[1],c))>=7])
})
orf3.stack<-stack(orf3.split)
length(names(orf3.split)[which(!names(orf3.split) %in% orf3.stack$ind)]) # 1,539 fusions do not generate sequences >=7aa; drop this levels
orf3.stack<-stack(orf3.split,drop=TRUE) # drop names that are not used
sum(duplicated(orf3.stack$ind)) # 0
orf3.stack$ID<-paste0(orf3.stack$ind,"_ORF3")
sum(duplicated(orf3.stack$values)) # 355 are duplicated sequences: aggregate
orf3.stack<-aggregate(ID~values,orf3.stack,paste,collapse=";")

fusions.fa<-rbind(orf1.stack,orf2.stack,orf3.stack)
fusions.fa<-aggregate(ID~values,fusions.fa,paste,collapse=";") # 15,351
fusions.fa<-data.frame(ID=fusions.fa$ID,sequence=fusions.fa$values)

##### Write .tsv with unique sequences #####
write.table(x=fusions.fa,file="Fragpipe/inputs/BeatAML_WholeCohort_arriba.tsv",row.names = FALSE,quote=FALSE,col.names=TRUE,sep="\t")

############################## APL rMATS #######################################
# Inputs are events that were annotated in APL vs AML and APL vs healthy
##### SE #####
SE<-read.delim("SE_common.txt",sep="\t")
# Select events
i1<-which(SE$FDR_APLvsAML<0.05 & SE$FDR_APLvsHealthy<0.05 &
             SE$dPSI_APLvsAML>0 & SE$dPSI_APLvsHealthy>0) # inc
i2<-which(SE$FDR_APLvsAML<0.05 & SE$FDR_APLvsHealthy<0.05 &
             SE$dPSI_APLvsAML<0 & SE$dPSI_APLvsHealthy<0) # exc

i3<-which(SE$all.0.healthy==19 & SE$all.0.aml==122 & SE$gt.0.apl>=1) # inc
i4<-which(SE$all.1.healthy==19 & SE$all.1.aml==122 & SE$lt.1.apl>=1) # exc

SE<-SE[unique(c(i1,i2,i3,i4)),]

UE.gr<-GRanges(seqnames=SE$chr,
               ranges=IRanges(start=SE$upstreamES+1,end=SE$upstreamEE),
               strand=SE$strand)
SE.gr<-GRanges(seqnames=SE$chr,
               ranges=IRanges(start=SE$exonStart_0base+1,end=SE$exonEnd),
               strand=SE$strand)
DE.gr<-GRanges(seqnames=SE$chr,
               ranges=IRanges(start=SE$downstreamES+1,end=SE$downstreamEE),
               strand=SE$strand)

UE.rna<-getSeq(BSgenome.Hsapiens.UCSC.hg38,UE.gr)
SE.rna<-getSeq(BSgenome.Hsapiens.UCSC.hg38,SE.gr)
DE.rna<-getSeq(BSgenome.Hsapiens.UCSC.hg38,DE.gr)

rna.inc<-ifelse(SE$strand[rownames(SE) %in% unique(c(i1,i3))]=="+",
                paste0(UE.rna,SE.rna,DE.rna)[rownames(SE) %in% unique(c(i1,i3))],
                paste0(DE.rna,SE.rna,UE.rna)[rownames(SE) %in% unique(c(i1,i3))])
rna.exc<-ifelse(SE$strand[rownames(SE) %in% unique(c(i2,i4))]=="+",
                paste0(UE.rna,DE.rna)[rownames(SE) %in% unique(c(i2,i4))],
                paste0(DE.rna,UE.rna)[rownames(SE) %in% unique(c(i2,i4))])

names(rna.inc)<-paste0("APL_SEinc_ID",SE$ID[rownames(SE) %in% unique(c(i1,i3))])
names(rna.exc)<-paste0("APL_SEexc_ID",SE$ID[rownames(SE) %in% unique(c(i2,i4))])


length(rna.inc)+length(rna.exc) # 10,549

inc.seq<-lapply(1:3,function(i){substring(rna.inc,i,nchar(rna.inc))})
exc.seq<-lapply(1:3,function(i){substring(rna.exc,i,nchar(rna.exc))})

# Translate. There are no events in chrM
orf1.inc<-as.character(Biostrings::translate(DNAStringSet(inc.seq[[1]]),
                                             if.fuzzy.codon = "solve",
                                             no.init.codon = TRUE))
orf2.inc<-as.character(Biostrings::translate(DNAStringSet(inc.seq[[2]]),
                                             if.fuzzy.codon = "solve",
                                             no.init.codon = TRUE))
orf3.inc<-as.character(Biostrings::translate(DNAStringSet(inc.seq[[3]]),
                                             if.fuzzy.codon = "solve",
                                             no.init.codon = TRUE))

orf1.exc<-as.character(Biostrings::translate(DNAStringSet(exc.seq[[1]]),
                                             if.fuzzy.codon = "solve",
                                             no.init.codon = TRUE))
orf2.exc<-as.character(Biostrings::translate(DNAStringSet(exc.seq[[2]]),
                                             if.fuzzy.codon = "solve",
                                             no.init.codon = TRUE))
orf3.exc<-as.character(Biostrings::translate(DNAStringSet(exc.seq[[3]]),
                                             if.fuzzy.codon = "solve",
                                             no.init.codon = TRUE))

##### Keep ORFs that are >=7aa #####
# ORF1
orf1.inc.split<-sapply(orf1.inc,function(x){
  a<-unlist(strsplit(x,"\\*"))
  nchar(a[1])>=7
  b<-a[-1]
  b<-b[grep("M",b)]
  pos<-sapply(b,function(x){unlist(gregexpr("M",b))[1]})
  c<-substring(b,pos,nchar(b))
  return(c(a[1],c)[nchar(c(a[1],c))>=7])
})
orf1.inc.stack<-stack(orf1.inc.split)
length(names(orf1.inc.split)[which(!names(orf1.inc.split) %in% orf1.inc.stack$ind)]) # 129 events do not generate sequences >=7aa; drop this levels
orf1.inc.stack<-stack(orf1.inc.split,drop=TRUE)
sum(duplicated(orf1.inc.stack$ind))
orf1.inc.stack$ID<-paste0(orf1.inc.stack$ind,"_ORF1.",unlist(sapply(table(orf1.inc.stack$ind),function(x){1:x})))
sum(duplicated(orf1.inc.stack$values)) # 1,027 are duplicated sequences: aggregate
orf1.inc.stack<-aggregate(ID~values,orf1.inc.stack,paste,collapse=";")

# ORF1
orf1.exc.split<-sapply(orf1.exc,function(x){
  a<-unlist(strsplit(x,"\\*"))
  nchar(a[1])>=7
  b<-a[-1]
  b<-b[grep("M",b)]
  pos<-sapply(b,function(x){unlist(gregexpr("M",b))[1]})
  c<-substring(b,pos,nchar(b))
  return(c(a[1],c)[nchar(c(a[1],c))>=7])
})
orf1.exc.stack<-stack(orf1.exc.split)
length(names(orf1.exc.split)[which(!names(orf1.exc.split) %in% orf1.exc.stack$ind)]) # 481 events do not generate sequences >=7aa; drop this levels
orf1.exc.stack<-stack(orf1.exc.split,drop=TRUE)
orf1.exc.stack$ID<-paste0(orf1.exc.stack$ind,"_ORF1.",unlist(sapply(table(orf1.exc.stack$ind),function(x){1:x})))
sum(duplicated(orf1.exc.stack$values)) # 3,682 are duplicated sequences: aggregate
orf1.exc.stack<-aggregate(ID~values,orf1.exc.stack,paste,collapse=";")

# ORF2
orf2.inc.split<-sapply(orf2.inc,function(x){
  a<-unlist(strsplit(x,"\\*"))
  nchar(a[1])>=7
  b<-a[-1]
  b<-b[grep("M",b)]
  pos<-sapply(b,function(x){unlist(gregexpr("M",b))[1]})
  c<-substring(b,pos,nchar(b))
  return(c(a[1],c)[nchar(c(a[1],c))>=7])
})
orf2.inc.stack<-stack(orf2.inc.split)
length(names(orf2.inc.split)[which(!names(orf2.inc.split) %in% orf2.inc.stack$ind)]) # 177 events do not generate sequences >=7aa; drop this levels
orf2.inc.stack<-stack(orf2.inc.split,drop=TRUE)
orf2.inc.stack$ID<-paste0(orf2.inc.stack$ind,"_ORF2.",unlist(sapply(table(orf2.inc.stack$ind),function(x){1:x})))
sum(duplicated(orf2.inc.stack$values)) # 1,064 are duplicated sequences: aggregate
orf2.inc.stack<-aggregate(ID~values,orf2.inc.stack,paste,collapse=";")


# ORF2
orf2.exc.split<-sapply(orf2.exc,function(x){
  a<-unlist(strsplit(x,"\\*"))
  nchar(a[1])>=7
  b<-a[-1]
  b<-b[grep("M",b)]
  pos<-sapply(b,function(x){unlist(gregexpr("M",b))[1]})
  c<-substring(b,pos,nchar(b))
  return(c(a[1],c)[nchar(c(a[1],c))>=7])
})
orf2.exc.stack<-stack(orf2.exc.split)
length(names(orf2.exc.split)[which(!names(orf2.exc.split) %in% orf2.exc.stack$ind)]) # 680 events do not generate sequences >=7aa; drop this levels
orf2.exc.stack<-stack(orf2.exc.split,drop=TRUE)
orf2.exc.stack$ID<-paste0(orf2.exc.stack$ind,"_ORF2.",unlist(sapply(table(orf2.exc.stack$ind),function(x){1:x})))
sum(duplicated(orf2.exc.stack$values)) # 3,567 are duplicated sequences: aggregate
orf2.exc.stack<-aggregate(ID~values,orf2.exc.stack,paste,collapse=";")

# ORF3
orf3.inc.split<-sapply(orf3.inc,function(x){
  a<-unlist(strsplit(x,"\\*"))
  nchar(a[1])>=7
  b<-a[-1]
  b<-b[grep("M",b)]
  pos<-sapply(b,function(x){unlist(gregexpr("M",b))[1]})
  c<-substring(b,pos,nchar(b))
  return(c(a[1],c)[nchar(c(a[1],c))>=7])
})
orf3.inc.stack<-stack(orf3.inc.split)
length(names(orf3.inc.split)[which(!names(orf3.inc.split) %in% orf3.inc.stack$ind)]) # 182 events do not generate sequences >=7aa; drop this levels
orf3.inc.stack<-stack(orf3.inc.split,drop=TRUE)
orf3.inc.stack$ID<-paste0(orf3.inc.stack$ind,"_ORF3.",unlist(sapply(table(orf3.inc.stack$ind),function(x){1:x})))
sum(duplicated(orf3.inc.stack$values)) # 991 are duplicated sequences: aggregate
orf3.inc.stack<-aggregate(ID~values,orf3.inc.stack,paste,collapse=";")

# ORF3
orf3.exc.split<-sapply(orf3.exc,function(x){
  a<-unlist(strsplit(x,"\\*"))
  nchar(a[1])>=7
  b<-a[-1]
  b<-b[grep("M",b)]
  pos<-sapply(b,function(x){unlist(gregexpr("M",b))[1]})
  c<-substring(b,pos,nchar(b))
  return(c(a[1],c)[nchar(c(a[1],c))>=7])
})
orf3.exc.stack<-stack(orf3.exc.split)
length(names(orf3.exc.split)[which(!names(orf3.exc.split) %in% orf3.exc.stack$ind)]) # 737 events do not generate sequences >=7aa; drop this levels
orf3.exc.stack<-stack(orf3.exc.split,drop=TRUE)
orf3.exc.stack$ID<-paste0(orf3.exc.stack$ind,"_ORF3.",unlist(sapply(table(orf3.exc.stack$ind),function(x){1:x})))
sum(duplicated(orf3.exc.stack$values)) # 3,339 are duplicated sequences: aggregate
orf3.exc.stack<-aggregate(ID~values,orf3.exc.stack,paste,collapse=";")

SE.fa<-rbind(orf1.inc.stack,orf1.exc.stack,orf2.inc.stack,orf2.exc.stack,orf3.inc.stack,orf3.exc.stack)
SE.fa<-aggregate(ID~values,SE.fa,paste,collapse=";") # 43,236
SE.fa<-data.frame(ID=SE.fa$ID,sequence=SE.fa$values)

##### Write .tsv with unique sequences #####
write.table(x=SE.fa,file="Fragpipe/inputs/BeatAPL_SE.tsv",row.names = FALSE,quote=FALSE,col.names=TRUE,sep="\t")

##### A5SS #####
A5SS<-read.delim("A5SS_common.txt",sep="\t")
# Select events
i1<-which(A5SS$FDR_APLvsAML<0.05 & A5SS$FDR_APLvsHealthy<0.05 &
               A5SS$dPSI_APLvsAML>0 & A5SS$dPSI_APLvsHealthy>0) # inc
i2<-which(A5SS$FDR_APLvsAML<0.05 & A5SS$FDR_APLvsHealthy<0.05 &
               A5SS$dPSI_APLvsAML<0 & A5SS$dPSI_APLvsHealthy<0) # exc

i3<-which(A5SS$all.0.healthy==19 & A5SS$all.0.aml==122 & A5SS$gt.0.apl>=1) # inc
i4<-which(A5SS$all.1.healthy==19 & A5SS$all.1.aml==122 & A5SS$lt.1.apl>=1) # exc

A5SS<-A5SS[unique(c(i1,i2,i3,i4)),] # rownames are row original positions
sum(unique(c(i1,i3)) %in% unique(c(i2,i4))) # 0 that could be inclusion or exclusion depending on rMATS PSI or my PSI


long.gr<-GRanges(seqnames=A5SS$chr,
                 ranges=IRanges(start=A5SS$longExonStart_0base+1,end=A5SS$longExonEnd),
                 strand=A5SS$strand)
short.gr<-GRanges(seqnames=A5SS$chr,
                  ranges=IRanges(start=A5SS$shortES+1,end=A5SS$shortEE),
                  strand=A5SS$strand)
flanking.gr<-GRanges(seqnames=A5SS$chr,
                     ranges=IRanges(start=A5SS$flankingES+1,end=A5SS$flankingEE),
                     strand=A5SS$strand)

long.rna<-getSeq(BSgenome.Hsapiens.UCSC.hg38,long.gr)
short.rna<-getSeq(BSgenome.Hsapiens.UCSC.hg38,short.gr)
flanking.rna<-getSeq(BSgenome.Hsapiens.UCSC.hg38,flanking.gr)

rna.inc<-paste0(long.rna,flanking.rna)[rownames(A5SS) %in% unique(c(i1,i3))]
rna.exc<-paste0(short.rna,flanking.rna)[rownames(A5SS) %in% unique(c(i2,i4))]

names(rna.inc)<-paste0("APL_A5SSinc_ID",A5SS$ID[rownames(A5SS) %in% unique(c(i1,i3))])
names(rna.exc)<-paste0("APL_A5SSexc_ID",A5SS$ID[rownames(A5SS) %in% unique(c(i2,i4))])


length(rna.inc)+length(rna.exc) # 1,618 OK

inc.seq<-lapply(1:3,function(i){substring(rna.inc,i,nchar(rna.inc))})
exc.seq<-lapply(1:3,function(i){substring(rna.exc,i,nchar(rna.exc))})

orf1.inc<-as.character(Biostrings::translate(DNAStringSet(inc.seq[[1]]),
                                             if.fuzzy.codon = "solve",
                                             no.init.codon = TRUE))
orf2.inc<-as.character(Biostrings::translate(DNAStringSet(inc.seq[[2]]),
                                             if.fuzzy.codon = "solve",
                                             no.init.codon = TRUE))
orf3.inc<-as.character(Biostrings::translate(DNAStringSet(inc.seq[[3]]),
                                             if.fuzzy.codon = "solve",
                                             no.init.codon = TRUE))

orf1.exc<-as.character(Biostrings::translate(DNAStringSet(exc.seq[[1]]),
                                             if.fuzzy.codon = "solve",
                                             no.init.codon = TRUE))
orf2.exc<-as.character(Biostrings::translate(DNAStringSet(exc.seq[[2]]),
                                             if.fuzzy.codon = "solve",
                                             no.init.codon = TRUE))
orf3.exc<-as.character(Biostrings::translate(DNAStringSet(exc.seq[[3]]),
                                             if.fuzzy.codon = "solve",
                                             no.init.codon = TRUE))

##### Keep ORFs that are >=7aa #####
# ORF1
orf1.inc.split<-sapply(orf1.inc,function(x){
  a<-unlist(strsplit(x,"\\*"))
  nchar(a[1])>=7
  b<-a[-1]
  b<-b[grep("M",b)]
  pos<-sapply(b,function(x){unlist(gregexpr("M",b))[1]})
  c<-substring(b,pos,nchar(b))
  return(c(a[1],c)[nchar(c(a[1],c))>=7])
})
orf1.inc.stack<-stack(orf1.inc.split)
length(names(orf1.inc.split)[which(!names(orf1.inc.split) %in% orf1.inc.stack$ind)]) # 13 event do not generate sequences >=7aa; drop this levels
orf1.inc.stack<-stack(orf1.inc.split,drop=TRUE)
orf1.inc.stack$ID<-paste0(orf1.inc.stack$ind,"_ORF1.",unlist(sapply(table(orf1.inc.stack$ind),function(x){1:x})))
sum(duplicated(orf1.inc.stack$values)) # 13 are duplicated sequences: aggregate
orf1.inc.stack<-aggregate(ID~values,orf1.inc.stack,paste,collapse=";")

# ORF1
orf1.exc.split<-sapply(orf1.exc,function(x){
  a<-unlist(strsplit(x,"\\*"))
  nchar(a[1])>=7
  b<-a[-1]
  b<-b[grep("M",b)]
  pos<-sapply(b,function(x){unlist(gregexpr("M",b))[1]})
  c<-substring(b,pos,nchar(b))
  return(c(a[1],c)[nchar(c(a[1],c))>=7])
})
orf1.exc.stack<-stack(orf1.exc.split)
length(names(orf1.exc.split)[which(!names(orf1.exc.split) %in% orf1.exc.stack$ind)]) # 104 events do not generate sequences >=7aa; drop this levels
orf1.exc.stack<-stack(orf1.exc.split,drop=TRUE)
orf1.exc.stack$ID<-paste0(orf1.exc.stack$ind,"_ORF1.",unlist(sapply(table(orf1.exc.stack$ind),function(x){1:x})))
sum(duplicated(orf1.exc.stack$values)) # 362 are duplicated sequences: aggregate
orf1.exc.stack<-aggregate(ID~values,orf1.exc.stack,paste,collapse=";")

# ORF2
orf2.inc.split<-sapply(orf2.inc,function(x){
  a<-unlist(strsplit(x,"\\*"))
  nchar(a[1])>=7
  b<-a[-1]
  b<-b[grep("M",b)]
  pos<-sapply(b,function(x){unlist(gregexpr("M",b))[1]})
  c<-substring(b,pos,nchar(b))
  return(c(a[1],c)[nchar(c(a[1],c))>=7])
})
orf2.inc.stack<-stack(orf2.inc.split)
length(names(orf2.inc.split)[which(!names(orf2.inc.split) %in% orf2.inc.stack$ind)]) # 26 events do not generate sequences >=7aa; drop this levels
orf2.inc.stack<-stack(orf2.inc.split,drop=TRUE)
orf2.inc.stack$ID<-paste0(orf2.inc.stack$ind,"_ORF2.",unlist(sapply(table(orf2.inc.stack$ind),function(x){1:x})))
sum(duplicated(orf2.inc.stack$values)) # 15 are duplicated sequences: aggregate
orf2.inc.stack<-aggregate(ID~values,orf2.inc.stack,paste,collapse=";")

# ORF2
orf2.exc.split<-sapply(orf2.exc,function(x){
  a<-unlist(strsplit(x,"\\*"))
  nchar(a[1])>=7
  b<-a[-1]
  b<-b[grep("M",b)]
  pos<-sapply(b,function(x){unlist(gregexpr("M",b))[1]})
  c<-substring(b,pos,nchar(b))
  return(c(a[1],c)[nchar(c(a[1],c))>=7])
})
orf2.exc.stack<-stack(orf2.exc.split)
length(names(orf2.exc.split)[which(!names(orf2.exc.split) %in% orf2.exc.stack$ind)]) # 146 events do not generate sequences >=7aa; drop this levels
orf2.exc.stack<-stack(orf2.exc.split,drop=TRUE)
orf2.exc.stack$ID<-paste0(orf2.exc.stack$ind,"_ORF2.",unlist(sapply(table(orf2.exc.stack$ind),function(x){1:x})))
sum(duplicated(orf2.exc.stack$values)) # 310 are duplicated sequences: aggregate
orf2.exc.stack<-aggregate(ID~values,orf2.exc.stack,paste,collapse=";")

# ORF3
orf3.inc.split<-sapply(orf3.inc,function(x){
  a<-unlist(strsplit(x,"\\*"))
  nchar(a[1])>=7
  b<-a[-1]
  b<-b[grep("M",b)]
  pos<-sapply(b,function(x){unlist(gregexpr("M",b))[1]})
  c<-substring(b,pos,nchar(b))
  return(c(a[1],c)[nchar(c(a[1],c))>=7])
})
orf3.inc.stack<-stack(orf3.inc.split)
length(names(orf3.inc.split)[which(!names(orf3.inc.split) %in% orf3.inc.stack$ind)]) # 15 events do not generate sequences >=7aa; drop this levels
orf3.inc.stack<-stack(orf3.inc.split,drop=TRUE)
orf3.inc.stack$ID<-paste0(orf3.inc.stack$ind,"_ORF3.",unlist(sapply(table(orf3.inc.stack$ind),function(x){1:x})))
sum(duplicated(orf3.inc.stack$values)) # 16 are duplicated sequences: aggregate
orf3.inc.stack<-aggregate(ID~values,orf3.inc.stack,paste,collapse=";")

# ORF3
orf3.exc.split<-sapply(orf3.exc,function(x){
  a<-unlist(strsplit(x,"\\*"))
  nchar(a[1])>=7
  b<-a[-1]
  b<-b[grep("M",b)]
  pos<-sapply(b,function(x){unlist(gregexpr("M",b))[1]})
  c<-substring(b,pos,nchar(b))
  return(c(a[1],c)[nchar(c(a[1],c))>=7])
})
orf3.exc.stack<-stack(orf3.exc.split)
length(names(orf3.exc.split)[which(!names(orf3.exc.split) %in% orf3.exc.stack$ind)]) # 130 events do not generate sequences >=7aa; drop this levels
orf3.exc.stack<-stack(orf3.exc.split,drop=TRUE)
orf3.exc.stack$ID<-paste0(orf3.exc.stack$ind,"_ORF3.",unlist(sapply(table(orf3.exc.stack$ind),function(x){1:x})))
sum(duplicated(orf3.exc.stack$values)) # 311 are duplicated sequences: aggregate
orf3.exc.stack<-aggregate(ID~values,orf3.exc.stack,paste,collapse=";")

# Rbind and aggregate duplicated sequences
A5SS.fa<-rbind(orf1.inc.stack,orf1.exc.stack,orf2.inc.stack,orf2.exc.stack,orf3.inc.stack,orf3.exc.stack)
A5SS.fa<-aggregate(ID~values,A5SS.fa,paste,collapse=";") # 7,013
A5SS.fa<-data.frame(ID=A5SS.fa$ID,sequence=A5SS.fa$values)

##### Write .tsv with unique sequences #####
write.table(x=A5SS.fa,file="Fragpipe/inputs/BeatAPL_A5SS.tsv",row.names = FALSE,quote=FALSE,col.names=TRUE,sep="\t")

##### A3SS #####
A3SS<-read.delim("A3SS_common.txt",sep="\t")
# Select events
i1<-which(A3SS$FDR_APLvsAML<0.05 & A3SS$FDR_APLvsHealthy<0.05 &
               A3SS$dPSI_APLvsAML>0 & A3SS$dPSI_APLvsHealthy>0) # inc
i2<-which(A3SS$FDR_APLvsAML<0.05 & A3SS$FDR_APLvsHealthy<0.05 &
               A3SS$dPSI_APLvsAML<0 & A3SS$dPSI_APLvsHealthy<0) # exc

i3<-which(A3SS$all.0.healthy==19 & A3SS$all.0.aml==122 & A3SS$gt.0.apl>=1) # inc
i4<-which(A3SS$all.1.healthy==19 & A3SS$all.1.aml==122 & A3SS$lt.1.apl>=1) # exc
A3SS<-A3SS[unique(c(i1,i2,i3,i4)),] # rownames are row original positions
sum(unique(c(i1,i3)) %in% unique(c(i2,i4))) # 0 that could be inclusion or exclusion depending on rMATS PSI or my PSI


long.gr<-GRanges(seqnames=A3SS$chr,
                 ranges=IRanges(start=A3SS$longExonStart_0base+1,end=A3SS$longExonEnd),
                 strand=A3SS$strand)
short.gr<-GRanges(seqnames=A3SS$chr,
                  ranges=IRanges(start=A3SS$shortES+1,end=A3SS$shortEE),
                  strand=A3SS$strand)
flanking.gr<-GRanges(seqnames=A3SS$chr,
                     ranges=IRanges(start=A3SS$flankingES+1,end=A3SS$flankingEE),
                     strand=A3SS$strand)

long.rna<-getSeq(BSgenome.Hsapiens.UCSC.hg38,long.gr)
short.rna<-getSeq(BSgenome.Hsapiens.UCSC.hg38,short.gr)
flanking.rna<-getSeq(BSgenome.Hsapiens.UCSC.hg38,flanking.gr)


rna.inc<-paste0(flanking.rna,long.rna)[rownames(A3SS) %in% unique(c(i1,i3))]
rna.exc<-paste0(flanking.rna,short.rna)[rownames(A3SS) %in% unique(c(i2,i4))]

names(rna.inc)<-paste0("APL_A3SSinc_ID",A3SS$ID[rownames(A3SS) %in% unique(c(i1,i3))])
names(rna.exc)<-paste0("APL_A3SSexc_ID",A3SS$ID[rownames(A3SS) %in% unique(c(i2,i4))])


length(rna.inc)+length(rna.exc) # 4,083 OK

inc.seq<-lapply(1:3,function(i){substring(rna.inc,i,nchar(rna.inc))})
exc.seq<-lapply(1:3,function(i){substring(rna.exc,i,nchar(rna.exc))})

orf1.inc<-as.character(Biostrings::translate(DNAStringSet(inc.seq[[1]]),
                                             if.fuzzy.codon = "solve",
                                             no.init.codon = TRUE))
orf2.inc<-as.character(Biostrings::translate(DNAStringSet(inc.seq[[2]]),
                                             if.fuzzy.codon = "solve",
                                             no.init.codon = TRUE))
orf3.inc<-as.character(Biostrings::translate(DNAStringSet(inc.seq[[3]]),
                                             if.fuzzy.codon = "solve",
                                             no.init.codon = TRUE))

orf1.exc<-as.character(Biostrings::translate(DNAStringSet(exc.seq[[1]]),
                                             if.fuzzy.codon = "solve",
                                             no.init.codon = TRUE))
orf2.exc<-as.character(Biostrings::translate(DNAStringSet(exc.seq[[2]]),
                                             if.fuzzy.codon = "solve",
                                             no.init.codon = TRUE))
orf3.exc<-as.character(Biostrings::translate(DNAStringSet(exc.seq[[3]]),
                                             if.fuzzy.codon = "solve",
                                             no.init.codon = TRUE))

##### Keep ORFs that are >=7aa #####
# ORF1
orf1.inc.split<-sapply(orf1.inc,function(x){
  a<-unlist(strsplit(x,"\\*"))
  nchar(a[1])>=7
  b<-a[-1]
  b<-b[grep("M",b)]
  pos<-sapply(b,function(x){unlist(gregexpr("M",b))[1]})
  c<-substring(b,pos,nchar(b))
  return(c(a[1],c)[nchar(c(a[1],c))>=7])
})
orf1.inc.stack<-stack(orf1.inc.split)
length(names(orf1.inc.split)[which(!names(orf1.inc.split) %in% orf1.inc.stack$ind)]) # 13 events do not generate sequences >=7aa; drop this levels
orf1.inc.stack<-stack(orf1.inc.split,drop=TRUE)
orf1.inc.stack$ID<-paste0(orf1.inc.stack$ind,"_ORF1.",unlist(sapply(table(orf1.inc.stack$ind),function(x){1:x})))
sum(duplicated(orf1.inc.stack$values)) # 3 are duplicated sequences: aggregate
orf1.inc.stack<-aggregate(ID~values,orf1.inc.stack,paste,collapse=";")

# ORF1
orf1.exc.split<-sapply(orf1.exc,function(x){
  a<-unlist(strsplit(x,"\\*"))
  nchar(a[1])>=7
  b<-a[-1]
  b<-b[grep("M",b)]
  pos<-sapply(b,function(x){unlist(gregexpr("M",b))[1]})
  c<-substring(b,pos,nchar(b))
  return(c(a[1],c)[nchar(c(a[1],c))>=7])
})
orf1.exc.stack<-stack(orf1.exc.split)
length(names(orf1.exc.split)[which(!names(orf1.exc.split) %in% orf1.exc.stack$ind)]) # 282 events do not generate sequences >=7aa; drop this levels
orf1.exc.stack<-stack(orf1.exc.split,drop=TRUE)
orf1.exc.stack$ID<-paste0(orf1.exc.stack$ind,"_ORF1.",unlist(sapply(table(orf1.exc.stack$ind),function(x){1:x})))
sum(duplicated(orf1.exc.stack$values)) # 1,079 are duplicated sequences: aggregate
orf1.exc.stack<-aggregate(ID~values,orf1.exc.stack,paste,collapse=";")



# ORF2
orf2.inc.split<-sapply(orf2.inc,function(x){
  a<-unlist(strsplit(x,"\\*"))
  nchar(a[1])>=7
  b<-a[-1]
  b<-b[grep("M",b)]
  pos<-sapply(b,function(x){unlist(gregexpr("M",b))[1]})
  c<-substring(b,pos,nchar(b))
  return(c(a[1],c)[nchar(c(a[1],c))>=7])
})
orf2.inc.stack<-stack(orf2.inc.split)
length(names(orf2.inc.split)[which(!names(orf2.inc.split) %in% orf2.inc.stack$ind)]) # 31 events do not generate sequences >=7aa; drop this levels
orf2.inc.stack<-stack(orf2.inc.split,drop=TRUE)
orf2.inc.stack$ID<-paste0(orf2.inc.stack$ind,"_ORF2.",unlist(sapply(table(orf2.inc.stack$ind),function(x){1:x})))
sum(duplicated(orf2.inc.stack$values)) # 12 are duplicated sequences: aggregate
orf2.inc.stack<-aggregate(ID~values,orf2.inc.stack,paste,collapse=";")


# ORF2
orf2.exc.split<-sapply(orf2.exc,function(x){
  a<-unlist(strsplit(x,"\\*"))
  nchar(a[1])>=7
  b<-a[-1]
  b<-b[grep("M",b)]
  pos<-sapply(b,function(x){unlist(gregexpr("M",b))[1]})
  c<-substring(b,pos,nchar(b))
  return(c(a[1],c)[nchar(c(a[1],c))>=7])
})
orf2.exc.stack<-stack(orf2.exc.split)
length(names(orf2.exc.split)[which(!names(orf2.exc.split) %in% orf2.exc.stack$ind)]) # 369 events do not generate sequences >=7aa; drop this levels
orf2.exc.stack<-stack(orf2.exc.split,drop=TRUE)
orf2.exc.stack$ID<-paste0(orf2.exc.stack$ind,"_ORF2.",unlist(sapply(table(orf2.exc.stack$ind),function(x){1:x})))
sum(duplicated(orf2.exc.stack$values)) # 1,111 are duplicated sequences: aggregate
orf2.exc.stack<-aggregate(ID~values,orf2.exc.stack,paste,collapse=";")

# ORF3
orf3.inc.split<-sapply(orf3.inc,function(x){
  a<-unlist(strsplit(x,"\\*"))
  nchar(a[1])>=7
  b<-a[-1]
  b<-b[grep("M",b)]
  pos<-sapply(b,function(x){unlist(gregexpr("M",b))[1]})
  c<-substring(b,pos,nchar(b))
  return(c(a[1],c)[nchar(c(a[1],c))>=7])
})
orf3.inc.stack<-stack(orf3.inc.split)
length(names(orf3.inc.split)[which(!names(orf3.inc.split) %in% orf3.inc.stack$ind)]) # 21 events do not generate sequences >=7aa; drop this levels
orf3.inc.stack<-stack(orf3.inc.split,drop=TRUE)
orf3.inc.stack$ID<-paste0(orf3.inc.stack$ind,"_ORF3.",unlist(sapply(table(orf3.inc.stack$ind),function(x){1:x})))
sum(duplicated(orf3.inc.stack$values)) # 9 are duplicated sequences: aggregate
orf3.inc.stack<-aggregate(ID~values,orf3.inc.stack,paste,collapse=";")


# ORF3
orf3.exc.split<-sapply(orf3.exc,function(x){
  a<-unlist(strsplit(x,"\\*"))
  nchar(a[1])>=7
  b<-a[-1]
  b<-b[grep("M",b)]
  pos<-sapply(b,function(x){unlist(gregexpr("M",b))[1]})
  c<-substring(b,pos,nchar(b))
  return(c(a[1],c)[nchar(c(a[1],c))>=7])
})
orf3.exc.stack<-stack(orf3.exc.split)
length(names(orf3.exc.split)[which(!names(orf3.exc.split) %in% orf3.exc.stack$ind)]) # 428 events do not generate sequences >=7aa; drop this levels
orf3.exc.stack<-stack(orf3.exc.split,drop=TRUE)
orf3.exc.stack$ID<-paste0(orf3.exc.stack$ind,"_ORF3.",unlist(sapply(table(orf3.exc.stack$ind),function(x){1:x})))
sum(duplicated(orf3.exc.stack$values)) # 1,062 are duplicated sequences: aggregate
orf3.exc.stack<-aggregate(ID~values,orf3.exc.stack,paste,collapse=";")

A3SS.fa<-rbind(orf1.inc.stack,orf1.exc.stack,orf2.inc.stack,orf2.exc.stack,orf3.inc.stack,orf3.exc.stack)
A3SS.fa<-aggregate(ID~values,A3SS.fa,paste,collapse=";") # 13,490
A3SS.fa<-data.frame(ID=A3SS.fa$ID,sequence=A3SS.fa$values)


##### Write .tsv with unique sequences #####
write.table(x=A3SS.fa,file="Fragpipe/inputs/BeatAPL_A3SS.tsv",row.names = FALSE,quote=FALSE,col.names=TRUE,sep="\t")


##### RI #####
RI<-read.delim("RI_common.txt",sep="\t")
# Select events
i1<-which(RI$FDR_APLvsAML<0.05 & RI$FDR_APLvsHealthy<0.05 &
             RI$dPSI_APLvsAML>0 & RI$dPSI_APLvsHealthy>0) # inc
i2<-which(RI$FDR_APLvsAML<0.05 & RI$FDR_APLvsHealthy<0.05 &
             RI$dPSI_APLvsAML<0 & RI$dPSI_APLvsHealthy<0) # exc

i3<-which(RI$all.0.healthy==19 & RI$all.0.aml==122 & RI$gt.0.apl>=1) # inc
i4<-which(RI$all.1.healthy==19 & RI$all.1.aml==122 & RI$lt.1.apl>=1) # exc


RI<-RI[unique(c(i1,i2,i3,i4)),] # rownames are row original positions
sum(unique(c(i1,i3)) %in% unique(c(i2,i4))) # 0 that could be inclusion or exclusion depending on rMATS PSI or my PSI


UE.gr<-GRanges(seqnames=RI$chr,
               ranges=IRanges(start=RI$upstreamES+1,end=RI$upstreamEE),
               strand=RI$strand)
DE.gr<-GRanges(seqnames=RI$chr,
               ranges=IRanges(start=RI$downstreamES+1,end=RI$downstreamEE),
               strand=RI$strand)
RI.gr<-GRanges(seqnames=RI$chr,
               ranges=IRanges(start=RI$riExonStart_0base+1,end=RI$riExonEnd),
               strand=RI$strand)

UE.rna<-getSeq(BSgenome.Hsapiens.UCSC.hg38,UE.gr)
DE.rna<-getSeq(BSgenome.Hsapiens.UCSC.hg38,DE.gr)
RI.rna<-getSeq(BSgenome.Hsapiens.UCSC.hg38,RI.gr)

# rna.inc # no event for inclusion
rna.exc<-ifelse(RI$strand=="+",
                paste0(UE.rna,DE.rna),
                paste0(DE.rna,UE.rna))

# names(rna.inc)<-paste0("APL_RIinc_ID",RI$ID[rownames(RI) %in% unique(c(i1,i3))])
names(rna.exc)<-paste0("APL_RIexc_ID",RI$ID)

length(rna.exc) # 448, OK
# length(rna.inc)+length(rna.exc)

# inc.seq<-lapply(1:3,function(i){substring(rna.inc,i,nchar(rna.inc))})
exc.seq<-lapply(1:3,function(i){substring(rna.exc,i,nchar(rna.exc))})

orf1.exc<-as.character(Biostrings::translate(DNAStringSet(exc.seq[[1]]),
                                             if.fuzzy.codon = "solve",
                                             no.init.codon = TRUE))
orf2.exc<-as.character(Biostrings::translate(DNAStringSet(exc.seq[[2]]),
                                             if.fuzzy.codon = "solve",
                                             no.init.codon = TRUE))
orf3.exc<-as.character(Biostrings::translate(DNAStringSet(exc.seq[[3]]),
                                             if.fuzzy.codon = "solve",
                                             no.init.codon = TRUE))

##### Keep ORFs that are >=7aa #####
# ORF1
orf1.exc.split<-sapply(orf1.exc,function(x){
  a<-unlist(strsplit(x,"\\*"))
  nchar(a[1])>=7
  b<-a[-1]
  b<-b[grep("M",b)]
  pos<-sapply(b,function(x){unlist(gregexpr("M",b))[1]})
  c<-substring(b,pos,nchar(b))
  return(c(a[1],c)[nchar(c(a[1],c))>=7])
})
orf1.exc.stack<-stack(orf1.exc.split)
length(names(orf1.exc.split)[which(!names(orf1.exc.split) %in% orf1.exc.stack$ind)]) # 37 events do not generate sequences >=7aa; drop this levels
orf1.exc.stack<-stack(orf1.exc.split,drop=TRUE)
orf1.exc.stack$ID<-paste0(orf1.exc.stack$ind,"_ORF1.",unlist(sapply(table(orf1.exc.stack$ind),function(x){1:x})))
sum(duplicated(orf1.exc.stack$values)) # 12 are duplicated sequences: aggregate
orf1.exc.stack<-aggregate(ID~values,orf1.exc.stack,paste,collapse=";")

# ORF2
orf2.exc.split<-sapply(orf2.exc,function(x){
  a<-unlist(strsplit(x,"\\*"))
  nchar(a[1])>=7
  b<-a[-1]
  b<-b[grep("M",b)]
  pos<-sapply(b,function(x){unlist(gregexpr("M",b))[1]})
  c<-substring(b,pos,nchar(b))
  return(c(a[1],c)[nchar(c(a[1],c))>=7])
})
orf2.exc.stack<-stack(orf2.exc.split)
length(names(orf2.exc.split)[which(!names(orf2.exc.split) %in% orf2.exc.stack$ind)]) # 43 event do not generate sequences >=7aa; drop this levels
orf2.exc.stack<-stack(orf2.exc.split,drop=TRUE)
orf2.exc.stack$ID<-paste0(orf2.exc.stack$ind,"_ORF2.",unlist(sapply(table(orf2.exc.stack$ind),function(x){1:x})))
sum(duplicated(orf2.exc.stack$values)) # 12 are duplicated sequences: aggregate
orf2.exc.stack<-aggregate(ID~values,orf2.exc.stack,paste,collapse=";")

# ORF3
orf3.exc.split<-sapply(orf3.exc,function(x){
  a<-unlist(strsplit(x,"\\*"))
  nchar(a[1])>=7
  b<-a[-1]
  b<-b[grep("M",b)]
  pos<-sapply(b,function(x){unlist(gregexpr("M",b))[1]})
  c<-substring(b,pos,nchar(b))
  return(c(a[1],c)[nchar(c(a[1],c))>=7])
})
orf3.exc.stack<-stack(orf3.exc.split)
length(names(orf3.exc.split)[which(!names(orf3.exc.split) %in% orf3.exc.stack$ind)]) # 46 event do not generate sequences >=7aa; drop this levels
orf3.exc.stack<-stack(orf3.exc.split,drop=TRUE)
orf3.exc.stack$ID<-paste0(orf3.exc.stack$ind,"_ORF3.",unlist(sapply(table(orf3.exc.stack$ind),function(x){1:x})))
sum(duplicated(orf3.exc.stack$values)) # 14 are duplicated sequences: aggregate
orf3.exc.stack<-aggregate(ID~values,orf3.exc.stack,paste,collapse=";")

RI.fa<-rbind(orf1.exc.stack,orf2.exc.stack,orf3.exc.stack)
RI.fa<-aggregate(ID~values,RI.fa,paste,collapse=";") # 1,762
RI.fa<-data.frame(ID=RI.fa$ID,sequence=RI.fa$values)


##### Write .tsv with unique sequences #####
write.table(x=RI.fa,file="Fragpipe/inputs/BeatAPL_RI.tsv",row.names = FALSE,quote=FALSE,col.names=TRUE,sep="\t")


############################### Beat AML GTF ###################################
# First we need to subset the GTF
beat.df<-as.data.frame(beat)
beat.df<-beat.df[grepl("MSTRG",beat.df$transcript_id) & # keep only novel transcripts
                   !beat.df$strand=="*" &  # we do not want transcripts without strand assigned
                   beat.df$type=="exon" & # if we do not subset, transcripts are translated from start to end
                   beat.df$seqnames %in% c(paste0("chr",1:22),"chrX","chrY","chrM"),] # Keep transcripts in these chromosomes
beat.df$strand<-droplevels(beat.df$strand) # drop * level because it raises an error in getTxsSeq
beat.df$transcript_id<-paste0("BEAT_",beat.df$transcript_id) # Add this prefix for further identification
length(unique(beat.df$transcript_id)) # 136,499 transcripts
table(beat.df$seqnames[!duplicated(beat.df$transcript_id)]) # Only 4 in the mitocondrial genome

# Get transcript sequence
beat.gr<-makeGRangesListFromDataFrame(beat.df,
                                      ignore.strand = FALSE,
                                      split.field="transcript_id",
                                      names.field="transcript_id") # gtf must be parsed as GRangesList to getTxsSeq

beat.tx.seq<-getTxsSeq(beat.gr,BSgenome.Hsapiens.UCSC.hg38) 


# Get 3 reading frames from each sequence
# Differentiate between std and mt genomes for next step
beat.seq.std<-lapply(1:3,function(i){substring(beat.tx.seq[names(beat.tx.seq) %in% unique(beat.df$transcript_id[beat.df$seqnames!="chrM"])],
                                               i,
                                               width(beat.tx.seq[names(beat.tx.seq) %in% unique(beat.df$transcript_id[beat.df$seqnames!="chrM"])]))})
beat.seq.mt<-lapply(1:3,function(i){substring(beat.tx.seq[names(beat.tx.seq) %in% unique(beat.df$transcript_id[beat.df$seqnames=="chrM"])],
                                              i,
                                              width(beat.tx.seq[names(beat.tx.seq) %in% unique(beat.df$transcript_id[beat.df$seqnames=="chrM"])]))})



# Check available genetic codes in Biostrings
SGC0<-getGeneticCode("SGC0")
SGC1<-getGeneticCode("SGC1")

# Translate (use no.init.codon=TRUE to avoid the use of alternative initiation codons)
# std: from ATG to TAA/TAG/TGA
orf1.std<-as.character(Biostrings::translate(DNAStringSet(beat.seq.std[[1]]),
                                             genetic.code = SGC0,
                                             if.fuzzy.codon = "solve",
                                             no.init.codon = TRUE))
orf2.std<-as.character(Biostrings::translate(DNAStringSet(beat.seq.std[[2]]),
                                             genetic.code = SGC0,
                                             if.fuzzy.codon = "solve",
                                             no.init.codon = TRUE))
orf3.std<-as.character(Biostrings::translate(DNAStringSet(beat.seq.std[[3]]),
                                             genetic.code = SGC0,
                                             if.fuzzy.codon = "solve",
                                             no.init.codon = TRUE))
# mt: from ATA/ATG to TAA/TAG/AGA/AGG)
orf1.mt<-as.character(Biostrings::translate(DNAStringSet(beat.seq.mt[[1]]),
                                            genetic.code = SGC1,
                                            if.fuzzy.codon = "solve",
                                            no.init.codon = TRUE))
orf2.mt<-as.character(Biostrings::translate(DNAStringSet(beat.seq.mt[[2]]),
                                            genetic.code = SGC1,
                                            if.fuzzy.codon = "solve",
                                            no.init.codon = TRUE))
orf3.mt<-as.character(Biostrings::translate(DNAStringSet(beat.seq.mt[[3]]),
                                            genetic.code = SGC1,
                                            if.fuzzy.codon = "solve",
                                            no.init.codon = TRUE))

##### Keep ORFs that are >=7aa #####
# 1. Split each transcript sequence by stop codon (*)
# 2. Substring each subsequence from first M (i.e. 1st ATG)
# 3. Report only sequences >=7aa
# ORF1 std
orf1.std.split<-sapply(orf1.std,function(x){
  a<-unlist(strsplit(x,"\\*"))
  a<-a[grep("M",a)]
  pos<-sapply(a,function(x){unlist(gregexpr("M",x))[1]})
  b<-substring(a,pos,nchar(a))
  return(b[nchar(b)>=7])
})
orf1.std.stack<-stack(orf1.std.split)
length(names(orf1.std.split)[which(!names(orf1.std.split) %in% orf1.std.stack$ind)]) # 3976 transcripts do not generate sequences >=7aa; drop this levels
orf1.std.stack<-stack(orf1.std.split,drop=TRUE)
orf1.std.stack$ID<-paste0(orf1.std.stack$ind,"_ORF1.",unlist(sapply(table(orf1.std.stack$ind),function(x){1:x})))
sum(duplicated(orf1.std.stack$values)) # 1,162,605 are duplicated sequences: aggregate
orf1.std.stack<-aggregate(ID~values,orf1.std.stack,paste,collapse=";") # 646,947

# ORF2 std
orf2.std.split<-sapply(orf2.std,function(x){
  a<-unlist(strsplit(x,"\\*"))
  a<-a[grep("M",a)]
  pos<-sapply(a,function(x){unlist(gregexpr("M",x))[1]})
  b<-substring(a,pos,nchar(a))
  return(b[nchar(b)>=7])
})
orf2.std.stack<-stack(orf2.std.split)
length(names(orf2.std.split)[which(!names(orf2.std.split) %in% orf2.std.stack$ind)]) # 4173 transcripts do not generate sequences >=7aa; drop this levels
orf2.std.stack<-stack(orf2.std.split,drop=TRUE)
orf2.std.stack$ID<-paste0(orf2.std.stack$ind,"_ORF2.",unlist(sapply(table(orf2.std.stack$ind),function(x){1:x})))
sum(duplicated(orf2.std.stack$values)) # 1,129,090 are duplicated sequences: aggregate
orf2.std.stack<-aggregate(ID~values,orf2.std.stack,paste,collapse=";") # 647,159

# ORF3 std
orf3.std.split<-sapply(orf3.std,function(x){
  a<-unlist(strsplit(x,"\\*"))
  a<-a[grep("M",a)]
  pos<-sapply(a,function(x){unlist(gregexpr("M",x))[1]})
  b<-substring(a,pos,nchar(a))
  return(b[nchar(b)>=7])
})
orf3.std.stack<-stack(orf3.std.split)
length(names(orf3.std.split)[which(!names(orf3.std.split) %in% orf3.std.stack$ind)]) # 4005 transcripts do not generate sequences >=7aa; drop this levels
orf3.std.stack<-stack(orf3.std.split,drop=TRUE)
orf3.std.stack$ID<-paste0(orf3.std.stack$ind,"_ORF3.",unlist(sapply(table(orf3.std.stack$ind),function(x){1:x})))
sum(duplicated(orf3.std.stack$values)) # 1,122,784 are duplicated sequences: aggregate
orf3.std.stack<-aggregate(ID~values,orf3.std.stack,paste,collapse=";") # 646,177



# ORF1 mt
orf1.mt.split<-sapply(orf1.mt,function(x){
  a<-unlist(strsplit(x,"\\*"))
  a<-a[grep("M",a)]
  pos<-sapply(a,function(x){unlist(gregexpr("M",x))[1]})
  b<-substring(a,pos,nchar(a))
  return(b[nchar(b)>=7])
})
orf1.mt.stack<-stack(orf1.mt.split)
length(names(orf1.mt.split)[which(!names(orf1.mt.split) %in% orf1.mt.stack$ind)]) # 0 transcripts do not generate sequences >=7aa; drop this levels
orf1.mt.stack<-stack(orf1.mt.split,drop=TRUE)
orf1.mt.stack$ID<-paste0(orf1.mt.stack$ind,"_ORF1.",unlist(sapply(table(orf1.mt.stack$ind),function(x){1:x})))
sum(duplicated(orf1.mt.stack$values)) # 0 are duplicated sequences: aggregate
orf1.mt.stack<-aggregate(ID~values,orf1.mt.stack,paste,collapse=";") # 25

# ORF2 mt
orf2.mt.split<-sapply(orf2.mt,function(x){
  a<-unlist(strsplit(x,"\\*"))
  a<-a[grep("M",a)]
  pos<-sapply(a,function(x){unlist(gregexpr("M",x))[1]})
  b<-substring(a,pos,nchar(a))
  return(b[nchar(b)>=7])
})
orf2.mt.stack<-stack(orf2.mt.split)
length(names(orf2.mt.split)[which(!names(orf2.mt.split) %in% orf2.mt.stack$ind)]) # 0 transcripts do not generate sequences >=7aa; drop this levels
orf2.mt.stack<-stack(orf2.mt.split,drop=TRUE)
orf2.mt.stack$ID<-paste0(orf2.mt.stack$ind,"_ORF2.",unlist(sapply(table(orf2.mt.stack$ind),function(x){1:x})))
sum(duplicated(orf2.mt.stack$values)) # 0 are duplicated sequences: aggregate
orf2.mt.stack<-aggregate(ID~values,orf2.mt.stack,paste,collapse=";") # 20

# ORF3 mt
orf3.mt.split<-sapply(orf3.mt,function(x){
  a<-unlist(strsplit(x,"\\*"))
  a<-a[grep("M",a)]
  pos<-sapply(a,function(x){unlist(gregexpr("M",x))[1]})
  b<-substring(a,pos,nchar(a))
  return(b[nchar(b)>=7])
})
orf3.mt.stack<-stack(orf3.mt.split)
length(names(orf3.mt.split)[which(!names(orf3.mt.split) %in% orf3.mt.stack$ind)]) # 0 transcripts do not generate sequences >=7aa; drop this levels
orf3.mt.stack<-stack(orf3.mt.split,drop=TRUE)
orf3.mt.stack$ID<-paste0(orf3.mt.stack$ind,"_ORF3.",unlist(sapply(table(orf3.mt.stack$ind),function(x){1:x})))
sum(duplicated(orf3.mt.stack$values)) # 0 are duplicated sequences: aggregate
orf3.mt.stack<-aggregate(ID~values,orf3.mt.stack,paste,collapse=";") # 23


beat.fa<-rbind(orf1.std.stack,orf2.std.stack,orf3.std.stack,
               orf1.mt.stack,orf2.mt.stack,orf3.mt.stack)
beat.fa<-aggregate(ID~values,beat.fa,paste,collapse=";") # 1,146,3311
beat.fa<-data.frame(ID=beat.fa$ID,sequence=beat.fa$values)

##### Write .tsv with unique sequences #####
write.table(x=beat.fa,file="Fragpipe/inputs/BeatAML_stringtie.tsv",row.names = FALSE,quote=FALSE,col.names=TRUE,sep="\t")


############################### TCGA-LAML GTF ###################################
# First we need to subset the GTF
tcga.df<-as.data.frame(tcga)
tcga.df<-tcga.df[grepl("MSTRG",tcga.df$transcript_id) & # keep only novel transcripts
                   !tcga.df$strand=="*" &  # we do not want transcripts without strand assigned
                   tcga.df$type=="exon" & # if we do not subset, transcripts are translated from start to end
                   tcga.df$seqnames %in% c(paste0("chr",1:22),"chrX","chrY","chrM"),] # Keep transcripts in these chromosomes
tcga.df$strand<-droplevels(tcga.df$strand) # drop * level because it raises an error in getTxsSeq
tcga.df$transcript_id<-paste0("TCGA_",tcga.df$transcript_id)
length(unique(tcga.df$transcript_id)) # 206,735 transcripts
table(tcga.df$seqnames[!duplicated(tcga.df$transcript_id)]) # Only 8 in the mitocondrial genome

# Get transcript sequence
tcga.gr<-makeGRangesListFromDataFrame(tcga.df,
                                      ignore.strand = FALSE,
                                      split.field="transcript_id",
                                      names.field="transcript_id") # gtf must be parsed as GRangesList to getTxsSeq

tcga.tx.seq<-getTxsSeq(tcga.gr,BSgenome.Hsapiens.UCSC.hg38) 


# Get 3 reading frames from each sequence
# Differentiate between std and mt genomes for next step
tcga.seq.std<-lapply(1:3,function(i){substring(tcga.tx.seq[names(tcga.tx.seq) %in% unique(tcga.df$transcript_id[tcga.df$seqnames!="chrM"])],
                                               i,
                                               width(tcga.tx.seq[names(tcga.tx.seq) %in% unique(tcga.df$transcript_id[tcga.df$seqnames!="chrM"])]))})
tcga.seq.mt<-lapply(1:3,function(i){substring(tcga.tx.seq[names(tcga.tx.seq) %in% unique(tcga.df$transcript_id[tcga.df$seqnames=="chrM"])],
                                              i,
                                              width(tcga.tx.seq[names(tcga.tx.seq) %in% unique(tcga.df$transcript_id[tcga.df$seqnames=="chrM"])]))})



# Check available genetic codes in Biostrings
SGC0<-getGeneticCode("SGC0")
SGC1<-getGeneticCode("SGC1")

# Translate (use no.init.codon=TRUE to avoid the use of alternative initiation codons)
# std: from ATG to TAA/TAG/TGA
orf1.std<-as.character(Biostrings::translate(DNAStringSet(tcga.seq.std[[1]]),
                                             genetic.code = SGC0,
                                             if.fuzzy.codon = "solve",
                                             no.init.codon = TRUE))
orf2.std<-as.character(Biostrings::translate(DNAStringSet(tcga.seq.std[[2]]),
                                             genetic.code = SGC0,
                                             if.fuzzy.codon = "solve",
                                             no.init.codon = TRUE))
orf3.std<-as.character(Biostrings::translate(DNAStringSet(tcga.seq.std[[3]]),
                                             genetic.code = SGC0,
                                             if.fuzzy.codon = "solve",
                                             no.init.codon = TRUE))
# mt: from ATA/ATG to TAA/TAG/AGA/AGG)
orf1.mt<-as.character(Biostrings::translate(DNAStringSet(tcga.seq.mt[[1]]),
                                            genetic.code = SGC1,
                                            if.fuzzy.codon = "solve",
                                            no.init.codon = TRUE))
orf2.mt<-as.character(Biostrings::translate(DNAStringSet(tcga.seq.mt[[2]]),
                                            genetic.code = SGC1,
                                            if.fuzzy.codon = "solve",
                                            no.init.codon = TRUE))
orf3.mt<-as.character(Biostrings::translate(DNAStringSet(tcga.seq.mt[[3]]),
                                            genetic.code = SGC1,
                                            if.fuzzy.codon = "solve",
                                            no.init.codon = TRUE))

##### Keep ORFs that are >=7aa #####
# 1. Split each transcript sequence by stop codon (*)
# 2. Substring each subsequence from first M (i.e. 1st ATG)
# 3. Report only sequences >=7aa
# ORF1 std
orf1.std.split<-sapply(orf1.std,function(x){
  a<-unlist(strsplit(x,"\\*"))
  a<-a[grep("M",a)]
  pos<-sapply(a,function(x){unlist(gregexpr("M",x))[1]})
  b<-substring(a,pos,nchar(a))
  return(b[nchar(b)>=7])
})
orf1.std.stack<-stack(orf1.std.split)
length(names(orf1.std.split)[which(!names(orf1.std.split) %in% orf1.std.stack$ind)]) # 1753 transcripts do not generate sequences >=7aa; drop this levels
orf1.std.stack<-stack(orf1.std.split,drop=TRUE)
orf1.std.stack$ID<-paste0(orf1.std.stack$ind,"_ORF1.",unlist(sapply(table(orf1.std.stack$ind),function(x){1:x})))
sum(duplicated(orf1.std.stack$values)) # 5,007,531 are duplicated sequences: aggregate
orf1.std.stack<-aggregate(ID~values,orf1.std.stack,paste,collapse=";") # 2,679,411

# ORF2 std
orf2.std.split<-sapply(orf2.std,function(x){
  a<-unlist(strsplit(x,"\\*"))
  a<-a[grep("M",a)]
  pos<-sapply(a,function(x){unlist(gregexpr("M",x))[1]})
  b<-substring(a,pos,nchar(a))
  return(b[nchar(b)>=7])
})
orf2.std.stack<-stack(orf2.std.split)
length(names(orf2.std.split)[which(!names(orf2.std.split) %in% orf2.std.stack$ind)]) # 1739 transcripts do not generate sequences >=7aa; drop this levels
orf2.std.stack<-stack(orf2.std.split,drop=TRUE)
orf2.std.stack$ID<-paste0(orf2.std.stack$ind,"_ORF2.",unlist(sapply(table(orf2.std.stack$ind),function(x){1:x})))
sum(duplicated(orf2.std.stack$values)) # 5,002,802 are duplicated sequences: aggregate
orf2.std.stack<-aggregate(ID~values,orf2.std.stack,paste,collapse=";") # 2,680,640

# ORF3 std
orf3.std.split<-sapply(orf3.std,function(x){
  a<-unlist(strsplit(x,"\\*"))
  a<-a[grep("M",a)]
  pos<-sapply(a,function(x){unlist(gregexpr("M",x))[1]})
  b<-substring(a,pos,nchar(a))
  return(b[nchar(b)>=7])
})
orf3.std.stack<-stack(orf3.std.split)
length(names(orf3.std.split)[which(!names(orf3.std.split) %in% orf3.std.stack$ind)]) # 1807 transcripts do not generate sequences >=7aa; drop this levels
orf3.std.stack<-stack(orf3.std.split,drop=TRUE)
orf3.std.stack$ID<-paste0(orf3.std.stack$ind,"_ORF3.",unlist(sapply(table(orf3.std.stack$ind),function(x){1:x})))
sum(duplicated(orf3.std.stack$values)) # 4,995,302,784 are duplicated sequences: aggregate
orf3.std.stack<-aggregate(ID~values,orf3.std.stack,paste,collapse=";") # 2,680,803



# ORF1 mt
orf1.mt.split<-sapply(orf1.mt,function(x){
  a<-unlist(strsplit(x,"\\*"))
  a<-a[grep("M",a)]
  pos<-sapply(a,function(x){unlist(gregexpr("M",x))[1]})
  b<-substring(a,pos,nchar(a))
  return(b[nchar(b)>=7])
})
orf1.mt.stack<-stack(orf1.mt.split)
length(names(orf1.mt.split)[which(!names(orf1.mt.split) %in% orf1.mt.stack$ind)]) # 0 transcripts do not generate sequences >=7aa; drop this levels
orf1.mt.stack<-stack(orf1.mt.split,drop=TRUE)
orf1.mt.stack$ID<-paste0(orf1.mt.stack$ind,"_ORF1.",unlist(sapply(table(orf1.mt.stack$ind),function(x){1:x})))
sum(duplicated(orf1.mt.stack$values)) # 24 are duplicated sequences: aggregate
orf1.mt.stack<-aggregate(ID~values,orf1.mt.stack,paste,collapse=";") # 78

# ORF2 mt
orf2.mt.split<-sapply(orf2.mt,function(x){
  a<-unlist(strsplit(x,"\\*"))
  a<-a[grep("M",a)]
  pos<-sapply(a,function(x){unlist(gregexpr("M",x))[1]})
  b<-substring(a,pos,nchar(a))
  return(b[nchar(b)>=7])
})
orf2.mt.stack<-stack(orf2.mt.split)
length(names(orf2.mt.split)[which(!names(orf2.mt.split) %in% orf2.mt.stack$ind)]) # 0 transcripts do not generate sequences >=7aa; drop this levels
orf2.mt.stack<-stack(orf2.mt.split,drop=TRUE)
orf2.mt.stack$ID<-paste0(orf2.mt.stack$ind,"_ORF2.",unlist(sapply(table(orf2.mt.stack$ind),function(x){1:x})))
sum(duplicated(orf2.mt.stack$values)) # 47 are duplicated sequences: aggregate
orf2.mt.stack<-aggregate(ID~values,orf2.mt.stack,paste,collapse=";") # 80

# ORF3 mt
orf3.mt.split<-sapply(orf3.mt,function(x){
  a<-unlist(strsplit(x,"\\*"))
  a<-a[grep("M",a)]
  pos<-sapply(a,function(x){unlist(gregexpr("M",x))[1]})
  b<-substring(a,pos,nchar(a))
  return(b[nchar(b)>=7])
})
orf3.mt.stack<-stack(orf3.mt.split)
length(names(orf3.mt.split)[which(!names(orf3.mt.split) %in% orf3.mt.stack$ind)]) # 0 transcripts do not generate sequences >=7aa; drop this levels
orf3.mt.stack<-stack(orf3.mt.split,drop=TRUE)
orf3.mt.stack$ID<-paste0(orf3.mt.stack$ind,"_ORF3.",unlist(sapply(table(orf3.mt.stack$ind),function(x){1:x})))
sum(duplicated(orf3.mt.stack$values)) # 33 are duplicated sequences: aggregate
orf3.mt.stack<-aggregate(ID~values,orf3.mt.stack,paste,collapse=";") # 59


tcga.fa<-rbind(orf1.std.stack,orf2.std.stack,orf3.std.stack,
               orf1.mt.stack,orf2.mt.stack,orf3.mt.stack)
tcga.fa<-aggregate(ID~values,tcga.fa,paste,collapse=";") # 4,556,029
tcga.fa<-data.frame(ID=tcga.fa$ID,sequence=tcga.fa$values)

##### Write .tsv with unique sequences #####
write.table(x=tcga.fa,file="Fragpipe/inputs/tcgaAML_stringtie.tsv",row.names = FALSE,quote=FALSE,col.names=TRUE,sep="\t")


########################### FLAIR GTF GTEx #####################################
table(flair$type) # transcript/exon
table(flair$transcript_id)
View(as.data.frame(flair))
sum(!grepl("ENST",flair$transcript_id[flair$type=="transcript"])) # 72,325
length(flair$transcript_id[flair$type=="transcript"]) # 93,718
table(strand(flair))
table(seqnames(flair))

# First we need to subset the GTF
flair.df<-as.data.frame(flair)
flair.df<-flair.df[!grepl("ENST",flair.df$transcript_id) & # keep only novel transcripts
                     flair.df$type=="exon",] # if we do not subset, transcripts are translated from start to end
flair.df$strand<-droplevels(flair.df$strand) # drop * level because it raises an error in getTxsSeq
flair.df$transcript_id<-paste0("FLAIR_",flair.df$transcript_id)
length(unique(flair.df$transcript_id)) # 72,325 transcripts
table(flair.df$seqnames[!duplicated(flair.df$transcript_id)]) # Only 9 in the mitocondrial genome

# Get transcript sequence
flair.gr<-makeGRangesListFromDataFrame(flair.df,
                                       ignore.strand = FALSE,
                                       split.field="transcript_id",
                                       names.field="transcript_id") # gtf must be parsed as GRangesList to getTxsSeq

flair.tx.seq<-getTxsSeq(flair.gr,BSgenome.Hsapiens.UCSC.hg38) 


# Get 3 reading frames from each sequence
# Differentiate between std and mt genomes for next step
flair.seq.std<-lapply(1:3,function(i){substring(flair.tx.seq[names(flair.tx.seq) %in% unique(flair.df$transcript_id[flair.df$seqnames!="chrM"])],
                                                i,
                                                width(flair.tx.seq[names(flair.tx.seq) %in% unique(flair.df$transcript_id[flair.df$seqnames!="chrM"])]))})
flair.seq.mt<-lapply(1:3,function(i){substring(flair.tx.seq[names(flair.tx.seq) %in% unique(flair.df$transcript_id[flair.df$seqnames=="chrM"])],
                                               i,
                                               width(flair.tx.seq[names(flair.tx.seq) %in% unique(flair.df$transcript_id[flair.df$seqnames=="chrM"])]))})



# Check available genetic codes in Biostrings
SGC0<-getGeneticCode("SGC0")
SGC1<-getGeneticCode("SGC1")

# Translate (use no.init.codon=TRUE to avoid the use of alternative initiation codons)
# std: from ATG to TAA/TAG/TGA
orf1.std<-as.character(Biostrings::translate(DNAStringSet(flair.seq.std[[1]]),
                                             genetic.code = SGC0,
                                             if.fuzzy.codon = "solve",
                                             no.init.codon = TRUE))
orf2.std<-as.character(Biostrings::translate(DNAStringSet(flair.seq.std[[2]]),
                                             genetic.code = SGC0,
                                             if.fuzzy.codon = "solve",
                                             no.init.codon = TRUE))
orf3.std<-as.character(Biostrings::translate(DNAStringSet(flair.seq.std[[3]]),
                                             genetic.code = SGC0,
                                             if.fuzzy.codon = "solve",
                                             no.init.codon = TRUE))
# mt: from ATA/ATG to TAA/TAG/AGA/AGG)
orf1.mt<-as.character(Biostrings::translate(DNAStringSet(flair.seq.mt[[1]]),
                                            genetic.code = SGC1,
                                            if.fuzzy.codon = "solve",
                                            no.init.codon = TRUE))
orf2.mt<-as.character(Biostrings::translate(DNAStringSet(flair.seq.mt[[2]]),
                                            genetic.code = SGC1,
                                            if.fuzzy.codon = "solve",
                                            no.init.codon = TRUE))
orf3.mt<-as.character(Biostrings::translate(DNAStringSet(flair.seq.mt[[3]]),
                                            genetic.code = SGC1,
                                            if.fuzzy.codon = "solve",
                                            no.init.codon = TRUE))

##### Keep ORFs that are >=7aa #####
# 1. Split each transcript sequence by stop codon (*)
# 2. Substring each subsequence from first M (i.e. 1st ATG)
# 3. Report only sequences >=7aa
# ORF1 std
orf1.std.split<-sapply(orf1.std,function(x){
  a<-unlist(strsplit(x,"\\*"))
  a<-a[grep("M",a)]
  pos<-sapply(a,function(x){unlist(gregexpr("M",x))[1]})
  b<-substring(a,pos,nchar(a))
  return(b[nchar(b)>=7])
})
orf1.std.stack<-stack(orf1.std.split)
length(names(orf1.std.split)[which(!names(orf1.std.split) %in% orf1.std.stack$ind)]) # 1589 transcripts do not generate sequences >=7aa; drop this levels
orf1.std.stack<-stack(orf1.std.split,drop=TRUE)
sum(duplicated(orf1.std.stack$ind))
orf1.std.stack$ID<-paste0(orf1.std.stack$ind,"_ORF1.",unlist(sapply(table(orf1.std.stack$ind),function(x){1:x})))
sum(duplicated(orf1.std.stack$values)) # 165,067 are duplicated sequences: aggregate
orf1.std.stack<-aggregate(ID~values,orf1.std.stack,paste,collapse=";") # 238,898

# ORF2 std
orf2.std.split<-sapply(orf2.std,function(x){
  a<-unlist(strsplit(x,"\\*"))
  a<-a[grep("M",a)]
  pos<-sapply(a,function(x){unlist(gregexpr("M",x))[1]})
  b<-substring(a,pos,nchar(a))
  return(b[nchar(b)>=7])
})
orf2.std.stack<-stack(orf2.std.split)
length(names(orf2.std.split)[which(!names(orf2.std.split) %in% orf2.std.stack$ind)]) # 1414 transcripts do not generate sequences >=7aa; drop this levels
orf2.std.stack<-stack(orf2.std.split,drop=TRUE)
sum(duplicated(orf2.std.stack$ind))
orf2.std.stack$ID<-paste0(orf2.std.stack$ind,"_ORF2.",unlist(sapply(table(orf2.std.stack$ind),function(x){1:x})))
sum(duplicated(orf2.std.stack$values)) # 160,929 are duplicated sequences: aggregate
orf2.std.stack<-aggregate(ID~values,orf2.std.stack,paste,collapse=";") # 239,554

# ORF3 std
orf3.std.split<-sapply(orf3.std,function(x){
  a<-unlist(strsplit(x,"\\*"))
  a<-a[grep("M",a)]
  pos<-sapply(a,function(x){unlist(gregexpr("M",x))[1]})
  b<-substring(a,pos,nchar(a))
  return(b[nchar(b)>=7])
})
orf3.std.stack<-stack(orf3.std.split)
length(names(orf3.std.split)[which(!names(orf3.std.split) %in% orf3.std.stack$ind)]) # 1547 transcripts do not generate sequences >=7aa; drop this levels
orf3.std.stack<-stack(orf3.std.split,drop=TRUE)
sum(duplicated(orf3.std.stack$ind))
orf3.std.stack$ID<-paste0(orf3.std.stack$ind,"_ORF3.",unlist(sapply(table(orf3.std.stack$ind),function(x){1:x})))
sum(duplicated(orf3.std.stack$values)) # 163,274 are duplicated sequences: aggregate
orf3.std.stack<-aggregate(ID~values,orf3.std.stack,paste,collapse=";") # 238,401



# ORF1 mt
orf1.mt.split<-sapply(orf1.mt,function(x){
  a<-unlist(strsplit(x,"\\*"))
  a<-a[grep("M",a)]
  pos<-sapply(a,function(x){unlist(gregexpr("M",x))[1]})
  b<-substring(a,pos,nchar(a))
  return(b[nchar(b)>=7])
})
orf1.mt.stack<-stack(orf1.mt.split)
length(names(orf1.mt.split)[which(!names(orf1.mt.split) %in% orf1.mt.stack$ind)]) # 0 transcripts do not generate sequences >=7aa; drop this levels
orf1.mt.stack<-stack(orf1.mt.split,drop=TRUE)
sum(duplicated(orf1.mt.stack$ind))
orf1.mt.stack$ID<-paste0(orf1.mt.stack$ind,"_ORF1.",unlist(sapply(table(orf1.mt.stack$ind),function(x){1:x})))
sum(duplicated(orf1.mt.stack$values)) # 37 are duplicated sequences: aggregate
orf1.mt.stack<-aggregate(ID~values,orf1.mt.stack,paste,collapse=";") # 42

# ORF2 mt
orf2.mt.split<-sapply(orf2.mt,function(x){
  a<-unlist(strsplit(x,"\\*"))
  a<-a[grep("M",a)]
  pos<-sapply(a,function(x){unlist(gregexpr("M",x))[1]})
  b<-substring(a,pos,nchar(a))
  return(b[nchar(b)>=7])
})
orf2.mt.stack<-stack(orf2.mt.split)
length(names(orf2.mt.split)[which(!names(orf2.mt.split) %in% orf2.mt.stack$ind)]) # 0 transcripts do not generate sequences >=7aa; drop this levels
orf2.mt.stack<-stack(orf2.mt.split,drop=TRUE)
sum(duplicated(orf2.mt.stack$ind))
orf2.mt.stack$ID<-paste0(orf2.mt.stack$ind,"_ORF2.",unlist(sapply(table(orf2.mt.stack$ind),function(x){1:x})))
sum(duplicated(orf2.mt.stack$values)) # 41 are duplicated sequences: aggregate
orf2.mt.stack<-aggregate(ID~values,orf2.mt.stack,paste,collapse=";") # 42

# ORF3 mt
orf3.mt.split<-sapply(orf3.mt,function(x){
  a<-unlist(strsplit(x,"\\*"))
  a<-a[grep("M",a)]
  pos<-sapply(a,function(x){unlist(gregexpr("M",x))[1]})
  b<-substring(a,pos,nchar(a))
  return(b[nchar(b)>=7])
})
orf3.mt.stack<-stack(orf3.mt.split)
length(names(orf3.mt.split)[which(!names(orf3.mt.split) %in% orf3.mt.stack$ind)]) # 0 transcripts do not generate sequences >=7aa; drop this levels
orf3.mt.stack<-stack(orf3.mt.split,drop=TRUE)
sum(duplicated(orf3.mt.stack$ind))
orf3.mt.stack$ID<-paste0(orf3.mt.stack$ind,"_ORF3.",unlist(sapply(table(orf3.mt.stack$ind),function(x){1:x})))
sum(duplicated(orf3.mt.stack$values)) # 57 are duplicated sequences: aggregate
orf3.mt.stack<-aggregate(ID~values,orf3.mt.stack,paste,collapse=";") # 48


flair.fa<-rbind(orf1.std.stack,orf2.std.stack,orf3.std.stack,
                orf1.mt.stack,orf2.mt.stack,orf3.mt.stack)
flair.fa<-aggregate(ID~values,flair.fa,paste,collapse=";") # 532,741
flair.fa<-data.frame(ID=flair.fa$ID,sequence=flair.fa$values)


##### Write .tsv with unique sequences #####
write.table(x=flair.fa,file="Fragpipe/inputs/GTEx_FLAIR_novel.tsv",row.names = FALSE,quote=FALSE,col.names=TRUE,sep="\t")

###################### HL60 ONT GSM5331246, GSM5331247 #########################
hl60.1.df<-as.data.frame(hl60.1)
hl60.1.df<-hl60.1.df[hl60.1.df$type=="exon" &
                       hl60.1.df$transcript_status=="NOVEL" &
                       hl60.1.df$seqnames %in% c(paste0("chr",1:22),"chrX","chrY","chrM"),]

hl60.1.df$strand<-droplevels(hl60.1.df$strand) # drop * level because it raises an error in getTxsSeq
length(unique(hl60.1.df$transcript_id)) # 197,364 transcripts
table(hl60.1.df$seqnames[!duplicated(hl60.1.df$transcript_id)]) # Only 2 in the mitochondrial genome

hl60.2.df<-as.data.frame(hl60.2)
hl60.2.df<-hl60.2.df[hl60.2.df$type=="exon" &
                       hl60.2.df$transcript_status=="NOVEL" &
                       hl60.2.df$seqnames %in% c(paste0("chr",1:22),"chrX","chrY","chrM"),]

hl60.2.df$strand<-droplevels(hl60.2.df$strand) # drop * level because it raises an error in getTxsSeq
length(unique(hl60.2.df$transcript_id)) # 189,191 transcripts
table(hl60.2.df$seqnames[!duplicated(hl60.2.df$transcript_id)])


# Duplicates between samples
# We can not compare transcript IDs directly
# Do not remove to keep transcript ID for each sample
exon.annot.1<-paste0(hl60.1.df$seqnames,":",hl60.1.df$start,"-",hl60.1.df$end,":",hl60.1.df$strand)
a<-aggregate(exon.annot.1~transcript_id,hl60.1.df,paste,collapse=";")

exon.annot.2<-paste0(hl60.2.df$seqnames,":",hl60.2.df$start,"-",hl60.2.df$end,":",hl60.2.df$strand)
b<-aggregate(exon.annot.2~transcript_id,hl60.2.df,paste,collapse=";")

sum(a$exon.annot.1 %in% b$exon.annot.2) # 21,756
sum(b$exon.annot.2 %in% a$exon.annot.1) # 21,756
length(b$transcript_id[b$exon.annot.2 %in% a$exon.annot.1]) # 21,756

hl60.merged.df<-rbind(hl60.1.df,hl60.2.df)
length(unique(hl60.1.df$transcript_id))+
  length(unique(hl60.2.df$transcript_id))-
  length(b$transcript_id[b$exon.annot.2 %in% a$exon.annot.1]) # 364,799

length(unique(hl60.merged.df$transcript_id)) # 386,555


# Get transcript sequence
hl60.merged.gr<-makeGRangesListFromDataFrame(hl60.merged.df,
                                             ignore.strand = FALSE,
                                             split.field="transcript_id",
                                             names.field="transcript_id") # gtf must be parsed as GRangesList to getTxsSeq

hl60.merged.tx.seq<-getTxsSeq(hl60.merged.gr,BSgenome.Hsapiens.UCSC.hg38) 


# Get 3 reading frames from each sequence
# Differentiate between std and mt genomes for next step
hl60.merged.seq.std<-lapply(1:3,function(i){substring(hl60.merged.tx.seq[names(hl60.merged.tx.seq) %in% unique(hl60.merged.df$transcript_id[hl60.merged.df$seqnames!="chrM"])],
                                                      i,
                                                      width(hl60.merged.tx.seq[names(hl60.merged.tx.seq) %in% unique(hl60.merged.df$transcript_id[hl60.merged.df$seqnames!="chrM"])]))})
hl60.merged.seq.mt<-lapply(1:3,function(i){substring(hl60.merged.tx.seq[names(hl60.merged.tx.seq) %in% unique(hl60.merged.df$transcript_id[hl60.merged.df$seqnames=="chrM"])],
                                                     i,
                                                     width(hl60.merged.tx.seq[names(hl60.merged.tx.seq) %in% unique(hl60.merged.df$transcript_id[hl60.merged.df$seqnames=="chrM"])]))})



# Check available genetic codes in Biostrings
SGC0<-getGeneticCode("SGC0")
SGC1<-getGeneticCode("SGC1")

# Translate (use no.init.codon=TRUE to avoid the use of alternative initiation codons)
# std: from ATG to TAA/TAG/TGA
orf1.std<-as.character(Biostrings::translate(DNAStringSet(hl60.merged.seq.std[[1]]),
                                             genetic.code = SGC0,
                                             if.fuzzy.codon = "solve",
                                             no.init.codon = TRUE))
orf2.std<-as.character(Biostrings::translate(DNAStringSet(hl60.merged.seq.std[[2]]),
                                             genetic.code = SGC0,
                                             if.fuzzy.codon = "solve",
                                             no.init.codon = TRUE))
orf3.std<-as.character(Biostrings::translate(DNAStringSet(hl60.merged.seq.std[[3]]),
                                             genetic.code = SGC0,
                                             if.fuzzy.codon = "solve",
                                             no.init.codon = TRUE))
# mt: from ATA/ATG to TAA/TAG/AGA/AGG)
orf1.mt<-as.character(Biostrings::translate(DNAStringSet(hl60.merged.seq.mt[[1]]),
                                            genetic.code = SGC1,
                                            if.fuzzy.codon = "solve",
                                            no.init.codon = TRUE))
orf2.mt<-as.character(Biostrings::translate(DNAStringSet(hl60.merged.seq.mt[[2]]),
                                            genetic.code = SGC1,
                                            if.fuzzy.codon = "solve",
                                            no.init.codon = TRUE))
orf3.mt<-as.character(Biostrings::translate(DNAStringSet(hl60.merged.seq.mt[[3]]),
                                            genetic.code = SGC1,
                                            if.fuzzy.codon = "solve",
                                            no.init.codon = TRUE))

##### Keep ORFs that are >=7aa #####
# 1. Split each transcript sequence by stop codon (*)
# 2. Substring each subsequence from first M (i.e. 1st ATG)
# 3. Report only sequences >=7aa
# ORF1 std
orf1.std.split<-sapply(orf1.std,function(x){
  a<-unlist(strsplit(x,"\\*"))
  a<-a[grep("M",a)]
  pos<-sapply(a,function(x){unlist(gregexpr("M",x))[1]})
  b<-substring(a,pos,nchar(a))
  return(b[nchar(b)>=7])
})
orf1.std.stack<-stack(orf1.std.split)
length(names(orf1.std.split)[which(!names(orf1.std.split) %in% orf1.std.stack$ind)]) # 7,034 transcripts do not generate sequences >=7aa; drop this levels
orf1.std.stack<-stack(orf1.std.split,drop=TRUE)
sum(duplicated(orf1.std.stack$ind))
orf1.std.stack$ID<-paste0(orf1.std.stack$ind,"_ORF1.",unlist(sapply(table(orf1.std.stack$ind),function(x){1:x})))
sum(duplicated(orf1.std.stack$values)) # 1,470,187 are duplicated sequences: aggregate
orf1.std.stack<-aggregate(ID~values,orf1.std.stack,paste,collapse=";") # 863,677

# ORF2 std
orf2.std.split<-sapply(orf2.std,function(x){
  a<-unlist(strsplit(x,"\\*"))
  a<-a[grep("M",a)]
  pos<-sapply(a,function(x){unlist(gregexpr("M",x))[1]})
  b<-substring(a,pos,nchar(a))
  return(b[nchar(b)>=7])
})
orf2.std.stack<-stack(orf2.std.split)
length(names(orf2.std.split)[which(!names(orf2.std.split) %in% orf2.std.stack$ind)]) # 7165 transcripts do not generate sequences >=7aa; drop this levels
orf2.std.stack<-stack(orf2.std.split,drop=TRUE)
sum(duplicated(orf2.std.stack$ind))
orf2.std.stack$ID<-paste0(orf2.std.stack$ind,"_ORF2.",unlist(sapply(table(orf2.std.stack$ind),function(x){1:x})))
sum(duplicated(orf2.std.stack$values)) # 1,501,962 are duplicated sequences: aggregate
orf2.std.stack<-aggregate(ID~values,orf2.std.stack,paste,collapse=";") # 858,215

# ORF3 std
orf3.std.split<-sapply(orf3.std,function(x){
  a<-unlist(strsplit(x,"\\*"))
  a<-a[grep("M",a)]
  pos<-sapply(a,function(x){unlist(gregexpr("M",x))[1]})
  b<-substring(a,pos,nchar(a))
  return(b[nchar(b)>=7])
})
orf3.std.stack<-stack(orf3.std.split)
length(names(orf3.std.split)[which(!names(orf3.std.split) %in% orf3.std.stack$ind)]) # 7714 transcripts do not generate sequences >=7aa; drop this levels
orf3.std.stack<-stack(orf3.std.split,drop=TRUE)
sum(duplicated(orf3.std.stack$ind))
orf3.std.stack$ID<-paste0(orf3.std.stack$ind,"_ORF3.",unlist(sapply(table(orf3.std.stack$ind),function(x){1:x})))
sum(duplicated(orf3.std.stack$values)) # 1,446,912 are duplicated sequences: aggregate
orf3.std.stack<-aggregate(ID~values,orf3.std.stack,paste,collapse=";") # 858,270



# ORF1 mt
orf1.mt.split<-sapply(orf1.mt,function(x){
  a<-unlist(strsplit(x,"\\*"))
  a<-a[grep("M",a)]
  pos<-sapply(a,function(x){unlist(gregexpr("M",x))[1]})
  b<-substring(a,pos,nchar(a))
  return(b[nchar(b)>=7])
})
orf1.mt.stack<-stack(orf1.mt.split)
length(names(orf1.mt.split)[which(!names(orf1.mt.split) %in% orf1.mt.stack$ind)]) # 0 transcripts do not generate sequences >=7aa; drop this levels
orf1.mt.stack<-stack(orf1.mt.split,drop=TRUE)
sum(duplicated(orf1.mt.stack$ind))
orf1.mt.stack$ID<-paste0(orf1.mt.stack$ind,"_ORF1.",unlist(sapply(table(orf1.mt.stack$ind),function(x){1:x})))
sum(duplicated(orf1.mt.stack$values)) # 7 are duplicated sequences: aggregate
orf1.mt.stack<-aggregate(ID~values,orf1.mt.stack,paste,collapse=";") # 14

# ORF2 mt
orf2.mt.split<-sapply(orf2.mt,function(x){
  a<-unlist(strsplit(x,"\\*"))
  a<-a[grep("M",a)]
  pos<-sapply(a,function(x){unlist(gregexpr("M",x))[1]})
  b<-substring(a,pos,nchar(a))
  return(b[nchar(b)>=7])
})
orf2.mt.stack<-stack(orf2.mt.split)
length(names(orf2.mt.split)[which(!names(orf2.mt.split) %in% orf2.mt.stack$ind)]) # 0 transcripts do not generate sequences >=7aa; drop this levels
orf2.mt.stack<-stack(orf2.mt.split,drop=TRUE)
sum(duplicated(orf2.mt.stack$ind))
orf2.mt.stack$ID<-paste0(orf2.mt.stack$ind,"_ORF2.",unlist(sapply(table(orf2.mt.stack$ind),function(x){1:x})))
sum(duplicated(orf2.mt.stack$values)) # 9 are duplicated sequences: aggregate
orf2.mt.stack<-aggregate(ID~values,orf2.mt.stack,paste,collapse=";") # 18

# ORF3 mt
orf3.mt.split<-sapply(orf3.mt,function(x){
  a<-unlist(strsplit(x,"\\*"))
  a<-a[grep("M",a)]
  pos<-sapply(a,function(x){unlist(gregexpr("M",x))[1]})
  b<-substring(a,pos,nchar(a))
  return(b[nchar(b)>=7])
})
orf3.mt.stack<-stack(orf3.mt.split)
length(names(orf3.mt.split)[which(!names(orf3.mt.split) %in% orf3.mt.stack$ind)]) # 0 transcripts do not generate sequences >=7aa; drop this levels
orf3.mt.stack<-stack(orf3.mt.split,drop=TRUE)
sum(duplicated(orf3.mt.stack$ind))
orf3.mt.stack$ID<-paste0(orf3.mt.stack$ind,"_ORF3.",unlist(sapply(table(orf3.mt.stack$ind),function(x){1:x})))
sum(duplicated(orf3.mt.stack$values)) # 4 are duplicated sequences: aggregate
orf3.mt.stack<-aggregate(ID~values,orf3.mt.stack,paste,collapse=";") # 11


hl60.merged.fa<-rbind(orf1.std.stack,orf2.std.stack,orf3.std.stack,
                      orf1.mt.stack,orf2.mt.stack,orf3.mt.stack)
hl60.merged.fa<-aggregate(ID~values,hl60.merged.fa,paste,collapse=";") # 1,911,822
hl60.merged.fa<-data.frame(ID=hl60.merged.fa$ID,sequence=hl60.merged.fa$values)

##### Write .tsv with unique sequences #####
write.table(x=hl60.merged.fa,file="Fragpipe/inputs/HL60_novel.tsv",row.names = FALSE,quote=FALSE,col.names=TRUE,sep="\t")

############################## GENCODE v36 GTF ##################################
table(strand(gencode))
table(seqnames(gencode)) # all chr are 1:22, X, Y, M
table(gencode$type) # keep exon, CD, 
# First we need to subset the GTF
gencode.df<-as.data.frame(gencode)
gencode.df<-gencode.df[gencode.df$type=="exon",] # if we do not subset, transcripts are translated from start to end
gencode.df$strand<-droplevels(gencode.df$strand) # drop * level because it raises an error in getTxsSeq
gencode.df$transcript_id<-paste0("GENCODEv36_",gencode.df$transcript_id) # Add this prefix for further identification
length(unique(gencode.df$transcript_id)) # 232,117 transcripts
table(gencode.df$seqnames[!duplicated(gencode.df$transcript_id)]) # Only 37 in the mitocondrial genome

# Get transcript sequence
gencode.gr<-makeGRangesListFromDataFrame(gencode.df,
                                         ignore.strand = FALSE,
                                         split.field="transcript_id",
                                         names.field="transcript_id") # gtf must be parsed as GRangesList to getTxsSeq

gencode.tx.seq<-getTxsSeq(gencode.gr,BSgenome.Hsapiens.UCSC.hg38) 


# Get 3 reading frames from each sequence
# Differentiate between std and mt genomes for next step
gencode.seq.std<-lapply(1:3,function(i){substring(gencode.tx.seq[names(gencode.tx.seq) %in% unique(gencode.df$transcript_id[gencode.df$seqnames!="chrM"])],
                                                  i,
                                                  width(gencode.tx.seq[names(gencode.tx.seq) %in% unique(gencode.df$transcript_id[gencode.df$seqnames!="chrM"])]))})
gencode.seq.mt<-lapply(1:3,function(i){substring(gencode.tx.seq[names(gencode.tx.seq) %in% unique(gencode.df$transcript_id[gencode.df$seqnames=="chrM"])],
                                                 i,
                                                 width(gencode.tx.seq[names(gencode.tx.seq) %in% unique(gencode.df$transcript_id[gencode.df$seqnames=="chrM"])]))})



# Check available genetic codes in Biostrings
SGC0<-getGeneticCode("SGC0")
SGC1<-getGeneticCode("SGC1")

# Translate (use no.init.codon=TRUE to avoid the use of alternative initiation codons)
# std: from ATG to TAA/TAG/TGA
orf1.std<-as.character(Biostrings::translate(DNAStringSet(gencode.seq.std[[1]]),
                                             genetic.code = SGC0,
                                             if.fuzzy.codon = "solve",
                                             no.init.codon = TRUE))
orf2.std<-as.character(Biostrings::translate(DNAStringSet(gencode.seq.std[[2]]),
                                             genetic.code = SGC0,
                                             if.fuzzy.codon = "solve",
                                             no.init.codon = TRUE))
orf3.std<-as.character(Biostrings::translate(DNAStringSet(gencode.seq.std[[3]]),
                                             genetic.code = SGC0,
                                             if.fuzzy.codon = "solve",
                                             no.init.codon = TRUE))
# mt: from ATA/ATG to TAA/TAG/AGA/AGG)
orf1.mt<-as.character(Biostrings::translate(DNAStringSet(gencode.seq.mt[[1]]),
                                            genetic.code = SGC1,
                                            if.fuzzy.codon = "solve",
                                            no.init.codon = TRUE))
orf2.mt<-as.character(Biostrings::translate(DNAStringSet(gencode.seq.mt[[2]]),
                                            genetic.code = SGC1,
                                            if.fuzzy.codon = "solve",
                                            no.init.codon = TRUE))
orf3.mt<-as.character(Biostrings::translate(DNAStringSet(gencode.seq.mt[[3]]),
                                            genetic.code = SGC1,
                                            if.fuzzy.codon = "solve",
                                            no.init.codon = TRUE))

##### Keep ORFs that are >=7aa #####
# 1. Split each transcript sequence by stop codon (*)
# 2. Substring each subsequence from first M (i.e. 1st ATG)
# 3. Report only sequences >=7aa
# ORF1 std
orf1.std.split<-sapply(orf1.std,function(x){
  a<-unlist(strsplit(x,"\\*"))
  a<-a[grep("M",a)]
  pos<-sapply(a,function(x){unlist(gregexpr("M",x))[1]})
  b<-substring(a,pos,nchar(a))
  return(b[nchar(b)>=7])
})
orf1.std.stack<-stack(orf1.std.split)
length(names(orf1.std.split)[which(!names(orf1.std.split) %in% orf1.std.stack$ind)]) # 22,444 transcripts do not generate sequences >=7aa; drop this levels
orf1.std.stack<-stack(orf1.std.split,drop=TRUE)
orf1.std.stack$ID<-paste0(orf1.std.stack$ind,"_ORF1.",unlist(sapply(table(orf1.std.stack$ind),function(x){1:x})))
sum(duplicated(orf1.std.stack$values)) # 300,953 are duplicated sequences: aggregate
orf1.std.stack<-aggregate(ID~values,orf1.std.stack,paste,collapse=";") # 643,329

# ORF2 std
orf2.std.split<-sapply(orf2.std,function(x){
  a<-unlist(strsplit(x,"\\*"))
  a<-a[grep("M",a)]
  pos<-sapply(a,function(x){unlist(gregexpr("M",x))[1]})
  b<-substring(a,pos,nchar(a))
  return(b[nchar(b)>=7])
})
orf2.std.stack<-stack(orf2.std.split)
length(names(orf2.std.split)[which(!names(orf2.std.split) %in% orf2.std.stack$ind)]) # 24,006 transcripts do not generate sequences >=7aa; drop this levels
orf2.std.stack<-stack(orf2.std.split,drop=TRUE)
orf2.std.stack$ID<-paste0(orf2.std.stack$ind,"_ORF2.",unlist(sapply(table(orf2.std.stack$ind),function(x){1:x})))
sum(duplicated(orf2.std.stack$values)) # 321,110 are duplicated sequences: aggregate
orf2.std.stack<-aggregate(ID~values,orf2.std.stack,paste,collapse=";") # 645,573

# ORF3 std
orf3.std.split<-sapply(orf3.std,function(x){
  a<-unlist(strsplit(x,"\\*"))
  a<-a[grep("M",a)]
  pos<-sapply(a,function(x){unlist(gregexpr("M",x))[1]})
  b<-substring(a,pos,nchar(a))
  return(b[nchar(b)>=7])
})
orf3.std.stack<-stack(orf3.std.split)
length(names(orf3.std.split)[which(!names(orf3.std.split) %in% orf3.std.stack$ind)]) # 26,125 transcripts do not generate sequences >=7aa; drop this levels
orf3.std.stack<-stack(orf3.std.split,drop=TRUE)
orf3.std.stack$ID<-paste0(orf3.std.stack$ind,"_ORF3.",unlist(sapply(table(orf3.std.stack$ind),function(x){1:x})))
sum(duplicated(orf3.std.stack$values)) # 305,416 are duplicated sequences: aggregate
orf3.std.stack<-aggregate(ID~values,orf3.std.stack,paste,collapse=";") # 635,639



# ORF1 mt
orf1.mt.split<-sapply(orf1.mt,function(x){
  a<-unlist(strsplit(x,"\\*"))
  a<-a[grep("M",a)]
  pos<-sapply(a,function(x){unlist(gregexpr("M",x))[1]})
  b<-substring(a,pos,nchar(a))
  return(b[nchar(b)>=7])
})
orf1.mt.stack<-stack(orf1.mt.split)
length(names(orf1.mt.split)[which(!names(orf1.mt.split) %in% orf1.mt.stack$ind)]) # 11 transcripts do not generate sequences >=7aa; drop this levels
orf1.mt.stack<-stack(orf1.mt.split,drop=TRUE)
orf1.mt.stack$ID<-paste0(orf1.mt.stack$ind,"_ORF1.",unlist(sapply(table(orf1.mt.stack$ind),function(x){1:x})))
sum(duplicated(orf1.mt.stack$values)) # 0 are duplicated sequences: aggregate
orf1.mt.stack<-aggregate(ID~values,orf1.mt.stack,paste,collapse=";") # 33

# ORF2 mt
orf2.mt.split<-sapply(orf2.mt,function(x){
  a<-unlist(strsplit(x,"\\*"))
  a<-a[grep("M",a)]
  pos<-sapply(a,function(x){unlist(gregexpr("M",x))[1]})
  b<-substring(a,pos,nchar(a))
  return(b[nchar(b)>=7])
})
orf2.mt.stack<-stack(orf2.mt.split)
length(names(orf2.mt.split)[which(!names(orf2.mt.split) %in% orf2.mt.stack$ind)]) # 18 transcripts do not generate sequences >=7aa; drop this levels
orf2.mt.stack<-stack(orf2.mt.split,drop=TRUE)
orf2.mt.stack$ID<-paste0(orf2.mt.stack$ind,"_ORF2.",unlist(sapply(table(orf2.mt.stack$ind),function(x){1:x})))
sum(duplicated(orf2.mt.stack$values)) # 0 are duplicated sequences: aggregate
orf2.mt.stack<-aggregate(ID~values,orf2.mt.stack,paste,collapse=";") # 46

# ORF3 mt
orf3.mt.split<-sapply(orf3.mt,function(x){
  a<-unlist(strsplit(x,"\\*"))
  a<-a[grep("M",a)]
  pos<-sapply(a,function(x){unlist(gregexpr("M",x))[1]})
  b<-substring(a,pos,nchar(a))
  return(b[nchar(b)>=7])
})
orf3.mt.stack<-stack(orf3.mt.split)
length(names(orf3.mt.split)[which(!names(orf3.mt.split) %in% orf3.mt.stack$ind)]) # 8 transcripts do not generate sequences >=7aa; drop this levels
orf3.mt.stack<-stack(orf3.mt.split,drop=TRUE)
orf3.mt.stack$ID<-paste0(orf3.mt.stack$ind,"_ORF3.",unlist(sapply(table(orf3.mt.stack$ind),function(x){1:x})))
sum(duplicated(orf3.mt.stack$values)) # 0 are duplicated sequences: aggregate
orf3.mt.stack<-aggregate(ID~values,orf3.mt.stack,paste,collapse=";") # 61


gencode.fa<-rbind(orf1.std.stack,orf2.std.stack,orf3.std.stack,
                  orf1.mt.stack,orf2.mt.stack,orf3.mt.stack)
gencode.fa<-aggregate(ID~values,gencode.fa,paste,collapse=";") # 1,427,112
gencode.fa<-data.frame(ID=gencode.fa$ID,sequence=gencode.fa$values)

##### Write .tsv with unique sequences #####
write.table(x=gencode.fa,file="Fragpipe/inputs/gencodev36_3RFs.tsv",row.names = FALSE,quote=FALSE,col.names=TRUE,sep="\t")

############################### IRIS TCGA ######################################
# File list (tsv files)
file.list<-list.files("C:/Users/bea_f/OneDrive/Escritorio/HematoLaFe/Neoantigens-RNA-splicing/IRIS/",
                      recursive = TRUE,
                      full.names = TRUE,
                      pattern = ".txt$")
SE.files<-file.list[grep("SE",basename(file.list))]
A5SS.files<-file.list[grep("A5SS",basename(file.list))]
A3SS.files<-file.list[grep("A3SS",basename(file.list))]
RI.files<-file.list[grep("RI",basename(file.list))]

#################################### SE ########################################
SE.read <- lapply(SE.files, fread)
# Collapse all files
SE.events<-lapply(SE.read,function(x){x[,1:8]})
SE.events<-rbindlist(SE.events,idcol="id")
SE.events$annotation<-with(SE.events,paste(chr,strand,exonStart,exonEnd,upstreamEE,downstreamES,sep = "|"))
sum(duplicated(SE.events$annotation)) # 6,524,014
SE<-SE.events[!duplicated(SE.events$annotation),] # 737,589

# Translate
# As we do not have UE start and DE end we will only translate 8 aa upstream and downstream
range(abs(SE$exonStart-SE$exonEnd)) # there are exons as short as 1bp
UE.gr<-GRanges(seqnames=SE$chr,
               ranges=IRanges(start=SE$upstreamEE-17,end=SE$upstreamEE),
               strand=SE$strand)
SE.gr<-GRanges(seqnames=SE$chr,
               ranges=IRanges(start=SE$exonStart+1,end=SE$exonEnd),
               strand=SE$strand)
DE.gr<-GRanges(seqnames=SE$chr,
               ranges=IRanges(start=SE$downstreamES+1,end=SE$downstreamES+18),
               strand=SE$strand)

UE.rna<-getSeq(BSgenome.Hsapiens.UCSC.hg19,UE.gr)
SE.rna<-getSeq(BSgenome.Hsapiens.UCSC.hg19,SE.gr)
DE.rna<-getSeq(BSgenome.Hsapiens.UCSC.hg19,DE.gr)

rna.inc<-ifelse(SE$strand=="+",
                paste0(UE.rna,SE.rna,DE.rna),
                paste0(DE.rna,SE.rna,UE.rna))
rna.exc<-ifelse(SE$strand=="+",
                paste0(UE.rna,DE.rna),
                paste0(DE.rna,UE.rna))

names(rna.inc)<-paste0("IRIS_",SE$annotation,"_SEinc")
names(rna.exc)<-paste0("IRIS_",SE$annotation,"_SEexc")


length(rna.inc)+length(rna.exc) # 1,475,178

inc.seq<-lapply(1:3,function(i){substring(rna.inc,i,nchar(rna.inc))})
exc.seq<-lapply(1:3,function(i){substring(rna.exc,i,nchar(rna.exc))})

orf1.inc<-as.character(Biostrings::translate(DNAStringSet(inc.seq[[1]]),
                                             if.fuzzy.codon = "solve",
                                             no.init.codon = TRUE))
orf2.inc<-as.character(Biostrings::translate(DNAStringSet(inc.seq[[2]]),
                                             if.fuzzy.codon = "solve",
                                             no.init.codon = TRUE))
orf3.inc<-as.character(Biostrings::translate(DNAStringSet(inc.seq[[3]]),
                                             if.fuzzy.codon = "solve",
                                             no.init.codon = TRUE))

orf1.exc<-as.character(Biostrings::translate(DNAStringSet(exc.seq[[1]]),
                                             if.fuzzy.codon = "solve",
                                             no.init.codon = TRUE))
orf2.exc<-as.character(Biostrings::translate(DNAStringSet(exc.seq[[2]]),
                                             if.fuzzy.codon = "solve",
                                             no.init.codon = TRUE))
orf3.exc<-as.character(Biostrings::translate(DNAStringSet(exc.seq[[3]]),
                                             if.fuzzy.codon = "solve",
                                             no.init.codon = TRUE))

##### Keep ORFs that are >=7aa #####
# ORF1
orf1.inc.split<-sapply(orf1.inc,function(x){
  a<-unlist(strsplit(x,"\\*"))
  nchar(a[1])>=7
  b<-a[-1]
  b<-b[grep("M",b)]
  pos<-sapply(b,function(x){unlist(gregexpr("M",b))[1]})
  c<-substring(b,pos,nchar(b))
  return(c(a[1],c)[nchar(c(a[1],c))>=7])
})
orf1.inc.stack<-stack(orf1.inc.split)
length(names(orf1.inc.split)[which(!names(orf1.inc.split) %in% orf1.inc.stack$ind)]) # 62,780 events do not generate sequences >=7aa; drop this levels
orf1.inc.stack<-stack(orf1.inc.split,drop=TRUE)
orf1.inc.stack$ID<-paste0(orf1.inc.stack$ind,"_ORF1.",unlist(sapply(table(orf1.inc.stack$ind),function(x){1:x})))
sum(duplicated(orf1.inc.stack$values)) # 212,487 are duplicated sequences: aggregate
orf1.inc.stack<-aggregate(ID~values,orf1.inc.stack,paste,collapse=";")

# ORF1
orf1.exc.split<-sapply(orf1.exc,function(x){
  a<-unlist(strsplit(x,"\\*"))
  nchar(a[1])>=7
  b<-a[-1]
  b<-b[grep("M",b)]
  pos<-sapply(b,function(x){unlist(gregexpr("M",b))[1]})
  c<-substring(b,pos,nchar(b))
  return(c(a[1],c)[nchar(c(a[1],c))>=7])
})
orf1.exc.stack<-stack(orf1.exc.split)
length(names(orf1.exc.split)[which(!names(orf1.exc.split) %in% orf1.exc.stack$ind)]) # 110,143 events do not generate sequences >=7aa; drop this levels
orf1.exc.stack<-stack(orf1.exc.split,drop=TRUE)
orf1.exc.stack$ID<-paste0(orf1.exc.stack$ind,"_ORF1.",unlist(sapply(table(orf1.exc.stack$ind),function(x){1:x})))
sum(duplicated(orf1.exc.stack$values)) # 241,685 are duplicated sequences: aggregate
orf1.exc.stack<-aggregate(ID~values,orf1.exc.stack,paste,collapse=";")

# ORF2
orf2.inc.split<-sapply(orf2.inc,function(x){
  a<-unlist(strsplit(x,"\\*"))
  nchar(a[1])>=7
  b<-a[-1]
  b<-b[grep("M",b)]
  pos<-sapply(b,function(x){unlist(gregexpr("M",b))[1]})
  c<-substring(b,pos,nchar(b))
  return(c(a[1],c)[nchar(c(a[1],c))>=7])
})
orf2.inc.stack<-stack(orf2.inc.split)
length(names(orf2.inc.split)[which(!names(orf2.inc.split) %in% orf2.inc.stack$ind)]) # 99,752 events do not generate sequences >=7aa; drop this levels
orf2.inc.stack<-stack(orf2.inc.split,drop=TRUE)
orf2.inc.stack$ID<-paste0(orf2.inc.stack$ind,"_ORF2.",unlist(sapply(table(orf2.inc.stack$ind),function(x){1:x})))
sum(duplicated(orf2.inc.stack$values)) # 265,478 are duplicated sequences: aggregate
orf2.inc.stack<-aggregate(ID~values,orf2.inc.stack,paste,collapse=";")

# ORF2
orf2.exc.split<-sapply(orf2.exc,function(x){
  a<-unlist(strsplit(x,"\\*"))
  nchar(a[1])>=7
  b<-a[-1]
  b<-b[grep("M",b)]
  pos<-sapply(b,function(x){unlist(gregexpr("M",b))[1]})
  c<-substring(b,pos,nchar(b))
  return(c(a[1],c)[nchar(c(a[1],c))>=7])
})
orf2.exc.stack<-stack(orf2.exc.split)
length(names(orf2.exc.split)[which(!names(orf2.exc.split) %in% orf2.exc.stack$ind)]) # 171,504 events do not generate sequences >=7aa; drop this levels
orf2.exc.stack<-stack(orf2.exc.split,drop=TRUE)
orf2.exc.stack$ID<-paste0(orf2.exc.stack$ind,"_ORF2.",unlist(sapply(table(orf2.exc.stack$ind),function(x){1:x})))
sum(duplicated(orf2.exc.stack$values)) # 217,858 are duplicated sequences: aggregate
orf2.exc.stack<-aggregate(ID~values,orf2.exc.stack,paste,collapse=";")


# ORF3
orf3.inc.split<-sapply(orf3.inc,function(x){
  a<-unlist(strsplit(x,"\\*"))
  nchar(a[1])>=7
  b<-a[-1]
  b<-b[grep("M",b)]
  pos<-sapply(b,function(x){unlist(gregexpr("M",b))[1]})
  c<-substring(b,pos,nchar(b))
  return(c(a[1],c)[nchar(c(a[1],c))>=7])
})
orf3.inc.stack<-stack(orf3.inc.split)
length(names(orf3.inc.split)[which(!names(orf3.inc.split) %in% orf3.inc.stack$ind)]) # 111,677 events do not generate sequences >=7aa; drop this levels
orf3.inc.stack<-stack(orf3.inc.split,drop=TRUE)
orf3.inc.stack$ID<-paste0(orf3.inc.stack$ind,"_ORF3.",unlist(sapply(table(orf3.inc.stack$ind),function(x){1:x})))
sum(duplicated(orf3.inc.stack$values)) # 212,211 are duplicated sequences: aggregate
orf3.inc.stack<-aggregate(ID~values,orf3.inc.stack,paste,collapse=";")

# ORF3
orf3.exc.split<-sapply(orf3.exc,function(x){
  a<-unlist(strsplit(x,"\\*"))
  nchar(a[1])>=7
  b<-a[-1]
  b<-b[grep("M",b)]
  pos<-sapply(b,function(x){unlist(gregexpr("M",b))[1]})
  c<-substring(b,pos,nchar(b))
  return(c(a[1],c)[nchar(c(a[1],c))>=7])
})
orf3.exc.stack<-stack(orf3.exc.split)
length(names(orf3.exc.split)[which(!names(orf3.exc.split) %in% orf3.exc.stack$ind)]) # 161,108 events do not generate sequences >=7aa; drop this levels
orf3.exc.stack<-stack(orf3.exc.split,drop=TRUE)
orf3.exc.stack$ID<-paste0(orf3.exc.stack$ind,"_ORF3.",unlist(sapply(table(orf3.exc.stack$ind),function(x){1:x})))
sum(duplicated(orf3.exc.stack$values)) # 219,915 are duplicated sequences: aggregate
orf3.exc.stack<-aggregate(ID~values,orf3.exc.stack,paste,collapse=";")

SE.fa<-rbind(orf1.inc.stack,orf1.exc.stack,orf2.inc.stack,orf2.exc.stack,orf3.inc.stack,orf3.exc.stack)
SE.fa<-aggregate(ID~values,SE.fa,paste,collapse=";") # 2,642,199
SE.fa<-data.frame(ID=SE.fa$ID,sequence=SE.fa$values)

##### Write .tsv with unique sequences #####
write.table(x=SE.fa,file="Fragpipe/inputs/IRIS_TCGA_SE.tsv",row.names = FALSE,quote=FALSE,col.names=TRUE,sep="\t")



#################################### A5SS ########################################
A5SS.read <- lapply(A5SS.files, fread)
# Collapse all files
A5SS.events<-lapply(A5SS.read,function(x){x[,1:10]})
A5SS.events<-rbindlist(A5SS.events,idcol="id")
A5SS.events$annotation<-with(A5SS.events,paste(chr,strand,longExonStart,longExonEnd,shortES,shortEE,flankingES,flankingEE,sep = "|"))
sum(duplicated(A5SS.events$annotation)) # 260,187
A5SS<-A5SS.events[!duplicated(A5SS.events$annotation),] # 17,050

# Translate (in this case we have complete coordinates)
range(abs(A5SS$longExonStart-A5SS$longExonEnd))
range(abs(A5SS$shortES-A5SS$shortEE))
range(abs(A5SS$flankingES-A5SS$flankingEE)) # there are exons as short as 1bp

long.gr<-GRanges(seqnames=A5SS$chr,
                 ranges=IRanges(start=A5SS$longExonStart+1,end=A5SS$longExonEnd),
                 strand=A5SS$strand)
short.gr<-GRanges(seqnames=A5SS$chr,
                  ranges=IRanges(start=A5SS$shortES+1,end=A5SS$shortEE),
                  strand=A5SS$strand)
flanking.gr<-GRanges(seqnames=A5SS$chr,
                     ranges=IRanges(start=A5SS$flankingES+1,end=A5SS$flankingEE),
                     strand=A5SS$strand)

long.rna<-getSeq(BSgenome.Hsapiens.UCSC.hg19,long.gr)
short.rna<-getSeq(BSgenome.Hsapiens.UCSC.hg19,short.gr)
flanking.rna<-getSeq(BSgenome.Hsapiens.UCSC.hg19,flanking.gr)

rna.inc<-paste0(long.rna,flanking.rna)
rna.exc<-paste0(short.rna,flanking.rna)

names(rna.inc)<-paste0("IRIS_",A5SS$annotation,"_A5Ssinc")
names(rna.exc)<-paste0("IRIS_",A5SS$annotation,"_A5SSexc")


length(rna.inc)+length(rna.exc) # 34,100

inc.seq<-lapply(1:3,function(i){substring(rna.inc,i,nchar(rna.inc))})
exc.seq<-lapply(1:3,function(i){substring(rna.exc,i,nchar(rna.exc))})

orf1.inc<-as.character(Biostrings::translate(DNAStringSet(inc.seq[[1]]),
                                             if.fuzzy.codon = "solve",
                                             no.init.codon = TRUE))
orf2.inc<-as.character(Biostrings::translate(DNAStringSet(inc.seq[[2]]),
                                             if.fuzzy.codon = "solve",
                                             no.init.codon = TRUE))
orf3.inc<-as.character(Biostrings::translate(DNAStringSet(inc.seq[[3]]),
                                             if.fuzzy.codon = "solve",
                                             no.init.codon = TRUE))

orf1.exc<-as.character(Biostrings::translate(DNAStringSet(exc.seq[[1]]),
                                             if.fuzzy.codon = "solve",
                                             no.init.codon = TRUE))
orf2.exc<-as.character(Biostrings::translate(DNAStringSet(exc.seq[[2]]),
                                             if.fuzzy.codon = "solve",
                                             no.init.codon = TRUE))
orf3.exc<-as.character(Biostrings::translate(DNAStringSet(exc.seq[[3]]),
                                             if.fuzzy.codon = "solve",
                                             no.init.codon = TRUE))

##### Keep ORFs that are >=7aa #####
# ORF1
orf1.inc.split<-sapply(orf1.inc,function(x){
  a<-unlist(strsplit(x,"\\*"))
  nchar(a[1])>=7
  b<-a[-1]
  b<-b[grep("M",b)]
  pos<-sapply(b,function(x){unlist(gregexpr("M",b))[1]})
  c<-substring(b,pos,nchar(b))
  return(c(a[1],c)[nchar(c(a[1],c))>=7])
})
orf1.inc.stack<-stack(orf1.inc.split)
length(names(orf1.inc.split)[which(!names(orf1.inc.split) %in% orf1.inc.stack$ind)]) # 745 events do not generate sequences >=7aa; drop this levels
orf1.inc.stack<-stack(orf1.inc.split,drop=TRUE)
orf1.inc.stack$ID<-paste0(orf1.inc.stack$ind,"_ORF1.",unlist(sapply(table(orf1.inc.stack$ind),function(x){1:x})))
sum(duplicated(orf1.inc.stack$values)) # 10,386 are duplicated sequences: aggregate
orf1.inc.stack<-aggregate(ID~values,orf1.inc.stack,paste,collapse=";")

# ORF1
orf1.exc.split<-sapply(orf1.exc,function(x){
  a<-unlist(strsplit(x,"\\*"))
  nchar(a[1])>=7
  b<-a[-1]
  b<-b[grep("M",b)]
  pos<-sapply(b,function(x){unlist(gregexpr("M",b))[1]})
  c<-substring(b,pos,nchar(b))
  return(c(a[1],c)[nchar(c(a[1],c))>=7])
})
orf1.exc.stack<-stack(orf1.exc.split)
length(names(orf1.exc.split)[which(!names(orf1.exc.split) %in% orf1.exc.stack$ind)]) # 1,087 event do not generate sequences >=7aa; drop this levels
orf1.exc.stack<-stack(orf1.exc.split,drop=TRUE)
orf1.exc.stack$ID<-paste0(orf1.exc.stack$ind,"_ORF1.",unlist(sapply(table(orf1.exc.stack$ind),function(x){1:x})))
sum(duplicated(orf1.exc.stack$values)) # 4,181 are duplicated sequences: aggregate
orf1.exc.stack<-aggregate(ID~values,orf1.exc.stack,paste,collapse=";")

# ORF2
orf2.inc.split<-sapply(orf2.inc,function(x){
  a<-unlist(strsplit(x,"\\*"))
  nchar(a[1])>=7
  b<-a[-1]
  b<-b[grep("M",b)]
  pos<-sapply(b,function(x){unlist(gregexpr("M",b))[1]})
  c<-substring(b,pos,nchar(b))
  return(c(a[1],c)[nchar(c(a[1],c))>=7])
})
orf2.inc.stack<-stack(orf2.inc.split)
length(names(orf2.inc.split)[which(!names(orf2.inc.split) %in% orf2.inc.stack$ind)]) # 1,199 events do not generate sequences >=7aa; drop this levels
orf2.inc.stack<-stack(orf2.inc.split,drop=TRUE)
orf2.inc.stack$ID<-paste0(orf2.inc.stack$ind,"_ORF2.",unlist(sapply(table(orf2.inc.stack$ind),function(x){1:x})))
sum(duplicated(orf2.inc.stack$values)) # 10,498 are duplicated sequences: aggregate
orf2.inc.stack<-aggregate(ID~values,orf2.inc.stack,paste,collapse=";")

# ORF2
orf2.exc.split<-sapply(orf2.exc,function(x){
  a<-unlist(strsplit(x,"\\*"))
  nchar(a[1])>=7
  b<-a[-1]
  b<-b[grep("M",b)]
  pos<-sapply(b,function(x){unlist(gregexpr("M",b))[1]})
  c<-substring(b,pos,nchar(b))
  return(c(a[1],c)[nchar(c(a[1],c))>=7])
})
orf2.exc.stack<-stack(orf2.exc.split)
length(names(orf2.exc.split)[which(!names(orf2.exc.split) %in% orf2.exc.stack$ind)]) # 1,692 events do not generate sequences >=7aa; drop this levels
orf2.exc.stack<-stack(orf2.exc.split,drop=TRUE)
orf2.exc.stack$ID<-paste0(orf2.exc.stack$ind,"_ORF2.",unlist(sapply(table(orf2.exc.stack$ind),function(x){1:x})))
sum(duplicated(orf2.exc.stack$values)) # 4,564 are duplicated sequences: aggregate
orf2.exc.stack<-aggregate(ID~values,orf2.exc.stack,paste,collapse=";")


# ORF3
orf3.inc.split<-sapply(orf3.inc,function(x){
  a<-unlist(strsplit(x,"\\*"))
  nchar(a[1])>=7
  b<-a[-1]
  b<-b[grep("M",b)]
  pos<-sapply(b,function(x){unlist(gregexpr("M",b))[1]})
  c<-substring(b,pos,nchar(b))
  return(c(a[1],c)[nchar(c(a[1],c))>=7])
})
orf3.inc.stack<-stack(orf3.inc.split)
length(names(orf3.inc.split)[which(!names(orf3.inc.split) %in% orf3.inc.stack$ind)]) # 1,219 events do not generate sequences >=7aa; drop this levels
orf3.inc.stack<-stack(orf3.inc.split,drop=TRUE)
orf3.inc.stack$ID<-paste0(orf3.inc.stack$ind,"_ORF3.",unlist(sapply(table(orf3.inc.stack$ind),function(x){1:x})))
sum(duplicated(orf3.inc.stack$values)) # 9,918 are duplicated sequences: aggregate
orf3.inc.stack<-aggregate(ID~values,orf3.inc.stack,paste,collapse=";")

# ORF3
orf3.exc.split<-sapply(orf3.exc,function(x){
  a<-unlist(strsplit(x,"\\*"))
  nchar(a[1])>=7
  b<-a[-1]
  b<-b[grep("M",b)]
  pos<-sapply(b,function(x){unlist(gregexpr("M",b))[1]})
  c<-substring(b,pos,nchar(b))
  return(c(a[1],c)[nchar(c(a[1],c))>=7])
})
orf3.exc.stack<-stack(orf3.exc.split)
length(names(orf3.exc.split)[which(!names(orf3.exc.split) %in% orf3.exc.stack$ind)]) # 1,769 events do not generate sequences >=7aa; drop this levels
orf3.exc.stack<-stack(orf3.exc.split,drop=TRUE)
orf3.exc.stack$ID<-paste0(orf3.exc.stack$ind,"_ORF3.",unlist(sapply(table(orf3.exc.stack$ind),function(x){1:x})))
sum(duplicated(orf3.exc.stack$values)) # 4,086 are duplicated sequences: aggregate
orf3.exc.stack<-aggregate(ID~values,orf3.exc.stack,paste,collapse=";")

A5SS.fa<-rbind(orf1.inc.stack,orf1.exc.stack,orf2.inc.stack,orf2.exc.stack,orf3.inc.stack,orf3.exc.stack)
A5SS.fa<-aggregate(ID~values,A5SS.fa,paste,collapse=";") # 105,536
A5SS.fa<-data.frame(ID=A5SS.fa$ID,sequence=A5SS.fa$values)

##### Write .tsv with unique sequences #####
write.table(x=A5SS.fa,file="Fragpipe/inputs/IRIS_TCGA_A5SS.tsv",row.names = FALSE,quote=FALSE,col.names=TRUE,sep="\t")

#################################### A3SS ########################################
A3SS.read <- lapply(A3SS.files, fread)
# Collapse all files
A3SS.events<-lapply(A3SS.read,function(x){x[,1:10]})
A3SS.events<-rbindlist(A3SS.events,idcol="id")
A3SS.events$annotation<-with(A3SS.events,paste(chr,strand,longExonStart,longExonEnd,shortES,shortEE,flankingES,flankingEE,sep = "|"))
sum(duplicated(A3SS.events$annotation)) # 391,615
A3SS<-A3SS.events[!duplicated(A3SS.events$annotation),] # 24,289

# Translate (in this case we have complete coordinates)
range(abs(A3SS$longExonStart-A3SS$longExonEnd))
range(abs(A3SS$shortES-A3SS$shortEE))
range(abs(A3SS$flankingES-A3SS$flankingEE)) # there are exons as short as 2bp

long.gr<-GRanges(seqnames=A3SS$chr,
                 ranges=IRanges(start=A3SS$longExonStart+1,end=A3SS$longExonEnd),
                 strand=A3SS$strand)
short.gr<-GRanges(seqnames=A3SS$chr,
                  ranges=IRanges(start=A3SS$shortES+1,end=A3SS$shortEE),
                  strand=A3SS$strand)
flanking.gr<-GRanges(seqnames=A3SS$chr,
                     ranges=IRanges(start=A3SS$flankingES+1,end=A3SS$flankingEE),
                     strand=A3SS$strand)

long.rna<-getSeq(BSgenome.Hsapiens.UCSC.hg19,long.gr)
short.rna<-getSeq(BSgenome.Hsapiens.UCSC.hg19,short.gr)
flanking.rna<-getSeq(BSgenome.Hsapiens.UCSC.hg19,flanking.gr)

rna.inc<-paste0(flanking.rna,long.rna)
rna.exc<-paste0(flanking.rna,short.rna)

names(rna.inc)<-paste0("IRIS_",A3SS$annotation,"_A3SSinc")
names(rna.exc)<-paste0("IRIS_",A3SS$annotation,"_A3SSexc")


length(rna.inc)+length(rna.exc) # 48,578

inc.seq<-lapply(1:3,function(i){substring(rna.inc,i,nchar(rna.inc))})
exc.seq<-lapply(1:3,function(i){substring(rna.exc,i,nchar(rna.exc))})

orf1.inc<-as.character(Biostrings::translate(DNAStringSet(inc.seq[[1]]),
                                             if.fuzzy.codon = "solve",
                                             no.init.codon = TRUE))
orf2.inc<-as.character(Biostrings::translate(DNAStringSet(inc.seq[[2]]),
                                             if.fuzzy.codon = "solve",
                                             no.init.codon = TRUE))
orf3.inc<-as.character(Biostrings::translate(DNAStringSet(inc.seq[[3]]),
                                             if.fuzzy.codon = "solve",
                                             no.init.codon = TRUE))

orf1.exc<-as.character(Biostrings::translate(DNAStringSet(exc.seq[[1]]),
                                             if.fuzzy.codon = "solve",
                                             no.init.codon = TRUE))
orf2.exc<-as.character(Biostrings::translate(DNAStringSet(exc.seq[[2]]),
                                             if.fuzzy.codon = "solve",
                                             no.init.codon = TRUE))
orf3.exc<-as.character(Biostrings::translate(DNAStringSet(exc.seq[[3]]),
                                             if.fuzzy.codon = "solve",
                                             no.init.codon = TRUE))

##### Keep ORFs that are >=7aa #####
# ORF1
orf1.inc.split<-sapply(orf1.inc,function(x){
  a<-unlist(strsplit(x,"\\*"))
  nchar(a[1])>=7
  b<-a[-1]
  b<-b[grep("M",b)]
  pos<-sapply(b,function(x){unlist(gregexpr("M",b))[1]})
  c<-substring(b,pos,nchar(b))
  return(c(a[1],c)[nchar(c(a[1],c))>=7])
})
orf1.inc.stack<-stack(orf1.inc.split)
length(names(orf1.inc.split)[which(!names(orf1.inc.split) %in% orf1.inc.stack$ind)]) # 1,327 events do not generate sequences >=7aa; drop this levels
orf1.inc.stack<-stack(orf1.inc.split,drop=TRUE)
orf1.inc.stack$ID<-paste0(orf1.inc.stack$ind,"_ORF1.",unlist(sapply(table(orf1.inc.stack$ind),function(x){1:x})))
sum(duplicated(orf1.inc.stack$values)) # 7,271 are duplicated sequences: aggregate
orf1.inc.stack<-aggregate(ID~values,orf1.inc.stack,paste,collapse=";")

# ORF1
orf1.exc.split<-sapply(orf1.exc,function(x){
  a<-unlist(strsplit(x,"\\*"))
  nchar(a[1])>=7
  b<-a[-1]
  b<-b[grep("M",b)]
  pos<-sapply(b,function(x){unlist(gregexpr("M",b))[1]})
  c<-substring(b,pos,nchar(b))
  return(c(a[1],c)[nchar(c(a[1],c))>=7])
})
orf1.exc.stack<-stack(orf1.exc.split)
length(names(orf1.exc.split)[which(!names(orf1.exc.split) %in% orf1.exc.stack$ind)]) # 1,740 event do not generate sequences >=7aa; drop this levels
orf1.exc.stack<-stack(orf1.exc.split,drop=TRUE)
orf1.exc.stack$ID<-paste0(orf1.exc.stack$ind,"_ORF1.",unlist(sapply(table(orf1.exc.stack$ind),function(x){1:x})))
sum(duplicated(orf1.exc.stack$values)) # 4,459 are duplicated sequences: aggregate
orf1.exc.stack<-aggregate(ID~values,orf1.exc.stack,paste,collapse=";")

# ORF2
orf2.inc.split<-sapply(orf2.inc,function(x){
  a<-unlist(strsplit(x,"\\*"))
  nchar(a[1])>=7
  b<-a[-1]
  b<-b[grep("M",b)]
  pos<-sapply(b,function(x){unlist(gregexpr("M",b))[1]})
  c<-substring(b,pos,nchar(b))
  return(c(a[1],c)[nchar(c(a[1],c))>=7])
})
orf2.inc.stack<-stack(orf2.inc.split)
length(names(orf2.inc.split)[which(!names(orf2.inc.split) %in% orf2.inc.stack$ind)]) # 1,838 events do not generate sequences >=7aa; drop this levels
orf2.inc.stack<-stack(orf2.inc.split,drop=TRUE)
orf2.inc.stack$ID<-paste0(orf2.inc.stack$ind,"_ORF2.",unlist(sapply(table(orf2.inc.stack$ind),function(x){1:x})))
sum(duplicated(orf2.inc.stack$values)) # 6,800 are duplicated sequences: aggregate
orf2.inc.stack<-aggregate(ID~values,orf2.inc.stack,paste,collapse=";")

# ORF2
orf2.exc.split<-sapply(orf2.exc,function(x){
  a<-unlist(strsplit(x,"\\*"))
  nchar(a[1])>=7
  b<-a[-1]
  b<-b[grep("M",b)]
  pos<-sapply(b,function(x){unlist(gregexpr("M",b))[1]})
  c<-substring(b,pos,nchar(b))
  return(c(a[1],c)[nchar(c(a[1],c))>=7])
})
orf2.exc.stack<-stack(orf2.exc.split)
length(names(orf2.exc.split)[which(!names(orf2.exc.split) %in% orf2.exc.stack$ind)]) # 2,457 events do not generate sequences >=7aa; drop this levels
orf2.exc.stack<-stack(orf2.exc.split,drop=TRUE)
orf2.exc.stack$ID<-paste0(orf2.exc.stack$ind,"_ORF2.",unlist(sapply(table(orf2.exc.stack$ind),function(x){1:x})))
sum(duplicated(orf2.exc.stack$values)) # 4,496 are duplicated sequences: aggregate
orf2.exc.stack<-aggregate(ID~values,orf2.exc.stack,paste,collapse=";")


# ORF3
orf3.inc.split<-sapply(orf3.inc,function(x){
  a<-unlist(strsplit(x,"\\*"))
  nchar(a[1])>=7
  b<-a[-1]
  b<-b[grep("M",b)]
  pos<-sapply(b,function(x){unlist(gregexpr("M",b))[1]})
  c<-substring(b,pos,nchar(b))
  return(c(a[1],c)[nchar(c(a[1],c))>=7])
})
orf3.inc.stack<-stack(orf3.inc.split)
length(names(orf3.inc.split)[which(!names(orf3.inc.split) %in% orf3.inc.stack$ind)]) # 1,895 events do not generate sequences >=7aa; drop this levels
orf3.inc.stack<-stack(orf3.inc.split,drop=TRUE)
orf3.inc.stack$ID<-paste0(orf3.inc.stack$ind,"_ORF3.",unlist(sapply(table(orf3.inc.stack$ind),function(x){1:x})))
sum(duplicated(orf3.inc.stack$values)) # 7,096 are duplicated sequences: aggregate
orf3.inc.stack<-aggregate(ID~values,orf3.inc.stack,paste,collapse=";")

# ORF3
orf3.exc.split<-sapply(orf3.exc,function(x){
  a<-unlist(strsplit(x,"\\*"))
  nchar(a[1])>=7
  b<-a[-1]
  b<-b[grep("M",b)]
  pos<-sapply(b,function(x){unlist(gregexpr("M",b))[1]})
  c<-substring(b,pos,nchar(b))
  return(c(a[1],c)[nchar(c(a[1],c))>=7])
})
orf3.exc.stack<-stack(orf3.exc.split)
length(names(orf3.exc.split)[which(!names(orf3.exc.split) %in% orf3.exc.stack$ind)]) # 2,454 events do not generate sequences >=7aa; drop this levels
orf3.exc.stack<-stack(orf3.exc.split,drop=TRUE)
orf3.exc.stack$ID<-paste0(orf3.exc.stack$ind,"_ORF3.",unlist(sapply(table(orf3.exc.stack$ind),function(x){1:x})))
sum(duplicated(orf3.exc.stack$values)) # 4,429 are duplicated sequences: aggregate
orf3.exc.stack<-aggregate(ID~values,orf3.exc.stack,paste,collapse=";")

A3SS.fa<-rbind(orf1.inc.stack,orf1.exc.stack,orf2.inc.stack,orf2.exc.stack,orf3.inc.stack,orf3.exc.stack)
A3SS.fa<-aggregate(ID~values,A3SS.fa,paste,collapse=";") # 148,487
A3SS.fa<-data.frame(ID=A3SS.fa$ID,sequence=A3SS.fa$values)

##### Write .tsv with unique sequences #####
write.table(x=A3SS.fa,file="Fragpipe/inputs/IRIS_TCGA_A3SS.tsv",row.names = FALSE,quote=FALSE,col.names=TRUE,sep="\t")


#################################### RI ########################################
RI.read <- lapply(RI.files, fread)
# CollapRI all files
RI.events<-lapply(RI.read,function(x){x[,1:10]})
RI.events<-rbindlist(RI.events,idcol="id")
RI.events$annotation<-with(RI.events,paste(chr,strand,riExonStart,riExonEnd,upstreamES,upstreamEE,downstreamES,downstreamEE, sep = "|"))
sum(duplicated(RI.events$annotation)) # 186,761
RI<-RI.events[!duplicated(RI.events$annotation),] # 6,793

# Translate
range(abs(RI$exonStart-RI$exonEnd)) # there are exons as short as 1bp

UE.gr<-GRanges(seqnames=RI$chr,
               ranges=IRanges(start=RI$upstreamES+1,end=RI$upstreamEE),
               strand=RI$strand)
DE.gr<-GRanges(seqnames=RI$chr,
               ranges=IRanges(start=RI$downstreamES+1,end=RI$downstreamEE),
               strand=RI$strand)
RI.gr<-GRanges(seqnames=RI$chr,
               ranges=IRanges(start=RI$riExonStart+1,end=RI$riExonEnd),
               strand=RI$strand)

UE.rna<-getSeq(BSgenome.Hsapiens.UCSC.hg19,UE.gr)
DE.rna<-getSeq(BSgenome.Hsapiens.UCSC.hg19,DE.gr)
RI.rna<-getSeq(BSgenome.Hsapiens.UCSC.hg19,RI.gr)

rna.inc<-RI.rna
rna.exc<-ifelse(RI$strand=="+",
                paste0(UE.rna,DE.rna),
                paste0(DE.rna,UE.rna))

names(rna.inc)<-paste0("IRIS_",RI$annotation,"_RIinc")
names(rna.exc)<-paste0("IRIS_",RI$annotation,"_RIexc")


length(rna.inc)+length(rna.exc) # 13,586

inc.seq<-lapply(1:3,function(i){substring(rna.inc,i,nchar(rna.inc))})
exc.seq<-lapply(1:3,function(i){substring(rna.exc,i,nchar(rna.exc))})

orf1.inc<-as.character(Biostrings::translate(DNAStringSet(inc.seq[[1]]),
                                             if.fuzzy.codon = "solve",
                                             no.init.codon = TRUE))
orf2.inc<-as.character(Biostrings::translate(DNAStringSet(inc.seq[[2]]),
                                             if.fuzzy.codon = "solve",
                                             no.init.codon = TRUE))
orf3.inc<-as.character(Biostrings::translate(DNAStringSet(inc.seq[[3]]),
                                             if.fuzzy.codon = "solve",
                                             no.init.codon = TRUE))

orf1.exc<-as.character(Biostrings::translate(DNAStringSet(exc.seq[[1]]),
                                             if.fuzzy.codon = "solve",
                                             no.init.codon = TRUE))
orf2.exc<-as.character(Biostrings::translate(DNAStringSet(exc.seq[[2]]),
                                             if.fuzzy.codon = "solve",
                                             no.init.codon = TRUE))
orf3.exc<-as.character(Biostrings::translate(DNAStringSet(exc.seq[[3]]),
                                             if.fuzzy.codon = "solve",
                                             no.init.codon = TRUE))

##### Keep ORFs that are >=7aa #####
# ORF1
orf1.inc.split<-sapply(orf1.inc,function(x){
  a<-unlist(strsplit(x,"\\*"))
  nchar(a[1])>=7
  b<-a[-1]
  b<-b[grep("M",b)]
  pos<-sapply(b,function(x){unlist(gregexpr("M",b))[1]})
  c<-substring(b,pos,nchar(b))
  return(c(a[1],c)[nchar(c(a[1],c))>=7])
})
orf1.inc.stack<-stack(orf1.inc.split)
length(names(orf1.inc.split)[which(!names(orf1.inc.split) %in% orf1.inc.stack$ind)]) # 217 event do not generate RIquences >=7aa; drop this levels
orf1.inc.stack<-stack(orf1.inc.split,drop=TRUE)
sum(duplicated(orf1.inc.stack$ind)) # 12,695
orf1.inc.stack$ID<-paste0(orf1.inc.stack$ind,"_ORF1.",unlist(sapply(table(orf1.inc.stack$ind),function(x){1:x})))
sum(duplicated(orf1.inc.stack$values)) # 4,390 are duplicated sequences: aggregate
orf1.inc.stack<-aggregate(ID~values,orf1.inc.stack,paste,collapse=";")

# ORF1
orf1.exc.split<-sapply(orf1.exc,function(x){
  a<-unlist(strsplit(x,"\\*"))
  nchar(a[1])>=7
  b<-a[-1]
  b<-b[grep("M",b)]
  pos<-sapply(b,function(x){unlist(gregexpr("M",b))[1]})
  c<-substring(b,pos,nchar(b))
  return(c(a[1],c)[nchar(c(a[1],c))>=7])
})
orf1.exc.stack<-stack(orf1.exc.split)
length(names(orf1.exc.split)[which(!names(orf1.exc.split) %in% orf1.exc.stack$ind)]) # 463 events do not generate RIquences >=7aa; drop this levels
orf1.exc.stack<-stack(orf1.exc.split,drop=TRUE)
sum(duplicated(orf1.exc.stack$ind))
orf1.exc.stack$ID<-paste0(orf1.exc.stack$ind,"_ORF1.",unlist(sapply(table(orf1.exc.stack$ind),function(x){1:x})))
sum(duplicated(orf1.exc.stack$values)) # 723 are duplicated sequences: aggregate
orf1.exc.stack<-aggregate(ID~values,orf1.exc.stack,paste,collapse=";")

# ORF2
orf2.inc.split<-sapply(orf2.inc,function(x){
  a<-unlist(strsplit(x,"\\*"))
  nchar(a[1])>=7
  b<-a[-1]
  b<-b[grep("M",b)]
  pos<-sapply(b,function(x){unlist(gregexpr("M",b))[1]})
  c<-substring(b,pos,nchar(b))
  return(c(a[1],c)[nchar(c(a[1],c))>=7])
})
orf2.inc.stack<-stack(orf2.inc.split)
length(names(orf2.inc.split)[which(!names(orf2.inc.split) %in% orf2.inc.stack$ind)]) # 326 events do not generate RIquences >=7aa; drop this levels
orf2.inc.stack<-stack(orf2.inc.split,drop=TRUE)
sum(duplicated(orf2.inc.stack$ind))
orf2.inc.stack$ID<-paste0(orf2.inc.stack$ind,"_ORF2.",unlist(sapply(table(orf2.inc.stack$ind),function(x){1:x})))
sum(duplicated(orf2.inc.stack$values)) # 4,154 are duplicated sequences: aggregate
orf2.inc.stack<-aggregate(ID~values,orf2.inc.stack,paste,collapse=";")

# ORF2
orf2.exc.split<-sapply(orf2.exc,function(x){
  a<-unlist(strsplit(x,"\\*"))
  nchar(a[1])>=7
  b<-a[-1]
  b<-b[grep("M",b)]
  pos<-sapply(b,function(x){unlist(gregexpr("M",b))[1]})
  c<-substring(b,pos,nchar(b))
  return(c(a[1],c)[nchar(c(a[1],c))>=7])
})
orf2.exc.stack<-stack(orf2.exc.split)
length(names(orf2.exc.split)[which(!names(orf2.exc.split) %in% orf2.exc.stack$ind)]) # 672 events do not generate RIquences >=7aa; drop this levels
orf2.exc.stack<-stack(orf2.exc.split,drop=TRUE)
sum(duplicated(orf2.exc.stack$ind))
orf2.exc.stack$ID<-paste0(orf2.exc.stack$ind,"_ORF2.",unlist(sapply(table(orf2.exc.stack$ind),function(x){1:x})))
sum(duplicated(orf2.exc.stack$values)) # 805 are duplicated sequences: aggregate
orf2.exc.stack<-aggregate(ID~values,orf2.exc.stack,paste,collapse=";")


# ORF3
orf3.inc.split<-sapply(orf3.inc,function(x){
  a<-unlist(strsplit(x,"\\*"))
  nchar(a[1])>=7
  b<-a[-1]
  b<-b[grep("M",b)]
  pos<-sapply(b,function(x){unlist(gregexpr("M",b))[1]})
  c<-substring(b,pos,nchar(b))
  return(c(a[1],c)[nchar(c(a[1],c))>=7])
})
orf3.inc.stack<-stack(orf3.inc.split)
length(names(orf3.inc.split)[which(!names(orf3.inc.split) %in% orf3.inc.stack$ind)]) # 306 events do not generate RIquences >=7aa; drop this levels
orf3.inc.stack<-stack(orf3.inc.split,drop=TRUE)
sum(duplicated(orf3.inc.stack$ind)) 
orf3.inc.stack$ID<-paste0(orf3.inc.stack$ind,"_ORF3.",unlist(sapply(table(orf3.inc.stack$ind),function(x){1:x})))
sum(duplicated(orf3.inc.stack$values)) # 4,357 are duplicated sequences: aggregate
orf3.inc.stack<-aggregate(ID~values,orf3.inc.stack,paste,collapse=";")

# ORF3
orf3.exc.split<-sapply(orf3.exc,function(x){
  a<-unlist(strsplit(x,"\\*"))
  nchar(a[1])>=7
  b<-a[-1]
  b<-b[grep("M",b)]
  pos<-sapply(b,function(x){unlist(gregexpr("M",b))[1]})
  c<-substring(b,pos,nchar(b))
  return(c(a[1],c)[nchar(c(a[1],c))>=7])
})
orf3.exc.stack<-stack(orf3.exc.split)
length(names(orf3.exc.split)[which(!names(orf3.exc.split) %in% orf3.exc.stack$ind)]) # 734 events do not generate RIquences >=7aa; drop this levels
orf3.exc.stack<-stack(orf3.exc.split,drop=TRUE)
sum(duplicated(orf3.exc.stack$ind)) 
orf3.exc.stack$ID<-paste0(orf3.exc.stack$ind,"_ORF3.",unlist(sapply(table(orf3.exc.stack$ind),function(x){1:x})))
sum(duplicated(orf3.exc.stack$values)) # 743 are duplicated sequences: aggregate
orf3.exc.stack<-aggregate(ID~values,orf3.exc.stack,paste,collapse=";")

RI.fa<-rbind(orf1.inc.stack,orf1.exc.stack,orf2.inc.stack,orf2.exc.stack,orf3.inc.stack,orf3.exc.stack)
RI.fa<-aggregate(ID~values,RI.fa,paste,collapse=";") # 60,274
RI.fa<-data.frame(ID=RI.fa$ID,sequence=RI.fa$values)

##### Write .tsv with unique sequences #####
write.table(x=RI.fa,file="Fragpipe/inputs/IRIS_TCGA_RI.tsv",row.names = FALSE,quote=FALSE,col.names=TRUE,sep="\t")





######################### Merge TSVs ###########################################
list.files("Fragpipe/inputs")
file.list<-list.files("Fragpipe/inputs",
                      recursive = TRUE,
                      full.names = TRUE,
                      pattern = ".tsv$")
tsv.read <- lapply(file.list, fread)
tsv.rbind <- rbindlist(tsv.read) 
# tsv.reduced<-aggregate(ID~sequence,tsv.rbind,paste,collapse=";"): takes a lot of time

e.g. !duplicated(c(1,1,1,2,3)) & !duplicated(c(1,1,1,2,3),fromLast=TRUE)
tsv.unique<-tsv.rbind[!duplicated(tsv.rbind$sequence) & !duplicated(tsv.rbind$sequence,fromLast=TRUE),]
tsv.duplicated<-tsv.rbind[duplicated(tsv.rbind$sequence) | duplicated(tsv.rbind$sequence,fromLast=TRUE),]
tsv.reduced<-aggregate(ID~sequence,tsv.duplicated,paste,collapse=";")

nrow(tsv.unique) + nrow(tsv.reduced) == length(unique(tsv.rbind$sequence)) # TRUE, ok

tsv.final<-rbind(tsv.unique,tsv.reduced)
write.table(x=tsv.final,file="Fragpipe/splicing.tsv",row.names = FALSE,quote=FALSE,col.names=TRUE,sep="\t")
