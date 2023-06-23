# Libraries
library(GenomicAlignments)

# Data
linear<-read.delim("C:/Users/bea_f/OneDrive/Escritorio/HematoLaFe/Neoantigens-RNA-splicing/RJunBase/detail_LS_annotation.txt")
ASannotations<-read.delim("C:/Users/bea_f/OneDrive/Escritorio/HematoLaFe/Neoantigens-RNA-splicing/RJunBase/Alternative_Splice_Junctions.tar.gz",sep=",",header=FALSE)
laml<-read.csv2("C:/Users/bea_f/OneDrive/Escritorio/HematoLaFe/Neoantigens-RNA-splicing/RJunBase/LAML_expr_specific.csv",sep=",")
gtex<-read.delim("C:/Users/bea_f/OneDrive/Escritorio/HematoLaFe/Neoantigens-RNA-splicing/GTEx_Analysis_2017-06-05_v8_STARv2.5.3a_junctions.gct.gz",sep="\t",skip=2)

cancers<-unique(unlist(strsplit(unique(linear$Tumor.specific.types),";")))
cancers.specific.count<-sapply(cancers,function(x){sum(linear$Tumor.specific.types==x)})
sum(cancers.specific.count[names(cancers.specific.count)!="non tumor-specific"]) # 25818
sum(table(linear$Tumor.specific.types[grep(";",linear$Tumor.specific.types)])) # 7019, OK

# Types of tumor-specific junctions
table(linear$Alternative.splicing.types[linear$Tumor.specific.types %in% cancers[cancers!="non tumor-specific"]]) # solo hay 7 que sean solo MXE
# 21819 junctions are not related to an AS event
25818-21819 # 3999 junctions are related to an AS event


# Hay junctions que no son espewcíficas de tumor pero que no se encuentran en tejido normal
# Ej. UN_SLCO1A2_LS0012
# Hay junctions que tienen una frecuencia 0(1), o sea que se encuentran en 1 muestra pero la frecuencia relativa al redondear es 0
pattern <- "\\(\\s*(.*?)\\s*\\)"
nt.freq<- regmatches(linear$NT.frequency.sample.number., regexec(pattern, linear$NT.frequency.sample.number.))
nt.freq.indiv<-unlist(lapply(nt.freq,function(x){x[[2]]}))

normal.freq<- regmatches(linear$Normal.frequency.sample.number., regexec(pattern, linear$Normal.frequency.sample.number.))
normal.freq.indiv<-unlist(lapply(normal.freq,function(x){x[[2]]}))

tumor.freq<- regmatches(linear$Tumor.frequency.sample.number., regexec(pattern, linear$Tumor.frequency.sample.number.))
tumor.freq.indiv<-unlist(lapply(tumor.freq,function(x){x[[2]]}))

# 133,646 junctions not found in normal tissue
sum(nt.freq.indiv==0 & normal.freq.indiv==0 & tumor.freq.indiv>1)
View(linear[nt.freq.indiv==0 & normal.freq.indiv==0 & tumor.freq.indiv>1,])

sum(linear$NT.median==0 & linear$Normal.median==0)
sum(linear$NT.median==0 & linear$Normal.median==0 & linear$Tumor.median>0)
sum(grepl("\\(0\\)",linear$NT.frequency.sample.number.) &
      grepl("\\(0\\)",linear$Normal.frequency.sample.number.) &
      grepl("\\(0\\)",linear$Normal.frequency.sample.number.))
sum(linear$Tumor.median[linear$Tumor.specific.types=="non tumor-specific"])

############## LAML: selection of 2.185 specific junctions #####################
laml.junctions<-laml$sample[laml$sample %in% linear$Junction.location[linear$Tumor.specific.types=="LAML"]]


############### CHOL: selection of 466 specific junctions ######################
chol<-read.csv2("C:/Users/bea_f/OneDrive/Escritorio/HematoLaFe/Neoantigens-RNA-splicing/CHOL_expr_specific.csv",sep=",")
chol.junctions<-chol$sample[chol$sample %in% linear$Junction.location[linear$Tumor.specific.types=="CHOL"]]


############################## BED file ########################################
chr.order<-c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY","chrM")

bed.file <- file("/RJunBase/RJunbase.txt", "wb")
df<-data.frame(chrom=sub(":.*","",linear$Junction.location),
               start=sub(".*: *(.*?) *\\|.*", "\\1", linear$Junction.location),
               end=sub(".*\\| *(.*?) *:.*", "\\1", linear$Junction.location),
               strand=sub(".*:","",linear$Junction.location))
df<-df[order(as.numeric(df$start)),]
df<-df[order(match(df$chrom,chr.order)),]
write.table(x=df,file=bed.file,row.names = FALSE,quote=FALSE,col.names=TRUE)
close(bed.file)


bed.annotated.file <- file("C:/Users/bea_f/OneDrive/Escritorio/HematoLaFe/Neoantigens-RNA-splicing/RJunBase/RJunbase_annotated.txt", "wb")
df<-data.frame(chrom=sub(":.*","",linear$Junction.location),
               start=sub(".*: *(.*?) *\\|.*", "\\1", linear$Junction.location),
               end=sub(".*\\| *(.*?) *:.*", "\\1", linear$Junction.location),
               strand=sub(".*:","",linear$Junction.location))
df<-cbind(df,linear)
df<-df[order(as.numeric(df$start)),]
df<-df[order(match(df$chrom,chr.order)),]
write.table(x=df,file=bed.annotated.file,row.names = FALSE,quote=FALSE,col.names=TRUE,sep="\t")
close(bed.annotated.file)


########################### Overlaping regions #################################
non.can.bam<-readGAlignments("Peptides/non-canonical.sorted.bam",use.names = TRUE)
junctions<-read.delim("RJunBase/RJunbase_annotated.txt",sep="\t")
junctions$chrom[junctions$chrom=="chrM"]<-"chrMT"

# Junction frequency
pattern <- "\\(\\s*(.*?)\\s*\\)"
nt.freq<- regmatches(junctions$NT.frequency.sample.number., regexec(pattern, junctions$NT.frequency.sample.number.))
nt.freq.indiv<-unlist(lapply(nt.freq,function(x){x[[2]]}))

normal.freq<- regmatches(junctions$Normal.frequency.sample.number., regexec(pattern, junctions$Normal.frequency.sample.number.))
normal.freq.indiv<-unlist(lapply(normal.freq,function(x){x[[2]]}))

tumor.freq<- regmatches(junctions$Tumor.frequency.sample.number., regexec(pattern, junctions$Tumor.frequency.sample.number.))
tumor.freq.indiv<-unlist(lapply(tumor.freq,function(x){x[[2]]}))


head(non.can.bam)
max(width(non.can.bam)) # 36

gr<-as.data.frame(non.can.bam)
gr<-makeGRangesFromDataFrame(gr)


### Search only in UE or DE
gr.UE<-GRanges(seqnames=junctions$chrom,
               ranges = IRanges(start=junctions$start-36,end=junctions$start),
               strand=junctions$strand)
gr.DE<-GRanges(seqnames=junctions$chrom,
               ranges = IRanges(start=junctions$end,end=junctions$end+36),
               strand=junctions$strand)
result.UE<-findOverlaps(gr,gr.UE,type="within")
result.DE<-findOverlaps(gr,gr.DE,type="within")

gr.df<-as.data.frame(gr)
colnames(gr.df)<-paste0("peptide_",colnames(gr.df))
View(cbind(gr.df[queryHits(result.UE),],junctions[subjectHits(result.UE),]))
write.table(x=cbind(gr.df[queryHits(result.UE),],junctions[subjectHits(result.UE),]),
            file="peptidesUE.txt",quote=FALSE)

write.table(x=cbind(gr.df[queryHits(result.DE),],junctions[subjectHits(result.DE),]),
            file="peptidesDE.txt",quote=FALSE)

sum(queryHits(result.UE) %in% queryHits(result.DE)) # 166

queries<-c(queryHits(result.UE),queryHits(result.DE))
unique.queries<-queries[!duplicated(queries) & !duplicated(queries,fromLast=TRUE)] # 2186
UE.DE.queries<-unique(queryHits(result.UE)[queryHits(result.UE) %in% queryHits(result.DE)]) # 103

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


# Search for more tumor-specific peptides
sum(nt.freq.indiv[UE.unique$subjectHits]==0 & normal.freq.indiv[UE.unique$subjectHits]==0) # 71
View(cbind(UE.unique,junctions[UE.unique$subjectHits,])[which(nt.freq.indiv[UE.unique$subjectHits]==0 & normal.freq.indiv[UE.unique$subjectHits]==0),]) # 71
View(cbind(DE.unique,junctions[DE.unique$subjectHits,])[which(nt.freq.indiv[DE.unique$subjectHits]==0 & normal.freq.indiv[DE.unique$subjectHits]==0),]) # 71

# 4 peptides with common hits in UE and DE from different junctions
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

UE.DE.unique[!is.na(UE.DE.unique)] # 4. ninguna junction clasificada como tumor-specific


UE.duplicated.unique<-vector()
for (i in 1:(length(unique(queryHits(UE.duplicated))))){
  index<-subjectHits(UE.duplicated)[queryHits(UE.duplicated)==unique(queryHits(UE.duplicated))[i]]
  if (all(c(normal.freq.indiv[index]==0 &
            nt.freq.indiv[index]==0))){
    UE.duplicated.unique[i]<-unique(queryHits(UE.duplicated))[i]
  }
}

UE.duplicated.unique[!is.na(UE.duplicated.unique)] # 1, ninguna junction clasificada como tumor-specific
View(cbind(as.data.frame(UE.duplicated[queryHits(UE.duplicated) %in% UE.duplicated.unique[!is.na(UE.duplicated.unique)]]),
           junctions[subjectHits(UE.duplicated[queryHits(UE.duplicated) %in% UE.duplicated.unique[!is.na(UE.duplicated.unique)]]),]))


DE.duplicated.unique<-vector()
for (i in 1:(length(unique(queryHits(DE.duplicated))))){
  index<-subjectHits(DE.duplicated)[queryHits(DE.duplicated)==unique(queryHits(DE.duplicated))[i]]
  if (all(c(normal.freq.indiv[index]==0 &
            nt.freq.indiv[index]==0))){
    DE.duplicated.unique[i]<-unique(queryHits(DE.duplicated))[i]
  }
}
DE.duplicated.unique[!is.na(DE.duplicated.unique)] # 52, ninguna junction clasificada como tumor-specific
View(cbind(as.data.frame(DE.duplicated[queryHits(DE.duplicated) %in% DE.duplicated.unique[!is.na(DE.duplicated.unique)]]),
           junctions[subjectHits(DE.duplicated[queryHits(DE.duplicated) %in% DE.duplicated.unique[!is.na(DE.duplicated.unique)]]),]))









# Overlaping hits: a hit may be related to more than one junction
unique(queryHits(result.UE)[queryHits(result.UE) %in% queryHits(result.DE)])
gr.UE[subjectHits(result.UE)[queryHits(result.UE)==749]]
gr.DE[subjectHits(result.DE)[queryHits(result.DE)==749]]


# Junction classification
# UE
# First remove related junctions
junctions.UE.NT<-junctions[rownames(junctions) %in% unique(subjectHits(result.UE)) &
  nt.freq.indiv==0 &
  normal.freq.indiv==0 &
  (duplicated(junctions$start) | duplicated(junctions$start,fromLast=TRUE)),]

View(junctions.UE.NT[junctions.UE.NT$start %in% names(table(junctions.UE.NT$start))[table(junctions.UE.NT$start)==1],]) # 82

junctions.UE.specific<-junctions[rownames(junctions) %in% unique(subjectHits(result.UE)) &
                                   junctions$Tumor.specific.types!="non tumor-specific",] # 22

View(junctions.UE.specific[junctions.UE.specific$start %in% names(table(junctions.UE.specific$start))[table(junctions.UE.specific$start)==1],])

# DE
junctions.DE<-junctions[unique(subjectHits(result.DE)),]
View(junctions.UE[junctions.DE$end %in% names(table(junctions.DE$end))[table(junctions.DE$end)==1],])


