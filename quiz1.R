#Quiz-week1
library(AnnotationHub)
ahub = AnnotationHub()
ahub = subset(ahub, species == "Homo sapiens")
qhs = query(ahub,"CpG Islands")
cpg1 = qhs[[1]] #AH5086


autosome = c(paste("chr",1:22,sep=""))
cpg1_splitted = split(cpg1, seqnames(cpg1))
cpg1_autosome = cpg1_splitted[autosome]


q1 = length(unlist(cpg1_autosome))


q2 = length(cpg1_splitted$chr4)


qhs = query(ahub,c("H3K4me3", "H1 Cells", "EpigenomeRoadMap"))
H3K4me3 = qhs[[2]] #AH29884 
H3K4me3_autosome = subset(H3K4me3, seqnames %in% autosome)
q3 = sum(width(H3K4me3_autosome))


qhs = query(ahub,c("H3K27me3", "H1 Cells", "EpigenomeRoadMap"))
H3K27me3 = qhs[[2]] #AH29892 
H3K27me3_autosome = subset(H3K27me3, seqnames %in% autosome)
q4 = mean(H3K27me3_autosome$signalValue)


bival_autosome = intersect(H3K4me3_autosome,H3K27me3_autosome, ignore.strand=TRUE)
q5 = sum(width(bival_autosome))


bival_ov_cpg1 = findOverlaps(bival_autosome, unlist(cpg1_autosome), ignore.strand=TRUE)
q6 = length(unique(queryHits(bival_ov_cpg1))) / length(bival_autosome)


cpg1_inters_bival = intersect(unlist(cpg1_autosome), bival_autosome , ignore.strand=TRUE)
q7 = sum(width(unique(cpg1_inters_bival))) / sum(width(unique(unlist(cpg1_autosome))))


cpg1a_10kb_centered = resize(cpg1_autosome, width = 20000+ width(cpg1_autosome),fix="center")
cpg1a10kbc_inters_bival = intersect(unlist(cpg1a_10kb_centered), bival_autosome , ignore.strand=TRUE)
q8=sum(width(unique(cpg1a10kbc_inters_bival)))


genome = ahub[["AH5018"]]
genome_autosome = keepSeqlevels(genome, autosome, pruning.mode = "coarse")                                               
q9 = sum(width(reduce(unlist(cpg1_autosome))))/sum(seqlengths(genome_autosome))



inOut = matrix(0,ncol=2,nrow = 2)
rownames(inOut) = c("in","out")
colnames(inOut) = c("in","out")

inOut[1,1] = sum(width(intersect(unlist(cpg1_autosome), bival_autosome , ignore.strand=TRUE)))
inOut[1,2] = sum(width(setdiff(unlist(cpg1_autosome), bival_autosome , ignore.strand=TRUE)))
inOut[2,1] = sum(width(setdiff(bival_autosome ,unlist(cpg1_autosome), ignore.strand=TRUE)))
genome_length = sum(seqlengths(genome_autosome))
inOut[2,2] = genome_length - sum(inOut)
OR = inOut[1,1] * inOut[2,2] / (inOut[2,1]*inOut[1,2])
q10 = OR