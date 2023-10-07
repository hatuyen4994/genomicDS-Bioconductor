library(BSgenome.Hsapiens.UCSC.hg19)
chr22 = Hsapiens$chr22
no_N = length(chr22) - unlist(letterFrequency(chr22,"N"))
GC_content = letterFrequency(chr22,"GC") / no_N
q1 = GC_content


library(AnnotationHub)
ahub = AnnotationHub()
ahub = subset(ahub, species == "Homo sapiens")
qhs = query(ahub, c("H3K27me3", "EpigenomeRoadMap", "H1 Cells"))
H3K27me3 = qhs[["AH29892"]]
H3K27me3_chr22 = subset(H3K27me3, seqnames=="chr22")
H3K27me3_chr22_vi = Views(Hsapiens,H3K27me3_chr22)
q2=mean(letterFrequency(H3K27me3_chr22_vi,"GC", as.prob = TRUE))

GC_chr22_prob = letterFrequency(H3K27me3_chr22_vi, "GC", as.prob=TRUE)
sigV = mcols(H3K27me3_chr22_vi)$signalValue
cor_test = cor.test(GC_chr22_prob, sigV)
q3 = cor_test$estimate

 
fs = qhs[["AH32033"]]
gr_chr22 = GRanges("chr22", ranges=IRanges(start = 1, end = 51304566))
fs_rle = import(fs, which=gr_chr22, as = "Rle")
fs_rle_chr22 = fs_rle$chr22
fs_chr22_view = Views(fs_rle_chr22, ranges(H3K27me3_chr22))
q4=cor(sigV, mean(fs_chr22_view))

q5= sum(fs_rle_chr22 > 1)


qhs_E055 = query(ahub, c("H3K27me3", "EpigenomeRoadMap", "E055"))
H3K27me3_fs_E055 = qhs_E055[["AH32470"]]
fs_rle_E055 = import(H3K27me3_fs_E055, which=gr_chr22, as = "Rle")
fs_rle_chr22_E055 = fs_rle_E055$chr22


E003_lower_0.5 = slice(fs_rle_chr22, upper = 0.5)
E055_upper_2 = slice(fs_rle_chr22_E055, lower = 2)

E003_lower_0.5 = as(E003_lower_0.5, "IRanges")
E055_upper_2 = as(E055_upper_2, "IRanges")

inter_region = intersect(E003_lower_0.5,E055_upper_2)
q6= sum(width(inter_region))


qhs = query(ahub, "CpG Islands")
cpg = qhs[["AH5086"]]
cpg_chr22 = keepSeqlevels(cpg,"chr22",pruning.mode = "coarse")
cpg_chr22_vi = Views(Hsapiens$chr22,ranges(cpg_chr22))

cpgs_chr22_lengths = elementNROWS(cpg_chr22_vi)
obs_CG = dinucleotideFrequency(cpg_chr22_vi)[,7]/cpgs_chr22_lengths
exp_CG = (letterFrequency(cpg_chr22_vi,"C")/cpgs_chr22_lengths)*(letterFrequency(cpg_chr22_vi,"G")/cpgs_chr22_lengths)
obs2exp_CG = obs_CG/exp_CG
q7=mean(obs2exp_CG)


TATA_boxes_count= countPattern("TATAAA", Hsapiens$chr22) + 
                  countPattern("TATAAA", reverseComplement(Hsapiens$chr22))
q8=TATA_boxes_count


library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb = TxDb.Hsapiens.UCSC.hg19.knownGene


transcripts_chr22 = subsetByOverlaps(transcripts(txdb),
                                     gr_chr22, ignore.strand=TRUE)
promoters_chr22 = promoters(transcripts_chr22,
                            upstream = 900,downstream = 100)
cdseq_chr22 = subsetByOverlaps(exons(txdb),gr_chr22,ignore.strand=TRUE)


proms_cds = subsetByOverlaps(promoters_chr22,cdseq_chr22)
proms_cds = unique(proms_cds)
proms_cds_vi = Views(chr22,ranges(proms_cds))

TATA_matched = matchPattern("TATAAA",proms_cds_vi)
TATA_ov_pcv = findOverlaps(as(TATA_matched,"IRanges"),
                           as(proms_cds_vi,"IRanges"))
q9 = length(unique(subjectHits(TATA_ov_pcv))) #193 is the right answer but how ?


tl_chr22 <- transcriptLengths(txdb,with.cds_len = TRUE)
tl_chr22  <- tl_chr22[tl_chr22$cds_len > 0,]
trans_eval <- promoters_chr22[mcols(promoters_chr22)$tx_id %in% tl.chr22$tx_id]
q10=sum(coverage(trans_eval)$chr22 > 1)



#
#H3K27me3_E055 = qhs[["AH30313"]]
#H3K27me3_chr22_E055 = subset(H3K27me3_E055, seqnames=="chr22")
#H3K27me3_chr22_vi_E055 = Views(Hsapiens,H3K27me3_chr22_E055)