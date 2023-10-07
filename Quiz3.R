packages <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg))
    install.packages(new.pkg, dependencies = TRUE,
                     repos='http://cran.rstudio.com/')
  sapply(pkg, require, character.only = TRUE)
}

packages(c("GEOquery", "Biobase", "ALL", "hgu95av2.db", "GenomicRanges", 
           "airway", "biomaRt", "dplyr", "minfiData", 
           "TxDb.Hsapiens.UCSC.hg19.knownGene", "AnnotationHub"))
library(ALL)
data(ALL)
q1 = mean(exprs(ALL)[,5])

# connect to Ensembl
mart <- useMart(host = 'feb2014.archive.ensembl.org', 
                biomart = "ENSEMBL_MART_ENSEMBL")
ensembl <- useDataset("hsapiens_gene_ensembl", mart)
names <- featureNames(ALL)
# find affymatrix attributes
attrs <- listAttributes(ensembl, page = "feature_page")
# return results
result <- getBM(attributes = c("affy_hg_u95av2", 
                               "ensembl_gene_id", 
                               "chromosome_name"),
                filters = "affy_hg_u95av2", 
                values = names,
                mart = ensembl)
tail(result)

prob_set <- result %>%
  group_by (affy_hg_u95av2) %>%
  summarise(prob_count = n())
tail(prob_set)

q2=sum(prob_set$prob_count > 1)

result_autosome <- subset(result, chromosome_name < 23)
prob_set_autosome <- result_autosome %>%
  group_by (affy_hg_u95av2) %>%
  summarise(prob_count = n())

tail(prob_set_autosome)

q3 = sum(prob_set_autosome$prob_count > 0)

data(MsetEx)
#pData(MsetEx)  # 5723646052_R04C01
sample_2 <- MsetEx[,2]  # return a MethylSet for sample #2

q4 = mean(getMeth(sample_2)) # produce mean of MethylSet for sample #2


query <- getGEO("GSE788")
data <- query[[1]]
#pData(data)
GSM9024 <- data[,2]

q5 = mean(exprs(GSM9024))

# library(airway)
data(airway)

q6 = mean(airway$avgLength)

SRR1039512 <- airway[,3]
counts <- assay(SRR1039512, "counts")

q7 = sum(counts >= 1) # summation(subset counts which are greater than or equal to 1)

# library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
exons <- exons(txdb) # extract exon info

autosome <- paste0("chr", c(1:22)) # autosomes = chr1:22
df <- data.frame(seqnames = autosome)
exons <- keepSeqlevels(exons, autosome) # use only autosomal exons
ncbiStyleLevels <- mapSeqlevels(seqlevels(exons),"NCBI")
exons <- renameSeqlevels(exons, ncbiStyleLevels)
subset <- subsetByOverlaps(airway, exons)

subset

#q8=26276

SRR1039508 <- airway[,1]
subset_SRR1039508 <- subsetByOverlaps(SRR1039508, exons)# use only autosomal exons
counts <- assay(SRR1039508, "counts")
subset_counts <- assay(subset_SRR1039508, "counts")

q9 = sum(subset_counts)/sum(counts) # probability (where: 0 < p < 1)


ah <- AnnotationHub()
qah_h1 <- query(ah, c("E096", "H3K4me3"))
h1 <- qah_h1[["AH30596"]] # AH30596 | E096-H3K4me3.narrowPeak.gz 
h1 <- keepSeqlevels(h1, autosome)
h1 <- renameSeqlevels(h1, ncbiStyleLevels)
t <- range(rowRanges(subset_SRR1039508))
ncbiByGroup <- extractSeqlevelsByGroup(species = "Homo sapiens", 
                                       style = "NCBI", 
                                       group = "auto")
t <- keepSeqlevels(t, ncbiByGroup)
p <- promoters(t)
ov <- subsetByOverlaps(p, h1)
t2 <- subsetByOverlaps(SRR1039508, ov)
counts <- assay(t2, "counts")

q10 = median(counts)