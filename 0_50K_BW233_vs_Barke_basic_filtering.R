library(tidyverse)

# 50K recombination data is a tab separated table with samples as rows and 
# markers as columns. Values in the table represent the genotype for each 
# marker in each line represented by base symbols, e.g. "GG"
fiftyk <- read.table("Bw233_vs_Barke.txt",
                     header = TRUE,
                     na.strings = "--")

# starting with 43800 markers

# Morex V3 marker positions
snp_postions <- read.table("SNPPositions_50k_on_MorexV3.txt",
                           header = T)
rownames(snp_postions) <- snp_postions$marker

# invert 50K data table
colnames(fiftyk)
fiftyk <- rename(fiftyk,
                 line = LINE.MARKER)

rownames(fiftyk) <- fiftyk$line

tfiftyk <- t(fiftyk)
colnames(tfiftyk)
rownames(tfiftyk)
tfiftyk <- tfiftyk[-1,]

# replace the dot(.) in imported marker names with a dash(-)
rownames(tfiftyk) <- gsub("\\.",
                          "-",
                          rownames(tfiftyk))

# merge morex v3 marker positions with 50K results 
merged <- merge(tfiftyk,
                snp_postions,
                by = 0)

# 42912 markers map to Morex v3

# get a list of the markers that couldn't be found in the V3 marker mapping
no_positions <- setdiff(rownames(tfiftyk),
                         rownames(snp_postions))

# remove monomorphic markers between parental genotypes
rownames(merged) <- merged$Row.names
merged <- merged[,-1]
colnames(merged)

nomono_parents <- merged[merged[,53] != merged[,90],]

# leaves 13414 het markers markers

# remove any markers with no call for any line 
no_blank <- nomono_parents[rowSums(is.na(nomono_parents)) != ncol(nomono_parents), ]

# leaves 12730 markers

# remove markers with any missing value
nafilter <- no_blank

nafilter <- na.omit(nafilter)
# leaves 12696 markers

# -----------------------------------------------------------------------------

# change calls to indication of match to parental genotype
# 0 == Barke,
# 2 == Bw233,
# 1 == Het
numbered_allele <- nafilter

colnames(numbered_allele)

for (x in c(1:52,54:89,91:192)) {
  numbered_allele[,x] <- ifelse(numbered_allele[,x] == numbered_allele[,53],
                        0,
                        ifelse(numbered_allele[,x] == numbered_allele[,90],
                               2,
                               1))
}

# put markers in order of chromosome and then physical position 
numbered_allele <- numbered_allele[order(numbered_allele$Chromosome_v3,
                                         numbered_allele$Position_v3),]

save(numbered_allele,
     file = "markers_sorted_post_filtering.Rdata")
