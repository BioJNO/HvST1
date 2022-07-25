# Plot HvST1 gene models using package "genemodel"
# Using Barley Reference Transcriptome v2.18 (BaRTv2.18)

library(genemodel)

# BaRTv2.18
st1 <- read.csv("ST1_gene_model_BaRT.csv")

genemodel.plot(model = st1,
      start = 530327838,
      bpstop = 530329122,
      orientation = "forward",
      xaxis = T)

mutation.plot(530328566,
              530328566,
              text="+G",
              col="#D81B60",
              drop=-.15)

mutation.plot(530328776,
              530328778,
              text="PTC",
              col="#D81B60",
              drop=-.15)

mutation.plot(530328647,
              530328770,
              text="RING domain",
              col="#004D40",
              drop=-0)

