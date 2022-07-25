# Make a heatmap of MaxQuant output from HvST1 autoubiquitination gel fraction
# MS/MS

library(reshape2)
intensity <- read.csv("ST1/MS_Intensity.csv")
rownames(intensity) <- intensity$Protein
intensity <- intensity[,-1]
t.intensity <- t(intensity)
t.intensity <- as.data.frame(t.intensity)
t.intensity$total <- rowSums(t.intensity[,1:292])
t.intensity.scaled <- t.intensity[1:292]/t.intensity$total
intensity.scaled <- t(t.intensity.scaled)
intensity.scaled <- as.data.frame(intensity.scaled)

intensity.scaled <- intensity.scaled[c(1,46,47,50,51,52), ]
intensity.scaled$Protein <- rownames(intensity.scaled)

long.data <- melt(intensity.scaled,
                  id="Protein",
                  variable.name="Sample",
                  value.name = "intensity")
hist(long.data$intensity, breaks = 100, col = "#d95f0e")
long.data$Sample <- gsub("Intensity.", "", long.data$Sample)

# Heatmap ---------------------------------------------------------------------
library(ggplot2)
library(viridis)

autoub.heat <- ggplot(long.data, aes(Sample, Protein)) +
  geom_tile(aes(fill = intensity), colour = "white") +
  scale_fill_viridis() +
  theme(axis.text.x = element_text(),
        axis.text.y = element_text(size = 6),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  labs(fill = "Scaled Intensity")
autoub.heat

autoub.heat + geom_text(aes(label = round(intensity, 2)), color = "white")

ggsave("AutoUbiqMaxQuant_scaled_labels.svg",
       width=17.4,
       units = "cm",
       dpi=600)
