# Make box plots of HvST1 expression under salt, heat, and Drought stress 
# conditinos using data from the EoRNA database

library(tidyverse)

eorna <- read.csv("HvST1_Eorna.csv")

colnames(eorna)

levels(eorna$Experiment.name)


salt <- eorna[grepl("E1-SaltStressMurdochUniv", eorna$Experiment.name), ]
drought <- eorna[grepl("Drought", eorna$Experiment.name), ]
heat <- eorna[grepl("Heat", eorna$Experiment.name), ]
cold <- eorna[grepl("Cold", eorna$Experiment.name), ]
morex_16 <- eorna[grepl("16samples", eorna$Experiment.name), ]

IBM <- c("#648FFF",
         "#785EF0",
         "#DC267F",
         "#FE6100",
         "#FFB000")

salt_box <- ggplot(data = salt,
                   aes(Expt.Condition, BART1_0.p54930.001)) +
  geom_boxplot(aes(fill=Tissue)) +
  scale_fill_manual(values = IBM) +
  theme(axis.text=element_text(size=8)) +
  ylab("HvST1 Expression (TPM)")

salt_box

ggsave("st1_salt_stress.svg",
       width = 183,
       units = "mm")

drought_box <- ggplot(data = drought,
                   aes(Expt.Condition, BART1_0.p54930.001)) +
  geom_boxplot(aes(fill=Tissue)) +
  scale_fill_manual(values = IBM) +
  theme(axis.text=element_text(size=8)) +
  ylab("HvST1 Expression (TPM)")

drought_box

ggsave("st1_drought_stress.svg",
       width = 183,
       units = "mm")

heat_box <- ggplot(data = heat,
                      aes(Expt.Condition, BART1_0.p54930.001)) +
  geom_boxplot(aes(fill=Tissue)) +
  scale_fill_manual(values = IBM) +
  theme(axis.text=element_text(size=8)) +
  ylab("HvST1 Expression (TPM)")

heat_box

ggsave("st1_heat_stress.svg",
       width = 183,
       units = "mm")

cold_box <- ggplot(data = cold,
                   aes(Expt.Condition, BART1_0.p54930.001)) +
  geom_boxplot(aes(fill=Tissue)) +
  scale_fill_manual(values = IBM) +
  theme(axis.text=element_text(size=8)) +
  ylab("HvST1 Expression (TPM)")

cold_box

ggsave("st1_cold_stress.svg",
       width = 183,
       units = "mm")

colnames(morex_16)

morex <- ggplot(data = morex_16,
                   aes(Tissue, BART1_0.p54930.001)) +
  geom_boxplot(aes()) +
  theme(axis.text=element_text(size=8),
        axis.text.x = element_text(angle = 90,
                                   hjust = 1,
                                   vjust = 0)) +
  ylab("HvST1 Expression (TPM)")

morex

ggsave("st1_morex_16_tissue.svg",
       width = 183,
       units = "mm")
