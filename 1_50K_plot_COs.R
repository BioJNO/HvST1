# Plot COs in BW233 vs Barke 50K iSelect recombination data
library(tidyverse)
library(multcomp)
library(multcompView)
library(ggpubr)
library(cowplot)

chr_file <- read.table("morex_v3.txt")
load("chromosome_CO_counts.Rdata")

chr_stats <- rename(chr_file,
                    chromosome = V1,
                    chr_start = V2,
                    chr_end = V3,
                    centromere = V4)

co_count_all <- rbind(chr1H_CO,
                      chr2H_CO,
                      chr3H_CO,
                      chr4H_CO,
                      chr5H_CO,
                      chr6H_CO,
                      chr7H_CO)

co_count_all <- merge(co_count_all,
                      chr_stats, 
                      by ="chromosome")

colnames(co_count_all)

# work out the per chromosome increase
chr_total <- co_count_all %>% group_by(chromosome) %>%
  summarise(Hvst1 = sum(aa_total_COs),
            HvST1 = sum(bb_total_COs))

chr_total$perc_increase <- (chr_total$Hvst1/chr_total$HvST1)-1

write.csv(chr_total,
          "chromosome_total_crossovers.csv",
          row.names = F)


# create percentage bins ------------------------------------------------------
co_count_all$percent <- (co_count_all$position/co_count_all$chr_end)*100

co_count_all$perc_bin <- cut(co_count_all$percent,
                             breaks = seq(0,
                               100,
                               by = 2),
                             labels = paste(seq(0,
                                                98,
                                                by = 2),
                                            seq(2,
                                                100,
                                                by = 2),
                                            sep = '-' ))

# Plot a summary of recombination in all chromosomes by bin --------------------
# get per position bin CO number
bin_total <- co_count_all %>% group_by(perc_bin) %>%
  summarise(Hvst1 = sum(aa_total_COs),
            HvST1 = sum(bb_total_COs))

IBM <- c("#648FFF",
         "#DC267F",
         "#FFB000")

# function to plot every nth x axis label
every_nth = function(n) {
  return(function(x) {x[c(TRUE, rep(FALSE, n - 1))]})
}

colnames(bin_total)

colors <- c("Hvst1" = IBM[1],
            "HvST1" = IBM[2])

all_chr_WT <- ggplot(bin_total, aes(perc_bin)) +
           geom_bar(aes(y=HvST1,
                        fill = "HvST1"),
                    stat = "identity") +
           scale_fill_manual(values = colors,
                             name = "Genotype") +
           theme(axis.text.x = element_text(angle = 90,
                                            hjust = 1,
                                            vjust = 0.5,
                                            size = 8),
                 axis.text.y = element_text(size = 8),
                 axis.title = element_text(size = 8),
                 legend.position = "bottom",
                 legend.title = element_text(size = 8),
                 legend.text = element_text(size = 8,
                                            face = "italic")) +
           scale_x_discrete(breaks = every_nth(n = 4)) +
           ylab("Total Crossovers") +
           xlab("Genomic positon percentile") +
  ylim(0,250)
all_chr_WT

all_chr_mt <- ggplot(bin_total, aes(perc_bin)) +
  geom_bar(aes(y=Hvst1,
               fill = "Hvst1"),
           stat = "identity") +
  scale_fill_manual(values = colors,
                    name = "Genotype") +
  theme(axis.text.x = element_text(angle = 90,
                                   hjust = 1,
                                   vjust = 0.5,
                                   size = 8),
        axis.text.y = element_text(size = 8),
        axis.title = element_text(size = 8),
        legend.position = "bottom",
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 8,
                                   face = "italic")) +
  scale_x_discrete(breaks = every_nth(n = 4)) +
  ylab("Total Crossovers") +
  xlab("Genomic positon percentile") +
  ylim(0,250)
all_chr_mt

plot_grid(all_chr_WT,
          all_chr_mt,
          labels = "auto")

ggsave("Fig4.png",
       width = 183,
       height = 80,
       units = "mm",
       dpi = 600)

ggsave("Fig4.svg",
       width = 183,
       height = 80,
       units = "mm")


# Plot recombination per chromosome -------------------------------------------
# get per position bin total CO number per chromosome
bin_total <- co_count_all %>% group_by(perc_bin, chromosome) %>%
             summarise(Hvst1 = sum(aa_total_COs),
                       HvST1 = sum(bb_total_COs))

# tidyr hates this column ID for some reason and simply will not work unless 
# we change it.
bin_total$x <- bin_total$perc_bin
colnames(bin_total)
bin_total <- bin_total[,-1]

bin_total <- complete(bin_total,
                      chromosome,
                      x,
                      fill = list(Hvst1 = 0,
                                  HvST1 = 0))
colnames(bin_total)

long <- pivot_longer(bin_total,
                     c(3:4),
                     names_to = "genotype",
                     values_to = "Crossovers")

colnames(long)

one_pos <- subset(bin_total,
                  chromosome == "chr1H")
two_pos <- subset(bin_total,
                  chromosome == "chr2H")
thr_pos <- subset(bin_total,
                  chromosome == "chr3H")
fou_pos <- subset(bin_total,
                  chromosome == "chr4H")
fiv_pos <- subset(bin_total,
                  chromosome == "chr5H")
six_pos <- subset(bin_total,
                  chromosome == "chr6H")
sev_pos <- subset(bin_total,
                  chromosome == "chr7H")

one_bin_mt <- ggplot(one_pos, aes(x)) +
  geom_bar(aes(y=Hvst1,
               fill = "Hvst1"),
           stat = "identity") +
  scale_fill_manual(values = colors,
                    name = "Genotype") +
  theme(axis.text.x = element_text(size = 8,
                                   angle = 90),
        axis.text.y = element_text(size = 8),
        axis.title = element_blank(),
        legend.position = "none") +
  ylim(0,60) +
  scale_x_discrete(breaks = every_nth(n = 4))
one_bin_mt

one_bin_WT <- ggplot(one_pos, aes(x)) +
  geom_bar(aes(y=HvST1,
               fill = "HvST1"),
           stat = "identity") +
  scale_fill_manual(values = colors,
                    name = "Genotype") +
  theme(axis.text.x = element_text(size = 8,
                                   angle = 90),
        axis.text.y = element_text(size = 8),
        axis.title = element_blank(),
        legend.position = "none") +
  ylim(0,60) +
  scale_x_discrete(breaks = every_nth(n = 4))
one_bin_WT

two_bin_mt <- ggplot(two_pos, aes(x)) +
  geom_bar(aes(y=Hvst1,
               fill = "Hvst1"),
           stat = "identity") +
  scale_fill_manual(values = colors,
                    name = "Genotype") +
  theme(axis.text.x = element_text(size = 8,
                                   angle = 90),
        axis.text.y = element_text(size = 8),
        axis.title = element_blank(),
        legend.position = "ntwo") +
  ylim(0,60) +
  scale_x_discrete(breaks = every_nth(n = 4))
two_bin_mt

two_bin_WT <- ggplot(two_pos, aes(x)) +
  geom_bar(aes(y=HvST1,
               fill = "HvST1"),
           stat = "identity") +
  scale_fill_manual(values = colors,
                    name = "Genotype") +
  theme(axis.text.x = element_text(size = 8,
                                   angle = 90),
        axis.text.y = element_text(size = 8),
        axis.title = element_blank(),
        legend.position = "ntwo") +
  ylim(0,60) +
  scale_x_discrete(breaks = every_nth(n = 4))
two_bin_WT

thr_bin_mt <- ggplot(thr_pos, aes(x)) +
  geom_bar(aes(y=Hvst1,
               fill = "Hvst1"),
           stat = "identity") +
  scale_fill_manual(values = colors,
                    name = "Genotype") +
  theme(axis.text.x = element_text(size = 8,
                                   angle = 90),
        axis.text.y = element_text(size = 8),
        axis.title = element_blank(),
        legend.position = "nthree") +
  ylim(0,60) +
  scale_x_discrete(breaks = every_nth(n = 4))
thr_bin_mt

thr_bin_WT <- ggplot(thr_pos, aes(x)) +
  geom_bar(aes(y=HvST1,
               fill = "HvST1"),
           stat = "identity") +
  scale_fill_manual(values = colors,
                    name = "Genotype") +
  theme(axis.text.x = element_text(size = 8,
                                   angle = 90),
        axis.text.y = element_text(size = 8),
        axis.title = element_blank(),
        legend.position = "nthree") +
  ylim(0,60) +
  scale_x_discrete(breaks = every_nth(n = 4))
thr_bin_WT

fou_bin_mt <- ggplot(fou_pos, aes(x)) +
  geom_bar(aes(y=Hvst1,
               fill = "Hvst1"),
           stat = "identity") +
  scale_fill_manual(values = colors,
                    name = "Genotype") +
  theme(axis.text.x = element_text(size = 8,
                                   angle = 90),
        axis.text.y = element_text(size = 8),
        axis.title = element_blank(),
        legend.position = "nfou") +
  ylim(0,60) +
  scale_x_discrete(breaks = every_nth(n = 4))
fou_bin_mt

fou_bin_WT <- ggplot(fou_pos, aes(x)) +
  geom_bar(aes(y=HvST1,
               fill = "HvST1"),
           stat = "identity") +
  scale_fill_manual(values = colors,
                    name = "Genotype") +
  theme(axis.text.x = element_text(size = 8,
                                   angle = 90),
        axis.text.y = element_text(size = 8),
        axis.title = element_blank(),
        legend.position = "nfou") +
  ylim(0,60) +
  scale_x_discrete(breaks = every_nth(n = 4))
fou_bin_WT

fiv_bin_mt <- ggplot(fiv_pos, aes(x)) +
  geom_bar(aes(y=Hvst1,
               fill = "Hvst1"),
           stat = "identity") +
  scale_fill_manual(values = colors,
                    name = "Genotype") +
  theme(axis.text.x = element_text(size = 8,
                                   angle = 90),
        axis.text.y = element_text(size = 8),
        axis.title = element_blank(),
        legend.position = "nfiv") +
  ylim(0,60) +
  scale_x_discrete(breaks = every_nth(n = 4))
fiv_bin_mt

fiv_bin_WT <- ggplot(fiv_pos, aes(x)) +
  geom_bar(aes(y=HvST1,
               fill = "HvST1"),
           stat = "identity") +
  scale_fill_manual(values = colors,
                    name = "Genotype") +
  theme(axis.text.x = element_text(size = 8,
                                   angle = 90),
        axis.text.y = element_text(size = 8),
        axis.title = element_blank(),
        legend.position = "nfiv") +
  ylim(0,60) +
  scale_x_discrete(breaks = every_nth(n = 4))
fiv_bin_WT

six_bin_mt <- ggplot(six_pos, aes(x)) +
  geom_bar(aes(y=Hvst1,
               fill = "Hvst1"),
           stat = "identity") +
  scale_fill_manual(values = colors,
                    name = "Genotype") +
  theme(axis.text.x = element_text(size = 8,
                                   angle = 90),
        axis.text.y = element_text(size = 8),
        axis.title = element_blank(),
        legend.position = "nsix") +
  ylim(0,60) +
  scale_x_discrete(breaks = every_nth(n = 4))
six_bin_mt

six_bin_WT <- ggplot(six_pos, aes(x)) +
  geom_bar(aes(y=HvST1,
               fill = "HvST1"),
           stat = "identity") +
  scale_fill_manual(values = colors,
                    name = "Genotype") +
  theme(axis.text.x = element_text(size = 8,
                                   angle = 90),
        axis.text.y = element_text(size = 8),
        axis.title = element_blank(),
        legend.position = "nsix") +
  ylim(0,60) +
  scale_x_discrete(breaks = every_nth(n = 4))
six_bin_WT


sev_bin_mt <- ggplot(sev_pos, aes(x)) +
  geom_bar(aes(y=Hvst1,
               fill = "Hvst1"),
           stat = "identity") +
  scale_fill_manual(values = colors,
                    name = "Genotype") +
  theme(axis.text.x = element_text(size = 8,
                                   angle = 90),
        axis.text.y = element_text(size = 8),
        axis.title = element_blank(),
        legend.position = "nsev") +
  ylim(0,60) +
  scale_x_discrete(breaks = every_nth(n = 4))
sev_bin_mt

sev_bin_WT <- ggplot(sev_pos, aes(x)) +
  geom_bar(aes(y=HvST1,
               fill = "HvST1"),
           stat = "identity") +
  scale_fill_manual(values = colors,
                    name = "Genotype") +
  theme(axis.text.x = element_text(size = 8,
                                   angle = 90),
        axis.text.y = element_text(size = 8),
        axis.title = element_blank(),
        legend.position = "nsev") +
  ylim(0,60) +
  scale_x_discrete(breaks = every_nth(n = 4))
sev_bin_WT

long$chromosome <- as.factor(long$chromosome)

chroms <- levels(long$chromosome)

plot_grid(one_bin_WT,
          one_bin_mt,
          two_bin_WT,
          two_bin_mt,
          thr_bin_WT,
          thr_bin_mt,
          fou_bin_WT,
          fou_bin_mt,
          fiv_bin_WT,
          fiv_bin_mt,
          six_bin_WT,
          six_bin_mt,
          sev_bin_WT,
          sev_bin_mt,
          nrow = 7,
          ncol = 2,
          labels = c("chr1H",
                     "",
                     "chr2H",
                     "",
                     "chr3H",
                     "",
                     "chr4H",
                     "",
                     "chr5H",
                     "",
                     "chr6H",
                     "",
                     "chr7H",
                     ""),
          label_size = 8)

ggsave("FigS14.svg",
       width = 183,
       height = 200,
       units = "mm")

# Box plots each chromosome all individual CO counts

colnames(co_count_all)

# want chromosome totals per individual
indiv_totals <- pivot_longer(co_count_all,
             8:197,
             values_to = "crossovers",
             names_to = "individuals")

colnames(indiv_totals)

indiv_totals <- indiv_totals[,c(1,13,14)]


indiv_grp <- indiv_totals %>% group_by(individuals, chromosome) %>%
  summarise(COs = sum(crossovers, na.rm = T))

mt_rows <- grep("AA", indiv_grp$individuals)
mt_rows <- indiv_grp$individuals[mt_rows]

indiv_grp$genotype <- ifelse(indiv_grp$individuals %in% mt_rows,
                             "Hvst1",
                             "HvST1")

indiv_grp$chromosome <- as.factor(indiv_grp$chromosome)
indiv_grp$genotype <- as.factor(indiv_grp$genotype)

# t test of difference between means
lapply(split(indiv_grp,
             factor(indiv_grp$chromosome)),
       function(x)t.test(data=x,
                         COs ~ genotype, paired=FALSE))

sigdat <- data.frame(x=c(0.875,
                         1.875,
                         2.875,
                         3.875,
                         4.875,
                         5.875,
                         6.875),
                     xend=c(1.125,
                            2.125,
                            3.125,
                            4.125,
                            5.125,
                            6.125,
                            7.125),
                     y=c(9,
                         8.5,
                         10,
                         9,
                         10,
                         7.5,
                         11),
                     annotation=c("***",
                                  "***",
                                  "**",
                                  "***",
                                  "**",
                                  "***",
                                  "***"))

ggplot(indiv_grp, aes(chromosome, COs)) +
  geom_boxplot(aes(fill = genotype,
                   colour = genotype),
               alpha = 0.5) +
  theme(legend.text = element_text(face = "italic",
                                   size = 12),
        axis.title = element_text(size = 12),
        legend.title = element_text(size = 12),
        axis.text = element_text(size = 12)) +
  scale_fill_manual(values = colors) +
  scale_color_manual(values = colors) +
  geom_signif(xmin = sigdat$x,
              xmax = sigdat$xend,
              y_position = sigdat$y,
              annotations = sigdat$annotation,
              tip_length = 0)

ggsave("crossover_boxplot.svg",
       width = 183,
       height = 120,
       units = "mm")

# sanity check sums
totals <- indiv_grp %>% group_by(genotype, chromosome) %>%
  summarise(COs = sum(COs, na.rm = T))
