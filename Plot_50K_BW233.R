# Make overlapping histograms of BW233 vs Barke 50K iSelect recombination data

library(tidyverse)

one <- read.csv("chr1H.csv")
two <- read.csv("chr2H.csv")
three <- read.csv("chr3H.csv")
four <- read.csv("chr4H.csv")
five <- read.csv("chr5H.csv")
six<- read.csv("chr6H.csv")
seven <- read.csv("chr7H.csv")

zones <- read.csv("zone_coords.csv")
chr_file <- read.table("morex_v3.txt")

morexv3_50K <- read.csv("SNPpositionsFixed.csv")
barlex_iselect <- read.csv("iSelectMarker.csv")

colnames(zones)
zones$chromosome <- gsub("Hv_chr_",
                         "chr",
                         zones$seq_id)
zones$H <- "H"

zones <- unite(zones,
               "chromosome",
               chromosome:H,
               remove = TRUE,
               sep = "")

colnames(zones)
zones <- zones[,c(7,6,2,3)]

one_uppers <- subset(zones,
                     chromosome=="chr1H")
one_uppers <- one_uppers$seq_upper_coor[1:4]

two_uppers <- subset(zones,
                     chromosome=="chr2H")
two_uppers <- two_uppers$seq_upper_coor[1:4]

three_uppers <- subset(zones,
                       chromosome=="chr3H")
three_uppers <- three_uppers$seq_upper_coor[1:4]

four_uppers <- subset(zones,
                      chromosome=="chr4H")
four_uppers <- four_uppers$seq_upper_coor[1:4]

five_uppers <- subset(zones,
                      chromosome=="chr5H")
five_uppers <- five_uppers$seq_upper_coor[1:4]

six_uppers <- subset(zones,
                     chromosome=="chr6H")
six_uppers <- six_uppers$seq_upper_coor[1:4]

seven_uppers <- subset(zones,
                       chromosome=="chr7H")
seven_uppers <- seven_uppers$seq_upper_coor[1:4]

all <- rbind(one,
             two,
             three,
             four,
             five,
             six,
             seven)

all$marker <- gsub("JHI_Hv50k_2016_",
                   "JHI-Hv50k-2016-",
                   all$marker)

all$marker <- gsub("12-",
                   "12_",
                   all$marker)

all$marker <- gsub("11-",
                   "11_",
                   all$marker)

v1tov3 <- merge(all,
                morexv3_50K,
                by="marker")


missing_markers <- setdiff(all$marker,
                           v1tov3$marker)

barlex_iselect_found <- barlex_iselect[barlex_iselect$iSelect.Marker %in% missing_markers,]

missing_markers

colnames(v1tov3)

found <- c("chr1H",
           24807848,
           "JHI-Hv50k-2016-16174")

morexv3_50K <- rbind(morexv3_50K,
                     found)

v1tov3 <- merge(all,
                morexv3_50K,
                by="marker")

setdiff(all$marker,
        v1tov3$marker)

chr_stats <- rename(chr_file,
                    Chromosome_v3 = V1,
                    chr_start = V2,
                    chr_end = V3,
                    centromere = V4)

v1tov3 <- merge(v1tov3,
                chr_stats, 
                by ="Chromosome_v3")

v1tov3$Position_v3 <- as.numeric(v1tov3$Position_v3)

# write.csv(v1tov3,
#           "markers_v3_coordinates.csv",
#           row.names = F)

# create percentage bins
v1tov3$percent <- (v1tov3$Position_v3/v1tov3$chr_end)*100

v1tov3$bin <- cut(as.numeric(v1tov3$percent),
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

colnames(v1tov3)

# get per position bin CO number
bin_total <- v1tov3 %>% group_by(bin) %>%
  summarise(Barke = sum(pop1_COs),
            Bw233 = sum(pop2_COs))

IBM <- c("#648FFF",
         "#DC267F",
         "#FFB000")

colnames(bin_total)

long <- pivot_longer(bin_total,
                     c(2:3),
                     names_to = "genotype",
                     values_to = "Crossovers")

colnames(long)

ggplot(long, aes(prcntg, Crossovers,
                 fill = genotype)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(labels = c("Bw233", "Barke"),
                    values = c(IBM[1], IBM[2])) +
  theme(axis.text.x = element_text(angle = 90,
                                   hjust = 1,
                                   vjust = 0.5)) +
  ylab("Average Crossovers") +
  xlab("Genomic position percentile") +
  theme(legend.position = "bottom")

ggplot(v1tov3, aes(bin)) +
  geom_bar(aes(y=pop2_COs),
           stat = "identity",
           fill = IBM[2]) +
  geom_bar(aes(y=pop1_COs),
           stat = "identity",
           fill = IBM[1],
           alpha = 0.6) +
  theme(axis.text.x = element_text(angle = 90,
                                   hjust = 1,
                                   vjust = 0.5)) +
  ylab("Crossovers") +
  xlab("Genomic posiiton percentile") +
  theme(legend.position = "bottom")

# ggsave("avg_COs_vs_percentile_bp_legend.svg",
#        width = 183,
#        units = "mm")

# get per position bin total CO number per chromosome
bin_total <- v1tov3 %>% group_by(bin, chromosome) %>%
  summarise(Barke = sum(pop1_COs),
            Bw233 = sum(pop2_COs))

bin_total <- complete(bin_total,
                      chromosome,
                      bin,
                      fill = 0)
colnames(bin_total)



colnames(bin_total)

long <- pivot_longer(bin_total,
                     c(3:4),
                     names_to = "genotype",
                     values_to = "Crossovers")

colnames(long)
colnames(bin_total)
colnames(v1tov3)
levels(long$bin)

# dplyr doesn't like "bin" for some reason (?)
long$prcntg <- long$bin

long <- long[,-1]

bin_total$percnt <- bin_total$bin

bin_total <- bin_total[,-1]

allbins <- complete(bin_total,
                    chromosome,
                    percnt,
                    fill = list(Barke = 0,
                                Bw233 = 0))

one_pos <- subset(allbins,
                  chromosome == "chr1H")
two_pos <- subset(allbins,
                  chromosome == "chr2H")
thr_pos <- subset(allbins,
                  chromosome == "chr3H")
fou_pos <- subset(allbins,
                  chromosome == "chr4H")
fiv_pos <- subset(allbins,
                  chromosome == "chr5H")
six_pos <- subset(allbins,
                  chromosome == "chr6H")
sev_pos <- subset(allbins,
                  chromosome == "chr7H")

every_nth = function(n) {
  return(function(x) {x[c(TRUE, rep(FALSE, n - 1))]})
}


one_bin <- ggplot(one_pos, aes(percnt)) +
  geom_bar(aes(y=Bw233),
           stat = "identity",
           fill = IBM[2]) +
  geom_bar(aes(y=Barke),
           stat = "identity",
           fill = IBM[1],
           alpha = 0.6) +
  theme(axis.text.x = element_text(angle = 90,
                                   hjust = 1,
                                   vjust = 0.5)) +
  scale_x_discrete(breaks = every_nth(n = 5)) +
  ylab("Crossovers") +
  xlab("Genomic posiiton percentile") +
  theme(legend.position = "bottom")
one_bin

two_bin <- ggplot(two_pos, aes(percnt)) +
  geom_bar(aes(y=Bw233),
           stat = "identity",
           fill = IBM[2]) +
  geom_bar(aes(y=Barke),
           stat = "identity",
           fill = IBM[1],
           alpha = 0.6) +
  theme(axis.text.x = element_text(angle = 90,
                                   hjust = 1,
                                   vjust = 0.5)) +
  scale_x_discrete(breaks = every_nth(n = 5)) +
  ylab("Crossovers") +
  xlab("Genomic posiiton percentile") +
  theme(legend.position = "bottom")
two_bin

thr_bin <- ggplot(thr_pos, aes(percnt)) +
  geom_bar(aes(y=Bw233),
           stat = "identity",
           fill = IBM[2]) +
  geom_bar(aes(y=Barke),
           stat = "identity",
           fill = IBM[1],
           alpha = 0.6) +
  theme(axis.text.x = element_text(angle = 90,
                                   hjust = 1,
                                   vjust = 0.5)) +
  scale_x_discrete(breaks = every_nth(n = 5)) +
  ylab("Crossovers") +
  xlab("Genomic posiiton percentile") +
  theme(legend.position = "bottom")
thr_bin

fou_bin <- ggplot(fou_pos, aes(percnt)) +
  geom_bar(aes(y=Bw233),
           stat = "identity",
           fill = IBM[2]) +
  geom_bar(aes(y=Barke),
           stat = "identity",
           fill = IBM[1],
           alpha = 0.6) +
  theme(axis.text.x = element_text(angle = 90,
                                   hjust = 1,
                                   vjust = 0.5)) +
  scale_x_discrete(breaks = every_nth(n = 5)) +
  ylab("Crossovers") +
  xlab("Genomic posiiton percentile") +
  theme(legend.position = "bottom")
fou_bin

fiv_bin <- ggplot(fiv_pos, aes(percnt)) +
  geom_bar(aes(y=Bw233),
           stat = "identity",
           fill = IBM[2]) +
  geom_bar(aes(y=Barke),
           stat = "identity",
           fill = IBM[1],
           alpha = 0.6) +
  theme(axis.text.x = element_text(angle = 90,
                                   hjust = 1,
                                   vjust = 0.5)) +
  scale_x_discrete(breaks = every_nth(n = 5)) +
  ylab("Crossovers") +
  xlab("Genomic posiiton percentile") +
  theme(legend.position = "bottom")
fiv_bin

six_bin <- ggplot(six_pos, aes(percnt)) +
  geom_bar(aes(y=Bw233),
           stat = "identity",
           fill = IBM[2]) +
  geom_bar(aes(y=Barke),
           stat = "identity",
           fill = IBM[1],
           alpha = 0.6) +
  theme(axis.text.x = element_text(angle = 90,
                                   hjust = 1,
                                   vjust = 0.5)) +
  scale_x_discrete(breaks = every_nth(n = 5)) +
  ylab("Crossovers") +
  xlab("Genomic posiiton percentile") +
  theme(legend.position = "bottom")
six_bin

sev_bin <- ggplot(sev_pos, aes(percnt)) +
  geom_bar(aes(y=Bw233),
           stat = "identity",
           fill = IBM[2]) +
  geom_bar(aes(y=Barke),
           stat = "identity",
           fill = IBM[1],
           alpha = 0.6) +
  theme(axis.text.x = element_text(angle = 90,
                                   hjust = 1,
                                   vjust = 0.5)) +
  scale_x_discrete(breaks = every_nth(n = 5)) +
  ylab("Crossovers") +
  xlab("Genomic posiiton percentile") +
  theme(legend.position = "bottom")
sev_bin

library(cowplot)

long$chromosome <- as.factor(long$chromosome)

chroms <- levels(long$chromosome)

plot_grid(one_bin,
          two_bin,
          thr_bin,
          fou_bin,
          fiv_bin,
          six_bin,
          sev_bin,
          ncol = 2,
          nrow = 4,
          labels = chroms,
          label_y = 1)

ggsave("individual_chroms_binned.svg",
       width = 180,
       height = 240,
       units = "mm")

