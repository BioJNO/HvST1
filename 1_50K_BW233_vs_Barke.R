# Filter 50K data for misplaced markers and F2 recombination events and plot
# the results 
library(runner)
library(tidyverse)
library(viridis)
library(cowplot)

load("markers_sorted_post_filtering.Rdata")

# pick a chromosome -----------------------------------------------------------
chr_sel <- "chr7H"

# Once run for all seven chromosomes save the results for further plotting
save(chr1H_CO,
     chr2H_CO,
     chr3H_CO,
     chr4H_CO,
     chr5H_CO,
     chr6H_CO,
     chr7H_CO,
     file = "chromosome_CO_counts.Rdata")

# Below here everything should run with no further changes --------------------
chroma <- subset(numbered_allele,
                 Chromosome_v3 == chr_sel)

# plot raw marker calls in each group in order for this chromosome ------------
colnames(chroma)
chroma[,53] <- 0
chroma[,90] <- 2
chroma_long <- pivot_longer(chroma,
                            1:192,
                            names_to = "line",
                            values_to = "call")

# keep markers in order of position on chromosome
chroma_long$marker <- factor(chroma_long$marker,
                             levels = chroma$marker)

raw_calls_aa <- ggplot(chroma_long[c(grep("AA", chroma_long$line)),],
                       aes(marker,
                           line)) +
  geom_tile(aes(fill = call)) +
  theme(axis.text = element_blank())
raw_calls_aa

raw_calls_bb <- ggplot(chroma_long[c(grep("BB", chroma_long$line)),],
                       aes(marker,
                           line)) +
  geom_tile(aes(fill = call)) +
  theme(axis.text = element_blank())
raw_calls_bb

# Deal with poorly mapped, or multi-mapped out of place markers ---------------
adjusted_call <- chroma

# remove parental genotypes from lines (53 = Barke, 90 = Bw233)
colnames(adjusted_call)
adjusted_call <- adjusted_call[,-c(53,90)]

# Using a sliding window function set the value of each call to the median 
# within a sliding window of 20 markers
adjusted_call <- apply(adjusted_call[,c(1:190)],
                       2,
                       function (x) runner(x,
                       k=20,
                       lag=1,
                       f = function(x) {median(x)}))

adjusted_call <- as.data.frame(adjusted_call)
adjusted_call[1,] <- adjusted_call[2,]

# plot adjusted call markers as above -----------------------------------------
adjusted_call$marker <- chroma$marker
adjusted_call$position <- chroma$Position_v3
adjusted_call$chromosome <- chroma$Chromosome_v3

adjusted_call_long <- pivot_longer(adjusted_call,
                                   1:190,
                                   names_to = "line",
                                   values_to = "call")

# keep markers in order
adjusted_call_long$marker <- factor(adjusted_call_long$marker,
                                    levels = chroma$marker)

adjusted_calls_aa <- ggplot(adjusted_call_long[c(grep("AA", adjusted_call_long$line)),],
                            aes(marker,
                                line)) +
  geom_tile(aes(fill = call)) +
  theme(axis.text = element_blank())
adjusted_calls_aa

adjusted_calls_bb <- ggplot(adjusted_call_long[c(grep("BB", adjusted_call_long$line)),],
                            aes(marker,
                                line)) +
  geom_tile(aes(fill = call)) +
  theme(axis.text = element_blank())
adjusted_calls_bb

# Assign crossover counts to adjusted call data -------------------------------
# 0.5 or 1.5 is a single CO
# if 1 comes in between a 0 and a 2 or vice versa that's a double CO
# NA for row one assume that marker matches second, wouldn't count a CO here anyway
co_count <- adjusted_call
colnames(co_count)
co_count$marker
for (x in 2:nrow(co_count)) {
  co_count[x,1:190] <- ifelse(adjusted_call[x,1:190] == 0.5 |
                           adjusted_call[x,1:190] == 1.5,
                             1,
                             ifelse(adjusted_call[x,1:190] == 1 &
                                      adjusted_call[x-1,1:190] == 0 &
                                      adjusted_call[x+1,1:190] == 2|
                                      adjusted_call[x,1:190] == 1 &
                                      adjusted_call[x-1,1:190] == 2 &
                                      adjusted_call[x+1,1:190] == 0,
                                      2,
                                      0))
}

co_count[1,1:190] <- 0

# plot crossover counts as they are now ---------------------------------------
colnames(co_count)
co_count$marker
co_count_long <- pivot_longer(co_count,
                              1:190,
                              names_to = "line",
                              values_to = "COs")

co_count_long$marker <- factor(co_count_long$marker,
                               levels = chroma$marker)

co_plot_unf_aa <- ggplot(co_count_long[c(grep("AA", co_count_long$line)), ],
                         aes(marker,
                             line)) +
  geom_tile(aes(fill = COs)) +
  theme(axis.text = element_blank()) +
  scale_fill_viridis(option = "magma")
co_plot_unf_aa

co_plot_unf_bb <- ggplot(co_count_long[c(grep("BB", co_count_long$line)), ],
                         aes(marker,
                             line)) +
  geom_tile(aes(fill = COs)) +
  theme(axis.text = element_blank()) +
  scale_fill_viridis(option = "magma")
co_plot_unf_bb

# crossovers at directly adjacent markers are probably not real, if any
# occur, get rid of them.
for (x in 1:nrow(co_count)) {
  co_count[x,1:190] <- ifelse(co_count[x,1:190] == 1 &
                         co_count[x+1,1:190] == 1 |
                         co_count[x,1:190] == 2 &
                         co_count[x+1,1:190] == 1 |
                         co_count[x,1:190] == 1 &
                         co_count[x+1,1:190] == 2,
                         0,
                         ifelse(co_count[x,1:190] == 0,
                                0,
                                co_count[x,1:190]))
}

# calculate the total number of crossovers per marker for each group ----------
colnames(co_count)
co_count$aa_total_COs <- rowSums(co_count[,c(grep("AA", names(co_count)))])
co_count$bb_total_COs <- rowSums(co_count[,c(grep("BB", names(co_count)))])

co_count$aa_recomb_freq <- co_count$aa_total_COs/95
co_count$bb_recomb_freq <- co_count$bb_total_COs/95

rownames(co_count) <- co_count$marker
rownames(adjusted_call) <- adjusted_call$marker

colnames(co_count)

# F2 recombination events show up as large monomorphic blocks which in turn
# results in individual markers with very high crossover counts due to a large
# number of individual lines with the same allelic state change in the same 
# place. Here we add a filter to total single marker COs to get rid of these
# markers.
remove_aa <- rownames(co_count[co_count$aa_total_COs > 5, ])
remove_bb <- rownames(co_count[co_count$bb_total_COs > 5, ])

colnames(co_count)

co_count_long_filt <- pivot_longer(co_count,
                                   1:190,
                                   names_to = "line",
                                   values_to = "COs")

# keep markers in order
co_count_long_filt$marker <- factor(co_count_long_filt$marker,
                                    levels = chroma$marker)

co_count_long_filt_AA <- co_count_long_filt[grep("AA",
                                                 co_count_long_filt$line),]
co_count_long_filt_AA <- co_count_long_filt_AA[!co_count_long_filt_AA$marker %in% remove_aa,]

co_count_long_filt_BB <- co_count_long_filt[grep("BB",
                                                      co_count_long_filt$line),]
co_count_long_filt_BB <- co_count_long_filt_BB[!co_count_long_filt_BB$marker %in% remove_bb,]

co_count_long_filt <- rbind(co_count_long_filt_AA,
                            co_count_long_filt_BB)

# plot filtered CO counts -----------------------------------------------------
co_plot_aa_filt <- ggplot(co_count_long_filt[grep("AA", co_count_long_filt$line),],
                      aes(marker,
                          line)) +
  geom_tile(aes(fill = COs)) +
  theme(axis.text = element_blank()) +
  scale_fill_viridis(option = "magma")
co_plot_aa_filt

co_plot_bb_filt <- ggplot(co_count_long_filt[grep("BB", co_count_long_filt$line),],
                          aes(marker,
                              line)) +
  geom_tile(aes(fill = COs)) +
  theme(axis.text = element_blank()) +
  scale_fill_viridis(option = "magma")
co_plot_bb_filt

# Make and save filtering plots -----------------------------------------------
plot_grid(raw_calls_aa,
          adjusted_calls_aa,
          co_plot_unf_aa,
          co_plot_aa_filt,
          ncol = 1,
          labels = "auto")

ggsave(paste(chr_sel,
             "_filtering_aa.png",
             sep = "_"),
       height = 183,
       width = 240,
       dpi = 1200,
       units = "mm")

plot_grid(raw_calls_bb,
          adjusted_calls_bb,
          co_plot_unf_bb,
          co_plot_bb_filt,
          ncol = 1,
          labels = "auto")

ggsave(paste(chr_sel,
             "_filtering_bb.png",
             sep = "_"),
       height = 183,
       width = 240,
       dpi = 1200,
       units = "mm")

# Save the final filtered dataset for each chromosome -------------------------
co_count_filt_final <- pivot_wider(co_count_long_filt,
                                   names_from = line,
                                   values_from = COs)

# re-calculate total COs
colnames(co_count)
co_count_filt_final$aa_total_COs <- rowSums(co_count_filt_final[,c(grep("AA", names(co_count_filt_final)))])
co_count_filt_final$bb_total_COs <- rowSums(co_count_filt_final[,c(grep("BB", names(co_count_filt_final)))])

co_count_filt_final$aa_recomb_freq <- co_count_filt_final$aa_total_COs/95
co_count_filt_final$bb_recomb_freq <- co_count_filt_final$bb_total_COs/95

# for markers removed in one group (now NA) set COs to 0
co_count_filt_final$aa_total_COs <- replace_na(co_count_filt_final$aa_total_COs,
                                        0)
co_count_filt_final$bb_total_COs <- replace_na(co_count_filt_final$bb_total_COs,
                                        0)

assign(paste(chr_sel, "CO", sep = "_"), co_count_filt_final)
