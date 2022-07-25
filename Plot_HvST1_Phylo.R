# Plot HvST1 IQ-TREE phylogeny

# BiocManager::install("ggtree")
# BiocManager::install("treeio")
# BiocManager::install("phytools")

library(treeio)
library(ggtree)
library(phytools)

tree <- read.tree("HvST1_tree.txt")

mrca <- findMRCA(tree)

# nice cb friendly colours
IBM <- c("#648FFF",
         "#DC267F",
         "#FFB000")

ggtree(tree,
       ladderize = FALSE) +
  geom_tiplab(align = TRUE,
              size = 4) +
  geom_nodelab(hjust = -0.2,
               size = 4) +
  scale_color_manual(values = c(IBM[1],IBM[2])) +
  geom_strip("T.aestivum_CLR-like_TraesCS7A02G391700",
             "B.distachyon_100844152",
             label = "BOP clade",
             offset = 3.5,
             offset.text = .1,
             color = IBM[1],
             barsize = 1,
             fontsize = 4) +
  geom_strip("Z.palustris_GUJ93_ZPchr0006g42405",
             "Z.mays_ZmCLR1_Zm00001d046663",
             label = "PACMAD clade",
             offset = 3.5,
             offset.text = .1,
             color = IBM[2],
             barsize = 1,
             fontsize = 4) +
  geom_treescale(linesize = 1,
                 fontsize = 4) +
  hexpand(.2, direction = 1) +
  geom_highlight(node = 4)

ggsave("HvST1_tree.svg")

