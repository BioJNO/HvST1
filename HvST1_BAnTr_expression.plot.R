# Plot expression of HvST1 in Barley anther and meiocyte transctiptome (BAnTr)

library(ggplot2)
library(multcomp)
library(ggpubr)

# Load and clean up data ---------------------------------------
exprdat <- read.csv("GeneExpressionR.csv")

hvst1 <- exprdat[exprdat$gene.id =="mikado.chr7_part2G2734", ]
hvst1 <- hvst1[,-c(1:2, 21)]
colnames(hvst1)

thvst1 <- t(hvst1)
thvst1 <- as.data.frame(thvst1)
thvst1$tissue <- c(rep("anther", 12),
                   rep("meiocyte", 6))
thvst1$stage <- c(rep("pre-meiosis", 3),
                  rep("leptotene/zygotene", 3),
                  rep("pachytene/diplotene",3),
                  rep("anaphase/telophase", 3),
                  rep("leptotene/zygotene", 3),
                  rep("pachytene/diplotene",3))

thvst1$cpm <- thvst1[,1]
thvst1$stage <- as.factor(thvst1$stage)
levels(thvst1$stage)
thvst1$stage <- factor(thvst1$stage, levels(thvst1$stage)[c(4,2,3,1)])
thvst1$sample <- gsub('.{1}$', '', rownames(thvst1))

thvst1$sample <- as.factor(thvst1$sample)
levels(thvst1$sample)
thvst1$sample <- factor(thvst1$sample, levels(thvst1$sample)[c(4,2,5,3,6,1)])

hvst1.bar <- ggplot(data=thvst1,
                    aes(sample, cpm)) +
             geom_boxplot(aes(fill=tissue)) +
             geom_point(aes(colour = tissue)) +
             stat_compare_means(method = "anova")
hvst1.bar

anther <- thvst1[thvst1$tissue == "anther", ]
meiocyte <- thvst1[thvst1$tissue == "meiocyte", ]

# Significant difference of expression at each stage in each sample ----------
expr.anova <- aov(cpm ~ sample, data = thvst1)
summary(expr.anova)

expr.tukey <- TukeyHSD(expr.anova)

generate_label_df <- function(tukey, variable){
  
  # Extract labels and factor levels from Tukey post-hoc 
  Tukey.levels <- tukey[[variable]][,4]
  Tukey.labels <- data.frame(multcompLetters(Tukey.levels)['Letters'])
  
  # Labels must be in the same order as in the boxplot:
  Tukey.labels$variable=rownames(Tukey.labels)
  Tukey.labels=Tukey.labels[order(Tukey.labels$variable) , ]
  return(Tukey.labels)
}

generate_label_df(expr.tukey, "sample")
