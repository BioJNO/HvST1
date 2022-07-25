# Make a boxplot of BW233 vs Bowman metaphase spreads and check for statistical
# significance using aov and tukey honest significant difference

library(tidyverse)
library(multcomp)
library(multcompView)
library(ggpubr)

my_data <- read.csv("metphase_spreads_Isabelle.csv")

colnames(my_data)

long.data <- pivot_longer(my_data, 
                          2:4,
                          names_to = "type",
                          values_to ="Count")

long.data$type <- factor(long.data$type,
                            levels = c("Uni",
                                       "Rod",
                                       "Ring"))

colnames(long.data)

# Significant difference of expression at each stage in each sample ----------
anova <- aov(Count ~ Genotype + type,
             data = long.data)

summary(anova)

hondiff <- TukeyHSD(anova)

print(hondiff)

generate_label_df <- function(tukey, variable){
  
  # Extract labels and factor levels from Tukey post-hoc 
  Tukey.levels <- tukey[[variable]][,4]
  Tukey.labels <- data.frame(multcompLetters(Tukey.levels)['Letters'])
  
  # Labels must be in the same order as in the boxplot:
  Tukey.labels$variable=rownames(Tukey.labels)
  Tukey.labels=Tukey.labels[order(Tukey.labels$variable) , ]
  return(Tukey.labels)
}

generate_label_df(hondiff,
                  "Genotype:type")

# plot a boxplot with ggplot2
plot <- ggplot(data = long.data, aes(type, Count)) +
        geom_boxplot(aes(fill = Genotype,
                         color = Genotype),
                     alpha = 0.4) +
        theme(legend.position = "bottom")
plot

# save the plot as an svg
ggsave("metaphase_boxplot.svg",
       width = 89,
       height = 89,
       units = "mm")


