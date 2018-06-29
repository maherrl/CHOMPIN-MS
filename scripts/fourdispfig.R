################################################
## Working script for beta diversity figures ##
## for Chompin manuscript.
## Created 6/28/18
## Created by Becca Maher
#-----------------------------------------------




# clear workspace--------------------------------
rm(list=ls())

# load libraries--------------------
library(ggplot2)
library(cowplot)
library(stringr)

# set working directory-------------------

# Load data
betdist <- read.csv(file = "../CHOMPIN-MS/data/betdist.csv")

# Data adjustments
betdist$interaction = factor(betdist$interaction, levels=c("Control", "High T", "Scarred", 
                                                             "NO3-", "NH4+", "High T, Scarred", 
                                                             "NO3-, Scarred", "NH4+, Scarred", 
                                                             "NO3-, High T", "NH4+, High T", 
                                                             "NO3-, High T, Scarred", 
                                                             "NH4+, High T, Scarred"))
betdist$stress <- factor(betdist$stress, levels = c("none", "single", "triple", "double"))


# Set color palette
enrichPal <- c('#00BFC4', '#7CAE00', '#F8766D', '#C77CFF')
# Significance codes are obtained from pairwise adonis and betadisper commands (beta_stats.R)
# Four distance measures: Weighted Unifrac (wu), Unweighted Unifrac (un), Bray Curtis (bc), Binary Jaccard (bj)
sigwuman <- c("A","B","ABC","AB","ABC","AC","ABC","ABC","AC","ABC","C","C") # 6/28/18
sigunman <- c("AB","A","C","A-D","A-D","ABC","A-D","AD","D","A-D","AB","BCD")# 6/28/18
sigbcman <- c("A","ABC","AB","B","BC","BC","BC","ABC","ABC","BC","C","C")# 6/28/18
sigbjman <- c() 

#sigwudisp <- c("A","B,D,E","E","B-E","B","D","B-E","C","B,E","E","B-E","B-E")
#sigwuman <- c("A", "B", "B", "A,B,C", "C", "B,C", "A,B,C", "A,C", "B,C", "A,B,C", "A,B,C", "B,C") #Updated

wu <- ggplot(betdist, aes(x = reorder(interaction, distwu, FUN = median), y = distwu, color = stress)) + 
  geom_boxplot(aes(color=stress)) + 
  scale_color_manual(values=enrichPal, name = "Stressors")  +
  theme_light() +
  theme(strip.background = element_blank(),
        axis.title.x = element_blank(), 
        axis.text.x = element_text(angle = 35, hjust = 1),
        legend.position = "none") +
  #ggtitle("Weighted Unifrac") + 
  ylim(0,1) +
  ylab("Weighted Unifrac Distance") +
  stat_summary(geom = 'text', label = sigwuman, fun.y = max, vjust = -1, size = 2.5) +
  annotate("text", label = "p<0.05", x=2, y =1)


un <- ggplot(betdist, aes(x = reorder(interaction, distun, FUN = median), y = distun, color = stress)) + 
  geom_boxplot(aes(color=stress)) + 
  scale_color_manual(values=enrichPal, name = "Stressors")  +
  theme_light() +
  theme(strip.background = element_blank(),
        axis.title.x = element_blank(), 
        axis.text.x = element_text(angle = 35, hjust = 1),
        legend.position = "none") +
  #ggtitle("Unweighted Unifrac") + 
  ylim(0.5,1) +
  ylab("Unweighted Unifrac Distance") +
  stat_summary(geom = 'text', label = sigunman, fun.y = max, vjust = -1, size = 2.5) +
  annotate("text", label = "p<0.01", x=2, y =1)

bc <- ggplot(betdist, aes(x = reorder(interaction, distbc, FUN = median), y = distbc, color = stress)) + 
  geom_boxplot(aes(color=stress)) + 
  scale_color_manual(values=enrichPal, name = "Stressors")  +
  theme_light() +
  theme(strip.background = element_blank(),
        axis.title.x = element_blank(), 
        axis.text.x = element_text(angle = 35, hjust = 1),
        legend.position = "none") +
  #ggtitle("Bray Curtis") + 
  ylim(0,1) +
  ylab("Bray Curtis Distance") +
  stat_summary(geom = 'text', label = sigbcman, fun.y = max, vjust = -1, size = 2.5) +
  annotate("text", label = "p<0.05", x=2, y =1)

bj <- ggplot(betdist, aes(x = reorder(interaction, distbj, FUN = median), y = distbj, color = stress)) + 
  geom_boxplot(aes(color=stress)) + 
  scale_color_manual(values=enrichPal, name = "Stressors")  +
  theme_light() +
  theme(strip.background = element_blank(),
        axis.title.x = element_blank(), 
        axis.text.x = element_text(angle = 35, hjust = 1),
        legend.position = "none") +
#  ggtitle("Binary Jaccard") + 
  ylim(0.5,1) +
  ylab("Binary Jaccard Distance") +
  annotate("text", label = "p<0.01", x=2, y =1) #+
#  stat_summary(geom = 'text', label = sigbjman, fun.y = max, vjust = -1, size = 2.5)
  
plot_grid(wu, bc, un, bj, ncol = 2, align = 'v', rel_heights = c(1,1,2,2))

#############################################
# Graphs with main effects
bj <- ggplot(betdist, aes(x = reorder(interaction, distbj, FUN = median), y = distbj, color = stress)) + 
  geom_boxplot(aes(color=stress)) + 
  scale_color_manual(values=enrichPal, name = "Stressors")  +
  theme_light() +
  theme(strip.background = element_blank(),
        text = element_text(size = 10), axis.title.x = element_blank(), 
        axis.text.x = element_text(angle = 45, hjust = 1),
        title = element_text(size = 9), legend.position = "none") +
  #ggtitle("Weighted Unifrac") + 
  ylim(0.5,1) +
  ylab("Beta Diversity")
#  stat_summary(geom = 'text', label = sigundisp, fun.y = max, vjust = -1, size = 2.5)

temp <- ggplot(betdist, aes(x = temp, y = distbj, color = temp)) +
  geom_boxplot(aes(color=temp)) +
  theme_light() +
  theme(legend.position = "none", axis.title.x = element_blank(),
        text = element_text(size = 10)) +
  scale_color_manual(values=enrichPal) +
  labs(y = "Beta Diversity") +
  ylim(0.5,1) + scale_x_discrete(labels=function(x) str_wrap(x, width = 10)) +
  geom_signif(comparisons = list(c("All Ambient","All High Temp")), annotation="*",
              color = "black", tip_length = 0, map_signif_level=TRUE, vjust = 0)

scar <- ggplot(betdist, aes(x = corallivory, y = distbj, color = corallivory)) +
  geom_boxplot(aes(color=corallivory)) +
  theme_light() +
  theme(legend.position = "none", axis.title.x = element_blank(), 
        axis.title.y = element_blank(), axis.text.y = element_blank(), 
        text = element_text(size = 10)) +
  scale_color_manual(values=enrichPal) +
  #labs(y = "Beta Diversity") +
  ylim(0.5,1) + scale_x_discrete(labels=function(x) str_wrap(x, width = 10)) 
# geom_signif(comparisons = list(c("All No Scar","All Scarred")), annotation="*", color = "black", tip_length = 0, map_signif_level=TRUE)

nut <- ggplot(betdist, aes(x = nutrient, y = distbj, color = nutrient)) +
  geom_boxplot(aes(color=nutrient)) +
  theme_light() +
  theme(legend.position = "none", axis.title.x = element_blank(), 
        axis.title.y = element_blank(), axis.text.y = element_blank(),
        text = element_text(size = 10)) +
  scale_color_manual(values=c('#00BFC4', '#7CAE00', '#7CAE00')) +
  #labs(y = "Beta Diversity") +
  ylim(0.5,1) + scale_x_discrete(labels=function(x) str_wrap(x, width = 5)) +
  geom_signif(comparisons = list(c("All Ambient","All NO3-")), annotation="*", color = "black", tip_length = 0, map_signif_level=TRUE)

all <- plot_grid(temp,scar, nut, nrow=1, align = 'w', rel_widths = c(1.4,1,1.4))
BJ  <- plot_grid(all,bj,nrow = 2, rel_heights = c(1,1.5))
ALL <- plot_grid(WU,UN,BC,BJ, nrow = 2)

temp <- ggplot(betdist, aes(x = temp, y = distwu, color = temp)) +
  geom_boxplot(aes(color=temp)) +
  theme_light() +
  theme(legend.position = "none", axis.title.x = element_blank(),
        text = element_text(size = 10)) +
  scale_color_manual(values=enrichPal) +
  labs(y = "Beta Diversity") +
  ylim(0,1) + scale_x_discrete(labels=function(x) str_wrap(x, width = 10)) +
  geom_signif(comparisons = list(c("All Ambient","All High Temp")), annotation="*", color = "black", tip_length = 0, map_signif_level=TRUE)

scar <- ggplot(betdist, aes(x = corallivory, y = distwu, color = corallivory)) +
  geom_boxplot(aes(color=corallivory)) +
  theme_light() +
  theme(legend.position = "none", axis.title.x = element_blank(), 
        axis.title.y = element_blank(), axis.text.y = element_blank(), 
        text = element_text(size = 10)) +
  scale_color_manual(values=enrichPal) +
  #labs(y = "Beta Diversity") +
  ylim(0,1) + scale_x_discrete(labels=function(x) str_wrap(x, width = 10))
#geom_signif(comparisons = list(c("Control","high")), annotation="*", color = "black", tip_length = 0, map_signif_level=TRUE)

nut <- ggplot(betdist, aes(x = nutrient, y = distwu, color = nutrient)) +
  geom_boxplot(aes(color=nutrient)) +
  theme_light() +
  theme(legend.position = "none", axis.title.x = element_blank(), 
        axis.title.y = element_blank(), axis.text.y = element_blank(),
        text = element_text(size = 10)) +
  scale_color_manual(values=c('#00BFC4', '#7CAE00', '#7CAE00')) +
  #labs(y = "Beta Diversity") +
  ylim(0,1) + scale_x_discrete(labels=function(x) str_wrap(x, width = 5))
#geom_signif(comparisons = list(c("Control","high")), annotation="*", color = "black", tip_length = 0, map_signif_level=TRUE)

all <- plot_grid(temp,scar, nut, nrow=1, align = 'w', rel_widths = c(1.4,1,1.4))
WU <- plot_grid(all,wu,nrow = 2, rel_heights = c(1,1.5))

