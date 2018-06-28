#################################################################################
# This script is the beta diversity analysis for the CHOMPIN Project:
# Here, I calculate ...
#
# Created by Rebecca Maher
# Created on 6/20/2018
# Edited on ...
#################################################################################

## clear workspace------------------------
rm(list=ls())

# load libraries
library('phyloseq'); packageVersion('phyloseq')
library('vegan'); packageVersion('vegan')

# set working directory-------------------
#setwd(~'CHOMPIN Ms/scripts')

## functions (if any)----------------------

## Data Analysis----------------------------

# First import qd_rarefied.RData object, this is the Phyloseq object with a rarefied OTU table

# Calculate all four distance matrices: weighted unifrac, unweighted unifrac, bray curtis, binary jaccard
qd_wu <- phyloseq::distance(qd, method = "wunifrac")
qd_un <- phyloseq::distance(qd, method = "unifrac")
qd_bc <- phyloseq::distance(qd, method = "bray")
qd_bj <- distance(qd, method = "jaccard", binary =TRUE)

# PERMANOVA's with Adonis -  Permutational Multivariate Analasis of Variance
# Make a data frame from the sample_data
sampledf <- data.frame(sample_data(qd))

# Adonis test for between group diversity, with full formula
adonis(qd_bc ~ temp*corallivory*nutrient, data = sampledf)
adonis(qd_bj ~ temp*nutrient*corallivory, data = sampledf) 
adonis(qd_wu ~ temp*nutrient*corallivory, data = sampledf) 
adonis(qd_un ~ temp*nutrient*corallivory, data = sampledf) 

# To test if the above differences betwee groups are due to differences in
# location (of group centroid) or dispersion (average distance to centroid)
# need to use PERMDISP. PERMDISP does not take a formula, so the
# treatment interaction is assessed as a categorical variable

# Adonis test for between group diversity by treatment
adonis(qd_bc ~ interaction, data = sampledf)
adonis(qd_bj ~ interaction, data = sampledf)
adonis(qd_wu ~ interaction, data = sampledf)
adonis(qd_un ~ interaction, data = sampledf)

# PERMDISP with betadisper - Multivariate Homogeneity of group dispersions
anova(betadisper(qd_bc, sampledf$interaction, bias.adjust = TRUE))
anova(betadisper(qd_bj, sampledf$interaction, bias.adjust = TRUE))
anova(betadisper(qd_wu, sampledf$interaction, bias.adjust = TRUE))
anova(betadisper(qd_un, sampledf$interaction, bias.adjust = TRUE))

# Add the pairwise=TRUE to get a pairwise comparision of all the results
permutest(betadisper(qd_bc, sampledf$interaction, bias.adjust = TRUE), pairwise=TRUE, nperm = 999)
permutest(betadisper(qd_bj, sampledf$interaction, bias.adjust = TRUE), pairwise=TRUE, nperm = 999)
permutest(betadisper(qd_wu, sampledf$interaction, bias.adjust = TRUE), pairwise=TRUE, nperm = 999)
permutest(betadisper(qd_un, sampledf$interaction, bias.adjust = TRUE), pairwise=TRUE, nperm = 999)
