########################################################
###           DESeq2 for manuscript               ###
########################################################

rm(list=ls())


library(phyloseq); packageVersion("phyloseq")
library("DESeq2"); packageVersion("DESeq2")
library(dplyr)
library(data.table)
library("BiocParallel")
library("purrr")
library("tibble")

# Functions
# To summarize otu table by a higher taxonomic rank
summarize_taxa = function(physeq, Rank, GroupBy = NULL){
  Rank <- Rank[1]
  if(!Rank %in% rank_names(physeq)){
    message("The argument to `Rank` was:\n", Rank,
            "\nBut it was not found among taxonomic ranks:\n",
            paste0(rank_names(physeq), collapse = ", "), "\n",
            "Please check the list shown above and try again.")
  }
  if(!is.null(GroupBy)){
    GroupBy <- GroupBy[1]
    if(!GroupBy %in% sample_variables(physeq)){
      message("The argument to `GroupBy` was:\n", GroupBy,
              "\nBut it was not found among sample variables:\n",
              paste0(sample_variables(physeq), collapse = ", "), "\n",
              "Please check the list shown above and try again.")
    }
  }
  # Start with fast melt
  mdt = fast_melt(physeq)
  if(!is.null(GroupBy)){
    # Add the variable indicated in `GroupBy`, if provided.
    sdt = data.table(SampleID = sample_names(physeq),
                     var1 = get_variable(physeq, GroupBy))
    setnames(sdt, "var1", GroupBy)
    # Join
    setkey(sdt, SampleID)
    setkey(mdt, SampleID)
    mdt <- sdt[mdt]
  }
  # Summarize
  summarydt = mdt[, list(meanRA = mean(RelativeAbundance),
                         sdRA = sd(RelativeAbundance),
                         minRA = min(RelativeAbundance),
                         maxRA = max(RelativeAbundance)),
                  by = c(Rank, GroupBy)]
  return(summarydt)
}

phyloseq_summarize_taxa <- function(psdata, taxonomic.ranks =
                                      rank_names(psdata)) {
  if(length(taxonomic.ranks) > 1) {
    names(taxonomic.ranks) <- taxonomic.ranks
    llply(taxonomic.ranks, phyloseq_summarize_taxa, psdata = psdata)
  } else {
    taxa <- as(tax_table(psdata)[, taxonomic.ranks], 'character')
    sum_tax_table <- summarize_taxa(as(otu_table(psdata), 'matrix'), taxa)
    phyloseq(otu_table(sum_tax_table, taxa_are_rows = TRUE),
             sample_data(psdata, FALSE))
  }
}

# To melt the data into long form
fast_melt = function(physeq){
  # supports "naked" otu_table as `physeq` input.
  otutab = as(otu_table(physeq), "matrix")
  if(!taxa_are_rows(physeq)){otutab <- t(otutab)}
  otudt = data.table(otutab, keep.rownames = TRUE)
  setnames(otudt, "rn", "taxaID")
  # Enforce character taxaID key
  otudt[, taxaIDchar := as.character(taxaID)]
  otudt[, taxaID := NULL]
  setnames(otudt, "taxaIDchar", "taxaID")
  # Melt count table
  mdt = melt.data.table(otudt, 
                        id.vars = "taxaID",
                        variable.name = "SampleID",
                        value.name = "count")
  # Remove zeroes, NAs
  mdt <- mdt[count > 0][!is.na(count)]
  # Calculate relative abundance
  mdt[, RelativeAbundance := count / sum(count), by = SampleID]
  if(!is.null(tax_table(physeq, errorIfNULL = FALSE))){
    # If there is a tax_table, join with it. Otherwise, skip this join.
    taxdt = data.table(as(tax_table(physeq, errorIfNULL = TRUE), "matrix"), keep.rownames = TRUE)
    setnames(taxdt, "rn", "taxaID")
    # Enforce character taxaID key
    taxdt[, taxaIDchar := as.character(taxaID)]
    taxdt[, taxaID := NULL]
    setnames(taxdt, "taxaIDchar", "taxaID")
    # Join with tax table
    setkey(taxdt, "taxaID")
    setkey(mdt, "taxaID")
    mdt <- taxdt[mdt]
  }
  return(mdt)
}

# Trouble shooting taxonomy table
validTaxonomyTable <- function(object){
  # Both dimensions must have non-zero length.
  if( any(dim(object)==0) ){
    return("\n Taxonomy Table must have non-zero dimensions.")
  }
  # Verify that it is character matrix
  if( !is.character(object@.Data[, 1]) ){
    text = "\n Non-character matrix provided as Taxonomy Table.\n"
    text = paste0(text, "Taxonomy is expected to be characters.")
    return(text)
  }
  return(TRUE)
}
## assign the function as the validity method for the sample_data class
setValidity("taxonomyTable", validTaxonomyTable)


# First load file = "/Users/Becca/Box Sync/CHOMPIN/Manuscript_mal16/R-info/qd_unrarefied.RData
qdt = fast_melt(qd)
prevdt = qdt[, list(Prevalence = sum(count > 0), 
                    TotalCounts = sum(count)),
             by = taxaID]
keepTaxa = prevdt[(Prevalence >=15 & TotalCounts >10), taxaID]
qd = prune_taxa(keepTaxa,qd)

# Summarize to genus
tax_table(qd_fam)[is.na(tax_table(qd_fam))] <- "NA"

des <- phyloseq_to_deseq2(qd_fam, ~1)
design(des) <- ~temp*nutrient*corallivory
cts = counts(des)
geoMeans = apply(cts, 1, function(row) if (all(row == 0)) 0 else exp(sum(log(row[row != 0]))))
des = estimateSizeFactors(des, geoMeans=geoMeans)
dess <- DESeq(des, betaPrior = F)

cts = counts(des)
geoMeans = apply(cts, 1, function(row) if (all(row == 0)) 0 else exp(sum(log(row[row != 0]))))
# Calculate geometric means prior to estimate size factors
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans = apply(counts(des), 1, gm_mean)

des = estimateSizeFactors(des, type = "iterate", geoMeans=geoMeans)


des_gen <- DESeq(des,betaPrior = F)

dess <- DESeq(des,betaPrior=F,full=mm1)


full_formula <- ~temp*nutrient*corallivory
des <- model_data(qd, full_formula)

dds_int <- phyloseq_to_deseq2(qd, ~temp*nutrient*corallivory)
# When I ran DESeq on this object, I had a problem with some rows not converging.
# So I used the following code to get from 6 to 1 non-converging rows
dds_int <- estimateSizeFactors(dds_int)
nc <- counts(dds_int, normalized=TRUE)
filter <- rowSums(nc >= 10) >= 3
dds_int <- dds_int[filter,]

dds <- DESeq(des, test = "Wald", fitType = "parametric", betaPrior = TRUE)
# Then I removed the one non-converging row because these are typically OTUs 
# with very small counts and little power. I also tried increasing the 
# iterations in the nbinomWaldTest but that didn't help.
dds_test_int_clean <- dds_test_int[which(mcols(dds_test_int)$betaConv),]

res_int = results(dds_test_int_clean, cooksCutoff = FALSE)
res_int = res_int[which(res_int$padj < alpha), ]
res_int_high_scarred = cbind(as(res_int, "data.frame"), as(tax_table(qd)[rownames(res_int), ], "matrix"))

alpha = 0.05
res1 <- results(dess, name = "temp_high_vs_Control", cooksCutoff = FALSE)
res1 = res1[which(res1$padj < alpha),]
res1 = cbind(as(res1, "data.frame"))
res1$comparison <- "Control vs. high"
res2 <- results(dess, name = "corallivory_Scarred_vs_Control", cooksCutoff = FALSE)
res2 = res2[which(res2$padj < alpha),]
res2 = cbind(as(res2, "data.frame"))
res2$comparison <- "Control vs. scarred"
res3 <- results(dess, name = "nutrient_NO3._vs_Control", cooksCutoff = FALSE)
res3 = res3[which(res3$padj < alpha),]
res3 = cbind(as(res3, "data.frame"))
res3$comparison <- "Control vs. NO3"
res4 <- results(dess, name = "nutrient_NH4._vs_Control", cooksCutoff = FALSE)
res4 = res4[which(res4$padj < alpha),]
res4 = cbind(as(res4, "data.frame"))
res4$comparison <- "Control vs. NH4"
res5 <- results(dess, name = "temphigh.corallivoryScarred", cooksCutoff = FALSE)
res5 = res5[which(res5$padj < alpha),]
res5 = cbind(as(res5, "data.frame"))
res5$comparison <- "High, scarred"
res6 <- results(dess, name = "nutrientNO3..corallivoryScarred", cooksCutoff = FALSE)
res6 = res6[which(res6$padj < alpha),]
res6 = cbind(as(res6, "data.frame"))
res6$comparison <- "NO3, scarred"
res7 <- results(dess, name = "nutrientNH4..corallivoryScarred", cooksCutoff = FALSE)
res7 = res7[which(res7$padj < alpha),]
res7 = cbind(as(res7, "data.frame"))
res7$comparison <- "NH4, scarred"
res8 <- results(dess, name = "temphigh.nutrientNO3.", cooksCutoff = FALSE)
res8 = res8[which(res8$padj < alpha),]
res8 = cbind(as(res8, "data.frame"))
res8$comparison <- "NO3, high"
res9 <- results(dess, name = "temphigh.nutrientNH4.", cooksCutoff = FALSE)
res9 = res9[which(res9$padj < alpha),]
res9 = cbind(as(res9, "data.frame"))
res9$comparison <- "NH4, high"
res10 <- results(dess, name = "temphigh.nutrientNO3..corallivoryScarred", cooksCutoff = FALSE)
res10 = res10[which(res10$padj < alpha),]
res10 = cbind(as(res10, "data.frame"))
res10$comparison <- "NO3, high, scarred"
res11 <- results(dess, name = "temphigh.nutrientNH4..corallivoryScarred", cooksCutoff = FALSE)
res11 = res11[which(res11$padj < alpha),]
res11 = cbind(as(res11, "data.frame"))
res11$comparison <- "NH4, high, scarred"

taxa <- read.csv(file="/Users/Becca/Box Sync/CHOMPIN/Manuscript_mal16/R-info/taxa_unrar.csv")
rownames(taxa) <- taxa[,1]
colnames(taxa)[1] <- "OTUID"

library(tidyverse)
library("pandas")
res <- list(res1,res2,res3,res4,res5,res6,res7,res8,res9,res10,res11) %>% map_df(rownames_to_column, 'OTUID')
resh <- merge(res,taxa)
write.csv(resh, file = "/Users/Becca/Box Sync/CHOMPIN/Manuscript_mal16/R-info/deseq_family_feb20.csv")

###################################################################
## Functions
summarize_taxa <- function(counts, taxonomy) {
  if(is.matrix(taxonomy)) {
    #message('multiple taxonomies')
    alply(taxonomy, 2, summarize_taxa, counts = counts, .dims = TRUE)
  } else if(is.matrix(counts)) {
    #message('multiple counts')
    require('plyr')
    apply(counts, 2, summarize_taxa, taxonomy = taxonomy)
  } else {
    #message('summarize')
    tapply(counts, taxonomy, sum)
  }
}

phyloseq_summarize_taxa <- function(psdata, taxonomic.ranks =
                                      rank_names(psdata)) {
  if(length(taxonomic.ranks) > 1) {
    names(taxonomic.ranks) <- taxonomic.ranks
    llply(taxonomic.ranks, phyloseq_summarize_taxa, psdata = psdata)
  } else {
    taxa <- as(tax_table(psdata)[, taxonomic.ranks], 'character')
    sum_tax_table <- summarize_taxa(as(otu_table(psdata), 'matrix'), taxa)
    phyloseq(otu_table(sum_tax_table, taxa_are_rows = TRUE),
             sample_data(psdata, FALSE))
  }
}

# Function fast_melt
fast_melt = function(physeq,
                     includeSampleVars = character(),
                     omitZero = FALSE){
  require("phyloseq")
  require("data.table")
  # supports "naked" otu_table as `physeq` input.
  otutab = as(otu_table(physeq), "matrix")
  if(!taxa_are_rows(physeq)){otutab <- t(otutab)}
  otudt = data.table(otutab, keep.rownames = TRUE)
  setnames(otudt, "rn", "TaxaID")
  # Enforce character TaxaID key
  otudt[, TaxaIDchar := as.character(TaxaID)]
  otudt[, TaxaID := NULL]
  setnames(otudt, "TaxaIDchar", "TaxaID")
  # Melt count table
  mdt = melt.data.table(otudt, 
                        id.vars = "TaxaID",
                        variable.name = "SampleID",
                        value.name = "count")
  if(omitZero){
    # Omit zeroes and negative numbers
    mdt <- mdt[count > 0]
  }
  # Omit NAs
  mdt <- mdt[!is.na(count)]
  # Calculate relative abundance
  mdt[, RelativeAbundance := count / sum(count), by = SampleID]
  if(!is.null(tax_table(physeq, errorIfNULL = FALSE))){
    # If there is a tax_table, join with it. Otherwise, skip this join.
    taxdt = data.table(as(tax_table(physeq, errorIfNULL = TRUE), "matrix"), keep.rownames = TRUE)
    setnames(taxdt, "rn", "TaxaID")
    # Enforce character TaxaID key
    taxdt[, TaxaIDchar := as.character(TaxaID)]
    taxdt[, TaxaID := NULL]
    setnames(taxdt, "TaxaIDchar", "TaxaID")
    # Join with tax table
    setkey(taxdt, "TaxaID")
    setkey(mdt, "TaxaID")
    mdt <- taxdt[mdt]
  }
  # includeSampleVars = c("DaysSinceExperimentStart", "SampleType")
  # includeSampleVars = character()
  # includeSampleVars = c()
  # includeSampleVars = c("aksjdflkas") 
  wh.svars = which(sample_variables(physeq) %in% includeSampleVars)
  if( length(wh.svars) > 0 ){
    # Only attempt to include sample variables if there is at least one present in object
    sdf = as(sample_data(physeq), "data.frame")[, wh.svars, drop = FALSE]
    sdt = data.table(sdf, keep.rownames = TRUE)
    setnames(sdt, "rn", "SampleID")
    # Join with long table
    setkey(sdt, "SampleID")
    setkey(mdt, "SampleID")
    mdt <- sdt[mdt]
  }
  setkey(mdt, "TaxaID")
  return(mdt)
}


df <- data.frame(sample_data(qd))

m2 <- model.matrix(~ interaction*tank, df)
colnames(m2)
unname(m2)
all.zero <- apply(m2, 2, function(x) all(x==0))
all.zero

idx <- which(all.zero)
m2 <- m2[,-idx]
unname(m2)

dds <- phyloseq_to_deseq2(qd, ~interaction+tank)

dds_test <- DESeq(dds, test = c("Wald"), full = m2, betaPrior = FALSE)

## Deseq results for interaction
des <- phyloseq_to_deseq2(qd, ~1)
des <- phyloseq_to_deseq2(qd, ~interaction)
dds <- DESeq(des, test = "Wald", fitType = "parametric", betaPrior = TRUE)
# Stopped with this analysis because the results were wonkey
alpha = 0.05
res1 <- results(dds, contrast = c("interaction", "High", "Control"))
res1 = res1[which(res1$padj < alpha),]
res1 = cbind(as(res1, "data.frame"))
res1$comparison <- "Control vs. high"
res2 <- results(dds, contrast = c("interaction", "Scarred", "Control"))
res2 = res2[which(res2$padj < alpha),]
res2 = cbind(as(res2, "data.frame"))
res2$comparison <- "Control vs. scarred"
res3 <- results(dds, contrast = c("interaction", "NO3.", "Control"))
res3 = res3[which(res3$padj < alpha),]
res3 = cbind(as(res3, "data.frame"))
res3$comparison <- "Control vs. NO3"
res4 <- results(dds, contrast = c("interaction", "NH4.", "Control"))
res4 = res4[which(res4$padj < alpha),]
res4 = cbind(as(res4, "data.frame"))
res4$comparison <- "Control vs. NH4"
res5 <- results(dds, name = "interactionHigh..scarred")
res5 = res5[which(res5$padj < alpha),]
res5 = cbind(as(res5, "data.frame"))
res5$comparison <- "High, scarred"
res6 <- results(dds_test_int_clean, name = "nutrientNO3..corallivoryScarred", cooksCutoff = FALSE)
res6 = res6[which(res6$padj < alpha),]
res6 = cbind(as(res6, "data.frame"))
res6$comparison <- "NO3, scarred"
res7 <- results(dds_test_int_clean, name = "nutrientNH4..corallivoryScarred", cooksCutoff = FALSE)
res7 = res7[which(res7$padj < alpha),]
res7 = cbind(as(res7, "data.frame"))
res7$comparison <- "NH4, scarred"
res8 <- results(dds_test_int_clean, name = "temphigh.nutrientNO3.", cooksCutoff = FALSE)
res8 = res8[which(res8$padj < alpha),]
res8 = cbind(as(res8, "data.frame"))
res8$comparison <- "NO3, high"
res9 <- results(dds_test_int_clean, name = "temphigh.nutrientNH4.", cooksCutoff = FALSE)
res9 = res9[which(res9$padj < alpha),]
res9 = cbind(as(res9, "data.frame"))
res9$comparison <- "NH4, high"
res10 <- results(dds_test_int_clean, name = "temphigh.nutrientNO3..corallivoryScarred", cooksCutoff = FALSE)
res10 = res10[which(res10$padj < alpha),]
res10 = cbind(as(res10, "data.frame"))
res10$comparison <- "NO3, high, scarred"
res11 <- results(dess, name = "temphigh.nutrientNH4..corallivoryScarred", cooksCutoff = FALSE)
res11 = res11[which(res11$padj < alpha),]
res11 = cbind(as(res11, "data.frame"))
res11$comparison <- "NH4, high, scarred"
