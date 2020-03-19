# ==============================================================================
# Author(s) : Heini M Natri, hnatri@asu.edu
# Date : July 2019
# Description: Integrating the D. melanogaste ATACseq and RNAseq data
# ==============================================================================

# ======================================
# Load libraries
# ======================================

library(ggplot2)

# ======================================
# Environment parameters
# ======================================

setwd("/Users/hnatri/Dropbox (ASU)/ICGC_HCC/")

# ======================================
# Reading data
# ======================================

metadata <- read.csv("metadata_for_de.csv", header = T)

metadata$sex_tissue <- paste(metadata$Gender, metadata$tissue, sep = "_") 
metadata$sex_tissue_infection <- paste(metadata$sex_tissue, metadata$Virus.infection, sep = "_")
rownames(metadata) <- metadata$sampleid

metadata$Edmondson.grade <- factor(metadata$Edmondson.grade)

# Adding in cell types
celltypes <- t(read.table("xCell_ICGC_HCC_FPKM_xCell_1911010120.txt", sep="\t", header=T, row.names=1, as.is=T))
rownames(celltypes) <- gsub("\\.", "\\-", rownames(celltypes))
metadata_celltypes <- merge(metadata, celltypes, by=0)

# Only keeping HCC cases
metadata_celltypes <- metadata_celltypes[metadata_celltypes$HistologicalType == "HCC" ,]

metadata_celltypes_nbnc <- metadata_celltypes[metadata_celltypes$Virus.infection == "NBNC" , ]
metadata_celltypes_hbv <- metadata_celltypes[metadata_celltypes$Virus.infection == "HBV" , ]
metadata_celltypes_hcv <- metadata_celltypes[metadata_celltypes$Virus.infection == "HCV" , ]
metadata_celltypes_both <- metadata_celltypes[metadata_celltypes$Virus.infection == "HBV, HCV" , ]

metadata_celltypes_nbnc_adjacent <- metadata_celltypes_nbnc[metadata_celltypes_nbnc$tissue == "adjacent" , ]
metadata_celltypes_hbv_adjacent <- metadata_celltypes_hbv[metadata_celltypes_hbv$tissue == "adjacent" , ]
metadata_celltypes_hcv_adjacent <- metadata_celltypes_hcv[metadata_celltypes_hcv$tissue == "adjacent" , ]
metadata_celltypes_both_adjacent <- metadata_celltypes_both[metadata_celltypes_both$tissue == "adjacent" , ]

metadata_celltypes_nbnc_tumor <- metadata_celltypes_nbnc[metadata_celltypes_nbnc$tissue == "tumor" , ]
metadata_celltypes_hbv_tumor <- metadata_celltypes_hbv[metadata_celltypes_hbv$tissue == "tumor" , ]
metadata_celltypes_hcv_tumor <- metadata_celltypes_hcv[metadata_celltypes_hcv$tissue == "tumor" , ]
metadata_celltypes_both_tumor <- metadata_celltypes_both[metadata_celltypes_both$tissue == "tumor" , ]

# ======================================
# t-test for different traists within groups
# ======================================

# For continuous variables
testvars <- c(colnames(celltypes), c("Age", "Tumor.size.mm", "Overall.survival.month"))
#colnames(celltypes), 
# Running t-test on all traits
p_vals <- data.frame(matrix(NA, ncol=1, nrow=length(testvars)))
colnames(p_vals) <- "pval"
rownames(p_vals) <- testvars
for (i in testvars) {
  
  formula <- as.formula(paste0(i, " ~ Gender"))
  
  ttest <- t.test(formula, metadata_celltypes_hcv_tumor)
  summary(fit)
  pval <- ttest$p.value
  
  p_vals[i, "pval"] = pval
}

#p_vals
p_vals$adjp <- p.adjust(p_vals$pval, n=length(rownames(p_vals)), method="fdr")
p_vals_sig <- p_vals[p_vals$adjp<0.1 , ]
p_vals_sig

p_vals

# Writing to a file
write.csv(p_vals, "celltypes_ttest_mvsf_hcv_tumor.csv")
# ======================================
# 2-way chi-square test for categorical variables
# ======================================

testvars <- c("T", "Portal.vein.invasion", "Hepatic.vein.invasion", "Bile.duct.invasion", "Liver.fibrosisc", "Alcohol.intake", "Smoking", "Prognosis")

# Running chi-square on all traits
p_vals <- data.frame(matrix(NA, ncol=1, nrow=length(testvars)))
colnames(p_vals) <- "pval"
rownames(p_vals) <- testvars
for (i in testvars) {
  
  formula <- as.formula(paste0(i, " ~ Gender"))
  
  fit <- chisq.test(xtabs(formula, metadata_celltypes_hcv_tumor))
  
  pval <- fit$p.value
  
  p_vals[i, "pval"] = pval
}

p_vals

