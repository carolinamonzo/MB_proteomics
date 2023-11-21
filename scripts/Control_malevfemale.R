---
  title: "analysis_for_Maarouf"
output: html_document
date: "2023-10-27"
author: "Carolina Monzo"
---
  
  Read in the files and format to remove negative values and mice that failed.
```{r}
library(dplyr)
library(stringr)
#library(mixOmics)
#library(impute)
library(tidyr)
#library(multcomp)
#library(clusterProfiler)
#library(enrichplot)
#library(ggupset)
#library(pathview)
#library(ggnewscale)

# Format removing the sequencing information from Olink
first <- read.csv("../data/Q-03195_Kroef_NPX.csv", sep = ";")
first <- first %>% 
  mutate(across(where(is.character), str_trim))
first <- first[-(1:2), , drop = FALSE]
names(first) <- as.character(unlist(first[1,]))
first <- first[-(1:3), , drop = FALSE]
first <- head(first,-4)
rownames(first) <- first$Assay
first$Assay <- NULL


# Read into a proper dataframe and remove negative values and mice that failed
df <- as.data.frame(first)
df$`QC Warning` <- NULL
df$`Plate ID` <- NULL
df <- mutate_all(df, function(x) as.numeric(as.character(as.numeric(gsub(",", ".", gsub("\\.", "", x))))))
# Change negative values into NA
df <- df %>% dplyr::mutate(across(everything(), function(x){replace(x, which(x<0), NA)}))
# Remove samples Olink selected as problematic
df <- df[!(row.names(df) %in% c("TNES-AC32", "TNES-AC86", "TNES-Y512", "TNES-AE58", "TNES-AH43")),]

## Read in the metadata
metadata <- read.csv("../metadata/metadata_Maarouf.csv", sep = ";")
metadata <- metadata[!(metadata$ID %in% c("TNES-AC32", "TNES-AC86", "TNES-Y512", "TNES-AE58", "TNES-AH43")),]
# Sort the rows in df according to metadata
df <- df[match(metadata$ID, rownames(df)),]
df <- df[(row.names(df) %in% metadata$ID), ]


# Extract only controls
metadata.controls <- metadata[metadata$Treatment == 'Control',]
df.control <- df[(row.names(df) %in% metadata.controls$ID), ]


get_stats <- function(combined_df){
  # Create a list to store results for each gene
  gene_results <- list()
  
  # Loop through each gene (assuming the gene names are the column names)
  for (gene_name in (2:length(colnames(combined_df))-1)) {
    print(colnames(combined_df)[gene_name])
    # Perform one-way ANOVA for "Continuous" vs "Control"
    anova_df <- aov(formula(paste(colnames(combined_df)[gene_name], " ~ Treatment")), 
                    data = combined_df)
    
    # Perform Tukey post-hoc tests for "Continuous" vs "Control"
    tukey_df <- TukeyHSD(anova_df)
    tuk <- as.data.frame(tukey_df$Treatment)
    
    
    # Calculate the fold change for each gene
    fold_change_continuous <- mean(combined_df[combined_df$Treatment == "f", colnames(combined_df)[gene_name]]) / mean(combined_df[combined_df$Treatment == "m", colnames(combined_df)[gene_name]])
    #fold_change_intermittent <- mean(combined_df[combined_df$Treatment == "Intermittent", gene_name]) / mean(combined_df[combined_df$Treatment == "Control", gene_name])
    
    # Store the results in a list
    gene_results[[colnames(combined_df)[gene_name]]] <- list(
      ANOVA = summary(anova_df)[[1]][[5]][1],
      Tukey_Continuous = tuk[row.names(tuk) %in% "F-M",]$`p adj`,
      #Tukey_Intermittent = tuk[row.names(tuk) %in% "Intermittent-Control",]$`p adj`,
      FoldChange_Continuous = fold_change_continuous
      #FoldChange_Intermittent = fold_change_intermittent
    )
  }
  
  gene_results <- do.call(rbind, lapply(gene_results, unlist))
  
  return(gene_results)
}





combined_df <- df.control
combined_df$Treatment <- metadata.controls$Sex
stats.control <- get_stats(combined_df)
stats.control <- as.data.frame(stats.control)
write.csv(stats.males, "../analysis/stats_femaleMale_controls_maarouf.csv")

View(stats.control)

vec_subset <- c("Ahr","Casp3","Ccl2","Ccl20","Ccl3","Ccl5","Csf2","Cxcl1","Cxcl9","Cyr61","Dll1","Eda2r","Epo","Erbb4","Fas"
                ,"Fstl3","Gcg","Ghrl","Il10","Il17a","Il17f","Il1a","Il1b","Il23r","Il5","Il6","Itgb6","Kitlg","Lgmn","Map2k6"
                ,"Parp1","Pla2g4a","Plin1","Prdx5","S100a4","Tgfb1","Tgfbr3","Tnf","Tnfrsf11b","Tnfrsf12a","Tnfsf12","Yes1")

immgenes <- stats.control[row.names(stats.control) %in% vec_subset, ]
write.csv(immgenes, "../analysis/stats_femaleMale_controls_Immune_maarouf.csv")
