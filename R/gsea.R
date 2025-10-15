source("renv/activate.R")
source(snakemake@input[["functions"]])

suppressPackageStartupMessages({
  library(tidyverse)
  library(clusterProfiler)
  library(org.Hs.eg.db)
})


# Enrichment analysis for each fiber type and intervention ---------------
print_section_header("Running GSEAs")
gsea_results <- list()

results_list <- list("res_i" = readRDS(snakemake@input[["results_res_i"]]),
                     "res_ii" = readRDS(snakemake@input[["results_res_ii"]]),
                    "ter_i" = readRDS(snakemake@input[["results_ter_i"]]),
                     "ter_ii" = readRDS(snakemake@input[["results_ter_ii"]])
                    )

for (comparison in names(results_list)){
  
  #Extract result
  res <- results_list[[comparison]]
  
  #Preparation of ranked protein list
  gsea_list <- as.numeric(res$logFC)
  names(gsea_list) = as.character(res$protein)
  gsea_list <- gsea_list[!is.na(gsea_list)]
  
  #GSEA analysis (GO:BP)
  gsea <- clusterProfiler::gseGO(
    geneList = gsea_list,
    OrgDb = org.Hs.eg.db,
    ont = "BP",
    pvalueCutoff = 0.05,
    keyType = "SYMBOL",
    eps=0,
    minGSSize=10
    )
  
  #Save to gsea_results
  gsea_results[[comparison]] <- gsea
}


# Extract and save to data folder ----------------------------------------
gsea_res_i <- gsea_results[["res_i"]]
gsea_res_ii <- gsea_results[["res_ii"]]
gsea_ter_i <- gsea_results[["ter_i"]]
gsea_ter_ii <- gsea_results[["ter_ii"]]

saveRDS(gsea_res_i, snakemake@output[["gsea_res_i"]])
saveRDS(gsea_res_ii, snakemake@output[["gsea_res_ii"]])
saveRDS(gsea_ter_i, snakemake@output[["gsea_ter_i"]])
saveRDS(gsea_ter_ii, snakemake@output[["gsea_ter_ii"]])

#Merge results
gsea_all <- clusterProfiler::merge_result(list("res_i" = gsea_res_i,
                              "res_ii" = gsea_res_ii,
                              "ter_i" = gsea_ter_i,
                              "ter_ii" = gsea_ter_ii))

#Save
saveRDS(gsea_all, snakemake@output[["gsea_all"]])

print_section_header("Running GSEA on interaction results")

#GSEAs on intervention x time comparison
#For fiber type independent, type 1, and type 2 comparisons
gsea_results <- list()
gsea_results_filtered <- list()

results_list <- list(
  "i_and_ii_interaction" = readRDS(snakemake@input[["results_i_and_ii_interaction"]]),
  "i_interaction" = readRDS(snakemake@input[["results_i_interaction"]]),
  "ii_interaction" = readRDS(snakemake@input[["results_ii_interaction"]])
)

for (comparison in names(results_list)){
  
  #Extract result
  res <- results_list[[comparison]]
  
  #Preparation of ranked protein list
  gsea_list <- as.numeric(res$logFC)
  names(gsea_list) = as.character(res$protein)
  gsea_list <- gsea_list[!is.na(gsea_list)]
  
  #GSEA analysis (GO:BP)
  gsea <- clusterProfiler::gseGO(
    geneList = gsea_list,
    OrgDb = org.Hs.eg.db,
    ont = "ALL",
    pvalueCutoff = 0.05,
    keyType = "SYMBOL",
    eps=0,
    minGSSize=10
    )
  
  #Save to gsea_results
  gsea_results[[comparison]] <- gsea
}

#Extract and save to data folder
gsea_i_and_ii_interaction <- gsea_results[["i_and_ii_interaction"]]
gsea_i_interaction <- gsea_results[["i_interaction"]]
gsea_ii_interaction <- gsea_results[["ii_interaction"]]

saveRDS(gsea_i_and_ii_interaction, snakemake@output[["gsea_i_and_ii_interaction"]])
saveRDS(gsea_i_interaction, snakemake@output[["gsea_i_interaction"]])
saveRDS(gsea_ii_interaction, snakemake@output[["gsea_ii_interaction"]])


#Export to supplementary files ------------------------------------------
export_list <- list(
  "Within-group (Fig.4f)" = gsea_all@compareClusterResult,
  "Interaction (I and II, Fig.5b)" = gsea_i_and_ii_interaction@result,
  "Interaction (type I, Fig.5c)" = gsea_i_interaction@result,
  "Interaction (type II, Fig.5d)" = gsea_ii_interaction@result
)

writexl::write_xlsx(export_list, snakemake@output[["supplementary_gsea"]])