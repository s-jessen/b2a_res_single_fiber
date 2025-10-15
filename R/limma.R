source("renv/activate.R")
source(snakemake@input[["functions"]])

suppressPackageStartupMessages({
  library(limma)
  library(SummarizedExperiment)
  library(tidyverse)
  library(qvalue)
  library(writexl)
})

se <- readRDS(snakemake@input[["se"]])

# Main effects, within-group effects, and interaction effects independent of group ----
print_section_header("Running Limma main effect analyses")

#Combine factors into a single term for fiber_type, time, and intervention
factors <- paste(se$time, se$fiber_type, sep=".")

#Set factor levels
factors <- factor(factors, levels=c("pre.I",
                                    "pre.II",
                                    "post.I",
                                    "post.II"))

#Create design matrix
design <- model.matrix(~0 + factors)

#Rename columns for easier interpretation
colnames(design) <- c("pre.I",
                      "pre.II",
                      "post.I",
                      "post.II")

#Estimate correlation
correlation <- limma::duplicateCorrelation(SummarizedExperiment::assay(se), design, block=se$id)

#Fit the linear model
fit <- limma::lmFit(SummarizedExperiment::assay(se),
                    design, 
                    block=se$id,
                    correlation=correlation$consensus)

#Set up contrasts
contrasts <- limma::makeContrasts(
  main = (post.I + post.II)/2 - (pre.I + pre.II)/2,
  type1 = post.I - pre.I,
  type2 = post.II - pre.II,
  interaction = (post.II - pre.II) - (post.I - pre.I),
  levels = design)

#Apply contrasts to the fit
fit_contrasts <- limma::contrasts.fit(fit, contrasts)

#Apply eBayes for moderated statistics
fit_final <- limma::eBayes(fit_contrasts)

#Loop over each contrast
contrast_list <- c(
  "main",
  "type1",
  "type2",
  "interaction"
)

results_list_res = list()

for (contrast in contrast_list) {

  results <- limma::topTable(fit_final, coef = contrast, number = Inf) %>% 
    dplyr::mutate(xiao = P.Value^abs(logFC),
                  protein = row.names(.),
                  regulated = ifelse(xiao < 0.05, "+", ""),
                  `-log10(p)` = -log10(P.Value),
                  q = qvalue::qvalue(.$P.Value)$qvalues
                  )%>% 
    dplyr::arrange(desc(logFC)) %>% 
    dplyr::mutate(contrast = contrast)
  
  results_list_res[[contrast]] <- results

}

#Extract results
results_main <- results_list_res[["main"]]
results_i <- results_list_res[["type1"]]
results_ii <- results_list_res[["type2"]]
results_interaction <- results_list_res[["interaction"]]

#Save to data folder
saveRDS(results_main, snakemake@output[["results_main"]])
saveRDS(results_i, snakemake@output[["results_i"]])
saveRDS(results_ii, snakemake@output[["results_ii"]])
saveRDS(results_interaction, snakemake@output[["results_interaction"]])


# Main effects, within-group effects, and interaction effects for terbutaline-group ----
print_section_header("Running Limma within-group & interaction effect (type II vs. type I) analyses for B2A group")

#Subset SE
se_ter <- se[,se$intervention == "terbutaline"]

#Combine factors into a single term for fiber_type, time, and intervention
factors <- paste(se_ter$time, se_ter$fiber_type, sep=".")

#Set factor levels
factors <- factor(factors, levels=c("pre.I",
                                    "pre.II",
                                    "post.I",
                                    "post.II"))

#Create design matrix
design <- model.matrix(~0 + factors)

# Rename columns for easier interpretation
colnames(design) <- c("pre.I",
                      "pre.II",
                      "post.I",
                      "post.II")

#Estimate correlation
correlation <- limma::duplicateCorrelation(SummarizedExperiment::assay(se_ter), design, block=se_ter$id)

#Fit the linear model
fit <- limma::lmFit(SummarizedExperiment::assay(se_ter),
                    design, block=se_ter$id,
                    correlation=correlation$consensus)

#Set up contrasts
contrasts <- limma::makeContrasts(
  terbutaline_main = (post.I + post.II)/2 - (pre.I + pre.II)/2,
  terbutaline_type1 = post.I - pre.I,
  terbutaline_type2 = post.II - pre.II,
  terbutaline_interaction = (post.II - pre.II) - (post.I - pre.I),
  levels = design)

#Apply contrasts to the fit
fit_contrasts <- limma::contrasts.fit(fit, contrasts)

#Apply eBayes for moderated statistics
fit_final <- limma::eBayes(fit_contrasts)

#Loop over each contrast
contrast_list <- c(
  "terbutaline_main",
  "terbutaline_type1",
  "terbutaline_type2",
  "terbutaline_interaction"
)

results_list_ter = list()

for (contrast in contrast_list) {

  results <- limma::topTable(fit_final, coef = contrast, number = Inf) %>% 
    dplyr::mutate(xiao = P.Value^abs(logFC),
                  protein = row.names(.),
                  regulated = ifelse(xiao < 0.05, "+", ""),
                  `-log10(p)` = -log10(P.Value),
                  q = qvalue::qvalue(.$P.Value)$qvalues
                  )%>% 
    dplyr::arrange(desc(logFC)) %>% 
    dplyr::mutate(contrast = contrast)
  
  results_list_ter[[contrast]] <- results

}

#Extract results
results_ter_main <- results_list_ter[["terbutaline_main"]]
results_ter_i <- results_list_ter[["terbutaline_type1"]]
results_ter_ii <- results_list_ter[["terbutaline_type2"]]
results_ter_interaction <- results_list_ter[["terbutaline_interaction"]]

#Save to data folder
saveRDS(results_ter_main, snakemake@output[["results_ter_main"]])
saveRDS(results_ter_i, snakemake@output[["results_ter_i"]])
saveRDS(results_ter_ii, snakemake@output[["results_ter_ii"]])
saveRDS(results_ter_interaction, snakemake@output[["results_ter_interaction"]])


# Main effects, within-group effects, and interaction effects for resistance-group ----
print_section_header("Running Limma within-group & interaction effect (type II vs. type I) analyses for RES group")

#Subset SE
se_res <- se[,se$intervention == "resistance"]

#Combine factors into a single term for fiber_type, time, and intervention
factors <- paste(se_res$time, se_res$fiber_type, sep=".")

#Set factor levels
factors <- factor(factors, levels=c("pre.I",
                                    "pre.II",
                                    "post.I",
                                    "post.II"))

#Create design matrix
design <- model.matrix(~0 + factors)

# Rename columns for easier interpretation
colnames(design) <- c("pre.I",
                      "pre.II",
                      "post.I",
                      "post.II")

#Estimate correlation
correlation <- limma::duplicateCorrelation(SummarizedExperiment::assay(se_res), design, block=se_res$id)

#Fit the linear model
fit <- limma::lmFit(SummarizedExperiment::assay(se_res),
                    design, block=se_res$id,
                    correlation=correlation$consensus)

#Set up contrasts
contrasts <- limma::makeContrasts(
  resistance_main = (post.I + post.II)/2 - (pre.I + pre.II)/2,
  resistance_type1 = post.I - pre.I,
  resistance_type2 = post.II - pre.II,
  resistance_interaction = (post.II - pre.II) - (post.I - pre.I),
  levels = design)

#Apply contrasts to the fit
fit_contrasts <- limma::contrasts.fit(fit, contrasts)

#Apply eBayes for moderated statistics
fit_final <- limma::eBayes(fit_contrasts)

#Loop over each contrast
contrast_list <- c(
  "resistance_main",
  "resistance_type1",
  "resistance_type2",
  "resistance_interaction"
)

results_list = list()

for (contrast in contrast_list) {

  results <- limma::topTable(fit_final, coef = contrast, number = Inf) %>% 
    dplyr::mutate(xiao = P.Value^abs(logFC),
                  protein = row.names(.),
                  regulated = ifelse(xiao < 0.05, "+", ""),
                  `-log10(p)` = -log10(P.Value),
                  q = qvalue::qvalue(.$P.Value)$qvalues
                  )%>% 
    dplyr::arrange(desc(logFC)) %>% 
    dplyr::mutate(contrast = contrast)
  
  results_list[[contrast]] <- results

}

#Extract results
results_res_main <- results_list[["resistance_main"]]
results_res_i <- results_list[["resistance_type1"]]
results_res_ii <- results_list[["resistance_type2"]]
results_res_interaction <- results_list[["resistance_interaction"]]

#Save to data folder
saveRDS(results_res_main, snakemake@output[["results_res_main"]])
saveRDS(results_res_i, snakemake@output[["results_res_i"]])
saveRDS(results_res_ii, snakemake@output[["results_res_ii"]])
saveRDS(results_res_interaction, snakemake@output[["results_res_interaction"]])


#Unified results
results <- dplyr::bind_rows(results_res_i,
                            results_res_ii,
                            results_ter_i,
                            results_ter_ii) %>% 
  tibble::remove_rownames()

#Save unified results
saveRDS(results, snakemake@output[["results"]])


# Interaction effects independent of fiber type --------------------------
print_section_header("Running Limma interaction effect (RES vs B2A, independent of fiber type) analyses")

#Subset SE
se_interaction <- se

#Combine factors
TS_i_and_ii_interaction <- paste(se$time, se$intervention, sep=".")

#Set factor levels
TS_i_and_ii_interaction <- factor(TS_i_and_ii_interaction, levels=c("pre.resistance",
                                                      "post.resistance",
                                                      "pre.terbutaline",
                                                      "post.terbutaline"))
#Create design matrix
design_i_and_ii_interaction <- model.matrix(~0+ TS_i_and_ii_interaction)

#Specify column names in design matrix
colnames(design_i_and_ii_interaction)=c("pre.resistance",
                                 "post.resistance",
                                 "pre.terbutaline",
                                 "post.terbutaline")

#Estimate correlation coefficients
correlation_i_and_ii_interaction <- limma::duplicateCorrelation(
  SummarizedExperiment::assay(se),
  design_i_and_ii_interaction,
  block=se$id)

#Fit the linear model
fit_i_and_ii_interaction <- limma::lmFit(SummarizedExperiment::assay(se),
                                  design_i_and_ii_interaction,
                                  block=se$id,
                                  correlation=correlation_i_and_ii_interaction$consensus)

#Set up contrast
contrast_i_and_ii_interaction <- limma::makeContrasts(
    interaction=(post.resistance - pre.resistance)-(post.terbutaline - pre.terbutaline),
                             levels = design_i_and_ii_interaction)

#Apply contrast
fit2_i_and_ii_interaction <- limma::eBayes(limma::contrasts.fit(fit_i_and_ii_interaction, contrast_i_and_ii_interaction))

#Extract result
results_i_and_ii_interaction <- limma::topTable(fit2_i_and_ii_interaction, coef = "interaction", number = Inf) %>% 
  dplyr::mutate(xiao = P.Value^abs(logFC),
                protein = row.names(.),
                regulated = ifelse(xiao < 0.05, "+", ""),
                `-log10(p)` = -log10(P.Value),
                q = qvalue::qvalue(.$P.Value)$qvalues
                )%>% 
  dplyr::arrange(desc(logFC))

#Save results
saveRDS(results_i_and_ii_interaction, snakemake@output[["results_i_and_ii_interaction"]])


#Interaction effects for type I fibers
print_section_header("Running Limma interaction effect (RES vs B2A, type I) analyses")

#Subset SE
se_i_interaction <- se[, se$fiber_type== 'I']

#Combine factors
TS_i_interaction <- paste(se_i_interaction$time, se_i_interaction$intervention, sep=".")

#Set factor levels
TS_i_interaction <- factor(TS_i_interaction, levels=c("pre.resistance",
                                                      "post.resistance",
                                                      "pre.terbutaline",
                                                      "post.terbutaline"))
#Create design matrix
design_i_interaction <- model.matrix(~0+ TS_i_interaction)

#Specify column names in design matrix
colnames(design_i_interaction)=c("pre.resistance",
                                 "post.resistance",
                                 "pre.terbutaline",
                                 "post.terbutaline")

#Estimate correlation coefficients
correlation_i_interaction <- limma::duplicateCorrelation(
  SummarizedExperiment::assay(se_i_interaction),
  design_i_interaction,
  block=se_i_interaction$id)

#Fit the linear model
fit_i_interaction <- limma::lmFit(SummarizedExperiment::assay(se_i_interaction),
                                  design_i_interaction,
                                  block=se_i_interaction$id,
                                  correlation=correlation_i_interaction$consensus)

#Set up contrast
contrast_i_interaction <- limma::makeContrasts(
    interaction=(post.resistance - pre.resistance)-(post.terbutaline - pre.terbutaline),
                             levels = design_i_interaction)

#Apply contrast
fit2_i_interaction <- limma::eBayes(limma::contrasts.fit(fit_i_interaction, contrast_i_interaction))

#Extract result
results_i_interaction <- limma::topTable(fit2_i_interaction, coef = "interaction", number = Inf) %>% 
  dplyr::mutate(xiao = P.Value^abs(logFC),
                protein = row.names(.),
                regulated = ifelse(xiao < 0.05, "+", ""),
                `-log10(p)` = -log10(P.Value),
                q = qvalue::qvalue(.$P.Value)$qvalues
                )%>% 
  dplyr::arrange(desc(logFC))

#Save results
saveRDS(results_i_interaction, snakemake@output[["results_i_interaction"]])


# Interaction effects for type II fibers ---------------------------------
print_section_header("Running Limma interaction effect (RES vs B2A, type II) analyses")

#Subset SE
se_ii_interaction <- se[, se$fiber_type== 'II']

#Combine factors
TS_ii_interaction <- paste(se_ii_interaction$time, se_ii_interaction$intervention, sep=".")

#Set factor levels
TS_ii_interaction <- factor(TS_ii_interaction, levels=c("pre.resistance",
                                                      "post.resistance",
                                                      "pre.terbutaline",
                                                      "post.terbutaline"))
#Create design matrix
design_ii_interaction <- model.matrix(~0+ TS_ii_interaction)

#Specify column names in design matrix
colnames(design_ii_interaction)=c("pre.resistance",
                                 "post.resistance",
                                 "pre.terbutaline",
                                 "post.terbutaline")

#Estimate correlation coefficients
correlation_ii_interaction <- limma::duplicateCorrelation(
  SummarizedExperiment::assay(se_ii_interaction),
  design_ii_interaction,
  block=se_ii_interaction$id)

#Fit the linear model
fit_ii_interaction <- limma::lmFit(SummarizedExperiment::assay(se_ii_interaction),
                                  design_ii_interaction,
                                  block=se_ii_interaction$id,
                                  correlation=correlation_ii_interaction$consensus)

#Set up contrast
contrast_ii_interaction <- limma::makeContrasts(
    interaction=(post.resistance - pre.resistance)-(post.terbutaline - pre.terbutaline),
                             levels = design_ii_interaction)

#Apply contrast
fit2_ii_interaction <- limma::eBayes(limma::contrasts.fit(fit_ii_interaction, contrast_ii_interaction))

#Extract result
results_ii_interaction <- limma::topTable(fit2_ii_interaction, coef = "interaction", number = Inf) %>% 
  dplyr::mutate(xiao = P.Value^abs(logFC),
                protein = row.names(.),
                regulated = ifelse(xiao < 0.05, "+", ""),
                `-log10(p)` = -log10(P.Value),
                q = qvalue::qvalue(.$P.Value)$qvalues
                )%>% 
  dplyr::arrange(desc(logFC))

#Save results
saveRDS(results_ii_interaction, snakemake@output[["results_ii_interaction"]])


# MYH7 vs. MYH2 fibers ---------------------------------------------------
print_section_header("Running Limma analysis for type II vs. type I fibers")

#Subset SE
se_myh_ii_vs_i <- se[, se$time == 'pre']

#Create design matrix
design_myh_ii_vs_i <- model.matrix(~0+ se_myh_ii_vs_i$fiber_type)

#Name design file columns
colnames(design_myh_ii_vs_i)=c("i", "ii")

#Estimate correlation coefficients
correlation_myh_ii_vs_i <- limma::duplicateCorrelation(
  SummarizedExperiment::assay(se_myh_ii_vs_i),
  design_myh_ii_vs_i, block=se_myh_ii_vs_i$id)

#Set up contrast
contrast_myh_ii_vs_i <- limma::makeContrasts(ii - i,
                             levels = design_myh_ii_vs_i)

#Fit the linear model
fit_myh_ii_vs_i <- limma::lmFit(SummarizedExperiment::assay(se_myh_ii_vs_i), 
  design_myh_ii_vs_i,
  block=se_myh_ii_vs_i$id, 
  correlation=correlation_myh_ii_vs_i$consensus)

#Apply contrast
fit2_myh_ii_vs_i <- limma::eBayes(limma::contrasts.fit(fit_myh_ii_vs_i, contrast_myh_ii_vs_i))

#Extract result
results_myh_ii_vs_i <- limma::topTable(fit2_myh_ii_vs_i, coef = 1, number = Inf, sort.by = "logFC") %>% 
  dplyr::mutate(xiao = P.Value^abs(logFC),
                protein = row.names(.),
                regulated = ifelse(xiao < 0.05, "+", ""),
                `-log10(p)` = -log10(P.Value),
                q = qvalue::qvalue(.$P.Value)$qvalues
                )%>% 
  dplyr::arrange(desc(logFC))
   
#Save results
saveRDS(results_myh_ii_vs_i, snakemake@output[["results_myh_ii_vs_i"]])


# Export to supplementary files ------------------------------------------
export_list <- list(
  "Baseline fiber diffs (Fig.2d)" = clean_results(results_myh_ii_vs_i),
  "RES main effect (Fig.3b)" = clean_results(results_res_main),
  "RES type I fibers (Fig.3b)" = clean_results(results_res_i),
  "RES type II fibers (Fig.3b)" = clean_results(results_res_ii),
  "RES interaction (Fig.3b)" = clean_results(results_res_interaction),
  "B2A main effect (Fig.3c)" = clean_results(results_ter_main),
  "B2A type I fibers (Fig.3c)" = clean_results(results_ter_i),
  "B2A type II fibers (Fig.3c)" = clean_results(results_ter_ii),
  "B2A interaction (Fig.3c)" = clean_results(results_ter_interaction),
  "Between-groups main (Fig.5b)" = clean_results(results_i_and_ii_interaction),
  "Between-groups type I (Fig.5c)" = clean_results(results_i_interaction),
  "Between-groups type II (Fig.5d)" = clean_results(results_ii_interaction)
)

writexl::write_xlsx(export_list, snakemake@output[["supplementary_limma"]])
