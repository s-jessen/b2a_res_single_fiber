source("renv/activate.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(snakecase)
  library(SummarizedExperiment)
  library(PhosR)
  library(tidyr)
  library(readxl)
  library(writexl)
})

#Columns to be excluded due to impurity
impure <- c("s_1", "s_2", "s_3", "s_4", "s_5", "s_6",
            "s_7", "s_8", "s_13", "s_15", "s_50", "s_52")

# Load and clean data ----------------------------------------------------

df <- read.csv2(snakemake@input[["raw_data"]]) %>%
  #Remove impure samples
  dplyr::select(-all_of(impure)) %>%
  #Make intensity data numeric
  dplyr::mutate_at(vars(-c("protein", "gene")), as.numeric) %>%
  #Replace NaN with NA
  dplyr::mutate_all(~replace(., is.nan(.), NA)) %>%
  #Remove genes without name
  dplyr::filter(!is.na(gene) & gene != "") %>%
  tibble::column_to_rownames(var="gene") %>%
  dplyr::select(-'protein') %>%
  dplyr::mutate_if(is.numeric, log2)

#Load metadata
metadata <- readxl::read_xlsx(snakemake@input[["design"]]) %>%
  tibble::column_to_rownames("sample") %>%
  dplyr::mutate(dplyr::across(dplyr::everything(), as.factor)) %>%
  #Remove impure samples
  dplyr::filter(!row.names(.) %in% impure) %>%
  dplyr::mutate(sample = row.names(.))

#Create Summarized Experiment (SE) of medianscaled data
se_raw <- PhosR::PhosphoExperiment(assay = PhosR::medianScaling(df), colData=metadata)

#Filter SE for 70% missing values in at least one fiber type
se <- PhosR::selectGrps(se_raw, colData(se_raw)$fiber_type, 0.7, n=1)

#Save files to data folder
saveRDS(metadata, snakemake@output[["metadata"]])
saveRDS(se_raw, snakemake@output[["se_raw"]])
saveRDS(se, snakemake@output[["se"]])


# Create long form data frame for all data -------------------------------

df_long <- SummarizedExperiment::assay(se) %>%
  tibble::rownames_to_column(var="protein") %>%
  tidyr::pivot_longer(
    cols = -protein,
    values_to = "abundance",
    names_to = "sample"
  ) %>%
  dplyr::left_join(metadata, by = "sample")

#Load annotations
#GO:terms and keywords
annotations <- readxl::read_xlsx(snakemake@input[["keywords"]]) %>%
  dplyr::rename_with(snakecase::to_snake_case) %>%
  dplyr::select(c("gene_names", "keywords",
                  "gene_ontology_biological_process",
                  "gene_ontology_cellular_component",
                  "gene_ontology_molecular_function")) %>%
  dplyr::rename(gobp = gene_ontology_biological_process,
                gocc = gene_ontology_cellular_component,
                gomf = gene_ontology_molecular_function,
                protein = gene_names) %>%
  dplyr::mutate(protein = gsub("\\ .*","", protein))%>%
  dplyr::mutate(protein = make.names(protein, unique=TRUE), protein)

#Mitocarta
mitocarta <- readxl::read_xls(snakemake@input[["mitocarta"]])%>%
  dplyr::select('symbol', 'pathways') %>%
  dplyr::rename(protein=symbol)

#Merge annotations
df_long <- df_long %>%
  merge(annotations, by="protein", all.x = T) %>%
  merge(mitocarta, by="protein", all.x = T) %>%
  dplyr::rename(mito = pathways)


# Create long form data frame with individual log2fold changes ------------

df_long_l2fc <- df_long %>%
  dplyr::group_by(protein, id, fiber_type) %>%
  dplyr::mutate(l2fc = abundance[time == "post"] - abundance[time == "pre"]) %>%
  dplyr::relocate(abundance, .after = id) %>%
  dplyr::relocate(l2fc, .after = abundance) %>%
  dplyr::ungroup() %>%
  dplyr::filter(time!="pre") %>%
  dplyr::select(-c("abundance", "time"))


# Create data frame of mean log2fold changes -----------------------------

df_long_l2fc_mean <- df_long_l2fc %>%
  dplyr::group_by(protein, fiber_type, intervention) %>%
  dplyr::summarize(mean_l2fc = mean(l2fc, na.rm = T),
                   n = sum(!is.na(l2fc)))


# Save to data folder ----------------------------------------------------

saveRDS(df_long, snakemake@output[["df_long"]])
saveRDS(df_long_l2fc, snakemake@output[["df_long_l2fc"]])
saveRDS(df_long_l2fc_mean, snakemake@output[["df_long_l2fc_mean"]])


# Export to supplementary files ------------------------------------------
#Original metadata
metadata <- readxl::read_xlsx(snakemake@input[["design"]])
writexl::write_xlsx(metadata, snakemake@output[["supplementary_metadata"]])


#Data files
supplementary_unfiltered <- read.csv2(snakemake@input[["raw_data"]]) %>% 
  #Make intensity data numeric
  dplyr::mutate_at(vars(-c("protein", "gene")), as.numeric) %>%
  #Replace NaN with NA
  dplyr::mutate_all(~replace(., is.nan(.), NA)) %>%
  #Remove genes without name
  dplyr::filter(!is.na(gene) & gene != "") %>%
  tibble::column_to_rownames(var="gene") %>%
  dplyr::select(-'protein') %>%
  dplyr::mutate_if(is.numeric, log2)

#Filter
supplementary_filtered <- PhosR::selectGrps(supplementary_unfiltered, metadata$fiber_type, 0.7, n=1)

#Retrieve "protein" column from rownames
supplementary_unfiltered <- supplementary_unfiltered %>% tibble::rownames_to_column("protein")
supplementary_filtered <- supplementary_filtered %>% tibble::rownames_to_column("protein")

#Export
export_list <- list(
  "Unfiltered" = supplementary_unfiltered,
  "Filtered for 70% valid values" = supplementary_filtered
)
writexl::write_xlsx(export_list, snakemake@output[["supplementary_data"]])
