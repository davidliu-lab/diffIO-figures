# load_hypoxia.R
# A helper script to load in all the data possibly required to plot the hypoxia
# related figures

immune_hi.v.lo_df <- readRDS(file.path(results_dir, "immune_hi.v.lo.rds"))
bulk_clin_metadata <- read_tsv(file.path(data_dir, "bulk", "bulk-clin-metadata.tsv")) %>% 
  filter(is.na(excluded_clin) & is.na(excluded_rna)) %>% 
  drop_na(sample_RNA) %>% 
  select(-immune_hi.v.lo, -immune_zscore, -immune_hi.v.lo_zscore) %>% 
  left_join(immune_hi.v.lo_df, by = "sample_RNA") 
immune_signatures <- read_tsv(file.path(results_dir, "immune_signatures.tsv"))$signatures

gs_vec <- readRDS(file.path(results_dir, "gs_vec.rds"))
ssgsea_df <- read_tsv(file.path(results_dir, "ssgsea_df.tsv"))

dge_immune_df <- read_tsv(file.path(results_dir, "dge_treatment_immune_df.tsv"))
gsea_immune_df <- read_tsv(file.path(results_dir, "gsea_R.v.NR_by_treatment_immune.tsv"))

dge_df <- read_tsv(file.path(results_dir, "dge_treatment.tsv"))
gsea_df <- read_tsv(file.path(results_dir, "gsea_treatment.tsv"))

immune_score_df <- readRDS(file.path(results_dir, "sc_immune_prop_df.rds"))
aucell_df <- readRDS(file.path(data_dir, "scRNAseq", "adata_aucell_scores_df.rds"))
obs_df <- read_csv(file.path(data_dir, "scRNAseq", "adata_pre_obs_2024-04-07.csv")) %>% 
  left_join(aucell_df, by = c("cell_label" = "cell"))

so_stromal <- readRDS(file.path(cache_dir, "sc-processing-sctransform-stromal", "so_stromal.rds"))
so_tumor <- readRDS(file.path(cache_dir, "sc-processing-sctransform-tumor", "so_tumor.rds"))
so_imm <- readRDS(file.path(cache_dir, "sc-processing-sctransform-immune", "so_immune.rds"))

deseq_obj_list <- readRDS(file.path(data_dir, "scRNAseq", "sc_deseq_obj_hypoxia_list.rds"))

mean_hypoxia_df <- readRDS(file.path(data_dir, "scRNAseq", "mean_hypoxia_df.rds"))
cpdb_tidy <- readRDS(file.path(results_dir, "cpdb_tidy.rds")) %>% 
  mutate(sample_ID = str_remove(id, "cpdb_")) %>%
  left_join(mean_hypoxia_df, by = "sample_ID") 
cpdb_tidy_sep <- cpdb_tidy %>% 
  separate(name, into = c("from", "to"), sep = "\\|")

new_obs_df <- obs_df %>% 
  group_by(compartment) %>% 
  mutate(HALLMARK_HYPOXIA_compartment_norm = as.vector(scale(HALLMARK_HYPOXIA))) 
scaled_mean_hypoxia_df <- new_obs_df %>% 
  group_by(sample_ID) %>% 
  summarise(mean_hypoxia = mean(HALLMARK_HYPOXIA_compartment_norm), 
            sd_hypoxia = sd(HALLMARK_HYPOXIA_compartment_norm)) %>% 
  left_join(immune_score_df, by = "sample_ID", suffix = c("", ".immune")) %>% 
  mutate(hypoxia_hi.v.lo = case_when(
    mean_hypoxia < quantile(mean_hypoxia, 0.333) ~ "low hypoxia", 
    mean_hypoxia < quantile(mean_hypoxia, 0.666) ~ "mid hypoxia", 
    TRUE ~ "high hypoxia"
  )) %>% 
  mutate(hypoxia_hi.v.lo = factor(hypoxia_hi.v.lo, levels = c("low hypoxia", "mid hypoxia", "high hypoxia"))) %>% 
  arrange(mean_hypoxia)
hypoxia_ordering <- scaled_mean_hypoxia_df$sample_ID

