

library(tidyverse)
library(here)


data_dir <- here("data")

list.files(data_dir)

tpm_df <- readRDS(file.path(data_dir, "bulk", "tpm_df.rds"))

takeda_sample_names <- 
  colnames(tpm_df) %>% str_subset("^D[:digit:][:digit:]-|^167...")

takeda_tpm_df <- tpm_df[,c("gene", takeda_sample_names)]

takeda_tpm_df2 <- takeda_tpm_df[complete.cases(takeda_tpm_df),]

write_tsv(takeda_tpm_df2, file.path(here("results", "seq_bulkCancer_tpm.tsv")))

