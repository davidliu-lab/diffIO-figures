library(Seurat)

list.files(here("data", "scRNAseq"))

so <- readRDS(here("data", "scRNAseq", "so_jerbyarnon_AHaddedclinmetadata.rds"))

so@meta.data$cell.type %>% unique()

expr_mtx <- so@assays$RNA@data %>% 
  as_tibble(rownames = "gene")

write_tsv(expr_mtx, file.path(here("data", "scRNAseq", "so_jerbyarnon_expression.txt")))

cell_type_metadata <- so@meta.data %>% 
  as_tibble(rownames = "cell") %>% 
  select(cell, cell.type)

write_tsv(cell_type_metadata, file.path(here("data", "scRNAseq", "so_jerbyarnon_metadata.txt")))
