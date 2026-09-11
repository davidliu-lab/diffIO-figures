# helper_Seurat.R

# Rsc (github.com/amyh25/Rsc) supplies the single-cell helpers that chapters
# call unqualified: get_gene_expr_from_so(), seurat_to_sce(), qc_sce(),
# make_pseudobulk(), plot_gene(), etc.
library(Rsc)


# Import from adata -------

#' seurat_from_adata
#'
#' @param adata_filepath path to adata
#' @param anndata specify an attached anndata
#' @param layer layer from adata to export. Default is X
#' @param project.name project name
#' @param ... Parameters passed to `CreateSeuratObject`
#'
#' @return a Seurat object

seurat_from_adata <-
  function(adata_filepath,
           anndata,
           layer = "X",
           project.name = "adata",
           ...) {
    adata <- anndata$read_h5ad(adata_filepath)
    if (layer == "X")
      counts <- Matrix::t(adata$X)
    else
      counts <- Matrix::t(adata$layers$get(layer))
    row.names(counts) <- row.names(adata$var)
    colnames(counts) <- row.names(adata$obs)
    so <- CreateSeuratObject(
      counts = counts,
      project = project.name,
      assay = "RNA",
      meta.data = adata$obs,
      row.names = adata$var_names,
      ...
    )
    return(so)
  }

# HTO helper functions ------

#'  get_hto_tidy
#'
#' @param so Seurat object
#' @param column_regex subset of string to select columns by using `str_subset`
#'
#' @return a tidy tibble with cellID, hashtag, and hashtag expr.

get_hto_tidy <- function(so, column_regex = "hashtag") {
  hto_tidy <- so@assays$HTO@data %>%
    t() %>%
    as_tibble(rownames = "cell") %>%
    pivot_longer(cols = str_subset(colnames(.), column_regex),
                 names_to = "hashtag",
                 values_to = "value")
}


#'  plot_hto_tidy
#'
#' @param hto_tidy a tidy tibble with cellID, hashtag, and hashtag expr.
#'
#' @return a ggplot2 object

plot_hto_ridges <- function(hto_tidy) {
  hto_tidy %>%
    ggplot(aes(value, hashtag, color = hashtag)) +
    ggridges::geom_density_ridges(fill = NA) +
    scale_y_discrete(expand = c(0.01, 0)) +
    scale_x_continuous(expand = c(0.01, 0)) +
    theme_ridges() +
    theme(panel.grid.major.x = element_blank(),
          legend.position = "none")
}


#'  plot_hto_density
#'
#' @param hto_df a wide tibble with cellID and hashtags as columns
#' @param x column var of hto_df, variable to plot on x-axis
#' @param y column var of hto_df, variable on y-axis
#' @param bins number of bins in the 2d density plots
#' @param limit limit of plot edge
#' @param group_str string specifying some clustering or other group to plot density by
#'
#' @return a ggplot2 object

plot_hto_density <- function(hto_df, x = hashtag1, y = hashtag2,
                             bins = 15,
                             limit = max(c(hto_df$hashtag1, hto_df$hashtag2))*0.5,
                             group_str = NA) {
  x <- enquo(x)
  y <- enquo(y)
  if (is.na(group_str)) {
    p <- hto_df %>%
      ggplot() +
      aes(!!x, !!y) +
      geom_density_2d(bins = bins)
  } else {
    if (!is.character(group_str))
      stop("Variable group_str must be a string")
    p <- hto_df %>%
      ggplot() +
      aes(!!x, !!y) +
      geom_density_2d(bins = bins, aes(color = get(group_str))) +
      labs(color = group_str)
  }
  p <- p +
    scale_y_continuous(expand = c(0, 0), limits = c(0, limit)) +
    scale_x_continuous(expand = c(0, 0), limits = c(0, limit)) +
    coord_fixed()
  return(p)
}


# Seurat helper functions -------

#'  run_seurat_rna_workflow
#'
#' @param so Seurat object
#' @param run_normalize bool specifying whether to run normalization
#' @param run_scale bool specifying whether to scale
#' @param scale_features char vector specifying which features to scale. Default to
#' @param run_SCT bool specifying whether to run SCTransform
#' @param run_clustering bool specifying whether to run dimensionality reduction
#' and clustering.
#' @param dims dimensions to use in finding neighbors and UMAP
#' @param resolutions clustering resolutions to calculate
#'
#' @return processed Seurat object

run_seurat_rna_workflow <- function(so,
                                    set_percentages = TRUE,
                                    run_normalize = TRUE,
                                    run_scale = TRUE,
                                    scale_features = NULL,
                                    run_SCT = FALSE,
                                    SCT_vst_flavor = "v2", 
                                    run_clustering = TRUE,
                                    dims = 1:30,
                                    resolutions = c(0.1, 0.5, 0.9)) {
  # set percentages
  if (set_percentages) {
    so <- PercentageFeatureSet(so, pattern = "^mt-", col.name = "percent.mt")
    so <- PercentageFeatureSet(so, pattern = "^Rpl|^Rps", col.name = "percent.ribo")
  }
  
  if (run_normalize) {
    so <- so %>%
      NormalizeData() %>%
      FindVariableFeatures()
  }
  if (run_scale) {
    so <- so %>% ScaleData(features = scale_features)
  }
  if (run_SCT) {
    so <- so %>% SCTransform(method = "glmGamPoi", vst.flavor = SCT_vst_flavor)
  }
  if (run_clustering) {
    so <- so %>%
      RunPCA(features = VariableFeatures(.)) %>%
      RunUMAP(dims = dims)
    for (res in resolutions) {
      so <- so %>%
        FindNeighbors(dims = dims) %>%
        FindClusters(resolution = res)
    }
  }
  return(so)
}

#'  run_seurat_hto_clustering
#'
#' This function runs HTO clustering using clara.
#'
#' @param so Seurat object
#' @param metric distance metric for clara clustering
#' May want to adjust if calculations are taking a long time.
#'
#' @return processed Seurat object

run_seurat_hto_clustering <- function(so, metric = "euclidean") {
  # calculate clusters
  k <- nrow(so@assays$HTO@data)
  clara_out <- cluster::clara(
    t(so@assays$HTO@data),
    metric = metric,
    k = k,
    samples = ncol(so@assays$HTO@data),
    cluster.only = TRUE
  )
  cluster_df <- as.data.frame(clara_out) %>%
    dplyr::rename(cluster = `clara_out`) %>%
    mutate(cluster = factor(cluster))
  
  # assign clusters to HTOs
  hto_df <- t(so@assays$HTO@data)
  cluster_hash.ID_key_df <- hto_df %>%
    bind_cols(cluster_df) %>%
    group_by(cluster) %>%
    summarise(across(.cols = starts_with("hashtag"),
                     .fns = median,
                     .names = "{.col}")) %>%
    rowwise() %>%
    # I know the below mutate is really gross to read but it works
    # and I'll think of a nicer way to write it later
    mutate(hash.ID = str_subset(colnames(.), "^hashtag")[which.max(c_across(starts_with("hashtag")))]) %>%
    select(cluster, hash.ID) %>%
    ungroup()
  cluster_df <- cluster_df %>%
    as_tibble(rownames = "cell") %>%
    left_join(cluster_hash.ID_key_df, by = "cluster") %>%
    column_to_rownames("cell")
  
  so <- AddMetaData(so, metadata = cluster_df)
  return(so)
}

#'  get_metadata_from_so
#'
#' @param so Seurat object
#' @param assay
#' @param slot
#' @param genes a character vector of genes or NULL
#' @param metadata a character vector of columns from metadata, "all", or NULL
#' @param reduction a string specifying the reduction.
#'
#' @return data table

get_metadata_from_so <- function(so,
                                 assay = "RNA",
                                 genes = NULL,
                                 metadata = "all",
                                 reduction = "umap") {
  if (is.null(metadata)) {
    out_df <- so@meta.data %>%
      as_tibble(rownames = "cell") %>%
      dplyr::select(cell)
  } else if (length(metadata) > 1) {
    out_df <- so@meta.data %>%
      as_tibble(rownames = "cell") %>%
      dplyr::select(cell, all_of(metadata))
  } else if (metadata == "all") {
    out_df <- so@meta.data %>%
      as_tibble(rownames = "cell")
  } else {
    out_df <- so@meta.data %>%
      as_tibble(rownames = "cell") %>%
      dplyr::select(cell, all_of(metadata))
  }
  
  if (!is.null(genes)) {
    # get index
    gene_i <- match(genes, rownames(so))
    
    # check that genes exists
    if (any(is.na(gene_i))) {
      stop(paste0(genes[which(is.na(gene_i))], " not found in data."))
    }
    
    # get gene expression
    if (length(genes) == 1) {
      gene_df <- so@assays[[assay]]@data[genes, ] %>%
        as_tibble(rownames = "cell")
      gene_df[[genes]] <- gene_df$value
      gene_df <- gene_df %>%
        dplyr::select(-value)
    } else {
      gene_df <- so@assays[[assay]]@data[genes, ] %>%
        as.matrix() %>% t() %>%
        as_tibble(rownames = "cell")
    }
    
    out_df <- out_df %>% left_join(gene_df, by = "cell")
  }
  
  if (!is.null(reduction)) {
    if (length(reduction) > 1) {
      stop(paste0("reduction must be of length 1"))
    } else if (!(reduction %in% names(so@reductions))) {
      stop(paste0(reduction, " not found in reductions. "))
    }
    
    reduc_df <- so@reductions[["umap"]]@cell.embeddings %>%
      as_tibble(rownames = "cell")
    # Seurat 5 lowercased the UMAP reduction key, so embeddings come back as
    # umap_1/umap_2. The book (and RFunctions::plot_umap) expect Seurat 4's
    # UMAP_1/UMAP_2.
    names(reduc_df) <- sub("^umap_", "UMAP_", names(reduc_df))
    out_df <- out_df %>% left_join(reduc_df, by = "cell")
  }
  
  return(out_df)
}

#'  set_idents
#' @param so Seurat object
#' @param ident string specifying identity
#' @return Seurat object


set_idents <- function(so, ident) {
  Idents(so) <- ident
  return(so)
}

#'  set_idents
#' @param so Seurat object
#' @param assay string specifying assay
#' @return Seurat object

set_assay <- function(so, assay) {
  DefaultAssay(so) <- assay
  return(so)
}

#'  get_wide_tbl_from_markers
#'
#' Makes wide table from markers and will sort columns such that
#' col1 (typically p-value) and col2 (typically fold change) are
#' next to each other.
#'
#' @param .tbl a tibble, typically the output to FindAllMarkers
#' @param gene_col column where gene information is found
#' @param cluster_col column where clustering information is found
#' @param col1 typically pval
#' @param col2 typically log2FC
#' @param sep separating character; an double underscore by default but if your
#' column names includes double underscores may want to change
#'
#' @return a tibble

get_wide_tbl_from_markers <- function(.tbl,
                                      gene_col = gene,
                                      cluster_col = cluster,
                                      col1 = p_val_adj,
                                      col2 = avg_log2FC,
                                      sep = "__") {
  gene_col <- enquo(gene_col)
  cluster_col <- enquo(cluster_col)
  col1 <- enquo(col1)
  col2 <- enquo(col2)
  
  markers_wide <- .tbl %>%
    pivot_wider(
      id_cols = !!gene_col,
      names_from = !!cluster_col,
      values_from = c(!!col1, !!col2),
      names_sep = sep
    )
  order_i <- colnames(markers_wide)[2:ncol(markers_wide)] %>%
    str_extract(paste0(sep,".*$")) %>%
    order()
  markers_wide %>%
    dplyr::select(1, order_i + 1)
}

# TODO implement other ways to get wide markers,
# Here's some sample code that I wrote that could be generalized:

# markers_wide <- markers_df %>% get_wide_tbl_from_markers()
# write_csv(markers_wide, file.path(results_dir, "so_filterRNA200mito10_subsetM_markers08wide.csv"))
#
# markers_readable <- markers_df %>%
#     filter(p_val_adj < 0.01) %>%
#     filter(avg_log2FC > 0) %>%
#     group_by(cluster) %>%
#     arrange(desc(avg_log2FC)) %>%
#     mutate(rank = row_number()) %>%
#     pivot_wider(id_cols = rank,
#                 names_from = cluster,
#                 values_from = gene) %>%
#     dplyr::select(rank, as.character(c(0:20)))
# write_csv(markers_readable, file.path(results_dir, "so_filterRNA200mito10_subsetM_markers08readable.csv"))

# SCP export functions -------

#'  scp_export_counts
#'
#' Exports counts from Seurat object and saves in Single Cell Portal format
#' as a compressed text file.
#'
#' @param so Seurat object
#' @param output_dir output directory
#' @param assay assay from which to export data
#' @param suffix suffix for filename
#'
#' @return NULL

scp_export_counts <- function(so, output_dir, assay = "RNA", suffix = "") {
  if (!dir.exists(output_dir)) {
    message(paste0(output_dir, " not found. Creating..."))
    dir.create(output_dir)
  }
  write_tsv(
    as_tibble(so@assays[["RNA"]]@counts, rownames = "GENE"),
    file.path(output_dir, paste0("scp_counts", suffix, ".txt.gz"))
  )
  message("Done!")
}

#'  scp_export_data
#'
#' Exports data from Seurat object and saves in Single Cell Portal format
#' as a compressed text file
#'
#' @param so Seurat object
#' @param output_dir output directory
#' @param assay assay from which to export data
#' @param suffix suffix for filename
#'
#' @return NULL

scp_export_data <- function(so, output_dir, assay = "RNA", suffix = "") {
  if (!dir.exists(output_dir)) {
    message(paste0(output_dir, " not found. Creating..."))
    dir.create(output_dir)
  }
  write_tsv(
    as_tibble(so@assays[["RNA"]]@data, rownames = "GENE"),
    file.path(output_dir, paste0("scp_lognormcounts", suffix, ".txt.gz"))
  )
}


#'  scp_export_metadata_file
#'
#' The purpose of this function is to export metadata from a Seurat object
#' to a single cell portal compatible text file. It will fill in default
#' metadata columns if none are specific, and will use `stringr` regular
#' expression to select columns when requested.
#'
#' See the help page from the Single Cell Portal
#' (https://singlecell.zendesk.com/hc/en-us/articles/360060609852-Required-Metadata)
#' for information on the proper format for metadata tables.
#'
#' @param so Seurat object
#' @param output_dir output directory
#' @param filename the output filename
#' @param biosample_id a symbol specifying a column in the Seurat object metadata
#' @param donor_id a symbol specifying a column in the Seurat object metadata
#' @param str_subset_regex_group a regular expression to select categorical columns
#' @param str_subset_regex_numeric a regular expression to select numerical columns
#' @param species code from the NCBI Taxon
#' @param species__ontology_label name from NCBI Taxon
#' @param disease code from MONDO
#' @param disease__ontology_label name from MONDO
#' @param organ code from UBERON
#' @param organ__ontology_label name from UBERON
#' @param library_preparation_protocol code from EFO
#' @param library_preparation_protocol__ontology_label name from EFO
#' @param sex a string, one of ["male", "female", "mixed", "unknown"]
#'
#' @return NULL

scp_export_metadata_file <- function(so,
                                     output_dir,
                                     filename = "scp_metadata",
                                     biosample_id = "orig.ident",
                                     donor_id = "orig.ident",
                                     str_subset_regex_group = NULL,
                                     str_subset_regex_numeric = NULL,
                                     species = "NCBITaxon_1",
                                     species__ontology_label = "root",
                                     disease = "MONDO_0000001",
                                     disease__ontology_label = "disease or disorder",
                                     organ = "UBERON_0001062",
                                     organ__ontology_label = "anatomical entity",
                                     library_preparation_protocol = "EFO_0010183",
                                     library_preparation_protocol__ontology_label = "single cell library construction",
                                     sex = "unknown") {
  
  metadata_df <- so@meta.data %>%
    as_tibble(rownames = "NAME")
  
  required_metadata_df <- metadata_df %>%
    mutate(biosample_id = get(biosample_id)) %>%
    mutate(donor_id = get(donor_id)) %>%
    mutate(species = species, species__ontology_label = species__ontology_label) %>%
    mutate(disease = disease, disease__ontology_label = disease__ontology_label) %>%
    mutate(organ = organ, organ__ontology_label = organ__ontology_label) %>%
    mutate(
      library_preparation_protocol = library_preparation_protocol,
      library_preparation_protocol__ontology_label = library_preparation_protocol__ontology_label
    ) %>%
    mutate(sex = sex) %>%
    select(
      NAME,
      biosample_id,
      donor_id,
      species,
      species__ontology_label,
      disease,
      disease__ontology_label,
      organ,
      organ__ontology_label,
      library_preparation_protocol,
      library_preparation_protocol__ontology_label,
      sex
    )
  
  type_vec <- rep("group", ncol(required_metadata_df) - 1)
  
  if (length(str_subset_regex_group) != 0) {
    message("Looking for matching group columns...")
    group_cols <-
      str_subset(colnames(metadata_df), str_subset_regex_group)
    group_metadata_df <- metadata_df %>%
      dplyr::select(NAME, all_of(group_cols))
    message(paste0("Found column: ", group_cols, "\n"))
    type_vec <-
      c(rep("group", ncol(group_metadata_df) - 1), type_vec)
  }
  if (length(str_subset_regex_numeric) != 0) {
    message("Looking for matching numeric columns...")
    numeric_cols <-
      str_subset(colnames(metadata_df), str_subset_regex_numeric)
    numeric_metadata_df <- metadata_df %>%
      dplyr::select(NAME, all_of(numeric_cols))
    message(paste0("Found column: ", numeric_cols, "\n"))
    type_vec <-
      c(rep("numeric", ncol(numeric_metadata_df) - 1), type_vec)
  }
  
  type_vec <- c("TYPE", type_vec)
  types_line <- paste0(type_vec, collapse = "\t")
  
  scp_metadata_df <- numeric_metadata_df %>%
    left_join(group_metadata_df, by = "NAME") %>%
    left_join(required_metadata_df, by = "NAME")
  colnames(scp_metadata_df) <-
    str_replace_all(colnames(scp_metadata_df), "\\.", "\\_")
  
  write_tsv(scp_metadata_df, file.path(output_dir, "tmp.txt"))
  f_raw <- file(file.path(output_dir, "tmp.txt"), open = "r")
  lines <- readLines(f_raw)
  close(f_raw)
  new_lines <- c(lines[1], types_line, lines[2:length(lines)])
  f_new <-
    file(file.path(output_dir, paste0(filename, ".txt")), open = "w")
  writeLines(new_lines, f_new)
  close(f_new)
  file.remove(file.path(output_dir, "tmp.txt"))
  message("done!")
  
}

#'  scp_export_clustering_file
#'
#' The purpose of this function is to export clustering from a Seurat object
#' to a single cell portal compatible text file. It will fill in default
#' metadata columns if none are specific, and will use `stringr` regular
#' expression to select columns when requested.
#'
#' @param so Seurat object
#' @param output_dir output directory
#' @param filename the output filename
#' @param str_subset_regex_group a regular expression to select categorical columns
#' @param str_subset_regex_numeric a regular expression to select numerical columns
#' @param reduction_str a regular expression to select the two columns specifying the reduction
#'
#' @return NULL


scp_export_clustering_file <- function (so,
                                        output_dir,
                                        filename = "scp_clustering",
                                        str_subset_regex_group = "RNA_snn_res.",
                                        str_subset_regex_numeric = NULL,
                                        reduction_str = "umap") {
  metadata_df <- so@meta.data %>%
    as_tibble(rownames = "NAME")
  
  umap_df <-
    so@reductions[[reduction_str]]@cell.embeddings %>%
    as_tibble(rownames = "NAME")
  colnames(umap_df) <- c("NAME", "X", "Y")
  
  clustering_cols <- c("NAME", "X", "Y")
  clustering_types_vec <- c("TYPE", "numeric", "numeric") %>%
    set_names(clustering_cols)
  
  if (length(str_subset_regex_group) != 0) {
    message("Looking for matching group columns...")
    group_cols <-
      str_subset(colnames(metadata_df), str_subset_regex_group)
    group_metadata_df <- metadata_df %>%
      dplyr::select(NAME, all_of(group_cols))
    message(paste0("Found column: ", group_cols, "\n"))
    clustering_types_vec <-
      c(clustering_types_vec, rep("group", ncol(group_metadata_df) - 1))
    umap_df <-
      umap_df %>% left_join(group_metadata_df, by = "NAME")
  }
  if (length(str_subset_regex_numeric) != 0) {
    message("Looking for matching numeric columns...")
    numeric_cols <-
      str_subset(colnames(metadata_df), str_subset_regex_numeric)
    numeric_metadata_df <- metadata_df %>%
      dplyr::select(NAME, all_of(numeric_cols))
    message(paste0("Found column: ", numeric_cols, "\n"))
    clustering_types_vec <-
      c(clustering_types_vec, rep("numeric", ncol(numeric_metadata_df) - 1))
    umap_df <-
      umap_df %>% left_join(numeric_metadata_df, by = "NAME")
  }
  
  clustering_types_line <-
    paste0(clustering_types_vec, collapse = "\t")
  clustering_scp <- umap_df
  
  write_tsv(clustering_scp, file.path(output_dir, "tmp.txt"))
  f_clustering_raw <- file(file.path(output_dir, "tmp.txt"),
                           open = "r")
  lines <- readLines(f_clustering_raw)
  close(f_clustering_raw)
  new_lines <-
    c(lines[1], clustering_types_line, lines[2:length(lines)])
  f_clustering <- file(file.path(output_dir, paste0(filename,
                                                    ".txt")), open = "w")
  writeLines(new_lines, f_clustering)
  close(f_clustering)
  file.remove(file.path(output_dir, "tmp.txt"))
  message("done!")
}


# Rsc::get_gene_expr_from_so reaches into so@assays[[assay]]@data, a slot that
# Seurat 5's Assay5 class does not have (it stores layers instead). Shadow it
# with a version that goes through the accessor, which works on Assay and
# Assay5 alike. Behaviour is otherwise identical.
get_gene_expr_from_so <- function(so, gene_vec, assay = "RNA", skip_missing = FALSE) {
  if (!sum(gene_vec %in% rownames(so)) == length(gene_vec)) {
    missing_genes <- gene_vec[!(gene_vec %in% rownames(so))]
    if (skip_missing) {
      message(paste0(missing_genes, " not found. Skipping..."))
      gene_vec <- gene_vec[gene_vec %in% rownames(so)]
    } else {
      stop(paste0(missing_genes, " not found. Please remove. "))
    }
  }
  SeuratObject::LayerData(so, layer = "data", assay = assay) %>%
    .[rownames(.) %in% gene_vec, ] %>%
    as.matrix() %>% t() %>% as_tibble(rownames = "cell")
}
