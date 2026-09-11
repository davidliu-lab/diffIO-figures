# helper.R


# dge functions -----

run_wilcox <- function(.data, var_str = "response") {
  .data %>% 
    group_by(gene) %>% 
    do(w = wilcox.test(value ~ get(var_str), exact = FALSE, data = .)) %>% 
    summarise(gene, p.value = w$p.value) %>% 
    ungroup() %>% 
    mutate(p_BH = p.adjust(p.value, method = "BH"), 
           p_bonferroni = p.adjust(p.value, method = "bonferroni")) 
}

calc_change <- function(.data, var, index_var) {
  var <- enquo(var)
  index_var <- enquo(index_var)
  
  # chance greater based analysis
  data_list <- .data %>% 
    select(!!index_var, gene, value, !!var) %>% 
    drop_na(!!var) %>% 
    named_group_split(!!var)
  
  chance_greater_df <- full_join(data_list[["1"]], data_list[["0"]], by = "gene") %>% 
    group_by(gene) %>% 
    summarise(chance_greater = sum(value.x > value.y) / n(), 
              chance_lower = sum(value.y > value.x) / n()) %>% 
    mutate(chance_greater_over_lower = chance_greater - chance_lower)
  
  # average based analysis
  avg_data <- .data %>% 
    group_by(gene, !!var) %>% 
    summarise(median = median(value), 
              average = mean(value)) 
  median_df <- avg_data %>% 
    pivot_wider(id_cols = gene, names_from = !!var, values_from = median) %>% 
    mutate(med_diff = `1` - `0`, 
           med_fc = `1` / `0`, 
           norm_med_diff = med_diff / (`1` + `0`)) %>% 
    rename(median_r = `1`) %>% 
    rename(median_nr = `0`)
  avg_df <- avg_data %>% 
    pivot_wider(id_cols = gene, names_from = !!var, values_from = average) %>% 
    mutate(avg_fc = `1` / `0`, 
           avg_diff = `1` / `0`) %>% 
    rename(mean_r = `1`) %>% 
    rename(mean_cr = `0`)
  full_join(median_df, avg_df, by = "gene") %>% 
    left_join(chance_greater_df, by = "gene")
}


run_dge <- function(.tidy_data, raw_metadata, var, index_var) {
  var <- enquo(var)
  index_var <- enquo(index_var)
  
  metadata <- raw_metadata %>% 
    select({{ index_var }}, {{ var }}) %>% 
    drop_na()
  
  count_df <- metadata %>% 
    group_by({{ var }}) %>% 
    dplyr::count() %>% 
    pivot_wider(names_from = {{ var }}, 
                values_from = n)
  
  .data <- left_join(metadata, .tidy_data, by = c(rlang::as_name(index_var)))
  wilcox_res <- run_wilcox(.data, rlang::as_name(var))
  change_res <- calc_change(.data, {{ var }})
  
  full_join(wilcox_res, change_res, by = "gene") %>% 
    mutate(neglog10_p_BH = -log10(p_BH), 
           signed_p = ifelse(med_diff > 0, p.value, -p.value), 
           signed_neglog10_p = ifelse(chance_greater_over_lower > 0, 
                                      -log10(p.value), log10(p.value)), 
           signed_neglog10_p_BH = ifelse(chance_greater_over_lower > 0, 
                                         neglog10_p_BH, -neglog10_p_BH)) %>% 
    mutate(`n_R` = count_df$`1`, 
           `n_NR` = count_df$`0`)
}

run_dge_split <- function(.tidy_data, metadata, var, index_var, split_var) {
  var <- enquo(var)
  index_var <- enquo(index_var)
  split_var <- enquo(split_var)
  
  metadata %>% 
    named_group_split(!!split_var) %>% 
    map_dfr(~{
      subset_df <- .tidy_data %>%
        dplyr::filter(!!index_var %in% ..1[[rlang::as_name(index_var)]])
      run_dge(subset_df, ..1, !!var, !!index_var)
    }, .id = rlang::as_name(split_var))
}

# plotting functions -------

plot_surv_fit <-
  function(surv_fit,
           title_str,
           legend.labs = NULL,
           palette = NULL) {
    gg <- ggsurvplot(
      surv_fit,
      pval = TRUE,
      conf.int = FALSE,
      ncensor.plot = FALSE,
      risk.table = TRUE,
      risk.table.col = "strata",
      palette = palette,
      legend.labs = legend.labs,
      pval.coord = c(5, 0.2),
      ggtheme = my_theme(),
      size = 0.5,
      pval.size = CAPTION_SIZE * GGPLOT_TEXT_SCALE_FACTOR,
      fontsize = CAPTION_SIZE * GGPLOT_TEXT_SCALE_FACTOR,
      censor.shape = "|",
      censor.size = 1
    )
    plot_grid(
      gg$plot +
        labs(subtitle = title_str,
             x = "Months") +
        scale_y_continuous(expand = c(0, 0)) +
        theme(
          legend.position = "none",
          plot.subtitle = element_text(hjust = 0.5),
          plot.margin = margin()
        ),
      gg$table +
        theme(
          plot.title = element_text(
            size = CAPTION_SIZE,
            vjust = -2,
            hjust = 0
          ),
          legend.position = "none",
          axis.line.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.y = element_blank(),
          axis.line.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.x = element_blank(),
          axis.text.x = element_blank()
        ),
      ncol = 1,
      align = "v",
      axis = "lr",
      rel_heights = c(2, 1)
    )
  }



plot_surv_fit_wpval <-
  function(surv_fit,
           pairwise, 
           data, 
           group_col, 
           title_str,
           legend.labs = NULL,
           palette = NULL) {
    # Main survival plot
    gg <- ggsurvplot(
      surv_fit,
      conf.int = FALSE,
      ncensor.plot = FALSE,
      risk.table = TRUE,
      risk.table.col = "strata",
      palette = palette,
      legend.labs = legend.labs,
      ggtheme = my_theme(),
      size = 0.5,
      fontsize = CAPTION_SIZE * GGPLOT_TEXT_SCALE_FACTOR,
      censor.shape = "|",
      censor.size = 1
    )
    
    # Calculate pairwise p-values
    pairwise_pvals <- pairwise
    
    # Extract significant pairwise comparisons and format them for annotation
    pairwise_pvals_df <- as.data.frame(pairwise_pvals$p.value)
    comparisons <- combn(levels(data[[group_col]]), 2, simplify = FALSE)
    
    # Prepare annotations for pairwise p-values
    annotations <- lapply(seq_along(comparisons), function(i) {
      group1 <- comparisons[[i]][1]
      group2 <- comparisons[[i]][2]
      pval <- pairwise_pvals_df[group2, group1]
      if (!is.na(pval)) {
        list(
          group1 = group1, 
          group2 = group2, 
          xmin = which(levels(data[[group_col]]) == group1),
          xmax = which(levels(data[[group_col]]) == group2),
          label = paste0("p = ", signif(pval, 3)), 
          pval = pval
        )
      } else {
        NULL
      }
    })
    annotations <- map_dfr(annotations, bind_rows) 
    new_row <- data.frame(group1 = "", 
                          group2 = "", 
                          xmin = 0, 
                          xmax = 0, 
                          label = "otherwise, n.s.", 
                          pval = 0)
    annotations <- rbind(annotations, new_row)

    # Add annotations to the main survival plot
    if (!is.null(annotations)) {
      max_time <- max(gg$plot$data$time)
      gg$plot <- gg$plot +
        geom_text(
          data = annotations %>% filter(pval < 0.05) %>% mutate(n = row_number()), 
          hjust = 1.2, 
          aes(y = 1 - n / nrow(annotations), x = max_time, label = paste(group1, " vs ", group2, label, sep = "  "))
        )
    }
    
    plot_grid(
      gg$plot +
        labs(subtitle = title_str,
             x = "Months") +
        scale_y_continuous(expand = c(0, 0)) +
        theme(
          legend.position = "none",
          plot.subtitle = element_text(hjust = 0.5), 
          plot.margin = margin()
        ),
      gg$table +
        theme(
          plot.title = element_text(
            size = CAPTION_SIZE,
            vjust = -2,
            hjust = 0
          ),
          legend.position = "none",
          axis.line.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.y = element_blank(),
          axis.line.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.x = element_blank(),
          axis.text.x = element_blank()
        ),
      ncol = 1,
      align = "v",
      axis = "lr",
      rel_heights = c(2, 1)
    )
  }
