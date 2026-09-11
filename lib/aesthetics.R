# aesthetics.R

# constants -------

# CAPTION_SIZE <- 7
# TEXT_SIZE <- 7
# TITLE_SIZE <- 9
# GGPLOT_TEXT_SCALE_FACTOR <- 0.3

CAPTION_SIZE <- 14
TEXT_SIZE <- 14
TITLE_SIZE <- 18
GGPLOT_TEXT_SCALE_FACTOR <- 0.3

# geom defaults -------

update_geom_defaults("hline", list(size = 0.2))
update_geom_defaults("vline", list(size = 0.2))
update_geom_defaults("line", list(size = 0.2))
update_geom_defaults("segment", list(size = 0.2))
update_geom_defaults("path", list(size = 0.2))
update_geom_defaults("abline", list(size = 0.2))
update_geom_defaults("errorbar", list(size = 0.2))
update_geom_defaults("boxplot", list(size = 0.2, 
                                     outlier.size = 0,
                                     outlier.stroke = 0))
update_geom_defaults("point", list(size = 0.2))
update_geom_defaults("text", list(size = TEXT_SIZE*GGPLOT_TEXT_SCALE_FACTOR))

# themes --------

my_theme <- function(..., 
                     title_size = TITLE_SIZE, 
                     text_size = TEXT_SIZE, 
                     caption_size = CAPTION_SIZE) {
  theme_classic() +
    theme(
      
      # lines
      line = element_line(size = 0.2),
      axis.ticks = element_line(size = 0.2),
      
      # text
      plot.caption = element_text(size = caption_size),
      text = element_text(size = text_size),
      axis.title = element_text(size = text_size),
      legend.title = element_text(size = text_size),
      plot.subtitle = element_text(size = text_size, hjust = 0.5),
      title = element_text(size = title_size, hjust = 0.5),
      plot.title = element_text(
        size = title_size,
        hjust = 0.5,
        face = "bold"
      ),
      
      # other elements
      plot.background = element_blank(),
      panel.background = element_blank(),
      strip.background = element_blank(),
      legend.background = element_blank(),
      legend.key.height = unit(caption_size, "pt"),
      ...
    )
}

theme_set(my_theme())

figure_theme <- function() {
  my_theme(title_size = 6.5, text_size = 5, caption_size = 5)
}

# palettes ----------

palette_response <- c("R" = "darkgreen", "NR" = "goldenrod")
palette_recist <-
  c(
    "CR" = "green4",
    "PR" = "green3",
    "SD" = "grey",
    "PD" = "goldenrod"
  )

palette_treatment <- c(
  "aCTLA4" = "#848C8E",
  "aPD1" = "#272932",
  "aPD1 (prior aCTLA4)" = "#5887FF",#, "#C0FDFB"
  "combo" = "#447604" #B5CA8D",
)

palette_cohort <- c(
  `2015 Van Allen` = "#696773",
  `2016 Weber` = "#FED766",
  `2017 Riaz` = "#272727",
  `2019 Gide` = "#6F7C12",
  `2019 Liu` = "#009FB7",
  `2022 Freeman` = "darkorange", 
  `2023 Campbell (067)` = "maroon",
  `2024 Yang` = "mediumpurple3"
)

palette_immune_state <- c(
  "immune low" = "#3A5FCD", 
  "immune high" = "#CD3700"
)

palette_sex <- c(
  "FEMALE" = "darkorange2", 
  "MALE" = "aquamarine3"
)

palette_treatment_expanded <- c(
  "ICI_combo" = palette_treatment[["combo"]],
  "ICI_PD1" = palette_treatment[["aPD1"]],
  "other" = "grey",
  "other_plus_ICI" = "grey50",
  "targeted" = "orange",
  "targeted_plus_ICI" = "darkorange4"
)

palette_subtype_expanded <- c(
  "cutaneous" = "black", 
  "muco-cutaneous" = "darkorange4", 
  "mucosal" = "orange", 
  "uveal" = "aquamarine4", 
  "unknown primary" = "grey"
)

palette_tissue_type <- c(
  "skin" = "black", 
  "colon" = "maroon", 
  "liver" = "orange", 
  "soft" = "pink", 
  "lymph" = "aquamarine4"
)

palette_compartment <- c(
  "Immune" = "#931621",
  "Tumor" = "grey60",
  "Stromal" = "grey20"
)

palette_general_cell_type <- c(
  "Tumor" = "grey60", 
  "Keratinocytes" = "#9DCBBA", 
  "Endothelial cells" = "#E3B23C", 
  "Fibroblasts" = "#36413E", 
  "Myeloid" = "#B80C09", 
  "T cell" = "darkgreen", 
  "B cell" = "#1C77C3", 
  "Contaminating" = "grey80"
)

palette_specific_immune <- c(
  "B cells" = "#1C77C3", 
  "Plasma" = "#192BC2", 
  "Naive-like T" = "#A6D3A0", #"#BBBE64", 
  "CD8+ T" = "darkgreen", 
  "Tregs" = "#618985", 
  "NK" = "#592941", 
  "DC" = "#B80C09", 
  "pDC" = "#9C27B0", 
  "TAM" = "#FF8C42", 
  "Neutrophil" = "#EED7C5",
  "Mast cell" = "#C19875" 
)

palette_cor <- c("#00688B", "white", "#FF8C00")

palette_hypoxia_state <- c(
  "low" = "pink3", 
  "mid" = "grey", 
  "high" = "grey20"
)

palette_bool <- c(
  "TRUE" = "red", 
  "FALSE" = "grey"
)
