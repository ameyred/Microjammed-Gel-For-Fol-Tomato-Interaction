library(ggplot2)
library(readr)
library(dplyr)
library(tidyr)
library(scales)
library(viridis)

# Step 1: Read CSV
df <- read_csv("~/Amey_Lab/Microgel_Project/Big_Experiment/Tomato_Analysis/Interpro_Term_Enrichment/Upregulated_Genes/IPR_With_Expression_Filtered.csv")

# Step 2: Pivot to long format
df_long <- df %>%
  pivot_longer(cols = c("Log2FC_Tip", "Log2FC_Root", "Log2FC_AR"),
               names_to = "Tissue_FC", values_to = "Log2FC") %>%
  pivot_longer(cols = c("padj_Tip", "padj_Root", "padj_AR"),
               names_to = "Tissue_padj", values_to = "padj") %>%
  filter(substr(Tissue_FC, 8, 10) == substr(Tissue_padj, 6, 8)) %>%
  mutate(
    Tissue = case_when(
      grepl("Tip", Tissue_FC) ~ "Tip",
      grepl("Root", Tissue_FC) ~ "Root",
      grepl("AR", Tissue_FC) ~ "AR",
      TRUE ~ NA_character_
    ),
    padj = ifelse(is.na(padj) | padj == 0, 1e-300, padj),
    neg_log10_padj = -log10(padj)
  )

# Step 3: Pad each IPR Term to 10 genes
df_padded <- df_long %>%
  group_by(`IPR Term`) %>%
  group_split() %>%
  lapply(function(group_df) {
    n_real <- length(unique(group_df$`Gene ID`))
    n_needed <- 10 - n_real
    if (n_needed > 0) {
      dummy_ids <- paste0("DUMMY_", seq_len(n_needed))
      dummy_rows <- expand.grid(
        `Gene ID` = dummy_ids,
        Tissue = c("Tip", "Root", "AR"),
        stringsAsFactors = FALSE
      ) %>%
        mutate(
          `IPR Term` = unique(group_df$`IPR Term`),
          Name = unique(group_df$Name),
          Log2FC = 0,
          padj = NA,
          neg_log10_padj = NA
        )
      bind_rows(group_df, dummy_rows)
    } else {
      group_df
    }
  }) %>%
  bind_rows()

# Step 4: Save each IPR plot
output_dir <- path.expand("~/Amey_Lab/Microgel_Project/Big_Experiment/Tomato_Analysis/Interpro_Term_Enrichment/Dot_Plots_With_AR")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

unique_iprs <- unique(df_padded$`IPR Term`)

for (ipr in unique_iprs) {
  ipr_data <- df_padded %>% filter(`IPR Term` == ipr)
  gene_order <- unique(ipr_data$`Gene ID`)
  ipr_data$`Gene ID` <- factor(ipr_data$`Gene ID`, levels = gene_order)
  ipr_name <- unique(ipr_data$Name)[1]
  ipr_safe <- gsub("[^A-Za-z0-9_]+", "_", ipr_name)
  
  p <- ggplot(ipr_data, aes(x = `Gene ID`, y = Tissue)) +
    geom_point(aes(size = Log2FC, color = neg_log10_padj), na.rm = TRUE) +
    scale_color_viridis(
      name = "-log10(padj)",
      option = "D",
      direction = -1,
      limits = c(0, 10),
      na.value = "grey90"
    ) +
    scale_size(
      name = "Log2FC",
      range = c(2, 10),
      limits = c(-6, 8),
      breaks = seq(-6, 8, by = 2),
      guide = guide_legend(nrow = 1)
    )+
    theme(
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank(),
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      aspect.ratio = 0.15
    ) +
    labs(
      title = paste("IPR:", ipr_name),
      x = "Gene ID",
      y = "Tissue"
    )
  
  ggsave(
    filename = file.path(output_dir, paste0(ipr_safe, ".png")),
    plot = p,
    width = 12,
    height = 5,
    dpi = 300
  )
}
