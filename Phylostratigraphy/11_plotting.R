# =========================================================
# LOAD LIBRARIES
# =========================================================

library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(scales)

# =========================================================
# CREATE OUTPUT DIRECTORY
# =========================================================

output_dir <- "C:/Users/..../Desktop/new_phylostratigraphy"

if(!dir.exists(output_dir)){
  dir.create(output_dir, recursive = TRUE)
}

# =========================================================
# BRAIN DATA
# =========================================================

brain_data <- data.frame(
  Tissue = "Brain",
  
  Stratum = c(
    "Cellular",
    "Eukaryota",
    "Bilateria",
    "Insecta",
    "Hymenoptera",
    "Apoidea",
    "Halictidae",
    "Lasioglossum"
  ),
  
  DEG = c(
    342,
    319,
    116,
    46,
    68,
    8,
    19,
    26
  ),
  
  Non_DEG = c(
    3283,
    4738,
    932,
    435,
    509,
    107,
    272,
    374
  )
)

# =========================================================
# FAT BODY DATA
# =========================================================

fatbody_data <- data.frame(
  Tissue = "Fat body",
  
  Stratum = c(
    "Cellular",
    "Eukaryota",
    "Bilateria",
    "Insecta",
    "Hymenoptera",
    "Apoidea",
    "Halictidae",
    "Lasioglossum"
  ),
  
  DEG = c(
    535,
    646,
    153,
    57,
    82,
    18,
    23,
    36
  ),
  
  Non_DEG = c(
    3090,
    4411,
    895,
    424,
    495,
    97,
    268,
    364
  )
)

# =========================================================
# COMBINE DATA
# =========================================================

all_data <- bind_rows(brain_data, fatbody_data)

# =========================================================
# LONG FORMAT
# =========================================================

plot_data <- all_data %>%
  pivot_longer(
    cols = c(DEG, Non_DEG),
    names_to = "GeneType",
    values_to = "Count"
  )

# =========================================================
# ORDER FACTORS
# =========================================================

plot_data$Stratum <- factor(
  plot_data$Stratum,
  levels = c(
    "Cellular",
    "Eukaryota",
    "Bilateria",
    "Insecta",
    "Hymenoptera",
    "Apoidea",
    "Halictidae",
    "Lasioglossum"
  )
)

plot_data$GeneType <- factor(
  plot_data$GeneType,
  levels = c("Non_DEG", "DEG")
)

# =========================================================
# COMMON THEME
# =========================================================

common_theme <- theme_classic(base_size = 14) +
  
  theme(
    
    text = element_text(
      family = "Arial",
      face = "bold"
    ),
    
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      size = 12
    ),
    
    axis.text.y = element_text(
      size = 12
    ),
    
    axis.title.y = element_text(
      size = 14
    ),
    
    plot.title = element_text(
      size = 18,
      hjust = 0.5
    ),
    
    legend.title = element_blank(),
    
    legend.text = element_text(
      size = 13
    ),
    
    legend.position = "right"
  )

# =========================================================
# BRAIN PLOT
# =========================================================

brain_plot <- ggplot(
  filter(plot_data, Tissue == "Brain"),
  aes(
    x = Stratum,
    y = Count,
    fill = GeneType
  )
) +
  
  geom_bar(
    stat = "identity",
    position = position_dodge(width = 0.8),
    width = 0.7,
    color = "black",
    linewidth = 0.3
  ) +
  
  scale_fill_manual(
    values = c(
      "Non_DEG" = "#8338EC",
      "DEG" = "#FFBE0B"
    ),
    
    labels = c(
      "Non-DEG",
      "DEG"
    )
  ) +
  
  labs(
    title = "Brain",
    x = "",
    y = "Number of genes"
  ) +
  
  common_theme

# =========================================================
# FAT BODY PLOT
# =========================================================

fatbody_plot <- ggplot(
  filter(plot_data, Tissue == "Fat body"),
  aes(
    x = Stratum,
    y = Count,
    fill = GeneType
  )
) +
  
  geom_bar(
    stat = "identity",
    position = position_dodge(width = 0.8),
    width = 0.7,
    color = "black",
    linewidth = 0.3
  ) +
  
  scale_fill_manual(
    values = c(
      "Non_DEG" = "#8338EC",
      "DEG" = "#FFBE0B"
    ),
    
    labels = c(
      "Non-DEG",
      "DEG"
    )
  ) +
  
  labs(
    title = "Fat body",
    x = "",
    y = "Number of genes"
  ) +
  
  common_theme

# =========================================================
# COMBINE PLOTS
# =========================================================

final_plot <- brain_plot + fatbody_plot

# =========================================================
# SHOW PLOT
# =========================================================

final_plot

# =========================================================
# SAVE HIGH-QUALITY TIFF
# REQUIRED FOR JOURNAL SUBMISSION
# =========================================================

ggsave(
  filename = file.path(
    output_dir,
    "Phylostratigraphic_distribution_DEGs.tiff"
  ),
  
  plot = final_plot,
  
  width = 14,
  height = 6.5,
  
  dpi = 1200,
  
  compression = "lzw",
  
  bg = "white"
)

# =========================================================
# SAVE VECTOR PDF
# =========================================================

ggsave(
  filename = file.path(
    output_dir,
    "Phylostratigraphic_distribution_DEGs.pdf"
  ),
  
  plot = final_plot,
  
  width = 14,
  height = 6.5,
  
  device = cairo_pdf,
  
  bg = "white"
)

# =========================================================
# DONE
# =========================================================

cat("\nPublication-quality figures saved in:\n")
cat(output_dir)
cat("\n")
