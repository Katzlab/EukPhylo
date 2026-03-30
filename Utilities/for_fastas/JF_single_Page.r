library(tidyverse)
library(ggforce)  


gc3 <- read_tsv('CUBOutput_all/SpreadSheets/CompTrans.ENc.Raw.tsv') %>%
  mutate(
    GC3S = as.numeric(GC3S),
    ObsWrightENc_6Fold = as.numeric(ObsWrightENc_6Fold),
    color_fixed = "#d9d9d9",  # fixed color for all points
    size_fixed = 1             # fixed size for all points
  )

enc_null <- read_tsv('CUBOutput_all/SpreadSheets/ENc.Null.tsv')


ncol_plot <- 1
nrow_plot <- 1

n_pages <- n_pages(
  ggplot(gc3) +
    facet_wrap_paginate(~File, ncol = ncol_plot, nrow = nrow_plot, page = 1)
)

pdf("gc3_plots_uncolored.pdf", width = 7, height = 5)

for (p in seq_len(n_pages)) {
  
  plot <- ggplot(gc3, aes(x = GC3S, y = ObsWrightENc_6Fold)) +
    geom_point(color = gc3$color_fixed, size = gc3$size_fixed) +
    geom_line(data = enc_null, aes(x = GC3S, y = ENc, group = 1), color = "black") +
    scale_x_continuous(breaks = seq(0, 100, by = 10)) +
    labs(
      x = "%GC at 3rd-pos 4-fold sites",
      y = "Observed Wright ENc (6Fold)"
    ) +
    ggtitle(paste("New Backbone - Page", p)) +
    theme(
      legend.position = "none",
      strip.text = element_text(size = 7),
      axis.text.x = element_text(colour = "black"),
      axis.text.y = element_text(colour = "black")
    ) +
    facet_wrap_paginate(~File, ncol = ncol_plot, nrow = nrow_plot, page = p)
  
  print(plot)
}

dev.off()