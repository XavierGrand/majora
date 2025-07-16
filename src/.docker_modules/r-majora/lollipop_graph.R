#!/usr/bin/Rscript

library(tidyverse, quietly = TRUE)
library(patchwork, quietly = TRUE)
library(conflicted, quietly = TRUE)
library(plotly)
library(gapminder)

conflict_prefer("filter", "dplyr")
conflict_prefer("lag", "dplyr")


# ajout de la possibilité de donner des arguments
 args <- commandArgs(trailingOnly = TRUE)

# vérification qu"un argument a bien été donné
if(length(args)==0){
   print("Please include differential expression results!")
   stop("Requires command line argument.")
 }

# :: pour utiliser la fonction read_tsv sans faire un library(readr) avant
vcf_table <- readr::read_tsv(args[1], col_names = FALSE, col_types = "icccddd")

header <- c("positions", "ref_allele", "alt_allele", "variation_type", "amb_reads", "ref_reads", "alt_reads")
colnames(vcf_table) <- header


# pour ajouter un axe Y
vcf_table <- vcf_table %>%
  mutate(percentage_of_reads_containing_variation = alt_reads)


if (nrow(vcf_table) == 1) { # création du lollipop linéaire
  main <- ggplot(vcf_table, aes(x = positions, y = percentage_of_reads_containing_variation, color = variation_type)) +
    # ajouter les barres
    geom_segment(
      aes(xend = positions, yend = 0),
      size = 0.6
    ) +
    # ajouter les ronds
    geom_point(aes(size = percentage_of_reads_containing_variation)) +
    scale_color_manual(values = c(
      "SNP" = "#D81B60",
      "INS" = "#1E88E5",
      "DEL" = "#004D40"
    ), guide = "legend") +
    # Assign a range for the size of the dots
    scale_size(
      range = c(1, 5), 
      guide = "none"
    ) +
    theme(
      legend.position = "right",
      legend.text = element_text(size = 12) 
      )
} else { # création du lollipop circulaire
  main <- ggplot(vcf_table, aes(x = positions, y = percentage_of_reads_containing_variation , color = variation_type)) +
    # ajouter les barres
    geom_segment(
      aes(xend = positions, yend = 0),
      linewidth = 0.6
    ) +
    # ajouter les ronds
    geom_point(aes(size = percentage_of_reads_containing_variation)) +
    # Use the Prism color scale for the categories
      scale_color_manual(values = c(
      "SNP" = "#D81B60",
      "insertion" = "#1E88E5",
      "deletion" = "#004D40"
    ), guide = "legend") +
    # Assign a range for the size of the dots
    scale_size(
      range = c(1, 5), 
      guide = "none"
    ) +
    # Make it circular!
    coord_polar() +
    theme(
      legend.position = "right",
      legend.text = element_text(size = 12) 
    )
}
ggsave("rplot.pdf", plot = main, device = "pdf")

interactive_plot <- ggplot(vcf_table, aes(x = positions, y = percentage_of_reads_containing_variation, color = variation_type, 
                                          text = paste("Position:", positions, "<br>Ref allele:", ref_allele, "<br>Alt allele:", alt_allele))) +
  geom_segment(
    aes(xend = positions, yend = 0),
    linewidth = 0.5
  ) +
  geom_point(aes(size = percentage_of_reads_containing_variation)) +
  scale_size(range = c(1, 2)) +
  scale_color_manual(values = c(
    "SNP" = "#D81B60",
    "INS" = "#1E88E5",
    "DEL" = "#004D40"
    ), guide = "legend") +
  theme(legend.position = "right") +
  guides(
    color = guide_legend(title = "variation_type"),  
    size = "none"                          
  )
  

interactive_plot <- ggplotly(interactive_plot, tooltip = "text")


htmlwidgets::saveWidget(interactive_plot, "interactive_plot.html")
