list.of.packages <- c("ggplot2", "ggtext", "patchwork", "conflicted", "rcartocolor", "plotly", "gapminder")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, dependencies = T)