# libreries
if (!require("pacman")) install.packages("pacman")
pacman::p_load("dplyr", "ggplot2", "treemap", "treemapify")

teemapPlot <- function(Stampa = data.frame(), level1 = NA, level2 = NA, Plot = "Richness", width = 5, height = 4, output = "."){
  
  table <- Stampa %>% dplyr::group_by_at(c(level1, level2)) %>% 
    dplyr::summarise(richness = n(), abundance = sum(abundance))
  colnames(table) <- c("level1", "level2", "richness", "abundance")
  
  division_names <- table$level1
  colour <-  unique(as.character(colorpalette(division_names)))
  
  if(Plot != "Richness"){
    g <- ggplot(table, aes(area = abundance, fill = level1, label = level1, subgroup = level2))+
      geom_treemap()+
      geom_treemap_subgroup_border(colour = "black", size = 1) +
      geom_treemap_subgroup_text(inherit.aes = T, size = 10, family = "serif", 
                                 place = "centre", grow = F, alpha = 1, 
                                 colour = "White",  min.size = 1)+
      scale_fill_manual(values = colour, name = "Taxonomy")+
      theme_bw()+
      theme(
        text = element_text(family="serif"),
        legend.title = element_text(size = 8, family="serif"),
        legend.background = element_blank(),
        legend.position = "right",
        legend.box = "vertical", legend.text = element_text(size = 8, family="serif"))
    }
  else{
    g <- ggplot(table, aes(area = richness, fill = level1, label = level1, subgroup = level2))+
      geom_treemap()+
      geom_treemap_subgroup_border(colour = "black", size = 1) +
      geom_treemap_subgroup_text(inherit.aes = T, size = 10, place = "centre", grow = F, alpha = 1, colour =
                                   "White",  min.size = 1)+
      scale_fill_manual(values = colour, name = "Taxonomy")+
      theme_bw()+
      theme(
        text = element_text(family="serif"),
        legend.title = element_text(size = 8, family="serif"),
        legend.background = element_blank(),
        legend.position = "right",
        legend.box = "vertical", legend.text = element_text(size = 8, family="serif"))
  }

  # treemap::treemap(table, index = c("level1", "level2"), vSize = "richness", 
  #                  type = "index", title = "", fontsize.labels = c(14,12),
  #                  fontcolor.labels = c("black", "white"), fontface.labels = c(2,1.5),
  #                  bg.labels = c("transparent"),
  #                  align.labels = list(c("left","top"),
  #                                      c("center","center")),
  #                  inflate.labels = F,
  #                  palette = colour
  # )
  return(g)
}

#function for barplot of the similarity between DB and Data (taxonomy, pident)
densityplot_stampa <- function(table = data.frame()){
  
  table$aux <- unlist(lapply(table$taxonomy, function(x) length(grep(T, table$taxonomy == x))))
  table <- table[order(table$aux, decreasing = F),]
  table$ord <- factor(table$taxonomy, levels=names(sort(table(table$taxonomy), de = F)))
  colour = unique(as.character(colorpalette(table$taxonomy)))
  
  g <- ggplot(table, aes(x = pident, fill = ord))+
    geom_density(alpha=0.7)+
    scale_colour_manual(values = colour)+
    scale_fill_manual(values = colour, name = "Taxonomy")+
    scale_x_continuous(labels = function(x) paste0(x, "%"))+
    labs(x = "% of Similarity", y = "Density", family="serif")+
    theme_bw()+
    theme(
      text = element_text(family="serif"),
      legend.position = c(0.2,0.8),
      legend.title = element_text(size = 12, family="serif"),
      legend.background = element_blank(),
      legend.box = "vertical", legend.text = element_text(size = 10, family="serif"),
      axis.text = element_text(size = 10, colour = "black", family="serif"))
  return(g)
}