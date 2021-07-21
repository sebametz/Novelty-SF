if (!require("pacman")) install.packages("pacman")
pacman::p_load("dplyr", "RColorBrewer")

#Palette function
colorpalette <- function(x){
  # palete generation
  qual_col_pals = brewer.pal.info[brewer.pal.info$category %in% c('qual'),]
  qual_col_pals = rbind(qual_col_pals, brewer.pal.info[brewer.pal.info$category %in% c("div"),])
  col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  col_vector = unique(col_vector)
  value <- as.numeric(x)
  color <- unlist(lapply(value, function(name) col_vector[name]))
  return(color)
}


#for clusters
colorpaletteC <- function(x){
  # palete generation
  qual_col_pals = brewer.pal.info[brewer.pal.info$category %in% c('qual'),]
  qual_col_pals = rbind(qual_col_pals, brewer.pal.info[brewer.pal.info$category %in% c("div"),])
  col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  col_vector = unique(col_vector)
  value <- ifelse(as.numeric(x)>length(col_vector),(length(col_vector)-as.numeric(x)+1), 
                  (as.numeric(x)+1))
  color <- unlist(lapply(value, function(name) col_vector[name]))
  return(color)
}
