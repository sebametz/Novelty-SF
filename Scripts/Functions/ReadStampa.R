# function to read stampa file

# libreries
if (!require("pacman")) install.packages("pacman")
pacman::p_load("dplyr", "stringr")

# functions

readStampa <- function(readfile="", taxnames = "", taxcorrect = F){

  stampa <- read.table(file = readfile, sep = "\t", stringsAsFactors = F)
  if(taxcorrect){stampa$V4 <- taxonomy_correction(stampa$V4)}
  
  stampa.n <- cbind(stampa, str_split(stampa$V4, "\\|", simplify = T))
  stampa.n <- stampa.n[,-4]
  if(length(taxnames) == 1){taxnames = rep(1, (ncol(stampa)-4))}
  colnames(stampa.n) <- c("ASV", "abundance", "pident", "ref_ids", taxnames)
  return(stampa.n)

}

#Taxonomy correction adl et. al.,  2019
taxonomy_correction <- function(list=""){
  list.n <- str_replace_all(list, "Hacrobia\\|Cryptophyta\\|", "Cryptista|Cryptophyta|")
  list.n <- str_replace_all(list.n, "Hacrobia\\|Katablepharidophyta\\|", "Cryptista|Katablepharidophyta|")
  list.n <- str_replace_all(list.n, "Hacrobia\\|Haptophyta\\|", "Haptista|Haptophyta|")
  list.n <- str_replace_all(list.n, "Hacrobia\\|Centroheliozoa\\|", "Haptista|Centroheliozoa|")
  return(list.n)
}

#plotStampa <-