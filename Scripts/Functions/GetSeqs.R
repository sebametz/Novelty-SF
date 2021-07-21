if (!require("pacman")) install.packages("pacman")
pacman::p_load("dplyr", "Biostrings", "stringi")

# get novel sequences from the fasta file
get_seqs <- function(tax_level = NA, group = "", fasta = FASTA, nodes = novels$NODES){
  if(is.na(tax_level)){
    aux <- nodes
  }
  else{
    aux <- nodes %>% dplyr::filter(!!as.symbol(tax_level) == group)
  }
  seq <- fasta[grep(T, stri_replace_all_regex(names(fasta), "(ASV_[0-9]+)_[0-9]+", "$1") %in% as.character(aux$name))]
  return(seq)
}


# get ref_ids from the novel ASVs
get_refs <- function(tax_level = "class", group = "", nodes = data.frame(), pr2 = NA, CIM = data.frame()){
  #STAMPA references
  ref_ids <- unique(unlist(lapply(nodes$ref_ids, function(x) str_split(x, "\\,"))))
  ref_cim <- unique(unlist(lapply(CIM$ref_ids, function(x) str_split(x, "\\,"))))
  # #PR2 references
  # pr2_refs <- pr2 %>% dplyr::filter(!!as.symbol(tax_level) == group & 
  #                                     pr2$reference_sequence == 1)
  
  # #merge list
  ref_ids <- unique(c(ref_ids, ref_cim))
  # 
  #create sequences
  aux <- pr2 %>% dplyr::filter(pr2_accession %in% ref_ids)
  seq <- data.frame(name = str_c(aux$pr2_accession, aux$class, aux$species, sep = "_"), 
                    sequence = aux$sequence)
  seqs <- Biostrings::DNAStringSet(seq$sequence)
  names(seqs) <- seq$name
  # taxonomy <- paste0(aux$pr2_accession, "\t" ,str_c(aux$kingdom, aux$supergroup, aux$division, aux$class,
  #                                                  aux$order, aux$family, aux$genus, sep = ";"))
  # return(list(REFERENCES=seqs, TAXONOMY=taxonomy))
  return(seqs)
}
