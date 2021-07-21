# libreries
if (!require("pacman")) install.packages("pacman")
pacman::p_load("dplyr", "ggplot2")

# Function to get network of potential novel ASVs from specific group and 
# all the nodes connected to it accoridng to a minimal identity
#Return a list with the NETWORK, NODES and RELATIONS

get_network <- function(level = NA, group = NA, edges = 95,
                        NODES = data.frame(), RELATIONS = data.frame()){
  
  #get nodes related to group
  if(is.na(level)){aux <- NODES}
  else{aux <- NODES %>% dplyr::filter(!!as.symbol(level) == group)}
          
  #get RELATIONS of AUX
  RELATIONS.n <- RELATIONS %>% dplyr::filter(from %in% aux$name | 
                                               to %in% aux$name)
  RELATIONS.n <- RELATIONS.n %>% dplyr::filter(identity >= edges)
  
  #get complete list of NODES
  to_select <- unique(c(as.character(aux$name), as.character(RELATIONS.n$from), 
                        as.character(RELATIONS.n$to)))
  NODES.n <- NODES %>% dplyr::filter(name %in% to_select)
  # print(to_select[grep(F, to_select %in% NODES.n$name)])
  # return network
  net <- graph_from_data_frame(RELATIONS.n, directed=FALSE, vertices=NODES.n)
  
  #get membership and add to NODES
  #names of membership
  cl <- cluster_label_prop(net)
  mem <- membership(cl)
  mem <- data.frame(name=names(mem), Cluster = as.numeric(unname(mem)))
  NODES.n <- merge(NODES.n, mem, all.x = T)
  
  #create a list of results and return
  list <- list(NETWORK = net, NODES = NODES.n, RELATIONS = RELATIONS.n)
  return(list)  
}


# Prepare data for plots
prepare_for_plot <- function(network = list(), color_level = "C"){
  #create file x1 y1 x2 y2 ASV level_color size
  relations <- network$RELATIONS
  nodes <- network$NODES
  
  relations.aux <- relations
  not_included <- relations[grep(F, relations$to %in% relations$from), ]
  not_included <- data.frame(from = not_included$to, to = not_included$from, 
                        identity = not_included$identity)
  relations <- rbind(relations.aux, not_included)
  
  # To itself
  relations.aux <- relations
  not_included <- nodes[grep(F, nodes$name %in% relations$from), ]
  not_included <- data.frame(from = not_included$name, to = not_included$name, 
                             identity = 100)
  relations <- rbind(relations.aux, not_included)
  
  
  #from
  x1 <- unlist(lapply(as.character(relations$from), function(x) nodes$CEM[nodes$name == x]))
  y1 <- unlist(lapply(as.character(relations$from), function(x) nodes$CIM[nodes$name == x]))
  
  #To
  x2 <- unlist(lapply(as.character(relations$to), function(x) nodes$CEM[nodes$name == x]))
  y2 <- unlist(lapply(as.character(relations$to), function(x) nodes$CIM[nodes$name == x]))
  
  abu <- unlist(lapply(as.character(relations$from), function(x) nodes$abundance[nodes$name == x]))
  
  division_names <- unlist(lapply(as.character(relations$from), 
                               function(x) nodes[nodes$name == x, 
                                                 which(colnames(nodes) == color_level)]))
  
  
  cluster <- unlist(lapply(as.character(relations$from), function(x) nodes$Cluster[nodes$name == x]))
  
  aux <- data.frame(X1 = x1, Y1 = y1, ASV = as.character(relations$from), 
                    X2 = x2, Y2 = y2, abundance = abu, group = division_names, 
                    color = colorpalette(division_names), Cluster = cluster)
  return(aux)
}

#plot clusters novelty
prepare_for_plot_c <- function(network = list(), color_level = "Cluster"){
  #create file x1 y1 x2 y2 ASV level_color size
  relations <- network$RELATIONS
  nodes <- network$NODES
  
  relations.aux <- relations
  not_included <- relations[grep(F, relations$to %in% relations$from), ]
  not_included <- data.frame(from = not_included$to, to = not_included$from, 
                             identity = not_included$identity)
  relations <- rbind(relations.aux, not_included)
  
  # To itself
  relations.aux <- relations
  not_included <- nodes[grep(F, nodes$name %in% relations$from), ]
  not_included <- data.frame(from = not_included$name, to = not_included$name, 
                             identity = 100)
  relations <- rbind(relations.aux, not_included)
  
  
  #from
  x1 <- unlist(lapply(as.character(relations$from), function(x) nodes$CEM[nodes$name == x]))
  y1 <- unlist(lapply(as.character(relations$from), function(x) nodes$CIM[nodes$name == x]))
  
  #To
  x2 <- unlist(lapply(as.character(relations$to), function(x) nodes$CEM[nodes$name == x]))
  y2 <- unlist(lapply(as.character(relations$to), function(x) nodes$CIM[nodes$name == x]))
  
  abu <- unlist(lapply(as.character(relations$from), function(x) nodes$abundance[nodes$name == x]))
  
  division_names <- unlist(lapply(as.character(relations$from), 
                                  function(x) nodes[nodes$name == x, 
                                                    which(colnames(nodes) == color_level)]))
  
  cluster <- unlist(lapply(as.character(relations$from), function(x) nodes$Cluster[nodes$name == x]))
  
  aux <- data.frame(X1 = x1, Y1 = y1, ASV = as.character(relations$from), 
                    X2 = x2, Y2 = y2, abundance = abu, group = as.character(division_names), 
                    color = colorpaletteC(division_names), Cluster = cluster)
  return(aux)
}


# Plot network
plot_network <- function(table = data.frame(), plot_lines = T, plot_mem = F, mcim = NA, mcem =NA, mcimg = NA, mcemg =NA, Title = NULL){
  aux <- table
  division_names <- aux$group
  color = unique(as.character(colorpalette(division_names)))
  #ploting
  g <- ggplot(aux, aes(x = X1, y = Y1, xend = X2, yend = Y2)) +
    geom_edges(color = "grey50", size = 0.3) + 
    geom_nodes(aes(fill = group), shape=21) + #, size = 0.1
    scale_colour_manual(values = color)+
    scale_fill_manual(values = color, name = "")+
    scale_x_continuous(labels = function(x) paste0(x, "%"), limits = c(floor(min(c(aux$X1,aux$X2))/5)*5,100))+
    scale_y_continuous(labels = function(x) paste0(x, "%"), limits = c(floor(min(c(aux$Y1,aux$Y2))/5)*5,100))

  #With size
  # g <- ggplot(aux, aes(x = X1, y = Y1, xend = X2, yend = Y2)) +
  #   geom_edges(color = "grey50") + 
  #   geom_nodes(aes(fill = group, size = abundance), shape=21) +
  #   scale_colour_manual(values = color)+
  #   scale_fill_manual(values = color)+
  #   scale_size_continuous(range = c(1, 10))
  
  #plot novelty lines?
  if(plot_lines){
    #mccm and mcem
    # mccm <- mean(nodes$CIM)
    # mcem <- mean(nodes$CEM)
    g <- g + 
      geom_hline(yintercept=mcim, linetype = "dashed", color = "red", cex = 0.5)+
      geom_vline(xintercept=mcem, linetype = "dashed", color = "red", cex = 0.5)
    if(!is.na(mcimg) & !is.na(mcemg)){
      g <- g + 
        geom_hline(yintercept=mcimg, linetype = "dashed", color = "gray70", cex = 0.5, alpha = 0.8)+
        geom_vline(xintercept=mcemg, linetype = "dashed", color = "gray70", cex = 0.5, alpha = 0.8)
    }
  }
  
  #plot cluster names?
  if(plot_mem){
    g <- g +
      geom_nodetext(aes(label = Cluster), family="serif")
  }
  
  #Format plot
  g <- g + 
      labs(x = NULL, y = NULL, title = Title)+ #x = "CEM", y = "CIM"
      guides(size="none", fill = guide_legend(override.aes = list(size = 1.5)))+
      theme_bw()+
      theme(
        panel.grid.minor = element_blank(),
        plot.title = element_text(family="serif", colour = "black", hjust=0.5),
        text=element_text(family="serif", colour = "black"), #, size = 12
        legend.position = c(0.25,0.8),
        legend.title = element_text( family="serif"), #size = 12,
        legend.background = element_blank(),
        legend.box.background = element_blank(),
        legend.box = "vertical", #legend.text = element_text(size = 14), 
        legend.spacing.x = unit(0.5, "cm"), legend.key.size = unit(0, "lines"),
        axis.text = element_text(colour = "black", family="serif")) #size = 12,
  
  return(g)
}


#function to filter potential novel ASVs and their related network
get_novels <- function(level = NA, group = NA, network = list(), table = ASVTABLE, 
                       stat = T, maxCIM = NA, maxCEM = NA){
  
  if(is.na(level)){aux <- network$NODES}
  else{aux <- network$NODES %>% dplyr::filter(!!as.symbol(level) == group)}
  
  #filter the novel ASVs
  potential_novels <- aux %>% dplyr::filter(CIM <= maxCIM & CEM <= maxCEM)
  
  relations.n <- network$RELATIONS %>% filter(from %in% potential_novels$name | 
                                              to %in% potential_novels$name)

  #get complete list of nodes
  to_select <- unique(c(as.character(potential_novels$name), as.character(relations.n$from), 
                        as.character(relations.n$to)))
  
  nodes.n <- network$NODES %>% dplyr::filter(name %in% to_select)
  
  #table of ASVs
  table <- data.frame(ASV=rownames(table), table)
  novel_table <- table %>% dplyr::filter(ASV %in% nodes.n$name)
  
  
  #novel network
  net <- graph_from_data_frame(relations.n, directed=FALSE, vertices=nodes.n)
  
  #STATISTIC (if the difference between CIM and CEM is significant i.e there is a gain of knowladge)
  if(!is.na(level) & stat == T){
    stat.table <- get_stat(nodes.n, level)
    print(stat.table)
  }
  
  #create a list of results and return
  list <- list(NETWORK = net, NODES = nodes.n, RELATIONS = relations.n, 
               NOVEL = novel_table)
  return(list)
}


get_stat <- function(x, level = NA){
  meanCIM <- mean(x$CIM)
  meanCEM <- mean(x$CEM)
  stat.table <- x %>% dplyr::group_by_at(level) %>%
    summarise(n = n(), mCIM = mean(CIM), seCIM =  sd(CIM)/sqrt(n), 
              mCEM = mean(CEM), seCEM = sd(CEM)/sqrt(n), 
              p.value = ifelse(n > 1 & mCIM != mCEM, round(t.test(CIM, CEM, var.equal = T, conf.level = 0.99)$p.value, digits = 4), NA),
              medCIM = median(CIM), medCEM = median(CEM), 
              Q1 = length(grep(T, CIM > meanCIM & CEM > meanCEM)),
              Q2 = length(grep(T, CIM > meanCIM & CEM <= meanCEM)),
              Q3 = length(grep(T, CIM <= meanCIM & CEM <= meanCEM)),
              Q4 = length(grep(T, CIM <= meanCIM & CEM > meanCEM)),
              NClusters = n_distinct(Cluster))
  return(stat.table)
}


# #function to filter potential novel ASVs and their related network
# get_env <- function(level = NA, group = NA, network = network, table = ASVTABLE){
#  
#   if(is.na(level)){aux <- network$NODES}
#   else{aux <- network$NODES %>% dplyr::filter(!!as.symbol(level) == group)}
#   
#   # mean ci (cultured and isolated) and env
#   minCEM <- mean(aux$CEM)
#   minCIM <- mean(aux$CIM)
#   
#   #filter the novel ASVs for the ASV table
#   potential_novels <- aux %>% dplyr::filter(CIM <= minCIM & CEM > minCEM)
#   table <- data.frame(ASV=rownames(table), table)
#   novel_table <- table %>% dplyr::filter(ASV %in% potential_novels$name)
#   
#   relations.n <- network$RELATIONS %>% filter(from %in% novel_table$ASV, 
#                                               to %in% novel_table$ASV)
#   
#   #get complete list of nodes
#   to_select <- c(as.character(novel_table$ASV), as.character(relations.n$from), as.character(relations.n$to))
#   nodes.n <- aux %>% dplyr::filter(name %in% to_select)
#   
#   #novel network
#   net <- graph_from_data_frame(relations.n, directed=FALSE, vertices=nodes.n)
#   
#   #create a list of results and return
#   list <- list(NETWORK = net, NODES = nodes.n, RELATIONS = relations.n, NOVEL = novel_table)
#   return(list)
# }