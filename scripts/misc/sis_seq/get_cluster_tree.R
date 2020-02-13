# Graphs
library(igraph)

# Plotting
library(ggraph)
library(viridis)

# Data manipulation
library(tidyverse)


getEdges <- function(clusterings) {
  
  # Loop over the different resolutions
  transitions <- lapply(1:(ncol(clusterings) - 1), function(i) {
    
    # Extract two neighbouring clusterings
    from.res <- sort(colnames(clusterings))[i]
    to.res <- sort(colnames(clusterings))[i + 1]
    
    # Get the cluster names
    from.clusters <- sort(unique(clusterings[, from.res]))
    to.clusters <- sort(unique(clusterings[, to.res]))
    
    # Get all possible combinations
    trans.df <- expand.grid(FromClust = from.clusters,
                            ToClust = to.clusters)
    
    # Loop over the possible transitions
    trans <- apply(trans.df, 1, function(x) {
      from.clust <- x[1]
      to.clust <- x[2]
      
      # Find the cells from those clusters
      is.from <- clusterings[, from.res] == from.clust
      is.to <- clusterings[, to.res] == to.clust

      # Count them up
      trans.count <- sum(is.from & is.to)
      
      # Get the sizes of the two clusters
      from.size <- sum(is.from)
      to.size <- sum(is.to)
      
      # Get the proportions of cells moving along this edge
      trans.prop.from <- trans.count / from.size
      trans.prop.to <- trans.count / to.size
      
      return(c(trans.count, trans.prop.from, trans.prop.to))
    })
    
    # Tidy up the results
    trans.df$FromRes <- as.numeric(gsub("res.", "", from.res))
    trans.df$ToRes <- as.numeric(gsub("res.", "", to.res))
    trans.df$TransCount <- trans[1, ]
    trans.df$TransPropFrom <- trans[2, ]
    trans.df$TransPropTo <- trans[3, ]
    
    return(trans.df)
  })
  
  # Bind the results from the different resolutions together
  transitions <- do.call("rbind", transitions)

  # Tidy everything up
  levs <- sort(as.numeric(levels(as.factor(transitions$ToClust))))

  transitions <- transitions %>%
    mutate(FromClust = factor(FromClust,
                              levels = levs))  %>%
    mutate(ToClust = factor(ToClust, levels = levs))
  
  return(transitions)
}

getNodes <- function(clusterings) {
  nodes <- clusterings %>%
    gather(key = Res, value = Cluster) %>%
    group_by(Res, Cluster) %>%
    summarise(Size = n()) %>%
    ungroup() %>%
    mutate(Res = stringr::str_replace(Res, "res.", "")) %>%
    mutate(Res = as.numeric(Res), Cluster = as.numeric(Cluster)) %>%
    mutate(Node = paste0("R", Res, "C", Cluster)) %>%
    select(Node, everything())
}





#clust = pData(sceset_cf)[,grep("sc3_(\\d+)_c",colnames(pData(sceset_cf)),value = TRUE)]
#colnames(clust) = paste0("res.",2:10)
clust = clust[,-which(colnames(clust) == "res.10")]
clust$res.1 = 1
clust = clust[,order(colnames(clust))]
edges <- getEdges(clust)
head(edges)

nodes <- getNodes(clust)
head(nodes)


graph <- edges %>%
  # Remove edges without any cell...
  filter(TransCount > 0) %>%
  # ...or making up only a small proportion of the new cluster
  filter(TransPropTo > 0.02) %>%
  # Rename the nodes
  mutate(FromNode = paste0("R", FromRes, "C", FromClust)) %>%
  mutate(ToNode = paste0("R", ToRes, "C", ToClust)) %>%
  # Reorder columns
  select(FromNode, ToNode, everything()) %>%
  # Build a graph using igraph
  graph_from_data_frame(vertices = nodes)

print(graph)

pdf("DC_clone_split_tree.pdf",width = 13,height = 15)
ggraph(graph, layout = "tree") +
  # Plot the edges, colour is the number of cells and transparency is the
  # proportion contribution to the new cluster
  geom_edge_link(arrow = arrow(length = unit(1, 'mm')),
                 end_cap = circle(3.5, "mm"), edge_width = 1,
                 aes(colour = log(TransCount), alpha = TransPropTo)) +
  # Plot the nodes, size is the number of cells
  geom_node_point(aes(colour = factor(Res),
                      size = Size)) +
  geom_node_text(aes(label = Cluster), size = 3) +
  # Adjust the scales
  scale_size(range = c(4, 15)) +
  scale_edge_colour_gradientn(colours = viridis(100)) +
  # Add legend labels
  guides(size = guide_legend(title = "Cluster Size", title.position = "top"),
         colour = guide_legend(title = "Clustering Number",
                               title.position = "top"),
         edge_colour = guide_edge_colorbar(title = "Cell Count (log)",
                                           title.position = "top"),
         edge_alpha = guide_legend(title = "Cluster Prop",
                                   title.position = "top", nrow = 2)) +
  # Remove the axes as they don't really mean anything
  theme_void() +
  theme(legend.position = "bottom")
dev.off()