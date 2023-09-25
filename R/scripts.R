#' Shows whether a gene is in a gene-set
#'
#' Takes a gene-set enrichment analysis result as input and creates a 
#' Data Frame that holds every gene-set provided in the input (columns) and 
#' every gene (rows) that is a part of these gene-sets. When a gene is part of a 
#' gene-set the corresponding cell contains 1, when the gene isn't part of the 
#' gene-set it contains 0.
#' @param gsea_result A Data Frame of a gene-set enrichment analysis result. Required columns: 
#' Description, NES, core_enrichment.
#' @export
make_set_gene = function(gsea_result) {
  #select the description and core_enrichment columns
  genes = gsea_result[, c("Description", "core_enrichment")]
  #split the core_enrichment column by / to get the separate gene ids and put them in different rows
  genes =  genes %>% mutate(core_enrichment = strsplit(as.character(core_enrichment), "/")) %>% 
    unnest(core_enrichment)
  
  #create an empty data frame
  #it will hold the gene ids as rows, pathways as columns
  set_gene = data.frame(matrix(ncol = 1, nrow = 0))
  names(set_gene) = "core_enrichment"
  
  #when a gene is present in the pathway put 1 into the cell, otherwise leave it empty (NA)
  for (i in gsea_result$Description) {
    col = genes[genes$Description == i, c("core_enrichment")]  #col holds the genes for the i-th description
    name = as.character(i)
    col$a = 1 # add a column with 1-s to col
    names(col) = c("core_enrichment", name) #change the column name from a to the i-th desciption
    set_gene = merge(set_gene, col, by = "core_enrichment",
                     all.x = TRUE, all.y = TRUE) #merge set_gene and col by core_enrichment
  }
  
  #change NA values to 0
  set_gene[is.na(set_gene)] = 0
  #set the core_enrichment column as row names
  row.names(set_gene) = set_gene$core_enrichment
  set_gene = set_gene[, 2:ncol(set_gene)]
  
  return(set_gene)
}



#' Calculates the number of intersections
#'
#' Takes a gene-set enrichment analysis result as input and calculates the 
#' number of intersections between each pair of gene-sets
#' @param set_gene A Data Frame object, output of the make_set_gene() function
#' @export
make_intersection_table = function(gsea_result) {
  
  #make set_gene table
  set_gene = make_set_gene(gsea_result)
  
  #make empty data frame
  intersection_table = data.frame(matrix(nrow = ncol(set_gene), ncol = ncol(set_gene)))
  
  #create a table where both cols and rows are the gene-sets and the values are 
  #the number of genes that are present in both gene-sets
  for (x in 1:ncol(set_gene)) {
    for (y in 1:ncol(set_gene)) {
      
      #add together y-th and x-th gene-sets
      #0 - that gene isn't in either of the gene-sets
      #1 - its in only one of the gene-sets
      #2 - it's in both of the gene-sets
      summed_col = set_gene[, y] + set_gene[, x] 
      
      #sum together the cells with 2 in them (the genes that are in both gene-sets)
      #with this I get the number of genes in each gene-sets too (diagonal of the matrix)
      intersection_table[y, x] = sum(summed_col == 2) 
    }
  }
  return(intersection_table)
}



#' Calculates Jaccard similarity between each pair of gene-sets
#'
#' Takes a gene-set enrichment analysis result as input and calculates the 
#' Jaccard similarity of every gene-set pair.
#' @param set_gene A Data Frame object, output of the make_set_gene() function
#' @return Returns a Data Frame object. Each column and row represents a gene-set
#' and values are the Jaccard similarity indexes of two gene-sets.
#' @export
make_jaccard_similarity_table = function(gsea_result){
  
  #create the intersection table
  intersection_table = make_intersection_table(gsea_result)
  
  #create empty data frame
  jaccard_table = data.frame(matrix(nrow = ncol(intersection_table), ncol = ncol(intersection_table)))
  
  #fill it up with jaccard similarity values
  for (x in 1:nrow(intersection_table)) {
    for (y in 1:ncol(intersection_table)) {
      #jaccard similarity: intersection / union
      jaccard_table[x,y] = intersection_table[x,y] / (intersection_table[x,x] + intersection_table[y,y] - intersection_table[x,y])
    }
  }
  return(jaccard_table)  
}



#' Make the edges of the graph
#' 
#' Creates a table that holds information about the edges of the graph. We can 
#' specify "intersect" or "jaccard similarity" as the value that will determine 
#' the width of the connecting lines on the graph.
#' @param df Output of make_jaccard_similarity_table() or make_intersection_table()
#' @param name "jaccard similarity" or "intersect"
#' @param cutoff Cutoff value of intersections/jaccard similarity index
#' @param gsea_result A Data Frame of a gene-set enrichment analysis result. Required columns: 
#' Description, NES, core_enrichment.
#' @return Returns a table that defines the edges of the graph. It holds the 
#' names of the two nodes an edge will connect and a value 
#' - intersection or jaccard similarity - that will determine the line_width
#' @export
make_edge = function(df, name, cutoff, gsea_result) {
  #name the columns according to description
  colnames(df) = gsea_result$Description
  
  #add a gene-sets column with the descriptions
  df$`gene-set` = gsea_result$Description
  

  edges = df %>%
    pivot_longer(!`gene-set`, names_to = "to", values_to = name)
  
  if (name == "intersection") {
    #delete gene-set pairs that have 0 genes in common - don't need edges 
    #between sets that don't have any gene sin common
    edges = edges[!(edges$intersection == 0),]
    #delete gene-set pairs under a specified cutoff value
    edges = edges[edges$intersection > cutoff,]
  }
  
  else {
    #delete gene-set pairs with jaccard sim. equal to 1, because these are identical sets
    edges = edges[!(edges$`jaccard similarity` == 1),]
    #delete gene-set pairs that are below the cutoff
    edges = edges[edges$`jaccard similarity` > cutoff,]
  }
  return(edges)
}



#' Draw the network graph
#' 
#' Creates a network of the gene-sets.
#' @param graph Graph object created by ggraph
#' @param layout Layouts types available in the ggraph package or created with create_layout()
#' @param low Color of the gene-sets with negative NES. Default color:"#1465AC"
#' @param mid Color of the gene-sets with zero NES. Default color: "white"
#' @param high Color of the gene-sets with positive NES. Default color: "#B31B21"
#' @param edges Table with information about the edges. Output of make_edge().
#' @param gsea_result A Data Frame of a gene-set enrichment analysis result. Required columns: 
#' Description, NES, core_enrichment.
#' @param name Can be either "jaccard similarity" or "intersect"
#' @return Returns a ggplot object onto which additional layers can be added.
#' @export
draw_graph = function(graph, layout, edges, name, gsea_result,
                      low = low, mid = mid, high = high) {
  
  #create table that holds info about the nodes
  #description, setSize(number of genes in a gene-set), NES(normalized enrichment score) and p adjusted value
  nodes = gsea_result[, c("Description", "setSize", "NES")]
  
  #create the graph object from Data Frame (ggraph function)
  graph = graph_from_data_frame(edges, vertices = nodes)
  
  plot = ggraph(graph, layout = layout) + 
         geom_edge_link(alpha = 0.2, aes(edge_width = unlist(edges[,3]))) +
         geom_node_point(aes(size = setSize, color = NES)) +
         scale_size(range = c(1, 10))  +
         geom_node_text(aes(label = nodes$Description), size = 3, repel = TRUE) +
         scale_edge_width(range = c(0.2, 2)) + 
         labs(edge_width = name, color = "NES", size = "Set size") +
         theme_void() + 
         scale_color_gradient2(low      = low,
                               mid      = mid,
                               high     = high,
                               midpoint = 0,
                               na.value = "grey80",
                               limits   = c(-4, 4),
                               guide    = "colourbar")
  
  return(plot)
}



#' Create networks from gene-sets
#' 
#' Expects a gene-set enrichment analysis result and creates a network graph from it.
#' The edges of the graph can be either based on the intersection of genes in each set
#' or Jaccard similarity.
#' Holds together all the functions from the package. 
#' 
#' @param gsea_result A Data Frame of a gene-set enrichment analysis result. Required columns: 
#' Description, NES, core_enrichment.
#' @param direction "none" (default) if we want gene-sets with both positive and negative normalized enrichment scores (NESs), "up" if we want gene-sets with only positive NESs, and "down" if we want gene-sets with only negative NESs
#' @param edges defaults to "intersection" otherwise Jaccard similarity
#' @param layout Layout types available in the ggraph package or created with create_layout()
#' @param low Color of the gene-sets with negative NES. Default color: "#1465AC"
#' @param mid Color of the gene-sets with zero NES. Default color: "white"
#' @param high Color of the gene-sets with positive NES. Default color: "#B31B21"
#' @export
network_function = function(gsea_result, direction = "none", edges = "intersection", 
                            layout, cutoff_intersection = 0, cutoff_jaccard = 0.7,
                            low = "#1465AC", mid = "white", high = "#B31B21"){
  
  #choose NES direction
  if (direction == "up") {
    gsea_result = gsea_result[gsea_result$NES > 0,]
  } 
  else if (direction == "down") {
    gsea_result = gsea_result[gsea_result$NES < 0,]
  } 
  else {gsea_result = gsea_result}
  
  #create edges based on intersection or jaccard similarity index
  if (edges == "intersection") {
    
    #variables that are needed for the functions
    name = "intersection"
    cutoff = cutoff_intersection
    
    #create table containing the number of genes that are intersecting in two gene-sets
    df = make_intersection_table(gsea_result)
    
    #create the edges table
    edges = make_edge(df, name, cutoff, gsea_result)
  }
  else {
    
    #variables that are needed for the functions
    name = "jaccard similarity"
    cutoff = cutoff_jaccard
    
    #create table containing the jaccard similarity of two gene-sets
    df = make_jaccard_similarity_table(gsea_result)
    
    #create the edges table
    edges = make_edge(df, name, cutoff, gsea_result)
  }
  

  
  #create the final graph
  plot = draw_graph(graph, layout, edges, name, gsea_result, 
                    low = low,
                    mid = mid,
                    high = high)
  return(plot)
}
