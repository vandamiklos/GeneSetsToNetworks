# GeneSetsToNetworks

A collection of functions to convert gene-set enrichment analysis results to a network graph.

## Installation

```
library(devtools)
install_github("vandamiklos/GeneSetsToNetworks")
```

## Usage

```
library(GeneSetsToNetworks)

#example gsea_result
libary(readr)
gsea_result = read_csv(url("https://raw.githubusercontent.com/vandamiklos/GeneSetsToNetworks/main/data/data.csv"))

#create the network graph
network_function(gsea_result, direction = "none", edges = "intersection", 
                 layout = "kk", cutoff_intersection = 1)
```

For help use:

```
?network_function
```

## Prerequisites:

- library(devtools)
- library(tidyverse)
- library(igraph)
- library(ggraph)

## Now you are ready to create networks such as this:
![](https://github.com/vandamiklos/GeneSetsToNetworks/blob/main/data/network.png)
