# GeneSetsToNetworks

A collection of functions to convert gene-set enrichment analysis results to a network graph.

## Install

```
library(devtools)
install_github("vandamiklos/GeneSetsToNetworks")
```

## Usage

```
library(GeneSetsToNetworks)

network_function(gsea_result, direction = "none", edges = "intersection", 
                 layout, cutoff_intersection = 0, cutoff_jaccard = 0.7,
                 low = "#1465AC", mid = "white", high = "#B31B21")
```

## Prerequisites:

- library(devtools)
- library(tidyverse)
- library(igraph)
- library(ggraph)

## Now you are ready to create networks such as this:
![](https://github.com/vandamiklos/GeneSetsToNetworks/blob/main/data/network.png)
