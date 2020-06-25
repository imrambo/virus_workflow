#install.packages("visNetwork")
#devtools::install_github("datastorm-open/visNetwork") #for development version
require(visNetwork)
require(tidyverse)
setwd("/Users/ian/Documents/phd_research/CRISPRCas_Sediment/GuaymasC/NetAn")
nodes <- read.table("asgard.tsv", sep = "\t", header=TRUE, comment.char = '!') 
 
# nodes_united <- nodes %>%
#   mutate(id = gsub("_G", "_Guay", id)) %>%
#   tidyr::unite(from, c("id", "Taxonomy"), sep = "_", remove = FALSE) %>%
#   select(-id) %>%
#   rename(id = from)
  
edges <- read.table("spacepharer_cat_guaymas_asgard_nocrispr_fdr_0.01_2020-06-17_23-19-46.edges", header = TRUE,sep="\t") #%>%
  #rename(id = "from") %>%
  #left_join(nodes_united %>% select(from, id)) %>%
  #select(from, to)

visNetwork(nodes,edges, height = "700px", width = "100%",
           main = "Guaymas Asgards CRISPRDetect Spacer Matches",
           submain = "VIBRANT mid, high quality :: SpacePHARER\nFDR cutoff = 0.01"
) %>%
  visEdges(arrows="to") %>%
  visOptions(selectedBy = "Taxonomy", 
             highlightNearest = TRUE, 
             nodesIdSelection = TRUE,
             manipulation = TRUE) %>%
  visPhysics(stabilization = FALSE)%>%
  visConfigure(enabled = TRUE)