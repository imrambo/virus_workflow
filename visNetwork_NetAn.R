#install.packages("visNetwork")
#devtools::install_github("datastorm-open/visNetwork") #for development version
require(visNetwork)
require(tidyverse)
require(RColorBrewer)
require(stringr)
require(ggplot2)

create.network <- function(n, e, tmain, tsub, select_by, h = "700px", w = "100%") { 
    ### Visualize a network using visNetwork
    visNetwork(n, e, height = h, width = w,
               main = tmain,
               submain = tsub) %>%
      visEdges(arrows="to") %>%
      visOptions(selectedBy = select_by, 
                 highlightNearest = TRUE, 
                 nodesIdSelection = TRUE,
                 manipulation = TRUE) %>%
      visPhysics(stabilization = FALSE) %>%
      visConfigure(enabled = TRUE) }

node_shapes <- c("square", "triangle", "box",
                 "dot", "star", "ellipse", "diamond")

#Qualitative palette colors
qual_col_pals <- brewer.pal.info[brewer.pal.info$category == 'qual',]

qual_col_vec <- unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

node_col_shape <- data.frame(qual_col_vec, rep(node_shapes, length.out = length(qual_col_vec)))
colnames(node_col_shape) <- c("color", "shape")
#=============================================================================
### Read in data

setwd("/Users/ian/Documents/phd_research/CRISPRCas_Sediment/GuaymasC")

#Generate master taxonomy mapping file
# archaea_taxonomy <- read.table("./GuaymasC_Archaea_Taxonomy_Master.txt", header = TRUE, sep = "\t") %>%
#   mutate(Domain = "Archaea")
# bacteria_taxonomy <- read.table("./GuaymasC_Bacteria_Taxonomy_Master.txt", header = TRUE, sep = "\t") %>%
#   mutate(Domain = "Bacteria")
# 
# gb19_missing_taxonomy <- read.csv("./Guaymas19_Bacteria_Archaea_Taxonomy.csv", header = TRUE) %>%
#   filter(!Bin %in% archaea_taxonomy$Bin & !Bin %in% bacteria_taxonomy$Bin) %>%
#   select(Bin, Taxonomy, Domain)
# 
# archaea_bacteria_taxonomy <- dplyr::bind_rows(archaea_taxonomy, bacteria_taxonomy, gb19_missing_taxonomy)
# 
# write.table(archaea_bacteria_taxonomy, file = "./GuaymasC_Archaea_Bacteria_Taxonomy_Master.txt",
#             quote = FALSE, col.names = TRUE, row.names = FALSE, fileEncoding = "UTF-8",
#             sep = "\t")

#Taxonomy mapping file 
archaea_bacteria_taxonomy <- read.table("/Users/ian/Documents/phd_research/CRISPRCas_Sediment/GuaymasC/GuaymasC_Archaea_Bacteria_Taxonomy_Master.txt",
                                        sep = "\t", header = TRUE)
#=============================================================================
#CRISPR Spacer :: VIBRANT virus BLASTN hits
blastn_tbl <- read.table("/Users/ian/Documents/phd_research/CRISPRCas_Sediment/GuaymasC/GuaymasC_virus_blastn_spacer_output/GuaymasC_CRISPRDetect_output_qt3_rl23_spacers___GuaymasC_VIBRANT_extracted_phages_all_high_medium.out",
                         header = FALSE, col.names = c("qseqid", "sseqid",
                                                       "pident", "alength",
                                                       "mismatch", "gapopen",
                                                       "qstart", "qend",
                                                       "sstart", "send",
                                                       "evalue", "bitscore",
                                                       "qlen", "slen"))

#Best BLASTN hits
blastn_bh <- blastn_tbl %>%
  #filter(qlen == alength & mismatch <= 1 & pident >= 95 & evalue <= 1e-5) %>%
  filter(qlen == alength & mismatch <= 1 & evalue <= 1e-4) %>%
  rename(crispr_spacer = qseqid, virus = sseqid) %>%
  mutate(Bin = gsub("_scaffold.*$|_redone.*$|_NODE.*$", "", crispr_spacer),
         Bin = gsub("Guay", "G", Bin),
         color = "skyblue") %>%
  left_join(archaea_bacteria_taxonomy) %>% 
  distinct()

blastn_asgard <- blastn_bh %>%
  dplyr::filter(grepl("Thor|Odin|Loki|Heimdall|Helarch|Sifarch|Freyarch", Taxonomy)) %>%
  rename(virus_len = slen) %>%
  select(-color)

asgard_blastn_output <- "/Users/ian/Documents/phd_research/CRISPRCas_Sediment/GuaymasC/GuaymasC_virus_blastn_spacer_output/Asgard_GuaymasC_virus_blastn_spacer_output.txt"

write.table(blastn_asgard, asgard_blastn_output, quote = FALSE,
            col.names = FALSE, row.names = FALSE,
            fileEncoding = "UTF-8", sep = "\t")

blastn_bh$color[which(blastn_bh$pident == 100)] <- "indianred"
blastn_bh$color[which(blastn_bh$pident >= 95 & blastn_bh$pident < 100)] <- "purple"
blastn_bh$color[which(blastn_bh$pident < 95)] <- "skyblue"

#=============================================================================
#CRISPR Spacer :: VIBRANT virus SpacePHARER hits
spacepharer_colnames <- c("crispr_spacer", "virus")

#Genbank 2018_09
spacepharer_allbin_genbank <- read.table("./spacepharer_output/GuaymasC_AllGenomes_CRISPRDetect_output_qt3_rl23_queryDB---GenbankPhage2018_09_target_fdr_0.01_2020-06-25_21-37-42_uniq_hits.tsv",
                                         header = FALSE, sep = "\t",
                                         col.names = spacepharer_colnames) %>%
  mutate(Bin = gsub("_scaffold.*$|_redone.*$|_NODE.*$", "", crispr_spacer),
         Bin = gsub("Guay", "G", Bin)) 

#VIBRANT GuaymasC virus draft genomes
spacepharer_allbin_vibrant <- read.table("./spacepharer_output/GuaymasC_AllGenomes_CRISPRDetect_output_qt3_rl23_queryDB---GuaymasC_VIBRANT_extracted_phages_high_medium_lytic_lysogenic_fdr_0.01_2020-06-25_23-53-19_uniq_hits.tsv",
                                         header = FALSE, sep = "\t",
                                         col.names = spacepharer_colnames) %>%
  mutate(Bin = gsub("_scaffold.*$|_redone.*$|_NODE.*$", "", crispr_spacer),
         Bin = gsub("Guay", "G", Bin))

#Eukaryotic viruses - false positive matches
spacepharer_allbin_eukvir <- read.table("./spacepharer_output/GuaymasC_AllGenomes_CRISPRDetect_output_qt3_rl23_queryDB---genbank_eukvir_2018_09_targetDB_fdr_0.01_2020-06-25_22-35-17_uniq_hits.tsv",
                                        header = FALSE, col.names = spacepharer_colnames)

#Asgard :: VIBRANT, Genbank hits from first analysis
spacepharer_asgard_genbank_vibrant <- read.table("./spacepharer_output/spacepharer_cat_guaymas_asgard_fdr_0.01_2020-06-17_23-19-46.tsv",
                                                 header = FALSE, sep = "\t",
                                                 col.names = spacepharer_colnames) %>%
  mutate(Bin = gsub("_scaffold.*$|_redone.*$|_NODE.*$", "", crispr_spacer),
         Bin = gsub("Guay", "G", Bin))
#-----------------------------------------------------------------------------
### MEBS output for viruses
#Non-normalized scores, just for viruses
guaymas_vibrant_mebs <- read.table("/Users/ian/Documents/phd_research/CRISPRCas_Sediment/GuaymasC/MEBS/IanVIRUS.tsv_completenes.tab",
                                   header = TRUE, sep = "\t")
#Long format non-normalized scores
guaymas_vibrant_mebs_lf <- guaymas_vibrant_mebs %>% tidyr::gather("pathway", "mebs_score", 2:ncol(guaymas_vibrant_mebs)) %>%
  rename(virus = X) %>%
  mutate()

#Normalized MEBS scores - virus scores normalized to reference microbial genomes
guaymas_mebs_norm <- read.table("/Users/ian/Documents/phd_research/CRISPRCas_Sediment/GuaymasC/MEBS/IanVIRUS.tsv_2_cluster_mebs.txt",
                                header = TRUE, sep = "\t")
#-----------------------------------------------------------------------------
#Dataframe containing all spacer hit and taxonomy information
spacepharer_allbin_genbank_vibrant_blastn <- dplyr::bind_rows(spacepharer_allbin_genbank,
                                                              spacepharer_allbin_vibrant,
                                                              spacepharer_asgard_genbank_vibrant) %>%
  dplyr::distinct() %>%
  dplyr::filter(!crispr_spacer %in% spacepharer_allbin_eukvir$crispr_spacer) %>%
  dplyr::left_join(archaea_bacteria_taxonomy) %>%
  dplyr::left_join(blastn_bh %>% select(crispr_spacer, virus, blastn_hit)) %>%
  dplyr::mutate(blastn_hit = tidyr::replace_na(blastn_hit, -1))


#Asgard-Bacteria-Virus links
#VIBRANT Virus hits to Asgards
asgard_vibrant_hits <- spacepharer_allbin_genbank_vibrant_blastn %>%
  dplyr::filter(grepl("Thor|Odin|Loki|Heimdall|Helarch|Sifarch|Freyarch", Taxonomy)
                & grepl("scaffold", virus))

# #VIBRANT Viruses that hit Asgards, but also other taxa
# asgard_vibrant_otherhits <- spacepharer_allbin_genbank_vibrant_blastn %>%
#   dplyr::filter(!grepl("Thor|Odin|Loki|Heimdall|Helarch|Sifarch|Freyarch", Taxonomy) &
#                   virus %in% asgard_vibrant_hits$virus)
# 
# asgard_vibrant_hits_alltax <- dplyr::bind_rows(asgard_vibrant_hits, asgard_vibrant_otherhits) %>%
#   dplyr::arrange(virus)
#   
# asgard_vibrant_hits_out <- "/Users/ian/Documents/phd_research/CRISPRCas_Sediment/GuaymasC/GuaymasC_Asgard_Virus_TaxaHits.txt"
# 
# write.table(asgard_vibrant_hits_alltax,
#             asgard_vibrant_hits_out,
#             quote = FALSE, col.names = TRUE,
#             row.names = FALSE, fileEncoding = "UTF-8",
#             sep = "\t")
# 
# #Geom_tile of viruses that hit an Asgard and other taxa
# #Rudimentary, intial look
# asgard_vibrant_hits_alltax %>% ggplot(aes(x = Taxonomy, y = virus,
#                                           width = 0.7, height = 0.7, fill = as.factor(blastn_hit))) + 
#   geom_tile() +
#   theme(panel.grid.major = element_line(color="gray"),
#                   panel.grid.minor = element_line(color="black"),
#                   #panel.background = element_blank(),
#                   axis.line = element_line(color="black"),
#                   axis.title=element_text(size=10),
#                   axis.text.x=element_text(color = "black", angle=70, hjust =1),
#                   axis.text.y=element_text(size = 7))
# 
# asgard_bacteria_virus_counts <- asgard_vibrant_hits_alltax %>% filter(Domain == "Bacteria") %>%
#   group_by(Taxonomy) %>% count()
# 
# #Viruses that infect an Asgard and a Bacteria
# for(btax in asgard_bacteria_virus_counts$Taxonomy){
#   avhit = asgard_vibrant_hits_alltax %>%
#     filter(grepl("Helarch|Lokiarch|Thorarch|Odinarch|Sifarch|Heimdall|Freyarch", Taxonomy) &
#              virus %in% asgard_vibrant_hits_alltax$virus[which(asgard_vibrant_hits_alltax$Taxonomy == btax)]) %>%
#     #select(Bin, Taxonomy) %>%
#     distinct() %>%
#     arrange(Taxonomy)
#   print(btax)
#   print.data.frame(avhit)
#     
# }

#=============================================================================
setwd("/Users/ian/Documents/phd_research/CRISPRCas_Sediment/GuaymasC/NetAn")

#Network setup and visualization for SpacePHARER and BLASTN hits

#Input file for NetAn undirected network
netan_input_undirected <- spacepharer_allbin_genbank_vibrant_blastn %>%
  dplyr::select(Bin, virus) %>%
  dplyr::distinct()

netan_iu <- paste("/Users/ian/Documents/phd_research/CRISPRCas_Sediment/GuaymasC/NetAn/GuaymasC_AllBin_Undirected_", Sys.Date(), ".txt", sep = '')
write.table(netan_input_undirected,
            netan_iu,
            quote = FALSE, col.names = FALSE, row.names = FALSE, fileEncoding = "UTF-8",
            sep = "\t")

#Output from virus_workflow/extract_netan_groups.py
cc_groups <- read.table("./GuaymasC_AllBin_Undirected_2020-07-01_connected-components.txt",
                        col.names = c("Group", "Bin")) %>%
  dplyr::select(Bin, Group)


nodes <- cc_groups %>%
  left_join(spacepharer_allbin_genbank_vibrant_blastn %>% 
              dplyr::select(Bin, Taxonomy)) %>%
  rename(id = Bin)


nodes$Taxonomy[which(is.na(nodes$Taxonomy) & grepl(".*?\\.[[:digit:]]", nodes$id))] <- "Genbank Phage"
nodes$Taxonomy[which(is.na(nodes$Taxonomy) & grepl("scaffold_2kb", nodes$id))] <- "Guaymas Phage"

#Create data frame of node Taxonomic id, and corresponding shape + color
node_taxonomy_uniq <- unique(nodes$Taxonomy)

node_taxonomy_col_shape <- data.frame(node_taxonomy_uniq,
                                      node_col_shape$color[1:length(node_taxonomy_uniq)],
                                      node_col_shape$shape[1:length(node_taxonomy_uniq)])
colnames(node_taxonomy_col_shape) <- c("Taxonomy", "color", "shape")

nodes <- nodes %>% left_join(node_taxonomy_col_shape) %>%
  dplyr::distinct() %>%
  tidyr::unite("id_tax", c("id", "Taxonomy"),
               sep = "; ", remove = FALSE)


edges <- spacepharer_allbin_genbank_vibrant_blastn %>%
  dplyr::select(Bin, virus, blastn_hit, Taxonomy) %>%
  dplyr::distinct() %>%
  rename(from = Bin, to = virus)

edges_from <- edges %>% select(from, blastn_hit) %>%
  rename(id = from) %>%
  left_join(nodes %>% select(id, id_tax)) %>%
  dplyr::select(id_tax) %>%
  dplyr::rename(from = id_tax)

edges_to <- edges %>% select(to, blastn_hit) %>%
  rename(id = to) %>%
  left_join(nodes %>% select(id, id_tax)) %>%
  dplyr::select(id_tax, blastn_hit) %>%
  dplyr::rename(to = id_tax)

edges <- dplyr::bind_cols(edges_from, edges_to)

edges$blastn_hit[which(edges$blastn_hit == 1)] <- "indianred"
edges$blastn_hit[which(edges$blastn_hit == -1)] <- "skyblue"

edges <- edges %>% dplyr::rename(color = blastn_hit) %>% distinct()

nodes <- nodes %>% dplyr::select(id_tax, Group:shape) %>%
  rename(id = id_tax) %>%
  dplyr::arrange(id)
  

network <- create.network(n = nodes, e = edges, tmain = "Guaymas CRISPR Spacer Matches",
               tsub = "CRISPRDetect - VIBRANT - Genbank - SpacePHARER - BLASTN",
               select_by = "Taxonomy")

network_file <- paste("/Users/ian/Documents/phd_research/CRISPRCas_Sediment/GuaymasC/NetAn/GuaymasC_AllBins_Spacer_Network_", Sys.Date(), ".html", sep = '')
visSave(network, file = network_file)

#Write edges and nodes for later use
edge_file <- paste("/Users/ian/Documents/phd_research/CRISPRCas_Sediment/GuaymasC/NetAn/GuaymasC_AllBin_Undirected_edges_", Sys.Date(), ".txt", sep = '')

node_file <- paste("/Users/ian/Documents/phd_research/CRISPRCas_Sediment/GuaymasC/NetAn/GuaymasC_AllBin_Undirected_nodes_", Sys.Date(), ".txt", sep = '')

write.table(edges,
            edge_file,
            quote = FALSE, col.names = FALSE,
            row.names = FALSE, fileEncoding = "UTF-8",
            sep = "\t")

write.table(nodes,
            node_file,
            quote = FALSE, col.names = FALSE,
            row.names = FALSE, fileEncoding = "UTF-8",
            sep = "\t")
#=============================================================================
setwd("/Users/ian/Documents/phd_research/CRISPRCas_Sediment/GuaymasC/NetAn")

#Network setup and visualization for BLASTN hits only
# 07/02/2020

# generate.nodes.edges <- function(hit_df, prefix, from, to, ccg_file, colshape_df){
#     ### Create nodes and edges dataframes for input into visNetwork
#     ### that contain taxonomic information.
#     ### Requires a dataframe of hits you want to generate the network from,
#     ### the path to a file
#   
#     netan_iu_df = hit_df %>%
#       dplyr::select(from, to) %>%
#       dplyr::distinct()
#     
#     netan_iu_infile = paste(prefix, "Undirected_NetAn_Input", Sys.Date(), ".txt", sep = '')
#     write.table(netan_iu_df,
#                 netan_iu_infile,
#                 quote = FALSE, col.names = FALSE, row.names = FALSE, fileEncoding = "UTF-8",
#                 sep = "\t")
#     
#     print(paste("Wrote NetAn input file:", netan_iu_infile, sep = " "))
#     print("Run NetAn, then run extract_netan_groups.py to extract connected components")
#     #Output from virus_workflow/extract_netan_groups.py
#     #ccg_file <- readline(prompt = "Please provide the path of the connected components file from extract_netan_groups.py")
#     
#     ccg = read.table(ccg_file, col.names = c("Group", `from`)) %>%
#       dplyr::select(`from`, Group)
#     
#     
#   
# }

#Input file for NetAn undirected network
netan_input_undirected_blastn <- blastn_bh %>%
  dplyr::select(Bin, virus) %>%
  dplyr::distinct()

netan_iu_blastn <- paste("/Users/ian/Documents/phd_research/CRISPRCas_Sediment/GuaymasC/NetAn/GuaymasC_AllBin_Undirected_BLASTN_", Sys.Date(), ".txt", sep = '')
write.table(netan_input_undirected_blastn,
            netan_iu_blastn,
            quote = FALSE, col.names = FALSE, row.names = FALSE, fileEncoding = "UTF-8",
            sep = "\t")

#Output from virus_workflow/extract_netan_groups.py
#cc_group_file <- paste("/Users/ian/Documents/phd_research/CRISPRCas_Sediment/GuaymasC/NetAn/GuaymasC_AllBin_Undirected_BLASTN_", Sys.Date(), "_undirected_connected-components.txt", sep = '')
cc_group_file <- "/Users/ian/Documents/phd_research/CRISPRCas_Sediment/GuaymasC/NetAn/GuaymasC_AllBin_Undirected_BLASTN_2020-07-07_undirected_connected-components.txt"
cc_groups_blastn <- read.table(cc_group_file,
                        col.names = c("Group", "Bin")) %>%
  dplyr::select(Bin, Group)


nodes_blastn <- cc_groups_blastn %>%
  left_join(archaea_bacteria_taxonomy) %>%
  rename(id = Bin)


nodes_blastn$Taxonomy[which(is.na(nodes_blastn$Taxonomy) & grepl(".*?\\.[[:digit:]]", nodes_blastn$id))] <- "Genbank Phage"
nodes_blastn$Taxonomy[which(is.na(nodes_blastn$Taxonomy) & grepl("scaffold_2kb", nodes_blastn$id))] <- "Guaymas Phage"

#Create data frame of node Taxonomic id, and corresponding shape + color
node_blastn_taxonomy_uniq <- unique(nodes_blastn$Taxonomy)

node_blastn_taxonomy_col_shape <- data.frame(node_blastn_taxonomy_uniq,
                                      node_col_shape$color[1:length(node_blastn_taxonomy_uniq)],
                                      node_col_shape$shape[1:length(node_blastn_taxonomy_uniq)])
colnames(node_blastn_taxonomy_col_shape) <- c("Taxonomy", "color", "shape")

nodes_blastn <- nodes_blastn %>% left_join(node_blastn_taxonomy_col_shape) %>%
  dplyr::distinct() %>%
  tidyr::unite("id_tax", c("id", "Taxonomy"),
               sep = "; ", remove = FALSE)


edges_blastn <- blastn_bh %>%
  mutate(spacer = gsub("^.*_____", "", crispr_spacer)) %>%
  tidyr::unite(title, c("spacer", "pident"), sep = "; ", remove = FALSE) %>%
  dplyr::select(Bin, virus, color, title, Taxonomy) %>%
  dplyr::distinct() %>%
  rename(from = Bin, to = virus)

edges_blastn_from <- edges_blastn %>% select(from, color, title) %>%
  rename(id = from) %>%
  left_join(nodes_blastn %>% select(id, id_tax)) %>%
  dplyr::select(id_tax) %>%
  dplyr::rename(from = id_tax)

edges_blastn_to <- edges_blastn %>% select(to, color, title) %>%
  rename(id = to) %>%
  left_join(nodes_blastn %>% select(id, id_tax)) %>%
  dplyr::select(id_tax, color, title) %>%
  dplyr::rename(to = id_tax)

edges_blastn <- dplyr::bind_cols(edges_blastn_from, edges_blastn_to)

#edges_blastn$blastn_hit[which(edges_blastn$blastn_hit == 1)] <- "indianred"
#edges_blastn$blastn_hit[which(edges_blastn$blastn_hit == -1)] <- "skyblue"

# edges_blastn <- edges_blastn %>% dplyr::rename(color = blastn_hit) %>%
#   distinct()

nodes_blastn <- nodes_blastn %>% dplyr::select(id_tax, Group:shape) %>%
  rename(id = id_tax) %>%
  dplyr::arrange(id) %>%
  select(-Domain)


network_blastn <- create.network(n = nodes_blastn, e = edges_blastn, tmain = "Guaymas CRISPR Spacer Matches",
                          tsub = "CRISPRDetect - VIBRANT - BLASTN",
                          select_by = "Taxonomy")

network_file <- paste("/Users/ian/Documents/phd_research/CRISPRCas_Sediment/GuaymasC/NetAn/GuaymasC_AllBins_Spacer_Network_BLASTN_", Sys.Date(), ".html", sep = '')
visSave(network_blastn, file = network_file)

#Write edges and nodes for later use
edges_blastn_file <- paste("/Users/ian/Documents/phd_research/CRISPRCas_Sediment/GuaymasC/NetAn/GuaymasC_AllBin_Undirected_edges_", Sys.Date(), ".txt", sep = '')

nodes_blastn_file <- paste("/Users/ian/Documents/phd_research/CRISPRCas_Sediment/GuaymasC/NetAn/GuaymasC_AllBin_Undirected_nodes_blastn_", Sys.Date(), ".txt", sep = '')

write.table(edges_blastn,
            edge_file,
            quote = FALSE, col.names = FALSE,
            row.names = FALSE, fileEncoding = "UTF-8",
            sep = "\t")

write.table(nodes_blastn,
            node_file,
            quote = FALSE, col.names = FALSE,
            row.names = FALSE, fileEncoding = "UTF-8",
            sep = "\t")
#=============================================================================
#all_genome_spacepharer_vibrant <- read.table("/Users/ian/Documents/phd_research/CRISPRCas_Sediment/GuaymasC/spacepharer_output/GuaymasC_AllGenomes_CRISPRDetect_output_qt3_rl23_queryDB---GuaymasC_VIBRANT_extracted_phages_targetDB_fdr_0.01_2020-06-25_23-53-19_HITS-NO-SPACER-INFO.txt")

#asgard_nodes <- read.table("/Users/ian/Documents/phd_research/CRISPRCas_Sediment/GuaymasC/NetAn/asgard.tsv", sep = "\t", header=TRUE, comment.char = '!') 
 
# nodes_united <- nodes %>%
#   mutate(id = gsub("_G", "_Guay", id)) %>%
#   tidyr::unite(from, c("id", "Taxonomy"), sep = "_", remove = FALSE) %>%
#   select(-id) %>%
#   rename(id = from)
  
#asgard_edges <- read.table("/Users/ian/Documents/phd_research/CRISPRCas_Sediment/GuaymasC/NetAn/spacepharer_cat_guaymas_asgard_nocrispr_fdr_0.01_2020-06-17_23-19-46.edges", header = TRUE,sep="\t") #%>%
  #rename(id = "from") %>%
  #left_join(nodes_united %>% select(from, id)) %>%
  #select(from, to)

# create.network(n = asgard_nodes, e = asgard_edges, tmain = "first asgard",
#                tsub = "vibrant - spacepharer",
#                select_by = "Taxonomy")
