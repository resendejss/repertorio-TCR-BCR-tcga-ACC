################################################################################
#
#
#
#
#
#
################################################################################

library(immunarch)

setwd("Documents/onedrive/Documentos/GitHub/repertorio-TCR-BCR-tcga-ACC/04_anlytics/00_clones")
dir_data <- "../../../../Projetos/Bigdata/BigData/BigData/repertorio_tcrbcr_acc/data/"

immdata <- repLoad(paste(dir_data, "outputTrust4_report", sep = ""))
load(paste(dir_data, "metadata.RData", sep = ""))

immdata$meta <- metadata[metadata$sample_id %in% gsub("_report", "", immdata$meta$Sample),]

# -- separar os grupos bcr e tcr

immdata$meta$Sample <- paste(immdata$meta$sample_id,"_report",sep = "")

immdata_bcr <- repFilter(immdata,
                         .method = "by.clonotype",
                         .query = list(V.name = include("IG")),
                         .match = "substring")

immdata_tcr <- repFilter(immdata,
                         .method = "by.clonotype",
                         .query = list(V.name = include("TR")),
                         .match = "substring")

# -- dist: group sequence by V.name, J.name and length CDR3
dist_BCR <- seqDist(immdata_bcr$data, 
                   .group_by = c("V.name", "J.name"),
                   .group_by_seqLength = TRUE,
                   .col = "CDR3.nt")
clust_BCR <- seqCluster(immdata_bcr$data, dist_BCR, .perc_similarity = 0.9)

dist_TCR <- seqDist(immdata_tcr$data[-12], 
                    .group_by = c("V.name", "J.name"),
                    .group_by_seqLength = TRUE,
                    .col = "CDR3.nt")
clust_TCR <- seqCluster(immdata_tcr$data[-12], dist_TCR, .perc_similarity = 0.95)

temp <- immdata_tcr$data[12]

temp$`130723_UNC9-SN296_0386_BC2E4WACXX_TAGCTT_L006_report`$J.name
nchar(temp$`130723_UNC9-SN296_0386_BC2E4WACXX_TAGCTT_L006_report`$CDR3.nt)

temp$`130723_UNC9-SN296_0386_BC2E4WACXX_TAGCTT_L006_report`$Cluster <- 
  c("TRAV9-2/TRAJ17_length_39","TRAV9-2/TRAJ17_length_39")

names(temp)

clust_TCR$`130723_UNC9-SN296_0386_BC2E4WACXX_TAGCTT_L006_report` <- 
  as.data.frame(temp$`130723_UNC9-SN296_0386_BC2E4WACXX_TAGCTT_L006_report`) 


clust_TCR$`130723_UNC9-SN296_0386_BC2E4WACXX_TAGCTT_L006_report`

clust_TCR$`130723_UNC9-SN296_0386_BC2E4WACXX_ACTTGA_L003_report`

library(dplyr)

# funcao para agrupar os dataframes
removendoi_NA <- function(df){
  df %>%
    filter(!rowSums(is.na(.)) == ncol(.))
}

agrupar_por_cluster <- function(df){
  df = df[!is.na(df$Clone),]
  df %>%
    group_by(Cluster) %>%
    mutate(Clones = sum(Clones),
           Proportion = sum(Proportion)) %>%
    ungroup() %>%
    distinct(Cluster, .keep_all = TRUE) %>%
    as.data.frame()
}


teste_TCR <- lapply(clust_TCR, removendoi_NA)
teste_BCR <- lapply(clust_BCR, removendoi_NA)

clust_TCR_agrupados <- lapply(teste_TCR, agrupar_por_cluster)
clust_BCR_agrupados <- lapply(teste_BCR, agrupar_por_cluster)


rownames(immdata$data$`130723_UNC9-SN296_0386_BC2E4WACXX_GATCAG_L005_report`)
table(is.na(immdata$data$`130723_UNC9-SN296_0386_BC2E4WACXX_GATCAG_L005_report`$Clones))

# aplicando a funcao na lista de dataframes
clust_TCR_agrupados <- lapply(clust_TCR, agrupar_por_cluster)
clust_BCR_agrupados <- lapply(clust_BCR, agrupar_por_cluster)

length(clust_BCR$`130723_UNC9-SN296_0386_BC2E4WACXX_GATCAG_L005_report`$Cluster)
length(unique(clust_BCR$`130723_UNC9-SN296_0386_BC2E4WACXX_GATCAG_L005_report`$Cluster))

length(clust_BCR_agrupados$`130723_UNC9-SN296_0386_BC2E4WACXX_GATCAG_L005_report`$Cluster)
length(unique(clust_BCR_agrupados$`130723_UNC9-SN296_0386_BC2E4WACXX_GATCAG_L005_report`$Cluster))

clust_BCR$`130723_UNC9-SN296_0386_BC2E4WACXX_GATCAG_L005_report`
clust_BCR_agrupados$`130723_UNC9-SN296_0386_BC2E4WACXX_GATCAG_L005_report`


table(clust_BCR$`130723_UNC9-SN296_0386_BC2E4WACXX_GATCAG_L005_report`$Cluster)

clust_BCR$`130723_UNC9-SN296_0386_BC2E4WACXX_GATCAG_L005_report`[
  clust_BCR$`130723_UNC9-SN296_0386_BC2E4WACXX_GATCAG_L005_report`$Cluster == "IGHV1-3/IGHJ4_length_24",]

clust_BCR_agrupados$`130723_UNC9-SN296_0386_BC2E4WACXX_GATCAG_L005_report`[
  clust_BCR_agrupados$`130723_UNC9-SN296_0386_BC2E4WACXX_GATCAG_L005_report`$Cluster == "IGHV1-3/IGHJ4_length_24",]


clust_TCR$`130723_UNC9-SN296_0386_BC2E4WACXX_ATCACG_L003_report`$CDR3.nt[1]
clust_TCR$`130723_UNC9-SN296_0386_BC2E4WACXX_ATCACG_L003_report`$CDR3.nt[1]

clust_BCR_agrupados$`130723_UNC9-SN296_0386_BC2E4WACXX_ACTTGA_L003_report`
clust_BCR$`130723_UNC9-SN296_0386_BC2E4WACXX_ATCACG_L003_report`





