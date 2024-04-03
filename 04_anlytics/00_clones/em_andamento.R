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

# separar os grupos high e low
immdata$meta$steroid

immdata_high <- immdata
immdata_high$meta <- immdata_high$meta[immdata_high$meta$steroid == "HSP",]
immdata_high$meta <- immdata_high$meta[!is.na(immdata_high$meta$steroid),]
immdata_high$data <- immdata_high$data[gsub("_report","",names(immdata_high$data)) %in%
                                         immdata_high$meta$sample_id]

immdata_low <- immdata
immdata_low$meta <- immdata_low$meta[immdata_low$meta$steroid == "LSP",]
immdata_low$meta <- immdata_low$meta[!is.na(immdata_low$meta$steroid),]
immdata_low$data <- immdata_low$data[gsub("_report","",names(immdata_low$data)) %in%
                                         immdata_low$meta$sample_id]

immdata$meta$Sample <- paste(immdata$meta$sample_id,"_report",sep = "")
immdata_low$meta$Sample <- paste(immdata_low$meta$sample_id,"_report",sep = "")
immdata_high$meta$Sample <- paste(immdata_high$meta$sample_id,"_report",sep = "")

immdata_bcr <- repFilter(immdata,
                         .method = "by.clonotype",
                         .query = list(V.name = include("IG")),
                         .match = "substring")

immdata_tcr <- repFilter(immdata,
                         .method = "by.clonotype",
                         .query = list(V.name = include("TR")),
                         .match = "substring")

immdata_low_tcr <- repFilter(immdata_low,
                             .method = "by.clonotype",
                             .query = list(V.name = include("TR")),
                             .match = "substring")

immdata_high_tcr <- repFilter(immdata_high,
                             .method = "by.clonotype",
                             .query = list(V.name = include("TR")),
                             .match = "substring")

immdata_low_bcr <- repFilter(immdata_low,
                             .method = "by.clonotype",
                             .query = list(V.name = include("IG")),
                             .match = "substring")

immdata_high_bcr <- repFilter(immdata_high,
                             .method = "by.clonotype",
                             .query = list(V.name = include("IG")),
                             .match = "substring")

immdata_low_tcr$meta <- immdata_low$meta[immdata_low$meta$Sample %in% 
                                             names(immdata_low_tcr$data),]

immdata_high_tcr$meta <- immdata_high$meta[immdata_high$meta$Sample %in% 
                                             names(immdata_high_tcr$data),]

immdata_low_bcr$meta <- immdata_low$meta[immdata_low$meta$Sample %in% 
                                           names(immdata_low_bcr$data),]

immdata_high_bcr$meta <- immdata_high$meta[immdata_high$meta$Sample %in% 
                                             names(immdata_high_bcr$data),]

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


dist_high_BCR <- seqDist(immdata_high_bcr$data, 
                   .group_by = c("V.name", "J.name"),
                   .group_by_seqLength = TRUE,
                   .col = "CDR3.nt")
clust_high_BCR <- seqCluster(immdata_high_bcr$data, dist_high_BCR, .perc_similarity = 0.9)

dist_low_BCR <- seqDist(immdata_low_bcr$data, 
                         .group_by = c("V.name", "J.name"),
                         .group_by_seqLength = TRUE,
                         .col = "CDR3.nt")
clust_low_BCR <- seqCluster(immdata_low_bcr$data, dist_low_BCR, .perc_similarity = 0.9)

immdata_high_tcr$data$`130723_UNC9-SN296_0386_BC2E4WACXX_TAGCTT_L006_report` <- 
  as.data.frame(immdata_high_tcr$data$`130723_UNC9-SN296_0386_BC2E4WACXX_TAGCTT_L006_report`)


dist_high_TCR <- seqDist(immdata_high_tcr$data, 
                        .group_by = c("V.name", "J.name"),
                        .group_by_seqLength = TRUE,
                        .col = "CDR3.nt")
clust_high_TCR <- seqCluster(immdata_high_tcr$data, dist_high_TCR, .perc_similarity = 0.95)


dist_high_TCR$`130723_UNC9-SN296_0386_BC2E4WACXX_TAGCTT_L006_report`
immdata_high_tcr$data$`130723_UNC9-SN296_0386_BC2E4WACXX_TAGCTT_L006_report`

dist_low_TCR <- seqDist(immdata_low_tcr$data, 
                        .group_by = c("V.name", "J.name"),
                        .group_by_seqLength = TRUE,
                        .col = "CDR3.aa")
clust_low_TCR <- seqCluster(immdata_low_tcr$data, dist_low_TCR, .perc_similarity = 0.95)



# distancia
distTCR_high_aa <- seqDist(immdata_high$data, .col = "CDR3.aa")
distTCR_high_nt <- seqDist(immdata_high$data, .col = "CDR3.nt")

distTCR_low_aa <- seqDist(immdata_low$data, .col = "CDR3.aa")
distTCR_low_nt <- seqDist(immdata_low$data, .col = "CDR3.nt")



# agrupamento de clonotipos
?pubRep
