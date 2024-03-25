################################################################################
#
#
#
#
#
#
################################################################################

cluster_trust4 <- function(name_report_trust4){
  
  require(dplyr)
  
  
  data <- read.delim(name_report_trust4)
  
  df_agrupado <- data %>%
    group_by(cid) %>%
    mutate(soma_count = sum(X.count)) %>%
    distinct(cid, .keep_all = TRUE)
  
  
  df_agrupado <- data %>%
    group_by(cid) %>%
    mutate(soma_count = sum(X.count)) %>%
    distinct(cid, .keep_all = TRUE)
  
  df_agrupado$X.count <- df_agrupado$soma_count
  df_agrupado <- df_agrupado[,-ncol(df_agrupado)]
  
  write.table(df_agrupado, file = name_report_trust4, sep = "\t", row.names = F)
}


files <- list.files("../data/outputTrust4_report_cluster/")

for (i in files) {
  cluster_trust4(paste("../data/outputTrust4_report_cluster/",i, sep = ""))
}


