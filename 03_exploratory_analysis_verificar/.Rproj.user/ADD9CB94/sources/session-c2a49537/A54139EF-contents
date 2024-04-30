################################################################################
#
#
#
#
#
#
################################################################################

setwd("~/GitHub/repertorio-TCR-BCR-tcga-ACC/03_exploratory_analysis")

load("tableGeneral.RData")
load("clinical.RData")
load("tcgaACC_pre_processed.RData")

sample_id <- read.csv("sampleidReads.csv")

barcodes <- sample_id$barcode[match(gsub("_report","",unique(tableGeneral$sample_id)),
                                    sample_id$sample_id)]

# -- gender
table(clinical$gender[clinical$submitter_id %in% substr(barcodes,1,12)])

# -- age
mean(clinical$age_at_index[clinical$submitter_id %in% substr(barcodes,1,12)])
range(clinical$age_at_index[clinical$submitter_id %in% substr(barcodes,1,12)])

# -- HSP and LSP
table(tcgaACC$Steroid[tcgaACC$patient %in% substr(barcodes,1,12)])

sum(tableGeneral$Abundance[tableGeneral$chain == "IGH" |
      tableGeneral$chain == "IGHA1" | tableGeneral$chain == "IGHA2"
    | tableGeneral$chain == "IGHD" | tableGeneral$chain == "IGHE"
    | tableGeneral$chain == "IGHG1"| tableGeneral$chain == "IGHG2"
    | tableGeneral$chain == "IGHG3" | tableGeneral$chain == "IGHG4"
    | tableGeneral$chain == "IGHM" | tableGeneral$chain == "IGK"
    | tableGeneral$chain == "IGL"])

sum(tableGeneral$Abundance[tableGeneral$chain == "TRA" |
                             tableGeneral$chain == "TRB" | tableGeneral$chain == "TRD"
                           | tableGeneral$chain == "TRG"])

sum(tableGeneral$Abundance)
200596 + 1399


table(tableGeneral$chain)

200596 / 76





barcode_tcga_rnaseq <- tcgaACC$barcode
idx <- match(gsub("_report","",unique(tableGeneral$sample_id)),
             coldataACC$sample_id)
barcode_tcga_tcrbcr <- coldataACC$barcode[idx]


table(clinical$gender[clinical$submitter_id %in% substr(barcode_tcga_tcrbcr,1,12)])

clinicalData <- data.frame(
  barcode = barcode_tcga_tcrbcr,
  sampleId = gsub("_report","",unique(tableGeneral$sample_id)),
  gender = clinical$gender[
    clinical$submitter_id %in% substr(barcode_tcga_tcrbcr,1,12)],
  age = clinical$age_at_index[
    clinical$submitter_id %in% substr(barcode_tcga_tcrbcr,1,12)])

#idx.barcode <- match(clinicalData$barcode,tcgaACC$barcode)
#idx.barcode <- idx.barcode[is.na(idx.barcode)==FALSE]

clinicalData$steroid <- NA


head(barcode_tcga_tcrbcr)
tcgaACC$barcode[idx.barcode[1:6]]

clinical$gender
