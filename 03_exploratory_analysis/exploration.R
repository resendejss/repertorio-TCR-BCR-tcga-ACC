################################################################################
#
#
#
#
#
#
################################################################################

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
