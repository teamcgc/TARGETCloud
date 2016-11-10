#' Get the AML Exome-based Somatic MAF file
#'
#' @importFrom readr read_tsv
#'
#' @export
getAMLExomeMutationMAF = function() {
  maf = read_tsv('ftp://caftpd.nci.nih.gov/pub/dcc_target/AML/WXS/L3/mutation/BCM/VerifiedSomatic/TARGET_AML_WXS_somatic_filtered_verified.mafplus.txt')
  return(maf)
}

#' Get the Neuroblastoms Exome-based Somatic MAF file
#'
#' @importFrom readr read_tsv
#'
#' @export
getNBLExomeMutationMAF = function() {
  maf = read_tsv(url('ftp://caftpd.nci.nih.gov/pub/dcc_target/NBL/WXS/L3/mutation/Broad/VerifiedSomatic/TARGET_NBL_WXS_somatic_verified.maf.txt'))
  return(maf)
}

#' Get the AML Exome-based Somatic MAF file
#'
#' @importFrom readxl read_excel
#'
#' @export
getALLExomeMutationMAF = function() {
  tmpfile = file.path(tempdir(),'abc.xlsx')
  download.file('ftp://caftpd.nci.nih.gov/pub/dcc_target/ALL/WXS/Phase2/L3/mutation/BCM/VerifiedSomatic/target-all-recurrent-somatic-verified-v2.4-mafplus.xlsx',
                            destfile=tmpfile)
  maf = read_excel(tmpfile)
}

#' Get the Wilms Tumor Exome-based Somatic MAF file
#'
#' @importFrom readr read_tsv
#'
#' @export
getWTExomeMutationMAF = function() {
  return(read_tsv(url('ftp://caftpd.nci.nih.gov/pub/dcc_target/WT/WXS/L3/mutation/NCI-Meerzaman/VerifiedSomatic/target-wt-primary-recurrent-NCI-somatic-exonic-verified.bcmmaf.txt')))
}


#' Fix the column names from a data.frame for insertion to BigQuery
#'
#' BigQuery requires that names be letters, underscores, or numbers.
#' This function returns a data.frame with column names that conform
#' to that standard.
#'
#' @param df a data.frame or data.frame-like object
#'
#' @export
fixColnamesForBigQuery = function(df) {
  colnames(df) = gsub('\\.','_',make.names(colnames(df)))
  return(df)
}

