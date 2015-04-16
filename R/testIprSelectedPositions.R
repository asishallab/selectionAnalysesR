project.file.path <- function( ..., dir.sep="/" ) {
    # Returns the full file path to a file in subdirectories given in arguments
    # '…'
    #
    # Args:
    #  ...     : The path of subdirectories and finally the file to generate the
    #            complete file path for.
    #  dir.sep : The character to divide dirs with, '/' by default.
    #
    # Returns: The full path
    #   
    paste( path.package( "pamlR" ), ..., sep=dir.sep )
}


#' Integration test to cut things short.
test.iprWithSelectedAA <- function() {
  iprscan.tbl <- read.table(project.file.path("group5792", "group5792_iprscan_result_tbl.tsv"), 
    stringsAsFactors = FALSE)
  aa.fasta <- readAAStringSet(project.file.path("group5792", "group5792_sanitized_macse_AA.fasta"))
  msa.fasta <- readAAMultipleAlignment(project.file.path("group5792", "group5792_AA_msa.fasta"))
  gene <- "PROT1"
  sel.aa <- 68
  iprs <- iprWithSelectedAA(sel.aa, gene, aa.fasta, msa.fasta, iprscan.tbl)
  checkEquals(iprs, c("IPR001223", "IPR013781", "IPR017853"))
} 
