#'################################################################################
#' Functions to identify which conserved protein domains (InterPro) overlap with #
#' selected homologous amino acids (codons).                                     #
#'################################################################################

readInterProScanResultTable <- function(path.2.iprscan.res) {
  #' Parses the result file of InterProScan into a data.frame.
  #'
  #' Args:
  #'  path.2.iprscan.res : The file path to the argument InterProScan result
  #'                       table.
  #'
  #' Returns: A data.frame with 14 columns.
  #'   
  read.table(path.2.iprscan.res, sep = "\t", stringsAsFactors = FALSE, fill = TRUE, 
    comment.char = "", quote = "", na.string = "")
} 

iprWithSelectedAA <- function(sel.aa, gene, aa.fasta, msa.fasta, iprscan.tbl) {
  # Identifies conserved protein domains (InterPro) overapping with the
  # conserved homologous amino acid (codon) 'sel.aa'.
  #
  # Args:
  #  sel.aa      : The index of the conserved amino acid under selection in the
  #                corresponding multiple sequence alignment (MSA)
  #  gene        : The gene accession / ID as used in both 'aa.fasta' and
  #                'msa.fasta'
  #  aa.fasta    : The result of readAAStringSet(path_2_AAs.fasta) [package
  #                Biostrings]
  #  msa.fasta   : The result of readAAMultipleAlignment(path_2_AAs_MSA.fasta)
  #                [package Biostrings]
  #  iprscan.tbl : The result of calling
  #                readInterProScanResultTable(path_2_interproscan_result_table.tsv)
  #
  # Returns: A character vector of matching InterPro entries.
  #   
  sel.aa.unaligned.pos <- unalignedAAforAlignedAA(gene, sel.aa, aa.fasta, msa.fasta)
  if (!is.na(sel.aa.unaligned.pos)) 
    domainsForPos(gene, sel.aa.unaligned.pos, iprscan.tbl) else NA
} 

unalignedAAforAlignedAA <- function(gene, sel.aa, aa.fasta, msa.fasta, min.radius = 3) {
  # Infers the corresponding un-aligned position of an selected amino acid
  # obtained from an multiple sequence alignment. This is done by iterativly
  # increasing a amino acid substring, without gap characters, and finding its
  # position on the original unaligned amino acid sequence.
  #
  # Args:
  #  gene       : The gene accession / ID as used in both 'aa.fasta' and
  #               'msa.fasta'
  #  sel.aa     : The index of the homologous amino acid subject to selection
  #               (integer coordinate)
  #  aa.fasta   : The result of readAAStringSet(path_2_AAs.fasta) [package
  #               Biostrings]
  #  msa.fasta  : The result of readAAMultipleAlignment(path_2_AAs_MSA.fasta)
  #               [package Biostrings]
  #  min.radius : The minimum radius with which to start selecting up- and
  #               downstream subsequences
  #
  # Returns: An integer; either NA, if no position could be inferred, or the
  # corresponding un-aligned position.
  #   
  aa.full <- aa.fasta[[gene]]
  radii   <- min.radius:(length(aa.full))
  aa.pos  <- as.integer(NA)
  i       <- 1
  while (i <= length(radii) && is.na(aa.pos)) {
    radius <- radii[[i]]
    x <- getAASubstring(msa.fasta, gene, sel.aa, radius = radius)
    aa.subseq <- removeNonAAChars(x$aa.subseq)
    up.strm.offset <- x$up.strm.offset - countNonAminoAcidCharsFromStart(x$aa.subseq)
    y <- subseqPos(aa.subseq, aa.full)
    if (length(y) == 1 && y > -1) 
      aa.pos <- y - up.strm.offset
    i <- i + 1
  }
  aa.pos
} 
 
getAASubstring <- function(msa.fasta, gene, sel.aa, radius = 3) {
  #' Selects the aligned amino acid subsequence around the given coordinate
  #' including radius up- and downstream AAs.
  #'
  #' Args:
  #'  msa.fasta : The amino acid multiple sequence alignment (MSA) as returned
  #'              by readAAMultipleAlignment(…) [package Biostrings]
  #'  gene      : The gene identifier / accession as found in the MSA
  #'  sel.aa    : The index of the amino acid coordinate, starting at 1
  #'  radius    : The number of up- and downstream alignment characters to
  #'              return
  #'
  #' Returns: A list with two entries. Entry 'aa.subseq' is an instance of
  #' AAString of max length 1+2*radius. If radius exceeds start or stop of the
  #' AA sequence the result is shorter accordingly. Entry 'up.strm.offset'
  #' informs about how many characters prepend the original selected amino acid
  #' according to 'sel.aa'.
  #'   
  algn.aa.seq <- attr(msa.fasta, "unmasked")[[gene]]
  up.strm <- if (sel.aa > radius) 
    sel.aa - radius else 1
  down.strm <- if (sel.aa + radius <= length(algn.aa.seq)) 
    sel.aa + radius else length(algn.aa.seq)
  list(aa.subseq = algn.aa.seq[up.strm:down.strm], up.strm.offset = (sel.aa - 
    up.strm))
} 

countNonAminoAcidCharsFromStart <- function(aa.seq, non.aa.regex = "^[-]+") {
  #' Counts the number of non amino acid characters at the start of the
  #' argument 'aa.seq'.
  #'
  #' Args:
  #'  aa.seq       : The aligned amino acid sequence possibly including gap
  #'                 characters
  #'  non.aa.regex : The PERL regular expression to identify non amino acid
  #'                 characters at the start of the argument 'aa.seq'
  #'
  #' Returns: An integer - the count of non amino acid chars at the beginning
  #' of the argument 'aa.seq'.
  #'   
  x <- attr(regexpr(non.aa.regex, toString(aa.seq), perl = TRUE), "match.length")
  if (x > 0) 
    x else 0
} 

removeNonAAChars <- function(aa.seq, non.aa.regex = "[-]") {
  #' Removes all non Amino Acid characters from an aligned amino acid sequence.
  #'
  #' Args:
  #'  aa.seq       : The aligned amino acid sequence, possibly including to be
  #'                 removed characters.
  #'  non.aa.regex : A PERL regular expression used to identify to be removed
  #'                 matches from 'aa.seq'
  #'
  #' Returns: A string in which all matches have been deleted.
  #'   
  gsub(non.aa.regex, "", toString(aa.seq), perl = TRUE)
} 

subseqPos <- function(aa.subseq, aa.full) {
  #' Infers the amino acid position(s) at which the argument 'aa.subseq' starts
  #' in the full 'aa.full' amino acid sequence.
  #'
  #' Args:
  #'  aa.subseq : The amino acid subsequence to be matched
  #'  aa.full   : The full amino acid sequence
  #'
  #' Returns: An integer vector with matching start positions.
  #'   
  as.integer(gregexpr(toString(aa.subseq), toString(aa.full), fixed = TRUE)[[1]])
} 

domainsForPos <- function(gene, aa.pos, iprscan.tbl, gene.col = "V1", start.col = "V7", 
  end.col = "V8", ipr.col = "V12") {
  #' Looks up the conserved protein domains (InterPro) that have been annotated
  #' to gene and overlap with amino acid position 'aa.pos'.
  #'
  #' Args:
  #'  gene      : The gene accession or ID as used in iprscan.tbl
  #'  aa.pos    : The amino acid position (non aligned) to be overlapped by any
  #'              InterPros
  #'  gene.col  : The column index or name of iprscan.tbl in which to look up
  #'              the gene accessions
  #'  start.col : The column index or name of iprscan.tbl in which to look up
  #'              the domain start positions
  #'  end.col   : The column index or name of iprscan.tbl in which to look up
  #'              the domain end positions
  #'  ipr.col   : The column index or name of iprscan.tbl in which to look up
  #'              the InterPro identifier
  #'
  #' Returns: A character vector of matching InterPro domains.
  #'   
  x <- iprscan.tbl[which(iprscan.tbl[, gene.col] == gene), ]
  sort(unique(x[which(x[, start.col] <= aa.pos & x[, end.col] >= aa.pos), ipr.col]), 
    na.last = NA)
} 

replaceSanitizedWithOriginalIDs <- function(xstring.set, name.maps, san.col = "sanitized", 
  orig.col = "original") {
  # Replaces sanitized protein identifiers in XStringSet 'xstring.set' with the
  # original ones held in 'name.maps' table.
  #
  # Args:
  #  xstring.set : An instance of class XStringSet [package: Biostrings]
  #  name.maps   : A data.frame with at least two columns one holding the
  #                sanitized and the other the original protein IDs
  #  san.col     : The column index or name of 'name.maps' holding the
  #                sanitized protein IDs
  #  orig.col    : The column index or name of 'name.maps' holding the original
  #                protein IDs
  #
  # Returns: A copy of argument 'xstring.set' with the sanitized protein IDs
  # replaced with their originals.
  #   
  names(xstring.set) <- as.character(lapply(names(xstring.set), function(x) {
    name.maps[which(name.maps[, san.col] == x), orig.col]
  }))
  xstring.set
} 
