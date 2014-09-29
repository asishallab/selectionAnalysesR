hyphy.branch.site.bf <- 'inputRedirect = {};
inputRedirect["01"]="Universal";
inputRedirect["02"]="Yes";
inputRedirect["03"]="Yes";
inputRedirect["04"]="<%= fam.cds.msa.path %>";
inputRedirect["05"]="<%= fam.tree.4.paml.path %>";
inputRedirect["06"]="<%= fam.hyphy.branch.site.output.path %>";

ExecuteAFile ("/biodata/dep_tsiantis/common/software/lib/hyphy/TemplateBatchFiles/BranchSiteREL.bf", inputRedirect);'

hyphy.meme.bf <- 'inputRedirect = {};
inputRedirect["01"]="Universal";
inputRedirect["02"]="New Analysis";
inputRedirect["03"]="<%= fam.cds.msa.path %>";
inputRedirect["04"]="Custom";
inputRedirect["05"]="110240";
inputRedirect["06"]="<%= fam.tree.4.paml.path %>";
inputRedirect["07"]="/biodata/dep_tsiantis/grp_gan/song/asis_test_hyphy/hyphy_log.txt";
inputRedirect["08"]="Estimate dN/dS only";
inputRedirect["09"]="MEME";
inputRedirect["10"]="0.1";
inputRedirect["11"]="N";
inputRedirect["12"]="<%= fam.hyphy.meme.output.path %>";

ExecuteAFile ("/biodata/dep_tsiantis/common/software/lib/hyphy/TemplateBatchFiles/QuickSelectionDetection.bf", inputRedirect);'

hyphy.fubar.bf <- 'inputRedirect = {};
inputRedirect["01"]="Universal";
inputRedirect["02"]="1";
inputRedirect["03"]="<%= fam.cds.msa.path %>";
inputRedirect["06"]="<%= fam.tree.4.paml.path %>";
inputRedirect["12"]="<%= fam.hyphy.fubar.output.path %>";

ExecuteAFile ("/biodata/dep_tsiantis/common/software/lib/hyphy/TemplateBatchFiles/FUBAR.bf", inputRedirect);'

removeStopCodon <- function( dna.str, st.cdns='(TAG|tag|TAA|taa|TGA|tga)$' ) {
  DNAString( gsub( st.cdns, '', toString( dna.str ), perl=TRUE ) )
}

validateCds <- function(cds) {
    if (class(cds) != "DNAString")
        stop("Argument 'cds' is not of class 'DNAString'.")
    grepl("\\S\\*\\S", translate(cds), perl = TRUE)
} 

validateDNAStringSet <- function(cds.set) {
    all(as.logical(lapply(names(cds.set), function(acc) {
        tryCatch({
            validateCds(cds.set[[acc]])
        }, error = function(e) {
            stop("CDS ", acc, " caused an error: ", e)
        })
    })))
} 

translate2AASeqs <- function(path.2.cds.fasta, macse.call = "java -Xmx600m -jar ~/bin/macse_v1.01b.jar -prog translateNT2AA") {
    # Translates a fasta file of coding sequences to amino acid sequences using
    # the MACSE program. Biostrings throws an error if ambiguous nucleotides,
    # e.g. 'N', are found. The latest version of Biostrings, which can
    # circumnavigate this problem, is not yet installed…
    #
    # Returns: TRUE, if and only if no error has occurred. The output file will
    # be named as the input file with a trailing '_macse_AA.fasta'.
    #   
    cmd <- paste(macse.call, "-seq", path.2.cds.fasta)
    system(cmd)
    TRUE
} 

validateAASeqs <- function(aa.seq) {
    if (class(aa.seq) != "AAString")
        stop("Argument 'aa.seq' is not of class 'AAString'.")
    ! grepl("\\S\\*\\S", aa.seq, perl = TRUE)
} 

validateAAStringSet <- function(aa.set) {
    as.logical(lapply(names(aa.set), function(acc) {
        tryCatch({
            validateAASeqs(aa.set[[acc]])
        }, error = function(e) {
            stop("Amino-Acid-Sequence ", acc, " caused an error: ", e)
        })
    }))
} 

alignCDSWithAlignedAASeq <- function(cds, aligned.aa.seq, aa.gap = "-", codon.gap = "---") {
    aa.chars <- strsplit(toString(aligned.aa.seq), split = NULL)[[1]]
    cds.chars <- strsplit(toString(cds), split = NULL)[[1]]
    i <- 1
    DNAString(paste(lapply(aa.chars, function(aa.char) {
        if (aa.char == aa.gap) {
            codon.gap
        } else {
            start.ind <- (i - 1) * 3 + 1
            stop.ind <- i * 3
            i <<- i + 1
            paste(cds.chars[start.ind:stop.ind], collapse = "")
        }
    }), collapse = ""))
} 

alignCDSSetWithAlignedAAsAsGuide <- function(unaligned.cds.set, aligned.aa.set) {
    cds.set <- DNAStringSet(lapply(names(aligned.aa.set), function(acc) {
        alignCDSWithAlignedAASeq(unaligned.cds.set[[acc]], aligned.aa.set[[acc]])
    }))
    names(cds.set) <- names(aligned.aa.set)
    DNAMultipleAlignment(cds.set)
} 

prepTreeForPAML <- function(phyl.tree, asgn.nd.lbls = length(phyl.tree$node.label) - 
    1) {
    tr.root <- setdiff(phyl.tree$edge[, 1], phyl.tree$edge[, 2])
    if (tr.root != length(phyl.tree$tip.label) + 1) 
        stop("Expecting argument phyl.tree's root to be number of leaves + 1")
    n.tips <- length(phyl.tree$tip.label)
    inner.nd.lab <- length(phyl.tree$tip.label) + 1
    for (i in 1:nrow(phyl.tree$edge)) {
        tr.nd <- phyl.tree$edge[[i, 2]]
        k <- tr.nd - length(phyl.tree$tip.label)
        if (tr.nd <= n.tips) {
            phyl.tree$tip.label[[tr.nd]] <- paste(phyl.tree$tip.label[[tr.nd]], 
                "#", tr.nd, sep = "")
        } else if (phyl.tree$edge[[i, 1]] != tr.root) {
            phyl.tree$node.label[[k]] <- paste("#", inner.nd.lab, sep = "")
            inner.nd.lab <- inner.nd.lab + 1
        } else {
            phyl.tree$node.label[[k]] <- ""
        }
    }
    phyl.tree
} 

removeNodeLabelsAndBranchLengths <- function( phyl.tree ) {
  phyl.tree$node.label <- NULL
  phyl.tree$edge.length <- NULL
  phyl.tree
}
