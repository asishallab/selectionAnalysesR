require(pamlR)

# Hail User:
message("Usage: Rscript runPamlR.R codingSequences.fasta [number_of_threads]")

# Read command line arguments:
c.args <- commandArgs(trailingOnly = TRUE)
cds.fst <- c.args[[1]]
no.threads <- if (length(c.args) == 2) {
  as.integer(c.args[[2]])
} else {
  1
} 

# Prepare output folder:
file.name <- sub("\\.\\S+$", "", str_match(cds.fst, "[^/]+$")[[1, 1]], perl = TRUE)
output.dir <- paste(dirname(normalizePath(cds.fst)), "/", file.name, "/", sep = "") 
dir.create( output.dir, showWarnings=FALSE )
fam.san.fasta.path <- paste(output.dir, file.name, "_sanitized.fasta", sep = "")
fam.aa.fasta.path <- sub("\\.fasta$", "_macse_AA.fasta", fam.san.fasta.path, perl = TRUE) 
fam.aa.san.fasta.path <- sub("\\.fasta$", "_macse_AA_sanitized.fasta", fam.san.fasta.path, perl = TRUE) 
fam.aa.msa.path <- paste( output.dir, file.name, "_AA_msa.fasta", sep="" )
fam.cds.msa.path <- paste( output.dir, file.name, "_CDS_msa.fasta", sep="" )
fam.name.maps.tbl.path <- paste(output.dir, file.name, "_name_mappings.txt", sep = "") 
fam.tree.path <- paste( output.dir, file.name, "_ml_tree.newick", sep="" )
fam.hyphy.branch.site.batch.file.path <- paste( output.dir, file.name, "_hyphy_branch_site_analysis.bf", sep="" )
fam.hyphy.branch.site.output.path <- paste( output.dir, file.name, "_hyphy_branch_site_output.txt", sep="" )
fam.hyphy.meme.batch.file.path <- paste( output.dir, file.name, "_hyphy_meme_analysis.bf", sep="" )
fam.hyphy.meme.output.path <- paste( output.dir, file.name, "_hyphy_meme_output.txt", sep="" )
fam.hyphy.fubar.batch.file.path <- paste( output.dir, file.name, "_hyphy_fubar_analysis.bf", sep="" )
fam.hyphy.fubar.output.path <- paste( output.dir, file.name, "_hyphy_fubar_output.txt", sep="" )
fam.tree.4.paml.path <- paste( output.dir, file.name, "_ml_tree_pure_topology.newick", sep="" )
system(paste("mkdir -p", file.name))
# Read Input and remove stop codons:
fam <- DNAStringSet(lapply(readDNAStringSet(cds.fst), removeStopCodon)) 
# Sanitize sequence names:
fam.name.maps <- data.frame(original = names(fam), sanitized = paste("PROT", 1:length(fam), 
  sep = ""), stringsAsFactors = FALSE)
write.table(fam.name.maps, fam.name.maps.tbl.path, row.names = FALSE)
names(fam) <- as.character(lapply(names(fam), function(x) fam.name.maps[which(fam.name.maps$original == 
  x), "sanitized"]))
writeXStringSet(fam, fam.san.fasta.path)
# Translate to AAs:
translate2AASeqs(fam.san.fasta.path)
# Remove invalid AA-Sequences, i.e. AA-Seqs with premature stop-codons:
fam.aas <- readAAStringSet(fam.aa.fasta.path)
fam.aas.san <- fam.aas[validateAAStringSet(fam.aas)]
# Warn about removed AA-Seqs:
if (length(fam.aas.san) < length(fam.aas)) {
  len.diff <- length(fam.aas) - length(fam.aas.san)
  warning(len.diff, " amino-acid-sequences had a premature stop codon and were removed from further analysis.")
}
# Write out the sanitized amino acid seqs:
writeXStringSet(fam.aas.san, fam.aa.san.fasta.path)
# Generate a multiple sequence alignment:
system(paste("mafft --thread", no.threads, "--auto", fam.aa.san.fasta.path, ">", 
  fam.aa.msa.path)) 
fam.aas.san.msa <- readAAMultipleAlignment(fam.aa.msa.path)
# Use the aligned AA-Seqs as quide to align the CDS Sequences:
fam.cds.msa <- alignCDSSetWithAlignedAAsAsGuide(fam, attr(fam.aas.san.msa, "unmasked"))
writeXStringSet(attr(fam.cds.msa, "unmasked"), fam.cds.msa.path)
# Generate Phylogenetic maximum likelihood Tree:
system(paste("OMP_NUM_THREADS=", no.threads, " FastTreeMP -nt -gtr -gamma < ", 
  fam.cds.msa.path, " > ", fam.tree.path, sep = "")) 
fam.tree <- read.tree(fam.tree.path)
fam.tree.4.paml <- removeNodeLabelsAndBranchLengths(fam.tree)
write.tree(fam.tree.4.paml, fam.tree.4.paml.path)
# Generate HYPHY batch files:
brew( text=hyphy.branch.site.bf, output=fam.hyphy.branch.site.batch.file.path )
brew( text=hyphy.meme.bf, output=fam.hyphy.meme.batch.file.path )
brew( text=hyphy.fubar.bf, output=fam.hyphy.fubar.batch.file.path )
