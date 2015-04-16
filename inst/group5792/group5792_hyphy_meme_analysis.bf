inputRedirect = {};
inputRedirect["01"]="Universal";
inputRedirect["02"]="New Analysis";
inputRedirect["03"]="/home/hallab/projects/pamlR/inst/group5792/group5792_CDS_msa.fasta";
inputRedirect["04"]="Custom";
inputRedirect["05"]="110240";
inputRedirect["06"]="/home/hallab/projects/pamlR/inst/group5792/group5792_ml_tree_pure_topology.newick";
inputRedirect["07"]="/biodata/dep_tsiantis/grp_gan/song/asis_test_hyphy/hyphy_log.txt";
inputRedirect["08"]="Estimate dN/dS only";
inputRedirect["09"]="MEME";
inputRedirect["10"]="0.1";
inputRedirect["11"]="N";
inputRedirect["12"]="/home/hallab/projects/pamlR/inst/group5792/group5792_hyphy_meme_output.txt";

ExecuteAFile ("/biodata/dep_tsiantis/common/software/lib/hyphy/TemplateBatchFiles/QuickSelectionDetection.bf", inputRedirect);
