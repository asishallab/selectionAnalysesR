inputRedirect = {};
inputRedirect["01"]="Universal";
inputRedirect["02"]="1";
inputRedirect["03"]="/home/hallab/projects/pamlR/inst/group5792/group5792_CDS_msa.fasta";
inputRedirect["04"]="/home/hallab/projects/pamlR/inst/group5792/group5792_ml_tree_pure_topology.newick";
inputRedirect["05"]="20";
inputRedirect["06"]="5";
inputRedirect["07"]="2000000";
inputRedirect["08"]="1000000";
inputRedirect["09"]="100";
inputRedirect["10"]="0.5";
inputRedirect["11"]="/home/hallab/projects/pamlR/inst/group5792/group5792_hyphy_fubar_output.txt";

ExecuteAFile ("/biodata/dep_tsiantis/common/software/lib/hyphy/TemplateBatchFiles/FUBAR.bf", inputRedirect);
