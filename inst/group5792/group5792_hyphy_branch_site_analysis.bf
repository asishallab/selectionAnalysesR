inputRedirect = {};
inputRedirect["01"]="Universal";
inputRedirect["02"]="Yes";
inputRedirect["03"]="Yes";
inputRedirect["04"]="/home/hallab/projects/pamlR/inst/group5792/group5792_CDS_msa.fasta";
inputRedirect["05"]="/home/hallab/projects/pamlR/inst/group5792/group5792_ml_tree_pure_topology.newick";
inputRedirect["06"]="/home/hallab/projects/pamlR/inst/group5792/group5792_hyphy_branch_site_output.txt";

ExecuteAFile ("/biodata/dep_tsiantis/common/software/lib/hyphy/TemplateBatchFiles/BranchSiteREL.bf", inputRedirect);
