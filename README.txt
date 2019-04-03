
#replace Pubchemid with Smiles and merge known Smiles list with imputed Smiles
python putSmiles4PubChemID.py ../data/input-data/pubchem2smiles.txt ../data/input-data/Drug_info_release.csv > ../data/input-data/drug_smiles.csv

# create fingerprint features
python createFingerprintFeatures.py ../data/input-data/drug_smiles.csv > ../data/features/Dream10-drugs-fingerprint-matrix.txt

#  create target feature (binary matrix)  If drugs in combination that have a target in Target Set, its corresponding entry is  1, otherwise 0
python createTargetFeatureMatrix.py ../data/input-data/Drug_info_release.csv ../data/input-data/targets-expansion.txt > ../data/features/Dream10-drugs-commontarget-matrix.txt 


# create common pathway feature (binary matrix)  if drugs target share common pathway  1, otherwise 0
python createCommonPathwayFeatureMatrix.py ../data/input-data/Drug_info_release.csv ../data/input-data/targets-expansion.txt ../data/input-data/kegg-pathway-mapping.txt > ../data/features/Dream10-drugs-commonpathway-matrix.txt 

# calculate drug-drug pair  topological features using synergistic drug network
python buildSynergticDrugCombNetwork.py ../data/input-data/ch1_train_combination_and_monoTherapy.csv > ../data/features/Dream10-drugs-topo-sim.txt 


# create protein domain features 
python createProteinDomainFeatures.py ../data/input-data/Drug_info_release.csv ../data/input-data/targets-expansion.txt ../data/input-data/uniprot-hugo-mapping.txt ../data/input-data/domainAnalysis_Pfam.csv ../data/input-data/domainAnalysis_SMART.csv ../data/input-data/domainAnalysis_ProSiteProfiles.csv ../data/input-data/domainAnalysis_ProSitePatterns.csv ../data/input-data/domainAnalysis_SUPERFAMILY.csv > ../data/features/Dream10-drugs-domains.txt

# creatre Drug combinations' features
python createCombFeatures.py ../data/input-data/ch1_train_combination_and_monoTherapy.csv  ../data/features/Dream10-drugs-commontarget-matrix.txt ../data/features/Dream10-drugs-commonpathway-matrix.txt  ../data/features/Dream10-drugs-fingerprint-matrix.txt ../data/features/Dream10-drugs-domains.txt ../data/features/Dream10-drugs-topo-sim.txt > ../data/features/ch1-train-drugfeatures-1A.txt


# mean gene expression value for each clustered module  
python createMeanExprFeatureForGeneModule.py ../data/input-data/gex.csv ../data/input-data/Genes\ in\ Modules  > ../data/features/Dream10-cell-gex-module-mean.txt


# mutations filtered by genes in cancer pathway and its expression varying significanty 
python filterMutationsByCancerPathway.py ../data/input-data/gex.csv ../data/input-data/mutations.csv ../data/input-data/kegg-cancerpathway-mapping.txt > ../data/features/Dream10-cell-mutations-cancerpathway-significant.txt 

# cnv filtered by genes correlated its gene expression with copy numbers significantly (above the mean correlation) and  involved in the KEGG cancer pathway
 python createCNVFeaturesByGeneLevel.py ../data/input-data/gex.csv ../data/input-data/cnv_gene.csv ../data/input-data/kegg-cancerpathway-mapping.txt  > ../data/features/Dream10-cell-cnv-genes-cancerpathway.txt 

# cell clinical features
python findCoverageOfTermsInOnto.py ../data/input-data/cellosaurus.txt ../data/input-data/cell_info.csv 


python createCombFeatures.py ../data/input-data/ch1_train_combination_and_monoTherapy.csv  ../data/features/Dream10-drugs-commontarget-matrix.txt ../data/features/Dream10-drugs-commonpathway-matrix.txt  ../data/features/Dream10-drugs-fingerprint-matrix.txt ../data/features/Dream10-drugs-domains.txt ../data/features/Dream10-drugs-topo-sim.txt > ../data/features/ch1-train-drugfeatures-1A.txt


python createCellLineFeatures.py ../data/features/ch1-train-drugfeatures-1A.txt ../data/features/Dream10-cell_phenotype_features.txt ../data/features/Dream10-cell-gex-module-mean.txt   ../data/features/Dream10-cell-mutations-cancerpathway-significant.txt ../data/features/Dream10-cell-cnv-genes-cancerpathway.txt > ../data/features/ch1-train-cell+drugfeatures-1A.txt


python createCombFeatures.py ../data/input-data/ch1_leaderBoard_monoTherapy.csv ../data/features/Dream10-drugs-commontarget-matrix.txt ../data/features/Dream10-drugs-commonpathway-matrix.txt  ../data/features/Dream10-drugs-fingerprint-matrix.txt ../data/features/Dream10-drugs-domains.txt ../data/features/Dream10-drugs-topo-sim.txt > ../data/features/ch1-leaderboard-drugfeatures-1A.txt

python createCellLineFeatures.py ../data/features/ch1-leaderboard-drugfeatures-1A.txt ../data/features/Dream10-cell_phenotype_features.txt ../data/features/Dream10-cell-gex-module-mean.txt   ../data/features/Dream10-cell-mutations-cancerpathway-significant.txt ../data/features/Dream10-cell-cnv-genes-cancerpathway.txt > ../data/features/ch1-leaderboard-cell+drugfeatures-1A.txt





