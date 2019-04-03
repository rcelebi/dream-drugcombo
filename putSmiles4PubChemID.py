import sys
import re
import pandas as pd

# replace Pubchemid with Smiles, return new dataframe with imputed Smiles
def replacePubChemidWithSmiles(pubchem_df, drug_df):
    drugPubchemDict = { str(x[0]):x[1] for x in pubchem_df[['pubchemid','smilesid']].values}
    
    def replaceIDS(pubchem):
        pubchem = str(pubchem)
        if pubchem.isdigit():        
            if  pubchem in drugPubchemDict:
                    return drugPubchemDict[pubchem]
        return pubchem
    
    drug_df['SMILES'] = drug_df['SMILES or PubChem ID'].apply(replaceIDS)
    return  drug_df

if __name__== '__main__':
	#if sys.argv is None or len(sys.argv) is not 2:
	#	print "Usage : python convertDrug.. in_file "
	#	exit()

    pubchem_df = pd.read_csv(sys.argv[1],delimiter='\t')
    #print(pubchem_df.head())
    drug_df = pd.read_csv(sys.argv[2],delimiter=',')
    drug_df =replacePubChemidWithSmiles(pubchem_df, drug_df)
    drug_df[ drug_df['SMILES'] != 'nan'][['ChallengeName','SMILES']].to_csv('drug_smiles.csv', index=False)
