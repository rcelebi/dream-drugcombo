from csv import reader
import sys
from oddt import toolkit
from oddt import fingerprints


def getFingerprint(smiles):
    mol = toolkit.readstring(format='smiles',string=smiles)
    fp =mol.calcfp(fptype='MACCS').raw
    return fp


if __name__== '__main__':
	#if sys.argv is None or len(sys.argv) is not 2:
	#	print "Usage : python createFingerprintFeatures.py .. drug_smiles.csv "
	#	exit()
	infile = open(sys.argv[1])
	smilesdict =dict()

	header= next(infile)
	for row in reader(infile):
		drugid = row[0]
		smiles = row[1].split(";")[0]
		if smiles == "": continue 
		smilesdict[drugid] =getFingerprint(smiles)

	
	count=0
	maccsLength =len(smilesdict[drugid])
	header="Drug1\tDrug2"
	for i in range(maccsLength):
		header+="\tFingerprint"+str(i)

	print (header)
	
	for drug1 in smilesdict:
		for drug2 in smilesdict:
			if drug1 == drug2: continue
			fingerprint1=smilesdict[drug1]
			fingerprint2=smilesdict[drug2]
			feature=""
			for i in range(maccsLength):
				if fingerprint1[i] and fingerprint2[i]:
					feature+="\t2"
				elif fingerprint1[i] or fingerprint2[i]:
					feature+="\t1"
				else:
					feature+="\t0"	
			print  (str(drug1)+"\t"+str(drug2)+feature)
