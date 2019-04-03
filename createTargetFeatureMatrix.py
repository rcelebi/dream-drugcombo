import sys
import re
from csv import reader

def expandTargetList(targetList, uncertainTargets):
	newTargetList=[]
	for t in targetList:
		if t in uncertainTargets:
			newTargetList.extend(uncertainTargets[t])
		else:
			newTargetList.append(t)
	return newTargetList

if __name__== '__main__':
	#if sys.argv is None or len(sys.argv) is not 2:
	#	print "Usage : python createTargetFeatureMatrix.py Drug_info_release.csv targets-extension.txt "
	#	exit()
	
    drugdict =dict()
	
    uncertainTargets=dict()
    with open(sys.argv[2]) as targetfile:
        for row in targetfile:
            row = row.strip().split("\t")	
            gene = row[0]
            extensions = row[1].split(",")
            uncertainTargets[gene]=extensions	


	
	#infile = file("../data/DREAM10/Drug_info_release.csv")
	
	# the protein targets are listed with '*' denoting any character 

    targetdict =dict()
    alltargets =[]
    targetFreq =dict()
    with open(sys.argv[1]) as infile:
        header= next(infile)
        for row in reader(infile):
            drugid = row[0]
            targets = row[1].strip("\"'").replace(" ","").split(",")
            expansions =expandTargetList(targets,uncertainTargets)
            targetdict[drugid]=expansions
            for t in set(expansions):
                alltargets.append(t)
       
    uniqueTargets=sorted(set(alltargets))
    header ="Drug1\tDrug2"
    for t in uniqueTargets:
        header+="\t"+t
    print (header)
	
    for drug1 in sorted(targetdict):
        for drug2 in sorted(targetdict):
            targetList1 = targetdict[drug1]
            targetList2 = targetdict[drug2]
            featureStr = drug1+"\t"+drug2
            for t in uniqueTargets:
                if t in targetList1 and t in targetList2:
                    featureStr+="\t2"
                elif t in targetList1 or t in targetList2:
                    featureStr+="\t1"
                else:
                    featureStr+="\t0"
            print (featureStr)

