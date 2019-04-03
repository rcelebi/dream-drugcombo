#!/usr/bin/env python
import os
import re  
import sys
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
    #    print "Usage : python convertDrug.. in_file "
    #    exit()
    
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

    allpathways =[]
    genePathwayDict ={}
    uniqueTargets=sorted(set(alltargets))
   
    with open(sys.argv[3], 'r') as searchfile:
        for l in searchfile:
            l= l.strip().split("\t")
            genesymbol = l[0]
            pathway = l[1].split(",")
            #print pathway,genesymbols
            symbol=genesymbol
            if symbol in uniqueTargets:
                if symbol in genePathwayDict:
                    genePathwayDict[symbol].extend(pathway)
                else:
                    genePathwayDict[symbol] =pathway
                allpathways.extend(pathway)   
            
    uniquePathways = set(allpathways) 
              
    header ="Drug1\tDrug2"
    for p in uniquePathways:
        #header+="\t"+p+"\t"+p
        header+="\t"+p
    print (header)
    uniquePathways=sorted(uniquePathways)
    for drug1 in sorted(targetdict):
        for drug2 in sorted(targetdict):
            targetList1 = targetdict[drug1]

            pathwayList1=[] 
            for t in targetList1:
                if t in genePathwayDict:
                    pathwayList1.extend(genePathwayDict[t])

            targetList2 = targetdict[drug2]
            pathwayList2=[] 
            for t in targetList2:
                if t in genePathwayDict:
                    pathwayList2.extend(genePathwayDict[t])

            featureStr = drug1+"\t"+drug2
            pathwayList1=set(pathwayList1)
            pathwayList2=set(pathwayList2)
            for t in uniquePathways:
                if t in pathwayList1 and t in pathwayList2:
                    featureStr+="\t1"
                else:
                    featureStr+="\t0"
                
            print (featureStr)

