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
                
    uniqueTargets=sorted(set(alltargets))
    uniprot2hugo={}
    with open(sys.argv[3], 'r') as mappingfile:
        for l in mappingfile:
             line=l.split("|")
             #print line
             uniprot=line[1]
             text= line[2].split(" ")
             for s in text:
                if s.startswith("GN="):
                    genesymbol=s[3:]
             #genesymbol =genesymbol[:-6]
             uniprot2hugo[uniprot]=genesymbol

    domainDict ={}
    pathFreq={}

    for f in sys.argv[4:]:
        with open(f,'r') as searchfile:
            next(searchfile)
            for l in searchfile:
                l= l.strip().split(",")
                uniprot = l[0]
                domains = l[1].split(" ")
              
                symbol=uniprot2hugo[uniprot]
                #print uniprot,domains,symbol
                if symbol in domainDict:
                    domainDict[symbol].extend(domains)
                else:
                    domainDict[symbol] =domains
   
    domainDict["BRAF_mut"]=domainDict["BRAF"]
    domainDict["BRAF_V600E"]=domainDict["BRAF"]
    #if domainDict.has_key("CD19"):
    domainDict["CD19antibody"]=domainDict["CD19"]
    domainDict["cMET"]=domainDict["MET"]


    drugDomainDict={}
    for drug1 in sorted(targetdict):
        if drug1 not in drugDomainDict:
            targetList1 = targetdict[drug1]
            pathwayList1=[] 
            for t in targetList1:
                if t in domainDict:
                    pathwayList1.extend(domainDict[t])
            drugDomainDict[drug1]=set(pathwayList1)
            #print drug1,pathwayList1

    alldomains =[]
    for drug in drugDomainDict:
        alldomains.extend(drugDomainDict[drug])  
   
    uniquedomains = set(alldomains)
    uniquedomains=sorted(uniquedomains) 
              
    header ="Drug1\tDrug2"
    for p in uniquedomains:
        #header+="\t"+p+"\t"+p
        header+="\t"+p
    print (header)
    
    for drug1 in sorted(targetdict):
        pathwayList1 =drugDomainDict[drug1]
        #print drug1,pathwayList1
        for drug2 in sorted(targetdict):
            if drug1 == drug2: continue
            
            featureStr = drug1+"\t"+drug2
            pathwayList2 =drugDomainDict[drug2]
            for t in uniquedomains:
                if t in pathwayList1 and t in pathwayList2:
                    featureStr+="\t1"
                else:
                    featureStr+="\t0"
            print (featureStr)

