import os
import numpy as np
from csv import reader
import pandas as pd 
import sys


if __name__== '__main__':

	#path="/Users/remzicelebi/Dropbox/Remzi & Semih Folder/DREAM CHALLENGE/outputFiles/Genes in Modules"
	geneModules={}
	path=sys.argv[2]
	for a in os.listdir(path):
		if not a.endswith(".txt"): continue
		if os.path.isdir(path+"/"+a): continue
		genes =[]
		f = open(path+"/"+a)
		next(f)
		module=a.split()[0]
		for l in f:
			l =l.strip().split("\t")
			genes.append(l[0])
		geneModules[module]=genes


	#pathGex="/Users/remzicelebi/projects/dumontierlab/drugcomb-challenge/data/input-data/gex.csv"
	pathGex=sys.argv[1]
	data = pd.read_csv(pathGex)
	data.set_index('Genes')
	
	temp=data.T
	datanew=temp.rename(columns= data.Genes)
	cellLines = set(temp.columns)
	datanew =datanew[1:]

	for module in geneModules:
		datanew[module]=datanew[geneModules[module]].mean(axis=1)

		cellLines=set(datanew.index)

	header="CellLines"
	for module in geneModules:
		header+="\t"+module
	print (header)

	for cell in cellLines:
		row=cell
		for module in geneModules:
			row+="\t"+str(datanew[module][cell])
		print (row)
				