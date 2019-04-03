import csv
import numpy as np
import pandas as pd
import math
import itertools
import networkx as nx
import time

import findspark
findspark.init()

from pyspark import SparkConf, SparkContext

    

def networkx_jaccard(G, u, v):
    union_size = len(set(G[u]) | set(G[v]))
    if union_size == 0:
        return 0
    return len(list(nx.common_neighbors(G, u, v))) / union_size

def netFeatureSet(G, u, v, cls):
    if u in G.nodes and v in G.node: 
        CN = len(list(nx.common_neighbors(G, u, v)))
        JC = networkx_jaccard(G, u, v)
        AA = sum(1 / math.log(G.degree(w)) for w in nx.common_neighbors(G, u, v))
    else:
        CN = 0.0
        JC = 0.0
        AA = 0.0
    return u, v, cls, JC, AA, CN

def linkPredSpark(bc_triples, pairs, classes):

    pairList = list(zip(pairs[:,0],pairs[:,1],classes))

    start_time =time.time()
    rdd = sc.parallelize(pairList).map(lambda x: netFeatureSet(bc_triples.value, x[0], x[1], x[2] ))
    pair_df = rdd.collect()
    elapsed_time = time.time() - start_time
    #print ('Time elapsed to generate walks:',time.strftime("%H:%M:%S", time.gmtime(elapsed_time)))
    return pd.DataFrame(pair_df, columns=['Drug1','Drug2','Class','JC','AA','CN'])



if __name__== '__main__':

    combo_df = pd.read_csv("data/input/ch1_train_combination_and_monoTherapy.csv", sep=',')
    combo_lb1_df = pd.read_csv("data/input/ch1_LB.csv", sep=',')
    combo_lb2_df = pd.read_csv("data/input/ch2_LB.csv", sep=',')
    
    combo_df = combo_df.append(combo_lb1_df)
    combo_df = combo_df.append(combo_lb2_df)
	
    ddiKnown = set()
    # synergistic combinations
    for g,df in combo_df.groupby('COMBINATION_ID'):
        n = len(df)
        n1 = len(df[df.SYNERGY_SCORE >20])
        if float(n1)/n > 0.5:
            #print (g,n,n1)
            ddiKnown.add(tuple(g.split('.')))
            
    drugs = set(combo_df.COMPOUND_A.unique())
    drugs = drugs.union(combo_df.COMPOUND_B.unique())
    
    print (len(drugs))
    
    pairs = []
    classes = []
    for comb in itertools.combinations(sorted(drugs), 2):
        dr1 = comb[0]
        dr2 = comb[1]
        
        
        if (dr1,dr2)  in ddiKnown or  (dr2,dr1)  in ddiKnown:
            cls=1
        else:
            cls=0
        pairs.append((dr1,dr2))
        classes.append(cls)

    pairs = np.array(pairs)        
    classes = np.array(classes)

    indices = np.where(classes == 1)
    positives = pd.DataFrame(list(zip(pairs[indices][:,0],pairs[indices][:,1],classes[indices])), columns=['Drug1','Drug2','Class'])

    indices = np.where(classes == 0)
    all_negatives = pd.DataFrame(list(zip(pairs[indices][:,0],pairs[indices][:,1],classes[indices])), columns=['Drug1','Drug2','Class'])

    print("DDI size: ", len(positives))
    print("non-DDI size: ",len(all_negatives))


    pair_df = pd.concat([positives,all_negatives], ignore_index=True)
    print("Train size: ", len(pair_df))
    
    
    weights=dict()
    G = nx.Graph()
    for index, row in pair_df[pair_df.Class == 1 ].iterrows():
        drug1=row['Drug1']
        drug2=row['Drug2']
        G.add_edge(drug1,drug2)
    
    config = SparkConf()
    config.setMaster("local[10]")
    #config.set("spark.executor.memory", "2g")
    sc = SparkContext(conf=config)
    print (sc)
    bc_net=sc.broadcast(G)
    
    
    all_pairs = linkPredSpark(bc_net, pair_df[['Drug1','Drug2']].values, pair_df.Class)
    all_pairs.drop(columns='Class', inplace=True)
    all_pairs1 = all_pairs.copy()
    
    all_pairs1['Drug1'] =all_pairs['Drug2']
    all_pairs1['Drug2'] =all_pairs['Drug1']
    all_pairs = all_pairs.append(all_pairs1)
    
    outfile ='data/features/drug-synergy-network-features.txt'
    print ('Output file was stored in %s'%(outfile))
    all_pairs.to_csv(outfile, sep='\t', index=False)