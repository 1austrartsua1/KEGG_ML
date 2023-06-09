from Bio.KEGG.KGML.KGML_parser import read
import pandas as pd 
from geneSetUtils import *
import copy
import numpy as np 

def getEntries(names):
    ids = []
    for ent in names:
        if ent == 'undefined':
            continue
        species,entid = ent.split(":")
        if species != 'hsa':
            continue
        ids.append(int(entid))
    return ids

def getGraph(pathway,filteredGenes):
    filtered_relations = 0
    graph = {}
    relationTypes = {}
    GSE_seen = set()

    for relation in pathway.relations:
        parentIds = getEntries(relation.entry1._names)
        childIds = getEntries(relation.entry2._names)
        
        filt_Parent = set(parentIds) & filteredGenes
        filt_Child = set(childIds) & filteredGenes
        
        if (len(filt_Parent) == 0) or (len(filt_Child) == 0):
            filtered_relations += 1
            continue
        
        for ind in filt_Parent.union(filt_Child):
            GSE_seen.add(ind)

        for parent in filt_Parent:
            if not graph.get(parent):
                graph[parent] = {}
            for child in filt_Child:
                if not graph[parent].get(child):
                    graph[parent][child] = set()
                graph[parent][child].add((relation.type,relation.subtypes[0][0]))
                
        # add edges between all parents and all children
        for parent1 in filt_Parent:
            for parent2 in filt_Parent:
                if parent1 == parent2:
                    continue 
                if not graph[parent1].get(parent2):
                    graph[parent1][parent2] = set()
                graph[parent1][parent2].add(("sameNode",None))

        for child1 in filt_Child:
            for child2 in filt_Child:
                if child1 == child2:
                    continue 
                if not graph.get(child1):
                    graph[child1] = {}
                if not graph[child1].get(child2):
                    graph[child1][child2] = set()
                graph[child1][child2].add(("sameNode",None))

        if not relationTypes.get((relation.type,relation.subtypes[0][0])):
            relationTypes[(relation.type,relation.subtypes[0][0])] = 0
        relationTypes[(relation.type,relation.subtypes[0][0])] += 1
        
    print('filtered relations / total relations:',filtered_relations,'/',len(pathway.relations))
    print('GSE genes seen / total',len(GSE_seen),'/',len(filteredGenes))

    return graph,relationTypes 



def cleanGraph(graph,acceptedRelations):
    numChildren = []

    
    for parent in graph:
        numChildren.append(len(graph[parent]))
        toDel = []
        for child in graph[parent]:
            if len(graph[parent][child])>1:
                print('Multiple relationships detected')
                print(graph[parent][child])
            toRemove = []
            for rel in graph[parent][child]:
                if (not rel in acceptedRelations):
                    toRemove.append(rel)
            
            for rel in toRemove:
                graph[parent][child].remove(rel)

        if len(graph[parent][child]) == 0:
            del graph[parent][child]
        
    numChildren.sort()
    return numChildren



class BN_Graph:
    def __init__(self,graph):
        self.node2children = {}
        self.node2parents = {}
        node_list = set()
        for parent in graph:
            node_list.add(parent)
            for child in graph[parent]:
                node_list.add(child)

        for node in node_list:
            self.node2parents[node] = []
            self.node2children[node] = []
        self.node_list = node_list

        for parent in graph:
            self.node2children[parent] = list(graph[parent].keys())
            for child in graph[parent]:
                self.node2parents[child].append(parent)
                
                
        
    def getParentList(self,order):
        adj = np.zeros((len(order),len(order)))
        for i,gene in enumerate(order):
            parentsi = self.node2parents[gene]
            for parent in parentsi:
                for j in range(len(order)):
                    if order[j] == parent:
                        adj[i,j] = 1.
                        break
        return adj 

    def topoSort(self):
        '''
        graph.node2parents: dic mapping node to a list of parents
        graph.node2children: dic mapping node to a list of children
        algorithm:
        start: loop through node2parents and find a node with no parents
        add parent to output list
        loop through its children
        for each child, remove parent from its parent list
        remove parent from node2parents
        loop to start
        '''
        node_list = copy.deepcopy(self.node_list)
        node2parents = copy.deepcopy(self.node2parents)
        node2children = copy.deepcopy(self.node2children)

        order = []
        while len(node_list)>0:
            cycleDet = True
            for a_root in node_list:
                if len(node2parents[a_root]) == 0:
                    cycleDet = False
                    break

            if cycleDet:
                raise Exception("cycle detected")
            order.append(a_root)
            for child in node2children[a_root]:
                node2parents[child].remove(a_root)

            node_list.remove(a_root)

        return order

if __name__ == "__main__":
    pathwayName = 'hsa04650'

    data = pd.read_csv('GSE/GSE43151_gs.csv')
    gse_genesets = pd.read_csv('GSE/geneSets.tsv')

    filteredGenes_all = getGeneLists(data,gse_genesets,'filtered')
    
    

    filteredGenes = set(filteredGenes_all[pathwayName])
    
    acceptedRelations = {('PPrel', 'activation'),('PPrel', 'inhibition')}

    
    pathway = read(open('KEGG_pathways/'+pathwayName+'.xml', 'r'))
    graph,relationTypes = getGraph(pathway,filteredGenes)
    #print(f"{relationTypes=}")     
    cleanGraph(graph,acceptedRelations)
    bn = BN_Graph(graph)
    order = bn.topoSort()
    print(order)
    adj = bn.getParentList(order)
    
    g2c = gene2col(data)
    order_inds = [g2c[str(gene)] for gene in order]
    labels = np.array(data['Dose'])
    data = data.iloc[:,:-1]
    data = data.to_numpy(np.float64)[:,order_inds]

    print(data.shape)

    means = np.mean(data,axis=0)
    print(means)
    stds = np.std(data,axis=0,ddof=1)
    print(stds)
    
    # check if edge-types matches pairwise covariance/correlation
    
'''
#plot pathway lengths
pwlens = [len(filteredGenes_all[key]) for key in  filteredGenes_all]
pwlens.sort(reverse=True)
import matplotlib.pyplot as plt 
plt.plot(pwlens)
plt.title('filtered pathway sizes')
plt.show()
exit()
'''