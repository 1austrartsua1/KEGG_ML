from Bio.KEGG.KGML.KGML_parser import read
import pandas as pd 
from geneSetUtils import getGeneLists
import graphviz
import os
from graph import *

def plotGraph(viz_dir,name,colors,graph):

    dot = graphviz.Digraph(name)  

    for parent in graph:
        dot.node(str(parent))
    usedRels = set()
    for parent in graph:
        for child in graph[parent]:
            for rel in graph[parent][child]:
                usedRels.add(rel)
                dot.edge(str(parent),str(child),color=colors[rel])
    
    legendStr = "<<TABLE BORDER=\"0\" CELLBORDER=\"1\" CELLSPACING=\"0\" CELLPADDING=\"4\">"
    
    for rel in usedRels:
        legendStr += "<TR><TD>"+str(rel)+"</TD><TD BGCOLOR=\""+str(colors[rel])+"\"></TD></TR>"
    legendStr += "</TABLE>>"
    
    dot.node('struct1', legendStr,shape='plaintext',fontsize='30')

    dot.render(directory=viz_dir)  

if __name__ == "__main__":
    lsdir = os.listdir('./KEGG_pathways/')
    pathwayNames = os.listdir('./KEGG_pathways/')
    pathwayNames = [p[:-4] for p in pathwayNames]
    #pathwayNames = ['hsa04630']
    data = pd.read_csv('GSE/GSE43151_gs.csv')
    gse_genesets = pd.read_csv('GSE/geneSets.tsv')
    
    filteredGenes_all = getGeneLists(data,gse_genesets,'filtered')
    

    #acceptedRelations = {('PPrel', 'activation'),('PPrel', 'inhibition'), ('PPrel', 'binding/association'),('sameNode',None)
    #                     ('ECrel', 'compound'),('PPrel', 'dissociation'),('GErel', 'expression'),('GErel', 'repression')}
    acceptedRelations = {('PPrel', 'activation'),('PPrel', 'inhibition'),('sameNode',None)}
    
    myColors = ['blue','red','green','black','purple','yellow','magenta','orange']
    colors = {}
    for i,relType in enumerate(acceptedRelations):
        colors[relType] = myColors[i]
    
    dags = []
    for pathwayName in pathwayNames:
        print(pathwayName)
        pathway = read(open('KEGG_pathways/'+pathwayName+'.xml', 'r'))
        filteredGenes = set(filteredGenes_all[pathwayName])
        graph,relationTypes = getGraph(pathway,filteredGenes)
        
        _,graph = cleanGraph(graph,acceptedRelations)
        if len(graph) == 0:
            continue 
        plotGraph('viz/',pathwayName+' graph',colors,graph)
        bn = BN_Graph(graph)
        order = bn.topoSort()
        if order is not None:
            #found DAG
            dags.append(pathwayName)
    print('dags:',dags)


'''
    #unfortunately it seems kegg does not allow you to programmatically
    #download xml files so we need to do it by hand
    # the below will not work
    except:
        print('kegg xml file not opened, trying to download...')

        url = 'https://www.kegg.jp/kegg-bin/download?entry='+pathwayName+'&format=kgml'
        r = requests.get(url, allow_redirects=True)
        open('KEGG_pathways/'+pathwayName+'.xml', 'wb').write(r.content)

        pathway = read(open('KEGG_pathways/'+pathwayName+'.xml', 'r'))
'''