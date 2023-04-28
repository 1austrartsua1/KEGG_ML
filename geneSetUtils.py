import pandas as pd
import copy

def gene2col(data):
    gene2col = {}
    for i,gene in enumerate(data.columns):
        gene2col[gene] = i
    return gene2col

def getGeneLists(data,pathways,filterType):
    colNames = set(data.columns)
    
    pathwayDic = {}
    for pathway in pathways.values:
        temp = copy.deepcopy(pathway)
        temp = temp[0]
        temp = temp.split('\t')
        key = temp[0]
        pathwayDic[key] = []
        for val in temp[1:]:
            pathwayDic[key].append(int(val))

    if filterType == 'unfiltered':
        return pathwayDic
    #####################################
    #filter genes
    #PJ: filter out genes that we don't have expression values for
    #####################################

    removePath = []
    for key in pathwayDic:
        geneLst = []
        for gene in pathwayDic[key]:
            if str(gene) in colNames:
                geneLst.append(gene)

        if geneLst:
            pathwayDic[key] = geneLst
        else:
            removePath.append(key)

    for key in removePath:
        del pathwayDic[key]

    return pathwayDic

if __name__ == "__main__":
    data = pd.read_csv('GSE/GSE43151_gs.csv')
    

    pathways = pd.read_csv('GSE/geneSets.tsv')
    pathwayDic = getFilteredGeneLists(data,pathways)
    print(pathwayDic['hsa04650'])
