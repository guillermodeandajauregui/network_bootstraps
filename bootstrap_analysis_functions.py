# -*- coding: utf-8 -*-
"""
Created on Fri Jul 29 11:30:36 2016

@author: gdj
"""

def graphfromsif(sif):
    s = pd.read_table(sif,
                      header =None,
                      delim_whitespace=True)
    s = s[[0,2,1]]
    s = s.values.tolist()
    def formatfunction(lista):
        return "{} {} {}".format(lista[0], lista[1], lista[2])
    reformat = []
    for elem in s:
        reformat.append(formatfunction(elem))
    g = nx.parse_edgelist(reformat, 
                           nodetype = str, 
                           data=(('weight',float),)
                           )
    return(g)
#MI p-value relationship
N = ""

def pvalue(mi, n):
    alfa = 1.062
    beta = -48.7
    gamma = -0.634
    #MI = (alfa - logP) / (-beta + (-gamma * n))
    p = math.exp(alfa -mi*(-beta + (-gamma * n)))
    return(p)

def threshold_subset(g, threshold):
    
    edges= [(from_node,to_node,edge_attributes) for from_node,to_node,edge_attributes in g.edges(data=True) if edge_attributes['weight'] > threshold]
    gg = nx.Graph()
    gg.add_edges_from(edges)
    return(gg)
