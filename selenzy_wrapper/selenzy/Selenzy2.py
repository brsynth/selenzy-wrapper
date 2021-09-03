import csv
import pandas as pd
import Selenzy 
import os
def superTax2(tax):
    taxNodes = {}
    for org in tax:
        parents = set(tax[org])
        rest = parents-set(taxNodes)
        if rest:
            for x in rest:
                taxNodes[x]=tax[org][tax[org].index(x):]
    return taxNodes
def taxDistsance2021(taxNodes, host, target):
    if host in taxNodes and target in taxNodes: 
        hostLineage = set(taxNodes[host])
        targetLineage = set(taxNodes[target])
        listaTarget= taxNodes[target]
        if  host in targetLineage:
            distance1 = -(listaTarget.index(host))
            vector1=[distance1,0,distance1]
            return vector1
        else:      
            distance2 = len(hostLineage ^ targetLineage)
            conjunto = hostLineage & targetLineage  
            conjuntoHostAnc = hostLineage - conjunto
            distanciaHostAnc = len(conjuntoHostAnc) 
            conjuntoTarget = targetLineage - conjunto
            distanciaTargetAnc = len(conjuntoTarget)
            vector2 = [distance2, distanciaHostAnc,distanciaTargetAnc]
            return vector2
    else:
        vector3=['-','-','-']
        return vector3
def analyse2(df,host,taxNodes,datadir,outfile = None):
    fileSeqOrg = "seq_org.tsv"
    seqorg = Selenzy.seqOrganism(datadir, fileSeqOrg)
    for rowi in df.index: 
          seqid = df.loc[rowi, 'Seq. ID'] 
          taxid = seqorg[seqid][0]  
          for org in host:
              distancia = taxDistsance2021(taxNodes, org, taxid)
              print(org,taxid,distancia) 
              df.loc[rowi,org+'_target_host_dist'] = distancia[0]
              df.loc[rowi,org+ ' ' '_host_ca_dist'] = distancia[1]
              df.loc[rowi,org+ ' ' '_target_ca_dist'] = distancia[2]
          [targetID,nameID] = seqorg.get(seqid)
          df.loc[rowi, 'target_ID']= targetID
          #df.loc[rowi, 'Source'] = score
          if outfile is not None:
              
              df.to_csv(outfile)
    return df