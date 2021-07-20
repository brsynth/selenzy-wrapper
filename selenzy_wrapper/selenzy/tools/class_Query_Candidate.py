#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 12 16:05:21 2017

@author: jerrywzy
"""
class Candidate:
    def __init__(self, up, desc, os, rxnid, cn, repid, cons, rxnsim, du, dp, h, s, t, c, mw, ip, pol):
        self.uniprot = up
        self.description = desc
        self.organism = os 
        self.reactionid = rxnid
        self.clusterno = cn
        self.representativeid = repid
        self.conservation = cons
        self.reactionsim = rxnsim
        self.directionused = du
        self.directionpref = dp
        self.helices = h
        self.sheets = s
        self.turns = t
        self.coils = c
        self.molweight = mw
        self.isoelectricpoint = ip
        self.polarity = pol
    
    def printCand(self):
        print(self)
    
class Query:
    def __init__(self, reaction, ntargets, header, candidates):
        self.reaction = 'SMILES'           #SMILES format
        self.ntargets = ntargets          # length of table -1 (header)
        self.header = header
        self.candidates = candidates             #array?
    
    def addrow(self, row):
        self.candidates.append(row)
        
    def printQuery(self):
        print (self)
        
if __name__ == '__main__':
    
    header = ['Seq. ID','Description', 'Organism Source', 'Rxn. ID', 'Clust. No.', 'Rep. ID.', 'Consv. Score', 'Rxn Sim.', "Direction Used", "Direction Preferred", '% helices', '% sheets', '% turns', '% coils', 'Mol. Weight', 'Isoelec. Point', 'Polar %']
    
    testquery = Query("myrxn", "20", "UNIPROT")
    
    testcand = Candidate('Q9ZKQ', 'D-aminoacid', 'Helicobacter', 'MNXR75706', "39903", "A3KEZ1", "41", "0.366", "1", "-1", "51.3", "24.6", "13.2", "15", "46101.62	", "8.4254",	"44.39")
    
    testquery.addrow(testcand)
    
    
    print (testquery.candidates)