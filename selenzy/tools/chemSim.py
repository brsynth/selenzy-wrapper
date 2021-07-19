from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import GetMorganFingerprint, GetConnectivityInvariants
from rdkit import DataStructs
from os import path
import cPickle as pickle

def getStructs(dbfile):
    structs = {}
    for l in open(dbfile):
        if l.startswith('#'):
            continue
        m = l.rstrip().split('\t')
        cid = m[0]
        smiles = m[6]
        if len(smiles) > 0:
            structs[cid] = smiles
    return structs

RECOMPUTE = False
if RECOMPUTE:
    radius = 5
    mol = {}
    fp = []
    fpNames = []
    cstr = getStructs(path.join('data', 'chem_prop.tsv'))
    for c in cstr:                  
        try:
            mol[c] = Chem.MolFromSmiles(cstr[c])
        except:
            continue
    pickle.dump(mol, open(path.join('data', 'mnxMol.pk'), 'w'))


    for c in mol:
        try:
            fp.append( GetMorganFingerprint(mol[c], radius) )
        except:
            continue
        fpNames.append(c)
    f = open(path.join('data', 'mnxFp.pk'), 'w')
    pickle.dump(fp, f)
    pickle.dump(fpNames, f)
    f.close()

else:
    print('Reading fingerprints...')
    f = open(path.join('data', 'mnxFp.pk'))
    fp = pickle.load(f)
    fpNames = pickle.load(f)
    f.close()
print('Computing similarity...')
tn = DataStructs.BulkTanimotoSimilarity(fp[2], fp)
print('Sorting...')
k = 0
for i in sorted(range(0, len(tn)), key = lambda j: -tn[j] ):
    print(fpNames[i], tn[i])
    k += 1
    if k > 10:
        break
