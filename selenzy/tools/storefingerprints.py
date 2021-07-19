

from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import GetMorganFingerprint, GetAtomPairFingerprint, GetTopologicalTorsionFingerprint
from rdkit.Chem.rdmolops import PatternFingerprint, RDKFingerprint
from rdkit import DataStructs
from os import path
import numpy as np
import csv, sys


def reactSMILES2FP(smi, smiles, fps, ffun, param):
    """ Return left / right reaction fingerprint """
    """ The reaction could be replaced by rdkit """
    left, right = smi.split('>>')
    rleft = left.split('.')
    rright = right.split('.')
    ok = True
    mleft = []
    mright = []
    for c in rleft:
        if c not in smiles:
            try:
                smiles[c] = Chem.MolFromSmiles(c)
            except:
                ok = False
                break
        mleft.append((c, smiles[c]))
    if not ok:
        return None
    for c in rright:
        if c not in smiles:
            try:
                smiles[c] = Chem.MolFromSmiles(c)
            except:
                ok = False
                break
        mright.append((c, smiles[c]))
    if not ok:
        return None
    rfp = [None, None]
    for c in mright:
        if c[0] not in fps:
            try:
                if param is not None:
                    if ffun == GetMorganFingerprint:
                        fps[c[0]] = ffun(c[1], radius=param)
                    elif ffun == RDKFingerprint:
                        fps[c[0]] = ffun(c[1], maxPath=param)
                    else:
                        fps[c[0]] = ffun(c[1])
                else:
                    fps[c[0]] = ffun(c[1])
            except:
                ok = False
                break
        if rfp[1] is None:
            rfp[1] = fps[c[0]]
        else:
            rfp[1] = rfp[1] | fps[c[0]]
    if not ok or rfp[1] is None:
        return None
    for c in mleft:
        if c[0] not in fps:
            try:
                if param is not None:
                    if ffun == GetMorganFingerprint:
                        fps[c[0]] = ffun(c[1], radius=param)
                    elif ffun == RDKFingerprint:
                        fps[c[0]] = ffun(c[1], maxPath=param)
                    else:
                        fps[c[0]] = ffun(c[1])
                else:
                    fps[c[0]] = ffun(c[1])
            except:
                ok = False
                break
        if rfp[0] is None:
            rfp[0] = fps[c[0]]
        else:
            rfp[0] = rfp[0] | fps[c[0]]
    if not ok or rfp[0] is None:
        return None
    return rfp



def reactionFingerprint(ffun, fname, param=None, bit=False):
    """ Reaction binary fingerprint based on prod-subs fingerprint logic difference """
    """ rsmifile: precomputed reaction SMILES from METANETX2 """
    csv.field_size_limit(sys.maxsize) # To avoid error with long csv fields
    rsmiFile = path.join('../data', 'reac_smi.csv')
    smiles = {}
    fps = {}
    rfp = {}
    with open(rsmiFile) as f:
        for row in csv.DictReader(f):
            rid = row['RID']
            smi = row['SMILES']
            fp = reactSMILES2FP(smi, smiles, fps, ffun, param)
            if fp is not None:
                rfp[rid] = fp

    fpNames = sorted(rfp)
    if bit:
        fp = [rfp[x].ToBitString() for x in fpNames]
    else:
        fp = [rfp[x] for x in fpNames]
    f = np.savez_compressed(fname, x=fp, y=fpNames)
    return rfp


def getReactants(dbfile):
    clist = set()
    for line in open(dbfile):
        if line.startswith('#'):
            continue
        row= line.rstrip().split('\t')
        for x in row[1].split(' '):
            if x.startswith('MNXM'):
                clist.add(x)
    return clist

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

def getMols():
    mol = {}
    clist = getReactants(path.join('../data', 'reac_prop.tsv'))
    cstr = getStructs(path.join('../data', 'chem_prop.tsv'))
    for c in set(cstr) & clist:                  
        try:
            mol[c] = Chem.MolFromSmiles(cstr[c])
        except:
            continue
    return mol

def storeFingerprint(mol, ffun, fname, param=None, bit=False):
    fp = []
    fpNames = []
    print('Computing fingerprints...')
    for c in mol:
        try:
            if param is not None:
                fp.append( ffun(mol[c], param) )
            else:
                if bit:
                    fp.append( ffun(mol[c]).ToBitString()  )
                else:
                    fp.append( ffun(mol[c]) )
        except:
            continue
        fpNames.append(c)
    print('Saving...') 
    f = np.savez_compressed(fname, x=fp, y=fpNames)

def testPattern(ptfile, bit=False):
    """ Test how to reload PatternFingerprint """
    print('Validating fingerprint...')
    data = np.load(ptfile)
    fps = data['x']
    fpNames = data['y']
    if bit:
        fp = [DataStructs.CreateFromBitString(z) for z in fps]
    else:
        fp = fps
    sim = DataStructs.BulkTanimotoSimilarity(fp[0], list(fp))
    return fp


print('Pattern fingerprint....')
reactionFingerprint(PatternFingerprint, 'ptrfp.npz', bit=True)
print('RDK fingerprint....')
for radius in range(1,11):
    reactionFingerprint(RDKFingerprint, 'rdkrfp%d.npz' % (radius,), param=radius, bit=True)
print('Morgan fingerprint....')
for radius in range(1,11):
    reactionFingerprint(GetMorganFingerprint, 'mgrfp%d.npz' % (radius,), param=radius)

sys.exit()

mol = getMols()
print('Pattern fingerprint....')
storeFingerprint(mol, PatternFingerprint, 'ptfp.npz', bit=True)
fp = testPattern('ptfp.npz', bit=True)          
print('RDK fingerprint....')
storeFingerprint(mol, RDKFingerprint, 'rdkfp.npz', bit=True)
fp = testPattern('rdkfp.npz', bit=True)
print('Atom pair fingerprint....')          
storeFingerprint(mol, GetAtomPairFingerprint, 'apfp.npz')
fp = testPattern('apfp.npz')
print('Topological torsion fingerprint....')          
storeFingerprint(mol, GetTopologicalTorsionFingerprint, 'ttfp.npz')
fp = testPattern('ttfp.npz')
print('Morgan fingerprint....')          
for radius in range(1,11):
     storeFingerprint(mol, GetMorganFingerprint, 'mgfp%d.npz' % (radius,) , radius)
     fp = testPattern('mgfp%d.npz' % (radius,))
