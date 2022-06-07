import sys, os, re
import numpy as np
from scipy.spatial import distance
from scipy.spatial.distance import cdist
import pandas as pd
from rdkit import Chem
from rdkit.ML.Descriptors.MoleculeDescriptors import MolecularDescriptorCalculator
import fileinput
from pharma import pharma
import featureSASA
import xgboost as xgb
import pickle
from openbabel import openbabel as ob
from openbabel import pybel
import alphaspace2 as al
import mdtraj


import calc_bridge_wat
import calc_ligCover_betaScore
import calc_rdkit
import calc_sasa
import calc_vina_features
import prepare_betaAtoms

Vina = '/home/cyang/paper_XGB/delta_LinF9_XGB/software/smina_feature'
Smina = '/home/cyang/paper_XGB/delta_LinF9_XGB/software/smina.static'
SF = '/home/cyang/paper_XGB/delta_LinF9_XGB/software/sf_vina.txt'
ADT = '/home/cyang/MGLTools-1.5.6/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_receptor4.py'
model_dir = '/home/cyang/paper_XGB/delta_LinF9_XGB/saved_model'

def run_XGB(pro, lig):

    if lig.endswith('.mol2'):
        mol = Chem.MolFromMol2File(lig, removeHs=False)
        lig = lig[:-5]+'.pdb'
        Chem.MolToPDBFile(mol, lig)
        
    elif lig.endswith('.sdf'):
        mol = Chem.MolFromMolFile(lig, removeHs=False)
        lig = lig[:-4]+'.pdb'
        Chem.MolToPDBFile(mol, lig)
        
    elif lig.endswith('.pdb'):
        mol = Chem.MolFromPDBFile(lig, removeHs=False)

    ## 1. prepare_betaAtoms
    beta = os.path.join(os.path.dirname(pro), 'betaAtoms.pdb')
    pro_pdbqt = prepare_betaAtoms.Prepare_beta(pro, beta, ADT)

    ## 2. Vina_features
    v = calc_vina_features.vina(pro_pdbqt, lig, Vina, Smina)
    vinaF = [v.LinF9]+v.features(48)

    ## 3. Beta_features
    betaScore, ligCover = calc_ligCover_betaScore.calc_betaScore_and_ligCover(lig, beta)

    ## 4. sasa_features
    datadir = os.path.dirname(os.path.abspath(pro))
    pro_ = os.path.abspath(pro)
    lig_ = os.path.abspath(lig)
    sasa_features = calc_sasa.sasa(datadir,pro_,lig_)
    sasaF = sasa_features.sasa+sasa_features.sasa_lig+sasa_features.sasa_pro

    ## 5. ligand_features
    ligF = list(calc_rdkit.GetRDKitDescriptors(mol))

    ## 6. water_features
    df = calc_bridge_wat.Check_bridge_water(pro, lig)
    if len(df) == 0:
        watF = [0,0,0]
    else:
        Nbw, Epw, Elw = calc_bridge_wat.Sum_score(pro, lig, df, Smina)
        watF = [Nbw, Epw, Elw]

    ## calculate XGB
    LinF9 = vinaF[0]*(-0.73349)
    LE = LinF9/vinaF[-4]
    sasa = sasaF[:18]+sasaF[19:28]+sasaF[29:]
    metal = vinaF[1:7]
    X = vinaF[7:]+[ligCover,betaScore,LE]+sasa+metal+ligF+watF
    X = np.array([X]).astype(np.float64)

    y_predict_ = []
    for i in range(1,11):
        xgb_model = pickle.load(open("%s/mod_%d.pickle.dat"%(model_dir,i),"rb"))
        y_i_predict = xgb_model.predict(X, ntree_limit=xgb_model.best_ntree_limit)
        y_predict_.append(y_i_predict)

    y_predict = np.average(y_predict_, axis=0)
    XGB = round(y_predict[0]+LinF9,3)
    
    return XGB

def main():
    args = sys.argv[1:]
    if not args:
        print ('usage: python runXGB.py pro lig')

        sys.exit(1)

    elif sys.argv[1] == '--help':
        print ('usage: python runXGB.py pro lig')

        sys.exit(1)

    elif len(args) == 2 and sys.argv[1].endswith('.pdb') and sys.argv[2].endswith(('.pdb','.mol2','sdf')):
        pro = sys.argv[1]
        lig = sys.argv[2]
        XGB = run_XGB(pro, lig)
        print ('XGB (in pK) : ', XGB)
        
    else:
        sys.exit(1)

if __name__ == '__main__':
    main()
