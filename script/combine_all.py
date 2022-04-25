import sys, os
import numpy as np
import pandas as pd
import xgboost as xgb
import pickle

vinaF_name = ['pdb','Lin_F9']+['metal%d'%i for i in range(1,7)]+['vina%d'%i for i in range(1,43)]

f_type = ["P","N","DA","D","A","AR","H","PL","HA","SA"]
SASA = ["P2." + i for i in f_type] + ["P2dl." + i for i in f_type] + ["P2dp." + i for i in f_type]
sasaF_name = ['pdb']+SASA

betaF_name = ['pdb','betaScore','ligCover']
watF_name = ["pdb","Nbw","Epw","Elw"]

ligF_name = ['pdb','HeavyAtomMolWt', 'NumValenceElectrons','FpDensityMorgan1',
             'FpDensityMorgan2', 'FpDensityMorgan3', 'LabuteASA',
             'TPSA', 'NHOHCount','MolLogP','MolMR'] 

def Combine_features(fea_dir):
    '''
    combine features from vina_features, sasa_features, beta_features, water_features, ligand_features
    '''
    types_dict = {'pdb':str}
    VinaF = pd.read_csv('%s/vina_features.csv'%fea_dir, header=None, names=vinaF_name, dtype = types_dict)
    sasaF = pd.read_csv('%s/sasa_features.csv'%fea_dir, header=None, names=sasaF_name, dtype = types_dict)
    betaF = pd.read_csv('%s/beta_features.csv'%fea_dir, header=None, names=betaF_name, dtype = types_dict)
    watF = pd.read_csv('%s/BW_features.csv'%fea_dir, header=None, names=watF_name, dtype = types_dict)
    ligF = pd.read_csv('%s/lig_features.csv'%fea_dir, header=None, names=ligF_name, dtype = types_dict)

    df = pd.concat([VinaF,sasaF,betaF,watF,ligF], axis=1)
    df = df.loc[:,~df.columns.duplicated()]
    df['Lin_F9'] = df['Lin_F9']*(-0.73349)
    df['LE'] = df['Lin_F9']/df['vina39']

    return df

def Calc_XGB(df):
    '''
    calculate XGB score (in pKd) based on the feature set.
    '''
    metal = ['metal%d'%i for i in range(1,7)]
    vina = ['vina%d'%i for i in range(1,43)]
    BW = ["Nbw","Epw","Elw"]
    ligandDescriptors = ['HeavyAtomMolWt', 'NumValenceElectrons','FpDensityMorgan1',
                         'FpDensityMorgan2', 'FpDensityMorgan3', 'LabuteASA',
                         'TPSA', 'NHOHCount','MolLogP','MolMR']
    f_type = ["P","N","DA","D","A","AR","H","PL","HA","SA"]
    SASA = ["P2." + i for i in f_type] + ["P2dl." + i for i in f_type] + ["P2dp." + i for i in f_type]
    SASA = [x for x in SASA if x not in ['P2dl.HA','P2dp.HA']]

    columns = vina+['ligCover','betaScore','LE']+SASA+metal+ligandDescriptors+BW

    X = np.c_[df[columns]]
    X = X.astype(np.float64)
    y_fix = np.r_[df['Lin_F9']]
    y_predict_ = []

    for i in range(1,11):
        xgb_model = pickle.load(open("/home/cyang/paper_XGB/delta_LinF9_XGB/saved_model/mod_%d.pickle.dat"%i,"rb"))
        y_i_predict = xgb_model.predict(X, ntree_limit=xgb_model.best_ntree_limit)
        y_predict_.append(y_i_predict)
        
    y_predict = np.average(y_predict_, axis=0)
    df['XGB'] = pd.Series(y_predict+y_fix, index=df.index)
    df = df[['pdb','vina39','ligCover','betaScore','Lin_F9','XGB']]

    return df

def main():
    args = sys.argv[1:]
    if not args:
        print ('usage: python combine_all.py feature_dir')

        sys.exit(1)

    elif sys.argv[1] == '--help':
        print ('usage: python combine_all.py feature_dir')

        sys.exit(1)

    elif len(args) == 1:
        feature_dir = sys.argv[1]
        df = Combine_features(feature_dir)
        df = Calc_XGB(df)
        df.to_csv('%s/final_score.csv'%feature_dir,index=False)
        print (df)
        
    else:
        sys.exit(1)

if __name__ == '__main__':
    main()




