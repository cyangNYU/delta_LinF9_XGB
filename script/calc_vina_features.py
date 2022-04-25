"""
Vina Features generation
"""

__author__ = "Chao Yang"
__copyright__ = "Copyright 2020, NYU"
__license__ = ""

#-----------------------------------------------------------------------------
# Imports
#-----------------------------------------------------------------------------

import os, sys
import numpy as np
import re

#-----------------------------------------------------------------------------
# Code
#-----------------------------------------------------------------------------

Vina = '/home/cyang/paper_XGB/delta_LinF9_XGB/software/smina_feature'
Smina = '/home/cyang/paper_XGB/delta_LinF9_XGB/software/smina.static'
SF = '/home/cyang/paper_XGB/delta_LinF9_XGB/software/sf_vina.txt'

def runVina(protpdbqt, ligpdbqt, Vina):
    """Run modified Smina program with Lin_F9 score and 48 features
    
    Parameters
    ----------
    protpdbqt : str
        PDBQT file name of protein
    ligpdbqt : str
        PDBQT file name of ligand
        
    Returns
    ----------
    vinalist : list[float]
        48 features by Smina
        
    """
    cmd = Vina + " -r" + protpdbqt + " -l " + ligpdbqt + \
          " --score_only --custom_scoring " + SF
    process = os.popen(cmd)
    
    List = process.read().strip('\n').split()[1:]
    vinalist = [float(i) for i in List]

    return vinalist

def Get_LinF9_from_pose(ligpdbqt):
    '''
    Get Lin_F9 score from docked pose

    '''
    f1 = open(ligpdbqt, 'r')
    score = 0
    
    for line in f1.readlines():
        if line.startswith('REMARK minimizedAffinity'):
            score = line.strip('\n').split()[2]
            score = float(score)
            break
        
    return score

def calc_LinF9(pro_file, lig_file, Smina):
    '''
    Get Lin_F9 score from Smina score_only
    '''
    cmd = Smina + " -r" + pro_file + " -l " + lig_file + \
          " --score_only --scoring Lin_F9"
    process = os.popen(cmd).read()

    match = re.search(r'Affinity:\s(\S+) ', process)

    if match:
        return float(match.group(1))
    else:
        return 0


class vina:
    """Vina score and vina features
    
    """
    
    def __init__(self, prot, lig, Vina, Smina):
        """Vina Socre and Vina Features
        
        Parameters
        ----------
        prot : str
            protein structure
        lig : str
            ligand structure
        
        """
        self.prot = prot
        self.lig = lig
        self.Vina = Vina
        self.Smina = Smina
        
        vinalist = runVina(self.prot, self.lig, self.Vina)
        score = Get_LinF9_from_pose(self.lig)
        if score == 0:
            score = calc_LinF9(self.prot, self.lig, self.Smina)
        
        self.LinF9 = score
        self.vinaFeatures = vinalist


    def features(self, num = 100):
        """Get subset of features
        
        Parameters
        ----------
        num : int (default 10)
            number of features to retrieve
        
        """
        idx6 = [10, 28, 37, 44, 49, 52]
        idx10 = [0, 2, 52, 54, 53, 55, 3, 51, 57, 47]
        
        if num == 10:
            return [self.vinaFeatures[i] for i in idx10]
        elif num == 6:
            return [self.vinaFeatures[i] for i in idx6]
        else:
            return self.vinaFeatures

def main():
    args = sys.argv[1:]
    if not args:
        print ('usage: python calc_vina_features.py pro lig outfile fn')

        sys.exit(1)
    elif len(args) == 4:
        pro = sys.argv[1]
        lig = sys.argv[2]
        outfile = sys.argv[3]
        fn = sys.argv[4]
        #fn = sys.argv[4].split('_')
        v = vina(pro, lig, Vina, Smina)
        out = open(outfile, "w")
        out.write(fn+","+str(round(v.LinF9, 5))+","+ ",".join([str(round(i,5)) for i in v.features(60)]) + "\n")
        out.close()


if __name__ == "__main__":
    main()


