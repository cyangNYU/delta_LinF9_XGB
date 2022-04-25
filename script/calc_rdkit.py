import sys, os
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.ML.Descriptors.MoleculeDescriptors import MolecularDescriptorCalculator

LigandDescriptors = ['HeavyAtomMolWt', 'NumValenceElectrons','FpDensityMorgan1',
                     'FpDensityMorgan2', 'FpDensityMorgan3', 'LabuteASA',
                     'TPSA', 'NHOHCount', 'MolLogP','MolMR'] #

DescCalc = MolecularDescriptorCalculator(LigandDescriptors)


def GetRDKitDescriptors(mol):
    # Function for the calculation of ligand descriptors
    mol.UpdatePropertyCache(strict=False)
    Chem.GetSymmSSSR(mol)
    return DescCalc.CalcDescriptors(mol)

def main():
    args = sys.argv[1:]
    if not args:
        print ('usage: python cal_rdkit.py lig outfile fn')

        sys.exit(1)

    elif sys.argv[1] == '--help':
        print ('usage: python calc_rdkit.py lig outfile fn')

        sys.exit(1)

    elif len(args) == 3:
        outfile = sys.argv[2]
        fn = sys.argv[3]
        
        if sys.argv[1].endswith('.pdb'):
            mol = Chem.MolFromPDBFile(sys.argv[1], removeHs=False)
            List = GetRDKitDescriptors(mol)
            out = open(outfile, 'w')
            out.write(fn+','+','.join([str(round(i,5)) for i in List]) + '\n')
            out.close()
            
        elif sys.argv[1].endswith('.mol2'):
            mol = Chem.MolFromMol2File(sys.argv[1], removeHs=False)
            List = GetRDKitDescriptors(mol)
            out = open(outfile, 'w')
            out.write(fn+','+','.join([str(round(i,5)) for i in List]) + '\n')
            out.close()

        elif sys.argv[1].endswith('.sdf'):
            mol = Chem.MolFromMolFile(sys.argv[1], removeHs=False)
            List = GetRDKitDescriptors(mol)
            out = open(outfile, 'w')
            out.write(fn+','+','.join([str(round(i,5)) for i in List]) + '\n')
            out.close()
            
        else:
            sys.exit(1)
    else:
        sys.exit(1)

if __name__ == '__main__':
    main()
