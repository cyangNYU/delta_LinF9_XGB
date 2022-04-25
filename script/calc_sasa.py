#-----------------------------------------------------------------------------
# SASA Feature
#-----------------------------------------------------------------------------
import os
import sys
import pandas as pd
import fileinput
from featureSASA import sasa


def cal_SASA(out,fn,lig,pro,datadir):
    # fn: pdbid
    # lig: ligand part
    # pro: protein part
    pro = os.path.join(datadir,pro)
    lig = os.path.join(datadir,lig)
    sasa_features = sasa(datadir,pro,lig)
    sasa_com = sasa_features.sasa
    sasa_pro = sasa_features.sasa_pro
    sasa_lig = sasa_features.sasa_lig
   
    out.write(fn + "," +  ",".join([str(round(i,2) )for i in sasa_com]) + "," + ",".join([str(round(i,2) )for i in sasa_lig]) + "," + ",".join([str(round(i,2) )for i in sasa_pro]))
    

def main():
    args = sys.argv[1:]
    if not args:
        print ('usage: python calc_sasa.py datadir lig_file pro_file outfile fn')
        
        sys.exit(1)
        
    elif len(args) == 5:
        datadir = sys.argv[1]
        lig = sys.argv[2]
        pro = sys.argv[3]
        outfile = sys.argv[4]
        fn = sys.argv[5]
        out = open(outfile, "w")
        cal_SASA(out,fn,lig,pro,datadir)
        out.close()
        #os.system('rm %s/%s'%(datadir, lig))
        #os.system('rm %s/%s'%(datadir, pro))

if __name__ == "__main__":
    main()
    
    
    


