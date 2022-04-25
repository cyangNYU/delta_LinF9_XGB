import sys, os
import numpy as np
import pandas as pd
import alphaspace2 as al
import mdtraj
 
ADT = '/home/cyang/MGLTools-1.5.6/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_receptor4.py'

def Protein_pdbqt(PDB, PDBQT, ADT):
    '''
    generate the protein pdbqt file.
    '''
    cmd = ADT + " -r " + PDB + " -o " + PDBQT + " -U nphs"
    os.system(cmd)

def Strip_h(input_file,output_file):
    '''
    input_file and output_file need to be in pdb or pdbqt format 
    '''
    inputlines = open(input_file,'r').readlines()
    output = open(output_file,'w')
    for line in inputlines:
        if not 'H' in line[12:14]:
            output.write(line)
    output.close()

def Write_betaAtoms(ss, outfile):
    '''
    ss is the input AlphaSpace object, outfile is the output pdb file. 
    '''
    betaAtoms = open(outfile,'w')
    count = 1
    c = 1
    for p in ss.pockets:
        for betaAtom in p.betas:
            count = count+1
            coord = betaAtom.centroid
            ASpace = '%.1f'%betaAtom.space
            Score = '%.1f'%betaAtom.score
            atomtype = betaAtom.best_probe_type
            x, y, z  = '%.3f'%coord[-3], '%.3f'%coord[-2], '%.3f'%coord[-1]
            line = 'ATOM  ' + str(count).rjust(5) + str(atomtype).upper().rjust(5) + ' BAC' + str(c).rjust(5) + '     ' + str(x).rjust(8) + str(y).rjust(8) + str(z).rjust(8) + ' ' + str(ASpace).rjust(5) + ' ' + str(Score).rjust(5) + '           %s\n'%atomtype
            betaAtoms.write(line)
    betaAtoms.close()
    
def Prepare_beta(pdb, outfile, ADT=ADT):
    pdbqt = pdb[:-4]+'.pdbqt'
    Protein_pdbqt(pdb, pdbqt, ADT)

    pdb_noh = pdb[:-4]+'_noh.pdb'
    pdbqt_noh = pdb[:-4]+'_noh.pdbqt'

    Strip_h(pdb, pdb_noh)
    Strip_h(pdbqt, pdbqt_noh)

    prot = mdtraj.load(pdb_noh)
    al.annotateVinaAtomTypes(pdbqt=pdbqt_noh, receptor=prot)
    ss = al.Snapshot()
    ss.run(prot)
    Write_betaAtoms(ss, outfile)
    os.system('rm %s %s'%(pdb_noh, pdbqt_noh))
    return pdbqt

def main():
    args = sys.argv[1:]
    if not args:
        print ('usage: python prepare_betaAtoms.py pro.pdb outfile')

        sys.exit(1)

    elif sys.argv[1] == '--help':
        print ('usage: python prepare_betaAtoms.py pro.pdb outfile')

        sys.exit(1)

    elif len(args) == 2 and sys.argv[1].endswith('.pdb'):
        pdb = sys.argv[1]
        outfile = sys.argv[2]
        pdbqt = Prepare_beta(pdb, outfile)
        
    else:
        sys.exit(1)

if __name__ == '__main__':
    main()
