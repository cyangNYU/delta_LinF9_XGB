import sys, os
import numpy as np
from scipy.spatial.distance import cdist

def get_lig_coord(filename):
    '''
    Get all heavy atom coordinates of the ligand and return a coordinate matrix

    '''
    f1 = open(filename, 'r')
    
    List = []
    for line in f1.readlines():
        if line.startswith(('HETATM', 'ATOM')) and line[77]!='H':
            coord_x, coord_y, coord_z = line[27:38], line[38:46], line[46:54]
            result = [float(coord_x), float(coord_y), float(coord_z)]
            List.append(result)
    f1.close()
    List = np.asarray(List)
    return List

def get_beta_info(filename):
    '''
    Get beta atom coordinates and the corresponding beta atom scores

    '''
    f1 = open(filename, 'r')
    
    List = []
    Score = []
    for line in f1.readlines():
        if line.startswith(('HETATM', 'ATOM')) and line[77]!='H':
            coord_x, coord_y, coord_z = line[27:38], line[38:46], line[46:54]
            result = [float(coord_x), float(coord_y), float(coord_z)]
            List.append(result)
            score = float(line[61:66])
            Score.append(score)
    f1.close()
    List = np.asarray(List)
    Score = np.asarray(Score)
    return List, Score

def calc_betaScore_and_ligCover(lig_file, beta_file):
    '''
    Given one ligand file and one betaAtoms file, calculate the betaScore and ligand Coverage
    '''
    lig_coords = get_lig_coord(lig_file)
    beta_coords, beta_scores = get_beta_info(beta_file)
    result1 = cdist(beta_coords, lig_coords)
    score = np.sum(beta_scores[np.min(result1, axis=1)<=1.6])
    result2 = cdist(lig_coords, beta_coords)
    lig_cover_coords = lig_coords[np.min(result2, axis=1)<=1.6]
    return round(score,3),round(len(lig_cover_coords)/len(lig_coords),3)
    

def main():
    args = sys.argv[1:]
    if not args:
        print ('usage: python calc_ligCover_betaScore.py betaAtom lig outfile fn')

        sys.exit(1)
    elif len(args) == 4:
        beta = sys.argv[1]
        lig = sys.argv[2]
        outfile = sys.argv[3]
        fn = sys.argv[4]
        betaScore, ligCover = calc_betaScore_and_ligCover(lig, beta)
        out = open(outfile, "w")
        out.write(fn+","+ str(round(betaScore, 2)) + "," + str(round(ligCover, 2)) + "\n")
        out.close()


if __name__ == "__main__":
    main()



