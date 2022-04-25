#! /bin/bash

datadir='../test/JG98'
lig='../test/JG98/JG98.pdb'
pro='../test/JG98/protein_ATP.pdb'

b=`basename $pro .pdb`
pro_pdbqt=$datadir/$b.pdbqt
betaAtoms=$datadir/betaAtoms.pdb
abs_dir=`realpath $datadir`
lig_=`basename $lig`
pro_=`basename $pro`


python prepare_betaAtoms.py $pro $betaAtoms

python calc_vina_features.py $pro_pdbqt $lig $datadir/vina_features.csv $datadir

python calc_ligCover_betaScore.py $betaAtoms $lig $datadir/beta_features.csv $datadir

python calc_sasa.py $abs_dir  $lig_ $pro_ $datadir/sasa_features.csv $datadir

python calc_rdkit.py $lig $datadir/lig_features.csv $datadir

python calc_bridge_wat.py $pro $lig $datadir/BW_features.csv $datadir

python combine_all.py $datadir

