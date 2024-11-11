import pandas as pd
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
%matplotlib inline
import math
import sys
from scipy import stats
import itertools
from itertools import combinations
from biopandas.pdb import PandasPdb
from joblib import Parallel, delayed

#PDB file to table
def pdb2table(file_path):
    col_names_pdb = ['type','n', 'atom', 'resName', 'chain', 'resN', 'x', 'y', 'z', 
                 'occup', 'tempfctr', 'element', 'empty']
    data = []
    with open(file_path, 'r') as f:
        for line in f:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                row_i = line.split()
                data.append(row_i)
    df = pd.DataFrame(data)
    if df.shape[1] == 12: 
        df.columns = col_names_pdb[:12]
        df = df.drop(columns = ['occup', 'tempfctr'])
    else:
        df.columns = col_names_pdb
        df = df.drop(columns = ['occup', 'tempfctr', 'empty'])
    return df

#PQR file to table
def pqr2table(file_path):
    dfpqr = pd.read_csv(file_path, sep='\s+', header=None)
    dfpqr.columns = ['type', 'n', 'atom', 'resName', 'resN', 'x', 'y', 'z', 'charge', 'colx']
    dfpqr['n'] = pd.to_numeric(dfpqr['n'], errors='coerce')
    dfpqr['resN'] = pd.to_numeric(dfpqr['resN'], errors='coerce')
    dfpqr = dfpqr[dfpqr['type'] == 'ATOM']
    dfpqr = dfpqr.drop(columns=['colx'])
    return dfpqr

# HEATM metals single or cluster 
def metal_dist(df, dist_max=5):
    df['x'] = pd.to_numeric(df['x'], errors='coerce')
    df['y'] = pd.to_numeric(df['y'], errors='coerce')
    df['z'] = pd.to_numeric(df['z'], errors='coerce')
    metpairs = pd.DataFrame(columns=['atoms', 'resName', 'resN', 'x', 'y', 'z'])
    clustered_metals = set()

    for i, j in combinations(range(len(df)), 2):
        metal1 = df.iloc[i]; metal2 = df.iloc[j]
        twometals = f"{metal1['atom']}{metal1['resN']}-{metal2['atom']}{metal2['resN']}"
        distance = np.linalg.norm(metal1[['x', 'y', 'z']] - metal2[['x', 'y', 'z']])

        if distance <= dist_max:
            if metal1['atom'] + str(metal1['resN']) not in clustered_metals and metal2['atom'] + str(metal2['resN']) not in clustered_metals:
                # Create a new cluster
                clustered_metals.add(metal1['atom'] + str(metal1['resN']))
                clustered_metals.add(metal2['atom'] + str(metal2['resN']))
                avg_coords = (metal1[['x', 'y', 'z']] + metal2[['x', 'y', 'z']]) / 2
                new_row = {'atoms': twometals, 'resName': metal1['resName'], 'resN': metal1['resN'],
                           'x': avg_coords['x'], 'y': avg_coords['y'], 'z': avg_coords['z']}
                metpairs = pd.concat([metpairs, pd.DataFrame([new_row])], ignore_index=True)
            else:
                # Add metals to existing cluster
                cluster_index = metpairs[(metpairs['atoms'].str.contains(metal1['atom'] + str(metal1['resN']))) | 
                                          (metpairs['atoms'].str.contains(metal2['atom'] + str(metal2['resN'])))]
                if cluster_index.shape[0] > 0:
                    cluster_index = cluster_index.index[0]
                    cluster_row = metpairs.loc[cluster_index]
                    new_coords = (cluster_row[['x', 'y', 'z']] + metal1[['x', 'y', 'z']] + metal2[['x', 'y', 'z']]) / 3
                    metpairs.loc[cluster_index, ['x', 'y', 'z']] = new_coords
                    metpairs.loc[cluster_index, 'atoms'] += f"-{metal2['atom']}{metal2['resN']}"
                    clustered_metals.add(metal2['atom'] + str(metal2['resN']))
    for index, row in df.iterrows():
        if row['atom'] + str(row['resN']) not in clustered_metals:
            row['atoms'] = row['atom']
            row = row[['atoms', 'resName', 'resN', 'x', 'y', 'z']]
            metpairs = pd.concat([metpairs, row.to_frame().T], ignore_index=True)
    return metpairs

#Define efield around metal sites
def efield_around_metal(pqrpath, pdbpath, cluster_maxdist=5, atomdist_shell=15):
    dfpdb = pdb2table(pdbpath)
    dfpdb = dfpdb[dfpdb['type'] == "HETATM"]
    dfmetals = metal_dist(dfpdb, cluster_maxdist)
    dfpqr = pqr2table(pqrpath)
    chain = dfpdb.iloc[0, 4]

    def compute_distance(atom1, atom2):
        pos1 = atom1[['x', 'y', 'z']].values
        pos2 = atom2[['x', 'y', 'z']].values
        distance = np.linalg.norm(pos1 - pos2)
        dir_vector = (pos1 - pos2) / distance
        ef = -abs(atom1['charge'] / distance) if atom1['charge'] < 0 else abs(atom1['charge'] / distance)
        ef_dir = abs(dir_vector) * ef
        return distance, dir_vector, ef, ef_dir

    def process_pair(i):
        atom1 = dfpqr.iloc[i[0]]
        atom2 = dfmetals.iloc[i[1]]
        distance, dir_vector, ef, ef_dir = compute_distance(atom1, atom2)
        if distance <= atomdist_shell:
            return {
                'n': atom1[1], 'atom': atom1[2], 'resName': atom1[3], 'resN': atom1[4], 'chain': chain, 'metal': atom2[0], 'metalN': atom2[2],
                'x': atom1[5], 'y': atom1[6], 'z': atom1[7], 'xm': atom2[3], 'ym': atom2[4], 'zm': atom2[5],
                'distance': distance, 'rx': dir_vector[0], 'ry': dir_vector[1], 'rz': dir_vector[2],
                'charge': atom1[8], 'E': ef, 'Ex': ef_dir[0], 'Ey': ef_dir[1], 'Ez': ef_dir[2]
            }
        return None
    i_atoms = range(len(dfpqr))
    i_metals = range(len(dfmetals))
    comb_atm_metal = list(itertools.product(i_atoms, i_metals))

    results = Parallel(n_jobs=-1)(delayed(process_pair)(i) for i in comb_atm_metal)
    atoms_near_metals = pd.DataFrame([res for res in results if res is not None])
    return atoms_near_metals

#efield for all data 
import os
def getall_efield(directory):
    all_data = []
    for root, _, files in os.walk(directory):
        for file in files:
            if file.endswith('.pdb'):
                pdbpath = os.path.join(root, file)
                pqrpath = pdbpath.replace("_Rlx.pdb", ".pqr")
                file_name = os.path.splitext(os.path.basename(file))[0].replace("_Rlx", "")
                try:
                    df = efield_around_metal(pqrpath, pdbpath, 5, 15) #obtain data for single PDB
                    df['structure'] = file_name
                    all_data.append(df)
                except:
                    print("Not enough data for " + file_name)    
    full_df = pd.concat(all_data, ignore_index=True)
    return full_df
