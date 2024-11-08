import pandas as pd
import math

# Splits each line in columns by position
def split_ATOM(raw):
    return(raw.str[6:11].str.replace(' ', ''),  
           raw.str[11:16].str.replace(' ', ''),
           raw.str[16].str.replace(' ', ''),
           raw.str[17:20].str.replace(' ', ''), 
           raw.str[21].str.replace(' ', ''),
           raw.str[22:26].str.replace(' ', ''),
           raw.str[27].str.replace(' ', ''),
           raw.str[30:37].str.replace(' ', ''),
           raw.str[38:45].str.replace(' ', ''),
           raw.str[46:53].str.replace(' ', ''),
           raw.str[54:59].str.replace(' ', ''),
           raw.str[60:65].str.replace(' ', ''),
           raw.str[72:75].str.replace(' ', ''),
           raw.str[76:78].str.replace(' ', ''),
           raw.str[79:].str.replace(' ', ''))

#Converts splitted lines into rows with column names
def get_ATOM_DF(ATOM, pdb=None):
    atom_data = split_ATOM(ATOM['raw'])
    ATOM = pd.DataFrame({
         'serial':atom_data[0],
         'atom_name':atom_data[1],
         'altLoc':atom_data[2],
         'resName':atom_data[3],
         'chainID':atom_data[4],
         'resNum':atom_data[5],
         'iCode':atom_data[6],
         'x':atom_data[7],
         'y':atom_data[8],
         'z':atom_data[9],
         'occupancy':atom_data[10],
         'tempFactor':atom_data[11],
         'segID':atom_data[12],
         'element':atom_data[13],
         'charge':atom_data[14] })
    ATOM['x']=ATOM['x'].astype(float)
    ATOM['y']=ATOM['y'].astype(float)
    ATOM['z']=ATOM['z'].astype(float)
    return(ATOM)

# Read PDB files and generates a table with ATOM and HETATM lines
def get_struc_atom_coords(fullpath2pdbfile):
    f = open(fullpath2pdbfile, "r")
    pdb_lines = f.readlines()
    pdb_file_df = pd.DataFrame(pdb_lines, columns = ['raw'])
    pdb_file_df['raw'] = pdb_file_df['raw'].str[:80]
    pdb_file_df['key'] = pdb_file_df['raw'].str[:6].str.replace(' ', '')
    atom_coords = get_ATOM_DF(pdb_file_df[pdb_file_df['key'] == 'ATOM'])
    if len(pdb_file_df[pdb_file_df['key'] == 'HETATM'])>0:
        hetatm_coords = get_ATOM_DF(pdb_file_df[pdb_file_df['key'] == 'HETATM'])
        atom_coords = pd.concat([atom_coords, hetatm_coords], ignore_index=True)
    return(atom_coords)

#read solvation energy files and 
def solvnrg(slvnrg_file):
     with open(slvnrg_file, 'r') as f:
          data = f.read()
     lines = data.split('\n')
     df = pd.DataFrame([line.split() for line in lines if line.startswith('SOLV')],
                    columns=['c1', 'c2', 'c3', 'atom', 'res', 'resNum', 'charge', 'solv_nrg'])
     df['charge'] = df['charge'].astype(float)
     return df[['atom', 'res', 'resNum', 'charge']]

def charge_res(nrg_file):
     nrg_res = nrg_file.groupby(['resNum', 'res'])[['charge']].sum().round(3)
     nrg_res.reset_index(inplace=True)
     return nrg_res

METAL_RES_CODES= [ "FE" ,"FE2","FES","FEO", "CU" ,"CU1","CUA","MG" ,"ZN" ,"MN","MO" ,"MOO","MOS","NI" ,"3CO","CO"]

# Get charge vectors and magnitudes for each metal sites in a single protein 
def get_charge_protein(charge_residues):
    charge_prot = charge_residues.groupby(['pdb', 'metal_num', 'metal_name'])[['charge','qx', 'qy', 'qz']].sum().round(5)
    charge_prot.reset_index(inplace=True)
    charge_prot["qt"] = 0.0
    df_charge = []
    for i, row in charge_prot.iterrows():
        row['qt'] = math.sqrt(row['qx']**2 + row['qy']**2 + row['qz']**2)
        if row['charge'] >= 0:
            row['qt'] = abs(row['qt'])
        else: 
            row['qt'] = -abs(row['qt'])
        irow = pd.DataFrame(row)
        df_charge.append(irow)
    df_charge = pd.concat(df_charge, axis=1).T
    return df_charge

# Generate a metal site or a cluster of metals


# Get vectors of charge for residues 
def get_charge_vectors(this_data_dir, this_struc_id, metal_maxdist):
    this_struc_file = "%s/%s_Rlx.pdb"%(this_data_dir, this_struc_id)
    slvnrg_file = "%s/%s_bluues.solv_nrg"%(this_data_dir, this_struc_id)

    # data from pdb and from solvation_energy
    df_pdb = get_struc_atom_coords(this_struc_file)
    df_slvnrg = solvnrg(slvnrg_file)
    df_slvnrg["charge"] = df_slvnrg["charge"].astype(float)
    df_slvnrg = charge_res(df_slvnrg)

    #merge solvation energy values
    atoms = df_pdb.loc[(~df_pdb["atom_name"].isin(METAL_RES_CODES)) & (df_pdb["atom_name"]=='C')]
    atoms = atoms.merge(df_slvnrg.set_index('resNum'), how='left', on='resNum', suffixes=('', '_y'))
    atoms = atoms.drop('charge', axis=1)
    atoms = atoms.rename(columns={'charge_y': 'charge'})
    
    #Select metal ions 
    metals = df_pdb.loc[df_pdb["atom_name"].isin(METAL_RES_CODES)]
    for mi, mr in mah2_sites.iterrows():
        metals_row = mr[['pdb_name', 'seqNum1', 'seqNum2', 'seqNum3', 'seqNum4']]
        metals_list = metals_row[['seqNum1', 'seqNum2', 'seqNum3', 'seqNum4']].values.tolist()
        mets_row = metals.loc[(metals["resNum"].isin(metals_list))&(metals["resName"]==metals_row["pdb_name"])]
        #print(metals_row, mets_row.head(1))


    # Create table of data for atoms and residues
    atoms['pdb'] = this_struc_id
    atoms['metal_name'] = ""
    atoms['metal_num'] = 0
    atoms['metal_dist'] = 10000
    atoms['ax'] = 0; atoms['ay'] = 0; atoms['az'] = 0
    atoms['ux'] = 0; atoms['uy'] = 0; atoms['uz'] = 0
    atoms['qx'] = 0; atoms['qy'] = 0; atoms['qz'] = 0
    df_atoms = []
    if len(metals)>0:
        for mindex, mrow in metals.iterrows():
            for index, row in atoms.iterrows():
                row['metal_dist'] = math.sqrt((mrow['x']-row['x'])**2 + (mrow['y']-row['y'])**2 + ((mrow['z']-row['z'])**2))
                row['ax'] = row['x']-mrow['x']; row['ux'] = row['ax']/row['metal_dist']
                row['ay'] = row['y']-mrow['y']; row['uy'] = row['ay']/row['metal_dist']
                row['az'] = row['z']-mrow['z']; row['uz'] = row['az']/row['metal_dist']
                row['metal_name']= mrow['resName']
                row['metal_num']= mrow['resNum']
                row['qx'] = row['charge']*abs(row['ux']/(row['metal_dist']))
                row['qy'] = row['charge']*abs(row['uy']/(row['metal_dist']))
                row['qz'] = row['charge']*abs(row['uz']/(row['metal_dist']))
                irow = pd.DataFrame(row)
                df_atoms.append(irow)
    df_atoms = pd.concat(df_atoms, axis=1).T
    df_atoms = df_atoms.loc[df_atoms['metal_dist'] <= metal_maxdist].copy().reset_index(drop=True)
    df_charge = get_charge_protein(df_atoms)
    return (df_charge, df_atoms)