'''
Copyright (c) 2021 Kristoffer Torbjørn Bæk, Kasper Planeta Kepp, Technical University of Denmark
'''

# Import modules
import os
import pandas as pd
import freesasa as fs
from Bio.PDB import PDBParser
from Bio import SeqIO
import pkg_resources
import json
from natsort import natsort_keygen

# Path to resource files
naccess_config = pkg_resources.resource_filename(__name__, 'naccess.config')
with open(pkg_resources.resource_filename(__name__, 'HVdiff_table.json')) as data_file:
        HVdiff_table = json.load(data_file)

# Define functions
def calc_RSA(filepath):
    '''Runs freesasa on PDB.'''
    classifier = fs.Classifier(naccess_config)
    structure = fs.Structure(filepath, classifier)
    result = fs.calc(structure)
    area_list = get_residueAreas(result)

    return pd.DataFrame(area_list, columns=('Chain', 'Number', 'Wild', 'RSA'))

def get_residueAreas(result):
    '''Extracts output values from freesasa and calculates RSA.'''
    standards = {'ALA': 107.95,
             'CYS': 134.28,
             'ASP': 140.39,
             'GLU': 172.25,
             'PHE': 199.48,
             'GLY': 80.1,
             'HIS': 182.88,
             'ILE': 175.12,
             'LYS': 200.81,
             'LEU': 178.63,
             'MET': 194.15,
             'ASN': 143.94,
             'PRO': 136.13,
             'GLN': 178.5,
             'ARG': 238.76,
             'SER': 116.5,
             'THR': 139.27,
             'VAL': 151.44,
             'TRP': 249.36,
             'TYR': 212.76}

    aa_codes = {'ALA': 'A',
             'CYS': 'C',
             'ASP': 'D',
             'GLU': 'E',
             'PHE': 'F',
             'GLY': 'G',
             'HIS': 'H',
             'ILE': 'I',
             'LYS': 'K',
             'LEU': 'L',
             'MET': 'M',
             'ASN': 'N',
             'PRO': 'P',
             'GLN': 'Q',
             'ARG': 'R',
             'SER': 'S',
             'THR': 'T',
             'VAL': 'V',
             'TRP': 'W',
             'TYR': 'Y'}

    l = []
    r = result.residueAreas()
    for chain, chainvalue in r.items():
        for residue, value in chainvalue.items():
            if len(value.residueType) == 3: # disregard nucleotides
                residue_code = aa_codes[value.residueType]
                total_rel = value.total / standards[value.residueType]
                l.append((chain,
                          value.residueNumber,
                          residue_code,
                          total_rel
                         ))

    return l

def calc_simba_IB(RSA, Vdiff, Hdiff):
    '''Calculates predicted ddG based on Simba-IB.'''
    ddG = (-0.692
    + 0.905 * RSA
    + 1.425 * Vdiff
    - 0.365 * Hdiff
    - 1.494 * RSA * Vdiff
    + 0.625 * RSA * Hdiff)

    return ddG

def calc_simba_SYM(RSA, Vdiff, Hdiff):
    '''Calculates predicted ddG based on Simba-SYM.'''
    ddG = (
    + 1.642 * Vdiff
    - 0.421 * Hdiff
    - 1.867 * RSA * Vdiff
    + 0.737 * RSA * Hdiff)

    return ddG

def check_chains(filepath):
    '''Checks if PDB is a multi_chain structure and homooligomer.'''
    chain_seq = []
    for chain in SeqIO.parse(filepath, "pdb-seqres"):
        chain_seq.append(chain.seq)

    chain_seq = [chain for chain in chain_seq if chain != '']

    multi_chain = len(chain_seq) > 1
    homo = len(set(chain_seq)) == 1

    return multi_chain, homo

def mean_RSA(df):
    '''Calculates mean RSA across chains.'''
    chain_mean = df.groupby('Number')['RSA'].mean().reset_index()
    chain_mean.columns = ['Number', 'RSA_mean']

    return df.join(chain_mean.set_index('Number'), on = 'Number')

def join_hvdiff(df):
    '''Imports Hdiff and Vdiff values and joins with dataframe.'''
    hvdiff = pd.DataFrame.from_dict(HVdiff_table)

    # Divide Vdiff with 100
    hvdiff['Vdiff'] = hvdiff.apply(lambda row : row['Vdiff'] / 100, axis=1)

    return df.join(hvdiff.set_index('Wild'), on='Wild')

def exists_pdb(name, pdb_dir):
    '''Checks if PDB was downloaded and returns path.'''
    pdb_filename = 'pdb' + name.lower() + '.ent'
    pdb_path = os.path.join(pdb_dir, pdb_filename)

    return pdb_path, os.path.exists(pdb_path)

def simba2_predict(name, pdb_path):
    '''Runs all calculations'''
    # Calculate RSA
    RSA_df = calc_RSA(pdb_path)

    # Test if multichain and homooligomer
    multi_chain, homo = check_chains(pdb_path)

    # If structure is a homooligomer, calculate mean RSA
    if multi_chain and homo:
        RSA_df = mean_RSA(RSA_df)

    # Join Hdiff and Vdiff
    df = join_hvdiff(RSA_df)

    # Calculate predicted ddG for each residue and each chain
    df['ddG_SimBa_IB'] = df.apply(lambda row : calc_simba_IB(row['RSA'], row['Vdiff'], row['Hdiff']), axis=1)
    df['ddG_SimBa_SYM'] = df.apply(lambda row : calc_simba_SYM(row['RSA'], row['Vdiff'], row['Hdiff']), axis=1)

    # If structure is a homooligomer, calculate predicted ddG for mean
    if multi_chain and homo:
        df['ddG_SimBa_IB_mean'] = df.apply(lambda row : calc_simba_IB(row['RSA_mean'], row['Vdiff'], row['Hdiff']), axis=1)
        df['ddG_SimBa_SYM_mean'] = df.apply(lambda row : calc_simba_SYM(row['RSA_mean'], row['Vdiff'], row['Hdiff']), axis=1)

    # Add column with PDB code
    df.insert(loc=0, column='PDB', value=name.upper())

    return df.sort_values(by = ['Chain', 'Number', 'Mutated'], key=natsort_keygen())
