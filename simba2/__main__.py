#!/usr/bin/python
# coding: utf-8

'''
Copyright (c) 2021 Kristoffer Torbjørn Bæk, Kasper Planeta Kepp, Technical University of Denmark
'''

import click
import os
from Bio.PDB import PDBList
from .methods import exists_pdb, simba2_predict, check_chains

@click.command()
@click.argument('pdbname', type=str)
@click.option('--dir', '-d', default=os.getcwd(), help='Alternative directory to save downloaded pdb file in.')
@click.option('--file', '-f', type=click.Path(exists=True), help='Path to existing pdb file (instead of downloading from RCSB).')
@click.option('--keep/--no-keep', '-k/ ', default=False, help='Keep downloaded PDB file (default: --no-keep).')
def main(pdbname, dir, file, keep):
    '''Predicts ddG with SimBa-IB and SimBa-SYM'''
    name = pdbname
    # Download PDB and check if PDB is downloaded
    if file is None:
        _, already_exists = exists_pdb(name, dir)

        pdbl = PDBList()
        pdbl.retrieve_pdb_file(name, pdir=dir, file_format='pdb')

        pdb_path, exists = exists_pdb(name, dir)

        if not exists:
            print('PDB', name.upper(), 'not downloaded')
            exit()

    # ... or use path to existing PDB
    else:
        pdb_path = file

    # Run calculations
    result = simba2_predict(name, pdb_path)
    multi_chain, homo = check_chains(pdb_path)

    # Keep or remove PDB
    if file is None and not already_exists:
        if keep:
            print('PDB', name.upper(), 'saved at', pdb_path)
        else:
            os.remove(pdb_path)

    # Print summary statements and export to csv
    print('The structure contains', str(len(result.Chain.unique())),
    'amino acid chain(s).')

    result['ddG_SimBa_IB'] = result['ddG_SimBa_IB'].round(decimals=1)
    result['ddG_SimBa_SYM'] = result['ddG_SimBa_SYM'].round(decimals=1)
    result['RSA'] = result['RSA'].round(decimals=3)
    result['Vdiff'] = result['Vdiff'].round(decimals=2)

    if multi_chain and homo:
        print('The structure is a homooligomer.')
        result['ddG_SimBa_IB'] = result['ddG_SimBa_IB_mean'].round(decimals=1)
        result['ddG_SimBa_SYM'] = result['ddG_SimBa_SYM_mean'].round(decimals=1)
        result['RSA_mean'] = result['RSA_mean'].round(decimals=3)
        result = result.drop(columns = ['ddG_SimBa_IB_mean', 'ddG_SimBa_SYM_mean'])

    elif multi_chain:
        print('The structure is a heterooligomer.')

    outputfile1 = 'SimBa-IB_' + name.upper()
    outputfile2 = 'SimBa-SYM_' + name.upper()
    result.drop(columns = 'ddG_SimBa_SYM').to_csv(os.path.join(os.getcwd(), outputfile1 + '.csv'), index=False)
    result.drop(columns = 'ddG_SimBa_IB').to_csv(os.path.join(os.getcwd(), outputfile2 + '.csv'), index=False)

if __name__ == "__main__": # if script is being invoked from command line
    main()
