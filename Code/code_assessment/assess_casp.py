# @Author:      Joana Pereira
# @Affiliation: MPI for Developmental Biology, Dept. Protein Evolution

import os
import matplotlib.cm
import matplotlib
import argparse
import math
import time
import sys
import requests
import urllib 

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import subprocess as sp

from matplotlib import rc

# define plots font to arial
rc('text', usetex=False)
rc('font', family='arial')

# Get inputs
parser = argparse.ArgumentParser('CASP_assessor: a python script for the high accurracy assession of CASP results after running analyse_casp.py')
requiredNamed = parser.add_argument_group('required arguments')
optionalNamed = parser.add_argument_group('optional arguments with defaults')

# required inputs
requiredNamed.add_argument('-casp', dest='which_casp', type=str, required=True, help='The casp version to analyse (e.g. CASP14, casp14 or 14')
# optional inputs
optionalNamed.add_argument('-tmp', dest='tmp_folder', type=str, default = '/tmp/', help='tmp folder (default: /tmp/')
optionalNamed.add_argument('-exclude_multidomain', dest='exclude_multidomain', type=str, default = 'False', help='Boolean statement to exclude multidomain targets. Will take only thos targets with "-D" (default: False')

# Define inputs
args = parser.parse_args()
which_casp = args.which_casp.upper()

if 'CASP' not in which_casp:
    which_casp = 'CASP{}'.format(which_casp)

tmp_folder = args.tmp_folder
tmp_folder = '{}{}'.format(tmp_folder, which_casp)

exclude_multidomain = args.exclude_multidomain
if exclude_multidomain == 'False':
    exclude_multidomain = False
else:
    exclude_multidomain = True

# HELPING ROUTINES

def parse_abstracts(names, which_casp = which_casp):

    print(' ... Parsing information from abstracts')

    # download the abstracts pdf
    abstracts_link = 'https://predictioncenter.org/{}/doc/{}_Abstracts.pdf'.format(which_casp.lower(), which_casp)
    abstracts = 'abstracts.pdf'
    urllib.request.urlretrieve(abstracts_link, abstracts)

    # convert pdf to txt
    txt_abstracts = '{}.txt'.format(abstracts.split('.')[0])
    convert_pdf = sp.Popen(['pdftotext', '-layout', abstracts, txt_abstracts], stdout=sp.PIPE, stderr=sp.PIPE, stdin=sp.PIPE)
    stdout, stderr = convert_pdf.communicate()

    # parse abstracts
    groups_data = {}
    found_name = False
    found_keywords  = False
    curr_name = None
    empty_key = None
    with open(txt_abstracts, 'r') as txt:
        for line in txt:
            line = line.strip()
            if len(line) > 0 and len(line.split()) == 1 and not line.startswith(' ') and line in names:
                curr_name = line
                if curr_name not in groups_data:
                    groups_data[curr_name] = {}
                found_name = True

            elif found_name:
                if 'Key:' in line and ';' in line:
                    line = line.replace('Key: ', '').replace(' ', '')
                    for key in line.split(';'):
                        key = key.split(':')
                        if len(key[-1]) > 0:
                            if 'MetaG' not in key[-1] and '.' in key[-1]:
                                groups_data[curr_name][key[0]] = key[-1].split('.')[0]
                            else:
                                groups_data[curr_name][key[0]] = key[-1]
                            empty_key = None
                        else:
                            empty_key = key[0]

                if '; MD:' in line:
                    line = line.replace(' ', '')
                    for key in line.split(';'):
                        key = key.split(':')

                        if empty_key is not None and len(empty_key)>0:
                            if 'MetaG' not in key[-1] and '.' in key[-1]:
                                groups_data[curr_name][empty_key] = key[-1].split('.')[0]
                            else:
                                groups_data[curr_name][empty_key] = key[-1]

                            empty_key = None
                        else:
                            if 'MetaG' not in key[-1] and '.' in key[-1]:
                                groups_data[curr_name][key[0]] = key[-1].split('.')[0]
                            else:
                                groups_data[curr_name][key[0]] = key[-1]
                    found_name = False


    return groups_data

def get_groups_names(which_casp = which_casp):

    print('Getting groups information')

    url = 'https://predictioncenter.org/{}/docs.cgi?view=groupsbyname'.format(which_casp.lower())

    html = requests.get(url).content
    df = pd.read_html(html)[-2]
    
    groups_names = {}
    for index, row in df.iterrows():
        gr = row['Group #']
        gr_name = row['Group Name']
        gr_type = row['Type']

        gr = f'{gr:03}'
        if 'server' in gr_type.lower():
            gr = '{}s'.format(gr)
            gr_type = 'Server'

        groups_names[gr] = {'name': gr_name, 'type': gr_type}

    groups_data = parse_abstracts(names=[groups_names[gr]['name'] for gr in groups_names])

    all_parameters = set([i for name in groups_data for i in groups_data[name].keys() if 'MetaG' not in i])
    
    for gr in groups_names:
        gr_name = groups_names[gr]['name']

        for parameter in all_parameters:
            if gr_name in groups_data and parameter in groups_data[gr_name]:
                groups_names[gr][parameter] = groups_data[gr_name][parameter]
            else:
                groups_names[gr][parameter] = np.nan

    return groups_names

def get_targets_individual_tables(data_folder):

    if os.path.isdir(data_folder):
        predictions_tables = sorted(['{}/{}'.format(data_folder, f) for f in os.listdir(data_folder) if '_predictions_updated_data' in f])
        loops_tables = sorted(['{}/{}'.format(data_folder, f) for f in os.listdir(data_folder) if '_templates_updated_data' in f])

        # check if there is any template table without a prediction table
        for table in loops_tables:
            curr_target = table.split('/')[-1].split('_')[0]
            if not os.path.isfile('{}/{}_predictions_updated_data.csv'.format(data_folder, curr_target)):
                print(' ---> Did not find predictions table for {}'.format(curr_target))
    else:
        predictions_tables = None
        loops_tables = None

    return predictions_tables, loops_tables

def get_targets_mr_data(data_folder):

    if os.path.isdir(data_folder):
        predictions_tables = sorted(['{}/{}'.format(data_folder, f) for f in os.listdir(data_folder) if '.csv' in f and 'templates' not in f])
        loops_tables = sorted(['{}/{}'.format(data_folder, f) for f in os.listdir(data_folder) if '-templates.csv' in f])

    else:
        predictions_tables = None
        loops_tables = None

    return predictions_tables, loops_tables

def compute_zscores(lst_values, zround = 0, z_scores = None, threshold = -2):

    # computes z-scores as described in Tristan's paper

    if z_scores != None:
        valid_lst_values = [lst_values[i] for i in range(len(lst_values)) if z_scores[i] > threshold]
    else:
        valid_lst_values = lst_values

    valid_lst_values = np.array(valid_lst_values)
    valid_lst_values = valid_lst_values[~np.isnan(valid_lst_values)]

    if len(valid_lst_values) > 0:

        avg = np.mean(valid_lst_values)
        sd  = np.std(valid_lst_values)

        z_scores = [(i-avg)/sd for i in lst_values]

        if zround == 1:
            z_scores = [0 if z < 0 or np.isnan(z) else z for z in z_scores]
        else:
            z_scores = compute_zscores(lst_values, zround = 1, z_scores = z_scores)

    else:
        z_scores = [np.nan for i in lst_values]

    return z_scores

def add_llgs(table, curr_target_llg, mode = None, new_columns = ['#Residues','LLG_Original','TFZeq_Original','Rms_Original','LLG_B_from_rms','TFZeq_B_from_rms','Rms_B_from_rms','LLG_const_B','TFZeq_const_B','Rms_const_B']):

    table = table.set_index('Model')

    if curr_target_llg is not None:

        lgg_table = pd.read_csv(curr_target_llg, sep=',')
        lgg_table['Model'] = ['' for i in range(len(lgg_table))]

        if mode is None:
            for i in range(len(lgg_table)):
                try:
                    val = lgg_table['LGA'][i].split('/')[1].replace('.sup.clean.pdb', '')
                except:
                    val = np.nan
                lgg_table['Model'][i] = val

        else:
            for i in range(len(lgg_table)):
                try:
                    val = lgg_table['LGA'][i].split('/')[1].split('_truncated')[0].replace('sculpt_','').replace('.pdb', '')
                    if 'null' in val:
                        val = 'null_model'
                except:
                    val = np.nan
                
                lgg_table['Model'][i] = val

        lgg_table = lgg_table.set_index('Model')
        table = pd.concat([table, lgg_table], axis=1, sort=False)  
        table = table.drop(columns=['LGA'])
        table.index.name = 'Model'

        return table
    else:
        for col in new_columns:
            table[col] = [np.nan for i in range(len(table))]

        return table

def compute_scores(predictions_tables, predictions_llg, groups_names):

    print('Computing models scores')

    major_table = pd.DataFrame()

    for index, table_file in enumerate(predictions_tables):

        print(index+1, len(predictions_tables), table_file)
        curr_target = table_file.split('/')[-1].split('_')[0]

        try:
            curr_target_llg = [i for i in predictions_llg if i.split('/')[-1].split('.')[0] == curr_target][0]
        except:
            curr_target_llg = None

        if not exclude_multidomain or '-D' in table_file:

            #try:
            table = pd.read_csv(table_file, sep='\t', index_col=0)
            table = add_llgs(table, curr_target_llg)
            
            if table_file.split('/')[-1].startswith('R'):
                table = table.drop(['Starting'])

            target_columns = [col for col in table.columns if 'Z' not in col and 'NP' not in col and 'LGA' not in col]

            z_table = pd.DataFrame({key: ['' for i in range(len(table))] for key in target_columns})
            z_table['Model'] = table.index
            z_table = z_table.set_index('Model')
            z_table['Target'] = [table_file.split('/')[-1].split('_')[0] for i in range(len(z_table))]

            for column in target_columns:
                z_table[column] = table[column]
                if not column in ['Model', 'GR#', 'Target_difficulty', 'Target_Neff', 'Target_classification', 'Target', 'Model_type', 'Model_template', 'nLoops', 'Target_title', 'Target_PDB_code', '#Residues', 'LGA']:
                    z_column = 'Z_{}'.format(column)
                    if 'BBscore' in column or 'SCscore' in column:
                        z_table[z_column] = compute_zscores([1-i for i in table[column]])
                    elif 'MolPrb' in column:
                        z_table[z_column] = compute_zscores([100-i for i in table[column]])
                    else:
                        z_table[z_column] = compute_zscores(table[column])

            if len(major_table) == 0:
                major_table = z_table.copy()
            else:
                major_table = pd.concat([major_table, z_table])

            # except:
            #     pass

    try:
        major_table['S_casp12'] = (1/3)*major_table['Z_GDT_HA'] + (1/9)*(major_table['Z_LDDT']+major_table['Z_CAD_AA']+major_table['Z_SphGr']) + (1/3)*major_table['Z_QSE']
    except:
        pass

    try:
        major_table['S_casp12-ASE'] = (1/2)*major_table['Z_GDT_HA'] + (1/6)*(major_table['Z_LDDT']+major_table['Z_CAD_AA']+major_table['Z_SphGr'])
    except:
        pass

    try:
        major_table['S_torsion_casp13'] = (2/3)*major_table['Z_Model_BBscore'] + (1/3)*major_table['Z_Model_SCscore']
    except:
        pass

    try:
        major_table['S_geom_casp13'] = (1/16)*(major_table['Z_LDDT']+major_table['Z_CAD_AA']+major_table['Z_SphGr']+major_table['Z_Model_SCscore']) + (1/8)*(major_table['Z_MolPrb_clash']+major_table['Z_Model_BBscore']) + (1/4)*(major_table['Z_GDT_HA']+major_table['Z_QSE'])
    except:
        pass

    try:
        major_table['S_geom_casp14'] = (1/16)*(major_table['Z_LDDT']+major_table['Z_CAD_AA']+major_table['Z_SphGr']+major_table['Z_Model_SCscore']) + (1/10)*(major_table['Z_MolPrb_clash']+major_table['Z_Model_BBscore']+major_table['Z_Model_DipDiff']) + (1/4)*(major_table['Z_GDT_HA']+major_table['Z_QSE'])
    except:
        pass

    all_parameters = set([i for name in groups_names for i in groups_names[name].keys() if 'MetaG' not in i])

    print(' ... Adding groups data')
    for i, parameter in enumerate(all_parameters):
        print(' ... ... {}.{} : {}'.format(i, len(all_parameters), parameter))
        major_table['GR_{}'.format(parameter)] = [groups_names[row['GR#']][parameter] if row['GR#'] in groups_names else np.nan for index, row in major_table.iterrows()]
    
    # save table to file
    major_table.to_csv('{}_Zscores_table.csv'.format(which_casp), sep='\t')
    
    return major_table

def combine_templates_data(templates_tables, templates_llg):

    print('Organizing template scores')

    llg_tables_analysed = []

    major_table = pd.DataFrame()

    for index, table_file in enumerate(templates_tables):

        print(index+1, len(templates_tables), table_file)
        curr_target = table_file.split('/')[-1].split('_')[0]

        try:
            curr_target_llg = [i for i in templates_llg if i.split('/')[-1].split('.')[0] == '{}-templates'.format(curr_target)][0]
            llg_tables_analysed.append(curr_target_llg)
        except:
            curr_target_llg = None

        if not exclude_multidomain or '-D' in table_file:

            table = pd.read_csv(table_file, sep='\t', index_col=0)
            table['Model'] = [i.split('/')[-1].split('_truncated')[0] for i in table['Truncated_pdb']]
            table = add_llgs(table, curr_target_llg, mode = 'template')

            target = table_file.split('/')[-1].split('_templates')[0]
            table['Target'] = [target for i in range(len(table))]

            if len(major_table) == 0:
                major_table = table.copy()
            else:
                major_table = pd.concat([major_table, table])

    # now add the LLG data of the null models and real structures for those cases for which there are no templates
    for index, table_file in enumerate(templates_llg):
        if table_file not in llg_tables_analysed:
            target = table_file.split('/')[-1].split('-templates')[0]

            table = pd.DataFrame({'Model': [target]})
            table = add_llgs(table, table_file, mode = 'template')
            table['Target'] = [target for i in range(len(table))]

            major_table = pd.concat([major_table, table])

    # save table to file
    major_table.to_csv('{}_templates_table.csv'.format(which_casp), sep='\t')

    return major_table

# MAIN CODE

data_folder = '{}/out_tables'.format(tmp_folder)
predictions_tables, templates_tables = get_targets_individual_tables(data_folder)

mr_folder = '../{}_MR'.format(which_casp)
predictions_llg, templates_llg = get_targets_mr_data(mr_folder)

groups_names = get_groups_names()

# 1. Compute major predictions table with z-scores for relevant CASP12 and CASP13 measures, 
#    including the scores used in those competitions and add LLGs

predictions_table = compute_scores(predictions_tables, predictions_llg, groups_names)

# # 2. Combine templates data into a unique table
if templates_tables is not None and len(templates_tables) > 0:
    templates_table = combine_templates_data(templates_tables, templates_llg)



