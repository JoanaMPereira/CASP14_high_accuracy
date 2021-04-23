# @Author:      Joana Pereira
# @Affiliation: MPI for Developmental Biology, Dept. Protein Evolution

import os
import urllib 
import shutil
import matplotlib.cm
import matplotlib
import argparse
import math
import time
import sys
import json
import re
import difflib

import subprocess as sp
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import multiprocessing as mp

from matplotlib import rc
from Bio.SeqUtils import seq3
from Bio.PDB import *
from Bio.SubsMat import MatrixInfo as matlist
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from scipy.stats.stats import pearsonr   

from pdbe import pyPDBeREST
p = pyPDBeREST()

# define plots font to arial
# rc('text', usetex=False)
# rc('font', family='arial')

# Define dependecies

muscle = '/ebio/abt1_share/toolkit_sync/bioprogs/tools/muscle/muscle'
dipcheck = '/ebio/xray/software/arp_warp_8.0/bin/bin-x86_64-Linux/dipcheck'
predictions_gdt = '/ebio/abt1_share/small_projects/jpereira/CASP_assessement_2020/code/predictions_gdt.sh'
templates_lga = '/ebio/abt1_share/small_projects/jpereira/CASP_assessement_2020/code/templates_lga.sh'
templates_gdt = '/ebio/abt1_share/small_projects/jpereira/CASP_assessement_2020/code/templates_gdt.sh'
dssp = '/ebio/xray/software/ccp4/bin/mkdssp'
molprobity_sidechains = '/ebio/xray/software/ccp4/bin/molprobity.rotalyze'
molprobity_ramangles = '/ebio/xray/software/ccp4/bin/molprobity.ramalyze'
molprobity_omegangles = '/ebio/xray/software/ccp4/bin/molprobity.omegalyze'
molprobity_cablam = '/ebio/xray/software/ccp4/bin/molprobity.cablam'
pdbtools = '/ebio/abt1_share/software/pdb-tools/pdbtools/'
phenix_pdbtools = '/ebio/xray/software/phenix-1.13-2998/build/bin/phenix.pdbtools'
sculptor = '/ebio/xray/software/ccp4/bin/phaser.sculptor'

# Get inputs
parser = argparse.ArgumentParser('CASP_analyser: a python script for the high accurracy analysis of CASP results. You have to set the memory limit to unlimited due to LGA.\nEXAMPLE USAGE: ulimit -s unlimited; python3 analyse_casp_v2.py -casp 14 -tmp /tmp/')
requiredNamed = parser.add_argument_group('required arguments')
optionalNamed = parser.add_argument_group('optional arguments with defaults')

# required inputs
requiredNamed.add_argument('-casp', dest='which_casp', type=str, required=True, help='The casp version to analyse (e.g. CASP14, casp14 or 14')
# optional inputs
optionalNamed.add_argument('-tmp', dest='tmp_folder', type=str, default = '/tmp/', help='tmp folder (default: /tmp/')
optionalNamed.add_argument('-n_threads', dest='n_threads', type=int, default = 1, help='Number of threads or the number of cpus to use (default: 1')
optionalNamed.add_argument('-update', dest='update', type=str, default = 'True', help='Boolean to check if there are new models (default: True')

# Define inputs
args = parser.parse_args()
which_casp = args.which_casp.upper()
n_threads = args.n_threads
update = args.update

if update == 'True':
    update = True
else:
    update = False

if 'CASP' not in which_casp:
    which_casp = 'CASP{}'.format(which_casp)

if which_casp == 'CASP14' or which_casp == 'CASP13':
    top_level_url = "https://predictioncenter.org/{}/assessors/".format(which_casp.lower())

    predictions_domains    = '{}TARBALLS/predictions_trimmed_to_domains'.format(top_level_url)
    predictions_fulleng    = '{}TARBALLS/predictions'.format(top_level_url)
    predictions_tables     = '{}TARBALLS/results/tables'.format(top_level_url)
    templates              = '{}TARBALLS/templates'.format(top_level_url)

    targets_trimmed   = '{}TARBALLS/targets/Targets_trimmed_to_domains.tar.gz'.format(top_level_url)
    targets_processed = '{}TARBALLS/targets/Targets_processed.tgz'.format(top_level_url)
    targets_submitted = '{}TARBALLS/targets/Targets_submitted.tgz'.format(top_level_url)
    all_targets       = '{}TARGETS/'.format(top_level_url)

    if which_casp == 'CASP14':
        sequence_deep = '{}/PSIBLAST-HHsearch_domains/PSI-HH_Templ-Neff_dom.csv'.format(templates)
    else:
        sequence_deep = '{}/PSI-HH__Templ-Neff.csv'.format(templates)

elif which_casp == 'CASP12':

    top_level_url = "https://predictioncenter.org/download_area/{}/".format(which_casp)

    predictions_domains    = '{}predictions_trimmed_to_domains'.format(top_level_url)
    predictions_fulleng    = '{}predictions'.format(top_level_url)
    predictions_tables     = '{}results/tables'.format(top_level_url)
    templates              = '{}templates'.format(top_level_url)

    targets_trimmed   = '{}targets/{}.domains_T0.releaseDec022016.tgz'.format(top_level_url, which_casp.lower())
    targets_processed = '{}targets/{}.targets_TR.releaseDec022016.tgz'.format(top_level_url, which_casp.lower())
    targets_submitted = None
    sequence_deep = None
    all_targets = None
    
    predictions_tables  = '{}SUMMARY_TABLES'.format(top_level_url)
else:
    print('NOT POSSIBLE TO ANALYSE {}'.format(which_casp))
    exit()

targets_summary = 'https://predictioncenter.org/casp14/targetlist.cgi?type=csv'.format(which_casp.lower())

tmp_folder = args.tmp_folder
tmp_folder = '{}{}'.format(tmp_folder, which_casp)
if not os.path.isdir(tmp_folder):
    os.mkdir(tmp_folder)

domains_classifications = '/ebio/abt1_share/small_projects/jpereira/CASP_assessement_2020/casp_domain_classifications/{}.csv'.format(which_casp)
mr_folder = '/ebio/abt1_share/small_projects/jpereira/CASP_assessement_2020/{}_MR'.format(which_casp)
if not os.path.isdir(mr_folder):
    os.mkdir(mr_folder)

# set parameters to choose templates (minimum coverage of the template and the maximum number of accepted templates)
min_cov = 50
top_n_templates = 10

# Define automatic logger
class Logger(object):
    def __init__(self, name):
        self.terminal = sys.stdout
        self.log = open(name, "a")

    def write(self, message):
        self.terminal.write(message)
        self.log.write(message)  

    def flush(self):
        #this flush method is needed for python 3 compatibility.
        #this handles the flush command by doing nothing.
        #you might want to specify some extra behavior here.
        pass 

logfile = '{}_analysis.log'.format(which_casp)
sys.stdout = Logger(logfile)

# Login into CASP ftp
password_mgr = urllib.request.HTTPPasswordMgrWithDefaultRealm()                   # create a password manager
password_mgr.add_password(None, top_level_url, 'USERNAME', 'PASSWORD')       # Add the username and password.
handler = urllib.request.HTTPBasicAuthHandler(password_mgr)
opener = urllib.request.build_opener(handler)

# DEFINE FUNCTIONS

# 0. General tools

def chunk_list(l, n, unique_targets):
    chunks = np.array_split(np.array(unique_targets), n)
    chunks = [[i for i in l for target in chunk if target in i] for chunk in chunks]
    return chunks

def add_chainid(pdb, chainid = 'A'):

    curr_temp = '{}_tmp'.format(pdb)

    with open(curr_temp, 'w') as tmp:
        with open(pdb, 'r') as inp:
            for line in inp:
                if 'ATOM' in line:
                    line = list(line)
                    line[21] = 'A'
                    line = ''.join(line)
                tmp.write(line)

    os.replace(curr_temp, pdb)

    return pdb

def fix_backbone_atom_order(pdb):

    atoms = {}

    with open(pdb, 'r') as inp:
        for line in inp:
            if 'ATOM' in line:
                resid = line[17:21].strip()
                resnum = line[22:27].strip()
                chain = line[21]
                atom = line[13:17].strip()

                if resnum not in atoms:
                    atoms[resnum] = {}
                atoms[resnum][atom] = line

    curr_temp = '{}_tmp'.format(pdb)
    with open(curr_temp, 'w') as tmp:
        for resnum in atoms:
            for atom in ['N', 'CA', 'C', 'O']:
                try:
                    tmp.write(atoms[resnum][atom])
                except:
                    print('   --> ATTENTION: atom {} not found for residue number {} in {}'.format(atom, resnum, pdb))
            for atom in atoms[resnum]:
                if atom not in ['N', 'CA', 'C', 'O']:
                    tmp.write(atoms[resnum][atom])

    os.replace(curr_temp, pdb)

    return pdb

def renumber_atoms(pdb):
    
    atom_numbers = []

    with open(pdb, 'r') as inp:
        for line in inp:
            if 'ATOM' in line and len(atom_numbers) < 30:
                atomnum = int(line[5:12].strip())
                atom_numbers.append(atomnum)

    run_pdbtools = sp.Popen(['python', '{}pdb_reatom.py'.format(pdbtools), '-{}'.format(min(atom_numbers)), pdb], stdout=sp.PIPE, stderr=sp.PIPE, stdin=sp.PIPE)
    stdout, stderr = run_pdbtools.communicate()

    with open(pdb, 'w') as outp: 
        outp.write(stdout.decode("utf-8"))

    return pdb

def set_bfactor(pdb, b = 10, out_folder = None):

    run_pdbtools = sp.Popen(['python', '{}pdb_b.py'.format(pdbtools), '-{}'.format(b), pdb], stdout=sp.PIPE, stderr=sp.PIPE, stdin=sp.PIPE)
    stdout, stderr = run_pdbtools.communicate()

    if '\nATOM' not in stdout.decode("utf-8"):
        stdout = stdout.decode("utf-8").replace('ATOM', '\nATOM')
    else:
        stdout = stdout.decode("utf-8")

    if out_folder == None:
        with open(pdb, 'w') as outp: 
            outp.write(stdout)
    else:
        if not os.path.isdir(out_folder):
            os.mkdir(out_folder)
        pdb = '{}/{}'.format(out_folder, pdb.split('/')[-1])
        with open(pdb, 'w') as outp: 
            outp.write(stdout)

    return pdb

def set_occ(pdb, occ = 1.0, out_folder = None):

    run_pdbtools = sp.Popen(['python', '{}pdb_occ.py'.format(pdbtools), '-{}'.format(occ), pdb], stdout=sp.PIPE, stderr=sp.PIPE, stdin=sp.PIPE)
    stdout, stderr = run_pdbtools.communicate()

    if '\nATOM' not in stdout.decode("utf-8"):
        stdout = stdout.decode("utf-8").replace('ATOM', '\nATOM')
    else:
        stdout = stdout.decode("utf-8")

    if out_folder == None:
        with open(pdb, 'w') as outp: 
            outp.write(stdout)
    else:
        if not os.path.isdir(out_folder):
            os.mkdir(out_folder)
        pdb = '{}/{}'.format(out_folder, pdb.split('/')[-1])
        with open(pdb, 'w') as outp: 
            outp.write(stdout)

    return pdb

def fix_pdbs(pdbs):

    for pdb in pdbs:
        #pdb = add_chainid(pdb)
        pdb = fix_backbone_atom_order(pdb)
        pdb = renumber_atoms(pdb)

    return pdbs

def download_pdb(pdb, chain, out_folder = tmp_folder):

    if not os.path.isdir(out_folder):
        os.mkdir(out_folder)

    out_pdb = '{}/{}_{}.pdb'.format(out_folder, pdb, chain)
    pdb_link = 'https://files.rcsb.org/download/{}.pdb'.format(pdb.upper())

    if not os.path.isfile(out_pdb):
        try:
            urllib.request.urlretrieve(pdb_link, out_pdb)
        except:
            out_pdb = None

    if out_pdb is not None and os.path.isfile(out_pdb) and chain is not None:
        out_pdb = extract_chain_from_pdb(out_pdb, chain)

    else:
        print('    ### PDB download error: Not possible to download {}_{}.pdb'.format(pdb, chain))
    
    return out_pdb

def extract_chain_from_pdb(pdb, chainID):

    out_pdb = '{}_chain{}'.format(pdb, chainID)

    collected_atoms = []
    collected_res_nums = []

    mod_pdb = remove_alternative_conformations(pdb)

    with open(out_pdb, 'w') as out:
        with open(mod_pdb, 'r') as inp:
            for line in inp:
                if 'ATOM' in line:
                    if line[21] == chainID:
                        line = list(line)
                        line[16] = ' ' # corrects for a VERY rare, but VERY STUPID event when VERY STUPID people put the chain name also before the resname, breaking everything! IF YOU DO THAT, I HATE YOU!
                        line = ''.join(line)

                        # now check that this is not a multiple conformation. if so, only take the first one
                        atom_name = line[13:16].strip()
                        resnum = line[22:27].strip()
                        if len(collected_atoms) == 0 or atom_name != collected_atoms[-1]:
                            out.write(line)
                            collected_atoms.append(atom_name)

    os.replace(out_pdb, pdb)

    return pdb    

def remove_alternative_conformations(pdb):

    out_pdb = '{}_modified.pdb'.format(pdb.split('.')[0])

    run_pdbtools = sp.Popen([phenix_pdbtools, pdb, "remove='altloc B'", "file_name='{}'".format(out_pdb)], stdout=sp.PIPE, stderr=sp.PIPE, stdin=sp.PIPE)
    stdout, stderr = run_pdbtools.communicate()

    if (os.path.isfile(out_pdb) and os.stat(out_pdb).st_size == 0) or not os.path.isfile(out_pdb):
        out_pdb = pdb

    return out_pdb    

def collect_atoms_coordinates(pdb):

    atoms = {}

    with open(pdb, 'r') as inp:
        for line in inp:
            if 'ATOM' in line:
                atom_name = line[13:17].strip()
                resid = line[17:21].strip()
                resnum = line[22:27].strip()

                res = '{}{}'.format(resid, resnum)

                if res not in atoms:
                    atoms[res] = {}

                atoms[res][atom_name] = [float(line[30:39].strip()), float(line[39:46].strip()), float(line[46:55].strip())]

    return atoms

def parse_muscle_out(clustal_out):
    
    alignment = []

    seq = ''
    for line in clustal_out.split('\n'):
        if line.startswith('>'):
            if seq != '':
                alignment.append(seq)
            seq = ''
        else:
            seq += line.strip()
    alignment.append(seq)
    
    return alignment

def get_overlap(s1, s2, label):

    # save sequences to file
    tmp_fasta = '{}alignment.muscle'.format(label)
    with open(tmp_fasta, 'w') as fp:
        fp.write('>s1\n{}\n>s2\n{}'.format(s1, s2))

    # align blades with muscle
    run_muscle = sp.Popen([muscle, '-in', tmp_fasta], stdout=sp.PIPE, stderr=sp.PIPE, stdin=sp.PIPE)
    stdout, stderr = run_muscle.communicate()
    os.remove(tmp_fasta)

    # parse clustal output
    alignment = parse_muscle_out(stdout.decode('ascii'))
    # for i in alignment:
    #     print(i)

    match = ''
    for i in range(len(alignment[0])):
        if alignment[0][i] != '-' and alignment[1][i] != '-':
            match += alignment[0][i]
        elif alignment[0][i] == 'X' and alignment[1][i] == '-' and len(match) > 5:
            break

    #match = ''.join([alignment[0][i] for i in range(len(alignment[0])) if alignment[0][i] != '-' and alignment[1][i] != '-'])

    return match

def truncate_pdb(pdb, template_seq, reset_resnum = False, label=None):

    pdb_data = run_dssp(pdb, out_label=pdb.split('/')[-1])

    # fill in unmodelled residues with X
    pdb_resnum = [int(i) for i in pdb_data['ResNum']]
    pdb_seq = []
    for i in range(pdb_resnum[0], pdb_resnum[-1]+1):
        if i in pdb_resnum:
            pdb_seq.append(pdb_data['AA'][pdb_resnum.index(i)])

        else:
            pdb_seq.append('X')

    pdb_seq = ''.join(pdb_seq)
    
    overlap_seq = get_overlap(pdb_seq, template_seq, label = '{}-{}'.format(pdb.split('/')[-1], label).replace('-', ''))

    # remove all X now
    overlap_seq = overlap_seq.replace('X','')
    pdb_seq = pdb_seq.replace('X','')

    # and now find the real interval in the pdb
    overlap_index = pdb_seq.index(overlap_seq)

    interval = []
    for i in range(len(overlap_seq)):
        if pdb_seq[overlap_index+i] == overlap_seq[i]:
            interval.append(pdb_resnum[overlap_index+i])

    interval = [interval[0], interval[-1]]
    interval = [int(i) for i in interval]

    if label is not None:
        truncated = '{}_{}-{}_truncated_to_{}.pdb'.format(pdb.split('.')[0], interval[0], interval[1], label)
    else:
        truncated = '{}_{}-{}_truncated.pdb'.format(pdb.split('.')[0], interval[0], interval[1])

    if not os.path.isfile(truncated):
        with open(truncated, 'w') as outp:
            found_first_res = False
            with open(pdb, 'r') as inp:
                for line in inp:
                    if line.startswith('ATOM'):
                        curr_resnum = int(line[22:26].strip())

                        if curr_resnum >= interval[0] and curr_resnum <= interval[1]:
                            outp.write(line)

    return truncated

# 1. Download stuff from CASP

def download_tables(ftp_tarballs_results = predictions_tables, output_folder = tmp_folder):

    print('\nDownloading {} predictions tables'.format(which_casp))

    if not os.path.isdir(output_folder):
        os.mkdir(output_folder)
        
    opener.open(ftp_tarballs_results)
    urllib.request.install_opener(opener)

    downloaded_tables = []

    with urllib.request.urlopen(ftp_tarballs_results) as f:
        lines = f.read().decode('utf-8')
        for line in lines.split('\n'):
            if 'alt="[TXT]"' in line:
                table_file = line.split('href=')[-1].split('>')[0].replace('"', '')
                if table_file.startswith('T') or table_file.startswith('R'):
                    print(' ... Downloading {}'.format(table_file))

                    table_link = '{}/{}'.format(ftp_tarballs_results, table_file)
                    table_outp = '{}/{}'.format(output_folder, table_file)

                    if not os.path.isfile(table_outp):
                        urllib.request.urlretrieve(table_link, table_outp)

                    downloaded_tables.append(table_outp)
    
    print(' --> Downloaded {} tables'.format(len(downloaded_tables)))

    return downloaded_tables

def download_targets(ftp_tarballs_targets = [targets_trimmed, targets_processed, targets_submitted, all_targets], output_folder = tmp_folder):

    print('\nDownloading {} targets and EUs'.format(which_casp))

    if not os.path.isdir(output_folder):
        os.mkdir(output_folder)
        
    for targets_link in ftp_tarballs_targets:

        print(' ... Checking {}'.format(targets_link))

        if targets_link != None:

            opener.open(targets_link)
            urllib.request.install_opener(opener)

            if 'gz' in targets_link:
                outp = '{}/{}'.format(output_folder, targets_link.split('/')[-1])
                urllib.request.urlretrieve(targets_link, outp)

                extract_structures = sp.Popen(['tar', '-xvzf', outp, '-C', output_folder], stdout=sp.PIPE, stderr=sp.PIPE, stdin=sp.PIPE)
                stdout, stderr = extract_structures.communicate()

                possible_outfolder = '{}/{}'.format(output_folder, targets_link.split('/')[-1].replace('.tar.gz', '').replace('.tgz', ''))

                if os.path.isdir(possible_outfolder) and 'submitted' not in possible_outfolder:
                    for pdb in os.listdir(possible_outfolder):
                        os.replace('{}/{}'.format(possible_outfolder, pdb), "{}/{}".format(output_folder, pdb))
                    shutil.rmtree(possible_outfolder)
                elif 'submitted' in possible_outfolder:
                    os.replace(possible_outfolder, "{}/Targets_submitted".format(output_folder))

                os.remove(outp)

            else:
                with urllib.request.urlopen(targets_link) as f:
                    lines = f.read().decode('utf-8')
                    for line in lines.split('\n'):
                        if 'alt="[   ]"' in line and '.pdb' in line:
                            file = line.split('href=')[-1].split('>')[0].replace('"', '')
                            outp = '{}/{}'.format(output_folder, file)

                            if (file.startswith('R') or file.startswith('T')) and not os.path.isfile(outp):

                                print(' ... ... Downloading {} ({})'.format(file, outp))

                                new_link = '{}{}'.format(targets_link, file)
                                opener.open(new_link)
                                urllib.request.install_opener(opener)

                                urllib.request.urlretrieve(new_link, outp)

    downloaded_targets = ['{}/{}'.format(output_folder, pdb) for pdb in os.listdir(output_folder) if (pdb.endswith('pdb') and (pdb.startswith('T') or pdb.startswith('R')))]
    
    print('Checking if all pdbs are fine')
    downloaded_targets = fix_pdbs(downloaded_targets)
    
    print(' --> Downloaded {} target structures (including the individual domains)'.format(len(downloaded_targets)))

    return output_folder

def download_templates_information(output_folder = tmp_folder, templates_link = templates):

    print('\nDownloading {} templates information'.format(which_casp))

    if not os.path.isdir(output_folder):
        os.mkdir(output_folder)

    lga_folder = '{}/LGA'.format(output_folder)
    if not os.path.isdir(lga_folder):
        os.mkdir(lga_folder)

    hhsearch_folder = '{}/HHsearches'.format(output_folder)
    if not os.path.isdir(hhsearch_folder):
        os.mkdir(hhsearch_folder)
        
    # check what is in the templates folder
    opener.open(templates_link)
    urllib.request.install_opener(opener)    

    with urllib.request.urlopen(templates_link) as f:
        lines = f.read().decode('utf-8')
        for line in lines.split('\n'):
            if 'alt="[DIR]"' in line or 'alt="[   ]"' in line or 'alt="[TXT]"' in line:
                file = line.split('href=')[-1].split('>')[0].replace('"', '')

                new_link = '{}/{}'.format(templates, file)
                opener.open(new_link)
                urllib.request.install_opener(opener)

                if 'HH' in file:
                    # we are dealing with hh searches files
                    if file.endswith('/'):
                        # this is a directory and we shall download its contents to the hhsearch_folder
                        with urllib.request.urlopen(new_link) as f:
                            lines = f.read().decode('utf-8')
                            for line in lines.split('\n'):
                                if 'alt="[   ]"' in line:
                                    curr_file = line.split('href=')[-1].split('>')[0].replace('"', '')
                                    if curr_file.endswith('.hhr'):
                                        curr_file_link = '{}/{}'.format(new_link, curr_file)
                                        curr_file_outp = '{}/{}'.format(hhsearch_folder, curr_file)

                                        if not os.path.isfile(curr_file_outp):
                                            urllib.request.urlretrieve(curr_file_link, curr_file_outp)

                    elif file.endswith('.tgz'):
                        # this is a tar file and we shall download it, extract it and save the files we want
                        outp = '{}/{}'.format(output_folder, file)
                        urllib.request.urlretrieve(new_link, outp)

                        extract_structures = sp.Popen(['tar', '-xvzf', outp, '-C', hhsearch_folder], stdout=sp.PIPE, stderr=sp.PIPE, stdin=sp.PIPE)
                        stdout, stderr = extract_structures.communicate()

                        possible_outfolders = ['{}/{}'.format(hhsearch_folder, f) for f in os.listdir(hhsearch_folder) if os.path.isdir('{}/{}'.format(hhsearch_folder, f))]

                        if len(possible_outfolders) > 0:
                            for possible_outfolder in possible_outfolders:
                                for curr_file in os.listdir(possible_outfolder):
                                    if curr_file.endswith('.hhr'):
                                        os.replace('{}/{}'.format(possible_outfolder, curr_file), "{}/{}".format(hhsearch_folder, curr_file))
                                shutil.rmtree(possible_outfolder)

                        os.remove(outp)

                elif '.csv' in file:
                    # we are dealing with a lga csv file and we shall download it to the LGA folder
                    outp = '{}/{}'.format(lga_folder, file)
                    urllib.request.urlretrieve(new_link, outp)

                elif 'LGA' in file and file.endswith('/'):
                    # we are dealing with a folder and shall download its contents to the LGA folder
                    with urllib.request.urlopen(new_link) as f:
                        lines = f.read().decode('utf-8')
                        for line in lines.split('\n'):
                            if 'alt="[TXT]"' in line:
                                curr_file = line.split('href=')[-1].split('>')[0].replace('"', '')
                                if '.csv' in curr_file:
                                    curr_file_link = '{}/{}'.format(new_link, curr_file)
                                    curr_file_outp = '{}/{}'.format(lga_folder, curr_file)

                                    if not os.path.isfile(curr_file_outp):
                                        urllib.request.urlretrieve(curr_file_link, curr_file_outp)

    return output_folder

def download_all_models(target, out_folder = tmp_folder):

    if '-D' in target and target.startswith('T'):
        ftp_tarballs_results = predictions_domains
    else:
        ftp_tarballs_results = predictions_fulleng

    in_out_folder = out_folder
    if not os.path.isdir(in_out_folder):
        os.mkdir(in_out_folder)

    out_folder = '{}/{}'.format(in_out_folder, target)
    if not os.path.isdir(out_folder):
        os.mkdir(out_folder)
        
    models_link = '{}/{}.tar.gz'.format(ftp_tarballs_results, target)
    model_outp = '{}/{}.tar.gz'.format(out_folder, target)
    models = '{}/{}'.format(out_folder, target)

    opener.open(ftp_tarballs_results)
    urllib.request.install_opener(opener)

    if not os.path.isdir(models):
        print(' ... Downloading models')
        try:
            if not os.path.isfile(model_outp):
                urllib.request.urlretrieve(models_link, model_outp)

            extract_structures = sp.Popen(['tar', '-xvzf', model_outp, '-C', out_folder], stdout=sp.PIPE, stderr=sp.PIPE, stdin=sp.PIPE)
            stdout, stderr = extract_structures.communicate()

            os.remove(model_outp)
        except:
            models = 'nan'

    elif update:
        print(' ... Updating models')
        if not os.path.isfile(model_outp):
            urllib.request.urlretrieve(models_link, model_outp)

        new_tmp_folder = '{}/{}_tmp'.format(out_folder, target)
        if not os.path.isdir(new_tmp_folder):
            os.mkdir(new_tmp_folder)

        extract_structures = sp.Popen(['tar', '-xvzf', model_outp, '-C', new_tmp_folder], stdout=sp.PIPE, stderr=sp.PIPE, stdin=sp.PIPE)
        stdout, stderr = extract_structures.communicate()

        os.remove(model_outp)

        for file in [f for f in os.listdir('{}/{}'.format(new_tmp_folder, target))]:
            os.replace('{}/{}/{}'.format(new_tmp_folder, target,file), '{}/{}'.format(models, file))
    
    if '-D' in target:
        tmp, full_models = download_all_models(target.split('-')[0], out_folder = in_out_folder)
    else:
        full_models = models

    return models, full_models

# 2. Parsers

def parse_targets_pdbs(targets_summary = targets_summary):

    print('\nChecking if there is any target with already deposited structures')

    targets_pdbs = {}
    targets_summary = pd.read_csv(targets_summary, error_bad_lines = False, sep = ';')

    for index, row in targets_summary.iterrows():
        target = row.Target
        description = row.Description

        if '<em>' in description:
            description = description.split('<em>')[0]

        if len(description.split()) > 1:
            if len(description.split()[-1]) == 4:
                pdb = description.split()[-1]
                if not pdb.isalpha() and not pdb.isnumeric(): # only acept mixes of numbers and letters
                    description = description.replace(pdb, '')
                    pdb = pdb.upper()
                else:
                    pdb = None
            else:
                pdb = None
        else:
            pdb = None

        if row['Cancellation Date'] == '-':
            status = 'Active'
        else:
            status = 'Cancelled'

        expiration_date = row['Human Exp.'].replace('-','')

        targets_pdbs[target] = {'Title': description, 'PDB_code': pdb, 'Status': status, 'Exp.Date': expiration_date}

    return targets_pdbs

def parse_targets_classification(classifications = domains_classifications):

    print('\nCollecting {} EUs classifications'.format(which_casp))

    targets_classification = {}

    with open(classifications, 'r') as inclassf:
        for line in inclassf:
            line = line.split(';')
            if len(line) > 1:
                input_target = line[0]
                domains = [d.replace(':', '') for d in line[3].split() if 'D' in d]
                if len(domains) > 0:
                    domains_class = line[4]
                    for domain in domains:
                        classification = domains_class.split(domain)[-1].split('D')[0].replace(':','').strip()
                        if len(classification) < 2:
                            classification = np.nan

                        if domain == 'D0':
                            target = input_target
                        else:
                            target = '{}-{}'.format(input_target, domain)

                        targets_classification[target] = classification

                    if len(domains) == 1:
                        targets_classification[input_target] = classification
                    elif input_target not in targets_classification:
                        targets_classification[input_target] = 'MultiDom'

    return targets_classification

def parse_templates_Neff(neffs = sequence_deep):

    print('\nCollecting {} EUs maximum Neff'.format(which_casp))

    templates_neff = {}

    if neffs is not None:
        neffs = pd.read_csv(neffs)

        for index, row in neffs.iterrows():
            target = row.domain
            neff_max = float(row['Neff/len_MAX'])

            templates_neff[target] = neff_max

    return templates_neff

def parse_table(table):
    
    print(' ... Parsing models table')

    all_data = {}
    
    with open(table, 'r') as intable:
        for line in intable:
            if line.startswith('#') and len(all_data) == 0:
                columns = [i for i in line.split(' ') if i != '' and i != '\n']
                all_data = {i: [] for i in columns}
            elif not line.startswith('#'):
                values = [i for i in line.split(' ') if i != '' and i != '\n']
                for j, value in enumerate(values):
                    if j > 2:
                        try:
                            value = float(value)
                        except:
                            value = value
                    try:
                        all_data[columns[j]].append(value)
                    except:
                        pass
    
    all_data = pd.DataFrame(all_data)
    all_data = pd.DataFrame(all_data).drop(columns=['#', 'RANK']).replace(to_replace='N/A', value=np.nan).replace(to_replace='-', value=np.nan)

    return all_data

# 3. For secondary structure

def parse_DSSPout(dsspout_file, mask_loops = False):

    found_table = False
    data_dict = {'ResNum':[], 'AA': [], 'AA3':[], 'SecStr': [], 'Phi': [], 'Psi': [], 'ACC': []}
    
    with open(dsspout_file, "r") as dsspout:
        for line in dsspout:
            if "#  RESIDUE AA STRUCTURE BP1" in line:
                found_table = True
            elif found_table:
                resnum = line[5:10].strip()
                secstr = line[16]
                res = line[13]
                phi = line[103:109].strip()
                psi = line[110:116].strip()
                acc = line[34:40].strip()
                
                if mask_loops and secstr in ['S', 'T', 'B', 'G', ' ']:
                    secstr = 'l'

                if res.isalpha():
                    data_dict['ResNum'].append(resnum)
                    data_dict['AA'].append(res)
                    data_dict['AA3'].append('{}{}'.format(seq3(res).upper(), resnum))
                    data_dict['SecStr'].append(secstr)
                    data_dict['Phi'].append(float(phi))
                    data_dict['Psi'].append(float(psi))
                    data_dict['ACC'].append(float(acc))
                     
    return data_dict

def run_dssp(inpdb, tmp_folder = '/tmp/', out_label = 'na'):
    
    print(' ... Running DSSP')

    dssp_outf = '{}{}_dssp.out'.format(tmp_folder, out_label)

    if not os.path.isfile(dssp_outf):
        run_dssp = sp.Popen([dssp, '-i', inpdb, '-o', dssp_outf], stdout=sp.PIPE, stderr=sp.PIPE, stdin=sp.PIPE)
        stdout, stderr = run_dssp.communicate()
                
    dssp_out = parse_DSSPout(dssp_outf, mask_loops = True)

    os.remove(dssp_outf)

    return dssp_out

# 4. For difficulty computations

def parse_hhfile(hhsearch_file, res_interval, exclude = None, max_date = None, min_cov = 0):

    prob = np.nan
    res_range = range(res_interval[0], res_interval[-1]+1)

    found_table = False
    prob = 0
    with open(hhsearch_file, 'r') as inhhr:
        for line in inhhr:
            if line.startswith(' No Hit'):
                index_prob = line.find('Prob')
                found_table = True
            elif found_table and len(line.strip()) > 0:
                curr_pdb = line.split()[1].split('_')[0]
                curr_prob = float(line[index_prob-1:].split()[0].strip())
                curr_interval = line[index_prob-1:].split()[-2].split('-')
                curr_interval = [int(i) for i in curr_interval]
                curr_range = range(curr_interval[0], curr_interval[-1]+1)

                curr_coverage = len(set(res_range).intersection(set(curr_range)))/len(res_range)

                pdb_data = json.loads(p.PDB.getSummary(pdbid=curr_pdb.lower()))
                pdb_release_date = pdb_data[curr_pdb.lower()][0]['release_date']

                if curr_prob > prob and curr_coverage > min_cov and (exclude is None or curr_pdb != exclude) and (max_date is None or pdb_release_date < max_date):
                    prob = curr_prob
            else:
                found_table = False

    return prob

def parse_lga_table(lga_file, exclude = None, max_date = None):

    lga = np.nan

    with open(lga_file, 'r') as inlga:
        for line in inlga:
            if np.isnan(lga):
                pdb = line.split(',')[0]
                pdb_data = json.loads(p.PDB.getSummary(pdbid=pdb.split('_')[0]))
                pdb_release_date = pdb_data[pdb.split('_')[0]][0]['release_date']

                if (exclude is None or exclude.lower() not in pdb) and (pdb_release_date < max_date or max_date is None):
                    lga = float(line.split(',')[-1])
                    break

    return lga

def compute_difficulty(target, templates_folder, deposited_pdb, max_date, res_interval, targets_classification = {}):

    print(' ... Computing {} difficulty'.format(target))

    hhsearch_file = ['{}/HHsearches/{}'.format(templates_folder, f) for f in os.listdir('{}/HHsearches'.format(templates_folder)) if target.split('-')[0] in f]
    lga_file = ['{}/LGA/{}'.format(templates_folder, f) for f in os.listdir('{}/LGA'.format(templates_folder)) if f.split('.csv')[0] == target]

    difficulty = np.nan

    if len(hhsearch_file) > 0:
        if len(lga_file) > 0:
            hhsearch_file = hhsearch_file[0]
            lga_file = lga_file[0]

            best_hhprob = parse_hhfile(hhsearch_file, res_interval, exclude = deposited_pdb, max_date = max_date)
            best_lga = parse_lga_table(lga_file, exclude = deposited_pdb, max_date = max_date)

            difficulty = (best_hhprob+best_lga)/2.0

        else:
            if '-' in target:
                difficulty, target_classification = compute_difficulty(target.split('-')[0], templates_folder, deposited_pdb, max_date, res_interval, targets_classification=targets_classification)
            else:
                pass

    if target in targets_classification:
        target_classification = targets_classification[target]
    else:
        target_classification = np.nan

    print(' ... ... Difficulty: {} (classified as "{}")'.format(difficulty, target_classification))

    return difficulty, target_classification

# 5. For loop analysis

def get_loops(dssp_out, min_size = 10, ss = 'l'):
    
    print(' ... Extracting loops with at least {} aa'.format(min_size))
    
    loops_found = []
    
    curr_loop = []
    curr_resi = []
    for i, aa in enumerate(dssp_out['AA3']):
        if dssp_out['SecStr'][i] == ss and (len(curr_resi) == 0 or int(dssp_out['ResNum'][i]) - curr_resi[-1] == 1):
            curr_loop.append(aa)
            try:
                curr_resi.append(int(dssp_out['ResNum'][i]))
            except:
                print(dssp_out['ResNum'][i])
                
        else:
            if len(curr_loop) >= min_size:
                loops_found.append(curr_loop)
            curr_loop = []
            curr_resi = []
    
    if len(curr_loop) >= min_size:
        loops_found.append(curr_loop)
    
    print(' ... ... Found {} loops'.format(len(loops_found)))
    
    return loops_found

def extract_loop(target_structure, loop, out_label = None, out_folder = tmp_folder):
    
    out_pdb = '{}/{}_loop{}.pdb'.format(out_folder, target_structure.split('/')[-1].replace('.pdb',''), out_label)
    
    if not os.path.isfile(out_pdb):
        with open(out_pdb, 'w') as outpdb:
            with open(target_structure, 'r') as inpdb:
                for line in inpdb:
                    if line.startswith('ATOM '):
                        resid = line[17:21].strip()
                        resnum = line[22:27].strip()
                        res = '{}{}'.format(resid, resnum)
                        if res in loop:
                            outpdb.write(line)
    return out_pdb


# 6. For DipDiff

def parse_dipcheck(dipcheck_out):
        
    results = {}
    chiscore = np.nan
    
    with open(dipcheck_out, 'r') as dipch:
        for line in dipch:
            if 'DipScore' in line and 'Percentile' in line:
                resid = line.split()[0]
                if not line.split()[1].isdigit():
                    chain = line.split()[1]
                else:
                    chain = ''
                resnum = line.split('DipScore')[0].split()[-1].strip()
                dipscore = line.split('DipScore')[1].split()[0]
                
                results['{}{}'.format(resid, resnum)] = float(dipscore)
            
            elif 'Overall Chi-score Percentile' in line:
                chiscore = float(line.split()[4].strip())
    
    return results, chiscore

def run_dipcheck(target):
    
    out_folder = '{}/dipcheck'.format('/'.join(target.split('/')[:-1]))
    #target = set_bfactor(target, out_folder = out_folder)
    target = set_occ(target, out_folder = out_folder)

    out_dipcheck = '{}.dipout'.format(target)

    if not os.path.isfile(out_dipcheck):

        run_dipcheck = sp.Popen([dipcheck, 'xyzin', target], stdout=sp.PIPE, stderr=sp.PIPE, stdin=sp.PIPE)
        stdout, stderr = run_dipcheck.communicate()

        with open(out_dipcheck, 'w') as outp:
            outp.write(stdout.decode("utf-8"))
        
    dipcheck_out, chiscore = parse_dipcheck(out_dipcheck)
    
    return dipcheck_out, chiscore

def run_dipdiff(target, model, target_interval = None, mode = None, make_plot = True):

    if model != None:
    
        target_dipcheck_out, target_chiscore = run_dipcheck(target)
        model_dipcheck_out, model_chiscore = run_dipcheck(model)
        
        if target_interval == None:
            residues_in_target = [key for key in target_dipcheck_out.keys()]
            residues_in_model = [key for key in model_dipcheck_out.keys()]

            chidiff = model_chiscore-target_chiscore

        else:
            residues_in_target = [key for key in target_dipcheck_out.keys() if key in target_interval]
            residues_in_model = [key for key in model_dipcheck_out.keys() if key in target_interval]

            chidiff = np.nan

        if mode == None:
            modelled_residues = [key for key in residues_in_target if key in residues_in_model]
        elif mode == 'template':
            modelled_residues = align_structures(target, model)

        dipdifferences = []
        accepted_residues = []
        target_scores = []
        model_scores = []
        for residue in modelled_residues:
            if type(modelled_residues) is dict:
                model_residue = modelled_residues[residue]
            else:
                model_residue = residue

            if residue in target_dipcheck_out and model_residue in model_dipcheck_out:
                dipdifferences.append(model_dipcheck_out[model_residue] - target_dipcheck_out[residue])
                accepted_residues.append(residue)
                target_scores.append(target_dipcheck_out[residue])
                model_scores.append(model_dipcheck_out[model_residue])

        if make_plot and len(accepted_residues) > 0:

            out_label = '{}_{}_DipDiff'.format(target.replace('.pdb', ''), model.split('/')[-1].replace('.pdb', ''))

            # plot dipscores
            plt.clf()
            fig = plt.figure(figsize = (10, 2))
            plt.plot(accepted_residues, target_scores, c = 'lightblue', label = '{}'.format(target.split('/')[-1]))
            plt.plot(accepted_residues, model_scores, c = 'steelblue', label = '{}'.format(model.split('/')[-1]))
            plt.xticks([], [])
            plt.ylabel('DipScore')
            plt.xlabel('Residues')
            plt.savefig('{}_DipCheck_scores.png'.format(out_label))
            plt.close()
            
            # plot dipscore difference
            plt.clf()
            fig = plt.figure(figsize = (10, 2))
            plt.plot(accepted_residues, dipdifferences, c = 'steelblue', label = '{}'.format(model.split('/')[-1]))
            plt.hlines(y=0, xmin=accepted_residues[0], xmax=accepted_residues[-1], linestyle = ':')
            plt.xticks([], [])
            plt.ylim(-1.1,1.1)
            plt.ylabel('DipScore difference')
            plt.xlabel('Residues')
            plt.savefig('{}_DipScore_diff.png'.format(out_label))
            plt.close()
            
            # plot histogram of dipscore differences
            plt.clf()
            fig = plt.figure(figsize = (4, 3))
            plt.hist(dipdifferences, color = 'steelblue', bins = 100 )
            plt.xlabel('DipScore difference')
            plt.ylabel('Count')
            plt.vlines(x=np.mean([i for i in dipdifferences if not np.isnan(i)]), ymin = 0, ymax = 50, color = 'red')
            plt.title('Mean/DipDiff: {}'.format(np.mean([i for i in dipdifferences if not np.isnan(i)])))
            plt.xlim(-1.1, 1.1)
            plt.savefig('{}_DipScore_diff_hist.png'.format(out_label))
            plt.close()

        dipdiffs = {'DipDiff': round(np.mean([i for i in dipdifferences if not np.isnan(i)]), 5), 'ChiDiff': chidiff}

    else:
        dipdiffs = {'DipDiff': np.nan, 'ChiDiff': np.nan}

    return dipdiffs

def run_loops_dipdiff(target_structure, curr_model, target_loops, target, out_folder = tmp_folder):

    if not os.path.isdir(out_folder):
        os.mkdir(out_folder)

    loops_folder = '{}/{}'.format(out_folder, target_structure.split('/')[-1].split('.pdb')[0])
    if not os.path.isdir(loops_folder):
        os.mkdir(loops_folder)

    loops_dipdiff = {}
    if len(target_loops) > 0:
        for i, loop in enumerate(target_loops):            
            loopID = '{}_{}-{}'.format(target.split('-')[0], loop[0], loop[-1])

            # extract target loop structure
            target_loop_structure = extract_loop(target_structure, loop, out_label = i, out_folder = loops_folder)
            # extract modelled loop structure
            model_loop_structure = extract_loop(curr_model, loop, out_label = i, out_folder = loops_folder)

            loop_dipdiff = run_dipdiff(target_loop_structure, model_loop_structure)
            loops_dipdiff[loopID] = loop_dipdiff

    return loops_dipdiff

# 7. For GDT_TS and LGA

def generate_lga_input(target, model):

    working_folder = '/'.join(target.split('/')[:-1])
    output_folder = '{}/MOL2'.format(working_folder)
    if not os.path.isdir(output_folder):
        os.mkdir(output_folder)

    targetID = target.split('/')[-1].replace('.pdb', '')
    modelID  = model.split('/')[-1].replace('.pdb', '')

    lga_input = '{}.{}'.format(targetID, modelID)
    lga_input_file = '{}/{}'.format(output_folder, lga_input)

    if not os.path.isfile(lga_input_file):
        with open(lga_input_file, 'w') as lgain:
            lgain.write('MOLECULE  {}\n'.format(targetID))
            with open(target, 'r') as intarget:
                for line in intarget:
                    lgain.write(line)
            lgain.write('END\n')

            lgain.write('MOLECULE  {}\n'.format(modelID))
            with open(model, 'r') as inmodel:
                for line in inmodel:
                    lgain.write(line)
            lgain.write('END\n')
    
    return lga_input, working_folder

def parse_predictions_gdt(gdtout):

    gdt_data = {'GDT_TS': np.nan, 'GDT_HA': np.nan, 'GDC_SC': np.nan}

    with open(gdtout, 'r') as gdt_out:
        for line in gdt_out:
            if '=' in line:
                if line.split('=')[0].strip() in gdt_data:
                    gdt_data[line.split('=')[0].strip()] = float(line.split('=')[1].strip())

    return gdt_data

def run_predictions_gdt(target, model):

    infile, wrkdir = generate_lga_input(target, model)

    outfile = '{}/{}_aligned.gdt'.format(wrkdir, infile)

    if not os.path.isfile(outfile):
    
        run_gdt = sp.Popen([predictions_gdt, wrkdir, infile], stdout=sp.PIPE, stderr=sp.PIPE, stdin=sp.PIPE)
        stdout, stderr = run_gdt.communicate()

        with open(outfile, 'w') as outp:
            outp.write(stdout.decode("utf-8"))
    
    gdt_out = parse_predictions_gdt(outfile)

    return gdt_out

def run_loops_gdt(target_structure, curr_model, target_loops, target, out_folder = tmp_folder):

    if not os.path.isdir(out_folder):
        os.mkdir(out_folder)

    loops_folder = '{}/{}'.format(out_folder, target_structure.split('/')[-1].split('.pdb')[0])
    if not os.path.isdir(loops_folder):
        os.mkdir(loops_folder)

    loops_gdt = {}

    if len(target_loops) > 0:
        for i, loop in enumerate(target_loops):            
            loopID = '{}_{}-{}'.format(target.split('-')[0], loop[0], loop[-1])

            # extract target loop structure
            target_loop_structure = extract_loop(target_structure, loop, out_label = i, out_folder = loops_folder)
            # extract modelled loop structure
            model_loop_structure = extract_loop(curr_model, loop, out_label = i, out_folder = loops_folder)
            loop_gdt = run_predictions_gdt(target_loop_structure, model_loop_structure)
            loops_gdt[loopID] = loop_gdt

    return loops_gdt

def parse_lga_residue_correspondence(lga_out):

    residue_correspondence = {}

    found_table = False
    with open(lga_out, 'r') as lgaout:
        for line in lgaout:
            if '#      Molecule1      Molecule2' in line:
                found_table = True

            elif found_table:
                if line.startswith('LGA') and line.split()[1].strip() != '-' and line.split()[3].strip() != '-':
                    mol1_res = '{}{}'.format(seq3(line.split()[1].strip()).upper(), line.split()[2].strip())
                    mol2_res = '{}{}'.format(seq3(line.split()[3].strip()).upper(), line.split()[4].strip()).split('_')[0]

                    residue_correspondence[mol1_res] = mol2_res

    return residue_correspondence

def align_structures(target, model):

    infile, wrkdir = generate_lga_input(target, model)
    outfile = '{}/{}_aligned.lga'.format(wrkdir, infile)

    if not os.path.isfile(outfile):
        run_lga = sp.Popen([templates_lga, wrkdir, infile], stdout=sp.PIPE, stderr=sp.PIPE, stdin=sp.PIPE)
        stdout, stderr = run_lga.communicate()

        with open(outfile, 'w') as out:
            out.write(stdout.decode("utf-8"))

    residue_correspondence = parse_lga_residue_correspondence(outfile)

    return residue_correspondence

# 8. For Tristan's and Randy's geometric scores

def prepare_pdb_for_molprobity(pdb):

    out_pdb = '{}_molp'.format(pdb)

    with open(out_pdb, 'w') as out:
        with open(pdb, 'r') as inp:
            for line in inp:
                if 'ATOM' in line:
                    out.write(line)

    return out_pdb


def parse_molprobity_rama(molpout):
    
    phi_psi = {}

    found_table = False

    with open(molpout, 'r') as molp:
        for line in molp:
            if 'residue:score' in line:
                found_table = True
            elif 'SUMMARY:' in line:
                found_table = False
            
            elif found_table:
                res = line.split(':')[0]
                res = '{}{}'.format(res.split()[-1].strip(), res.split()[-2].strip())
                phi = float(line.split(':')[2])
                psi = float(line.split(':')[3])

                phi_psi[res] = {'Phi': np.radians(phi), 'Psi': np.radians(psi)}
    
    return phi_psi

def parse_molprobity_omega(molpout, phi_psi):

    found_table = False

    with open(molpout, 'r') as molp:
        for line in molp:
            if 'residues:type' in line:
                found_table = True
            elif 'SUMMARY:' in line:
                found_table = False
            
            elif found_table:
                res = line.split(':')[0]
                res = '{}{}'.format(res.split()[-1].strip(), res.split()[-2].strip())
                omega = float(line.split(':')[2])

                if res in phi_psi:
                    phi_psi[res]['Omega'] = np.radians(omega)
                else:
                    phi_psi[res] = {'Omega': np.radians(omega)}

    return phi_psi

def parse_molprobity_chiangles(molpout):

    chi1_chi2 = {}

    found_table = False

    with open(molpout, 'r') as molp:
        for line in molp:
            if 'residue:occupancy' in line:
                found_table = True
            elif 'SUMMARY:' in line:
                found_table = False

            elif found_table:
                res = line.split(':')[0]
                res = '{}{}'.format(res.split()[-1].strip(), res.split()[-2].strip())

                try:
                    chi1 = np.radians(float(line.split(':')[3]))
                except:
                    chi1 = np.nan

                try:
                    chi2 = np.radians(float(line.split(':')[4]))
                except:
                    chi2 = np.nan

                n_torsions = len([value for value in line.split(':')[3:-2] if value != ''])

                chi1_chi2[res] = {'n_chis': n_torsions, 'Chi1': chi1, 'Chi2': chi2}

    return chi1_chi2

def calculate_backbone_dihedrals(pdb):

    prepared_target = prepare_pdb_for_molprobity(pdb)

    rama_out = '{}_molprobity_ramangles.txt'.format(prepared_target)
    omega_out = '{}_molprobity_omegangles.txt'.format(prepared_target)

    if not os.path.isfile(rama_out):
        run_molprobity = sp.Popen([molprobity_ramangles, prepared_target], stdout=sp.PIPE, stderr=sp.PIPE, stdin=sp.PIPE)
        stdout, stderr = run_molprobity.communicate()

        with open(rama_out, 'w') as ramaout:
            ramaout.write(stdout.decode("utf-8"))

    if not os.path.isfile(omega_out):
        run_molprobity = sp.Popen([molprobity_omegangles, prepared_target, 'nontrans_only=False'], stdout=sp.PIPE, stderr=sp.PIPE, stdin=sp.PIPE)
        stdout, stderr = run_molprobity.communicate()

        with open(omega_out, 'w') as omegaout:
            omegaout.write(stdout.decode("utf-8"))

    phi_psi = parse_molprobity_rama(rama_out)
    phi_psi_omega = parse_molprobity_omega(omega_out, phi_psi = phi_psi)

    return phi_psi_omega

def get_backbone_geometric_score(target, model, mode = None):

    bb_sc = {}

    target_bbangles = calculate_backbone_dihedrals(target)
    model_bbangles  = calculate_backbone_dihedrals(model)

    residues_in_target = [key for key in target_bbangles.keys()]
    residues_in_model = [key for key in model_bbangles.keys()]

    if mode == None:
        modelled_residues = [key for key in residues_in_target if key in residues_in_model]
    elif mode == 'template':
        modelled_residues = align_structures(target, model)

    for residue in modelled_residues:
        if type(modelled_residues) is dict:
            model_residue = modelled_residues[residue]
        else:
            model_residue = residue

        if residue in target_bbangles and model_residue in model_bbangles:
            if len(target_bbangles[residue]) == 3 and len(model_bbangles[model_residue]) == 3:
                Sbackbone = 0

                for angle in target_bbangles[residue]:
                    Sbackbone += (1-math.cos(model_bbangles[model_residue][angle] - target_bbangles[residue][angle]))/2
                Sbackbone = Sbackbone/3
                bb_sc[residue] = Sbackbone

    bb_sc = np.mean([bb_sc[residue] for residue in bb_sc if not np.isnan(bb_sc[residue])])

    return bb_sc

def calculate_sidechain_dihedrals(pdb):

    prepared_target = prepare_pdb_for_molprobity(pdb)

    chi_out = '{}_molprobity_chiangles.txt'.format(prepared_target)

    if not os.path.isfile(chi_out):
        run_molprobity = sp.Popen([molprobity_sidechains, prepared_target], stdout=sp.PIPE, stderr=sp.PIPE, stdin=sp.PIPE)
        stdout, stderr = run_molprobity.communicate()

        with open(chi_out, 'w') as chis:
            chis.write(stdout.decode("utf-8"))

    chi1_chi2 = parse_molprobity_chiangles(chi_out)

    return chi1_chi2

def calculate_sidechain_burial_score(pdb):

    betas = {}

    atoms_in_pdb = collect_atoms_coordinates(pdb)
    
    for residue in atoms_in_pdb:
        betas[residue] = 0
        no_sc_atoms = len([atom for atom in atoms_in_pdb[residue] if atom not in ['N', 'C', 'CA', 'O']])
        for atom in atoms_in_pdb[residue]:
            if atom not in ['N', 'C', 'CA', 'O']:
                for other_residue in atoms_in_pdb:
                    for other_atom in atoms_in_pdb[other_residue]:
                        if atoms_in_pdb[residue][atom][0] - atoms_in_pdb[other_residue][other_atom][0] <= 4:
                            if atoms_in_pdb[residue][atom][1] - atoms_in_pdb[other_residue][other_atom][1] <= 4:
                                if atoms_in_pdb[residue][atom][2] - atoms_in_pdb[other_residue][other_atom][2] <= 4:

                                    distance = np.sqrt(np.dot(np.array(atoms_in_pdb[residue][atom])-np.array(atoms_in_pdb[other_residue][other_atom]), 
                                                              np.array(atoms_in_pdb[residue][atom])-np.array(atoms_in_pdb[other_residue][other_atom])))
                                    if distance <= 4:
                                        betas[residue]+=1

        if no_sc_atoms > 0:
            betas[residue] = min([betas[residue]/(3*no_sc_atoms), 1])
        else:
            betas[residue] = 0

    return betas

def get_sidechains_geometric_score(target, model, mode = None, make_plot = True):

    sc_sc = {}

    target_scangles = calculate_sidechain_dihedrals(target)
    model_scangles  = calculate_sidechain_dihedrals(model)

    target_betas    = calculate_sidechain_burial_score(target)

    residues_in_target = [key for key in target_scangles.keys()]
    residues_in_model = [key for key in model_scangles.keys()]

    if mode == None:
        modelled_residues = [key for key in residues_in_target if key in residues_in_model]
    elif mode == 'template':
        modelled_residues = align_structures(target, model)

    for residue in modelled_residues:
        if type(modelled_residues) is dict:
            model_residue = modelled_residues[residue]
        else:
            model_residue = residue

        if residue in target_scangles and model_residue in model_scangles:
            Ssidechain = np.nan

            if target_scangles[residue]['n_chis'] == 0:
                Ssidechain = 0
            else:
                delta_chi1 = model_scangles[model_residue]['Chi1']-target_scangles[residue]['Chi1']

                if target_scangles[residue]['n_chis'] == 1 and model_scangles[model_residue]['n_chis'] == 1:
                    Sdelta_chi1 = (1-math.cos(delta_chi1))/2
                    Ssidechain = target_betas[residue]*Sdelta_chi1
                else:
                    Sdelta_chi1 = (1-math.cos(delta_chi1))/2

                    delta_chi2 = model_scangles[model_residue]['Chi2']-target_scangles[residue]['Chi2']
                    Sdelta_chi2 = (1-math.cos(delta_chi2))/2
                
                    Ssidechain = target_betas[residue]*((1-(2/3)*(1-Sdelta_chi1))-(1/3)*math.exp(-(delta_chi1/np.radians(30))**2)*(1-Sdelta_chi2))

            sc_sc[residue] = Ssidechain

    sc_sc = np.mean([sc_sc[residue] for residue in sc_sc if not np.isnan(sc_sc[residue])])

    return sc_sc

def run_tristan_geometry_scores(target, model, mode = None):

    # based on Croll, Sammito, Kryshtafovych and Read (2019) Proteins
    
    if model != None:
        model_bb_geometry = get_backbone_geometric_score(target, model, mode = mode)
        model_sc_geometry = get_sidechains_geometric_score(target, model, mode = mode)

    else:
        model_bb_geometry = np.nan
        model_sc_geometry = np.nan

    return {'BBscore': model_bb_geometry, 'SCscore': model_sc_geometry}

def run_loops_tristan_scores(target_structure, curr_model, target_loops, target, out_folder = tmp_folder):

    if not os.path.isdir(out_folder):
        os.mkdir(out_folder)

    loops_folder = '{}/{}'.format(out_folder, target_structure.split('/')[-1].split('.pdb')[0])
    if not os.path.isdir(loops_folder):
        os.mkdir(loops_folder)

    loops_scores = {}

    if len(target_loops) > 0:
        for i, loop in enumerate(target_loops):            
            loopID = '{}_{}-{}'.format(target.split('-')[0], loop[0], loop[-1])

            # extract target loop structure
            target_loop_structure = extract_loop(target_structure, loop, out_label = i, out_folder = loops_folder)
            # extract modelled loop structure
            model_loop_structure = extract_loop(curr_model, loop, out_label = i, out_folder = loops_folder)
            loop_scores = run_tristan_geometry_scores(target_loop_structure, model_loop_structure)
            loops_scores[loopID] = loop_scores

    return loops_scores

# 9. For template analysis

def parse_tophits_hhsearch(target_hhsearch, min_prob, exclude = None, exclude_date = None, min_cov = min_cov):

    templates = {'HitNo': [],
                 'pdb': [],
                 'pdb_release_date': [],
                 'pdb_method': [],
                 'prob': [],
                 'Eval': [],
                 'Query_interval': [],
                 'Template_interval': [],
                 'Target_coverage': [],
                 'Template_score': []}

    # get the table of hits first
    found_table = False
    len_target = np.nan
    if os.path.isfile(target_hhsearch):
        with open(target_hhsearch, 'r') as hhr:
            for line in hhr:
                if len(line) > 1:
                    if 'No' in line and 'Hit' in line and 'Prob' in line and not found_table:
                        found_table = True
                    elif found_table and 'No' in line:
                        found_table = False
                        break

                    elif not found_table:
                        if 'Match_columns' in line:
                            len_target = int(line.split()[1])

                    elif found_table:
                        hitno = line.split()[0].strip()
                        hit = line.split()[1].strip()
                        prob = float(line[35:41])
                        e_val = float(line[41:49])
                        query_interval = [int(i) for i in line[74:84].split('-')]
                        template_interval = [int(i) for i in line[84:93].split('-')]
                        target_coverage = (query_interval[1] - query_interval[0])*100/len_target

                        if prob >= min_prob and (exclude is None or exclude not in hit):

                            pdb_data = json.loads(p.PDB.getSummary(pdbid=hit.lower().split('_')[0]))
                            pdb_release_date = pdb_data[hit.lower().split('_')[0]][0]['release_date']
                            pdb_method = pdb_data[hit.lower().split('_')[0]][0]['experimental_method_class'][0]

                            if pdb_method == 'x-ray' and target_coverage >= min_cov and (exclude_date is None or exclude_date > pdb_release_date) and hit not in templates['pdb']:
                                templates['HitNo'].append(hitno)
                                templates['pdb'].append(hit)
                                templates['pdb_release_date'].append(pdb_release_date)
                                templates['pdb_method'].append(pdb_method)
                                templates['prob'].append(prob)
                                templates['Eval'].append(e_val)
                                templates['Query_interval'].append(query_interval)
                                templates['Template_interval'].append(template_interval)
                                templates['Target_coverage'].append(target_coverage)
                                templates['Template_score'].append((target_coverage+prob)/2.0)

    templates = pd.DataFrame(templates)
    templates['Alignment'] = ['' for i in range(len(templates))]

    if len(templates) > 0:
        # now get the sequence alignments
        found_alignment = False
        curr_hitno = None
        curr_hit = None
        query_seq = ''
        template_seq = ''
        template_dssp = ''

        with open(target_hhsearch, 'r') as hhr:
            for line in hhr:
                if len(line) > 1:
                    if 'No' in line and 'Hit' in line and 'Prob' in line and not found_table:
                        found_table = True
                    elif found_table and 'No ' in line and len(line.split())==2:
                        found_alignment = True
                        if query_seq != '':
                            # correct for absent modelled structure
                            if set(template_dssp) == set('-') or len(template_dssp) == 0:
                                templates.loc[templates['HitNo'] == curr_hitno, ['Alignment']] = np.nan
                            else:
                                templates.loc[templates['HitNo'] == curr_hitno, ['Alignment']] = '{}\n{}'.format(query_seq, template_seq)

                        curr_hitno = line.split()[1]

                    elif found_alignment and curr_hitno in list(templates['HitNo']):
                        if line.startswith('>'):
                            curr_hit = line.split()[0].replace('>','')
                            query_seq = ''
                            template_seq = ''
                            template_dssp = ''

                        elif line.startswith('Q') and 'Consensus' not in line:
                            query_seq += line[22:].split()[0].strip()

                        elif line.startswith('T') and curr_hit in line:
                            template_seq += line[22:].split()[0].strip()

                        elif line.startswith('T') and 'ss_dssp' in line:
                            template_dssp += line[22:].strip()

    templates = templates[templates['Alignment'].notna()]

    return templates

def get_target_templates(target, templates_folder, deposited_pdb, max_date, n = 5, min_cov = 25, min_prob = 70, copy_to = mr_folder):

    print(" ... Finding top {} templates at a minimum probability of {}%".format(n, min_prob))

    output_folder = '{}/structures'.format(templates_folder)
    if not os.path.isdir(output_folder):
        os.mkdir(output_folder)

    target_hhsearch = '{}/HHsearches/{}.hhr'.format(templates_folder, target)
    
    templates = parse_tophits_hhsearch(target_hhsearch, min_prob, exclude = deposited_pdb, exclude_date = max_date)

    if len(templates) > 0:
        templates.to_csv('{}/{}_templates_info.csv'.format(mr_folder, target))

        top_templates = templates.head(n).copy()
        top_truncated = []

        for index, row in top_templates.iterrows():
            curr_template = row.pdb
            curr_pdb = curr_template.split('_')[0]
            curr_chain = curr_template.split('_')[1]

            target_seq = row.Alignment.split('\n')[0]
            template_seq = row.Alignment.split('\n')[1]

            curr_template = download_pdb(curr_pdb, curr_chain, out_folder = output_folder)
            curr_template = fix_pdbs([curr_template])[0]
            
            if curr_template is not None:
            
                truncated_template = truncate_pdb(curr_template, template_seq.replace('-',''), reset_resnum = True, label = target)
                top_truncated.append(truncated_template)

                #now prepare the template for MR if we are dealing with CASP14
                if which_casp == 'CASP14':
                    alignment = '{}_{}_alignment.ali'.format(truncated_template, target)
                    with open(alignment, 'w') as alig:

                        if target_seq.replace('-', '') == template_seq.replace('-', ''):
                            sculpt_mode = 'same'
                        else:
                            sculpt_mode = 'diff'

                        alig.write('>{}\n{}\n>{}\n{}\n'.format(target, target_seq, truncated_template, template_seq))

                    sculpted_template = run_sculptor(truncated_template, alignment, mode = sculpt_mode, n = top_n_templates, min_cov = min_cov)

            else:
                top_truncated.append(np.nan)

        top_templates['Truncated_pdb'] = top_truncated
        top_templates = top_templates[top_templates['Truncated_pdb'].notna()]

        return top_templates

    else:
        return None

def run_sculptor(pdb, alignment, n, min_cov, mode='diff'):

    phil_file, sculptor_out = write_phil(pdb, alignment, n, min_cov)

    out_pdb = '{}/sculpt_{}'.format(sculptor_out, pdb.split('/')[-1])

    if not os.path.isfile(out_pdb):
        if mode == 'diff':
            run_sculptor = sp.Popen([sculptor, phil_file], stdout=sp.PIPE, stderr=sp.PIPE, stdin=sp.PIPE)
            stdout, stderr = run_sculptor.communicate()
        
        else:
            shutil.copyfile(pdb, out_pdb)


    return out_pdb

def write_phil(pdb, alignment, n, min_cov, mr_folder = mr_folder):

    template_phil = '{}/SCULPTOR_phil_template.txt'.format(mr_folder)

    sculptor_out = '{}/SCULPTOR_templates_minCov_{}_top{}'.format(mr_folder, min_cov, n)
    if not os.path.isdir(sculptor_out):
        os.mkdir(sculptor_out)

    out_phil = '{}_sculptor_phil'.format(pdb)

    if os.path.isfile(template_phil):
        with open(out_phil, 'w') as phil:
            found_model_field = False
            found_align_field = False
            found_outpt_foldr = False

            with open(template_phil, 'r') as tmp_phil:
                for line in tmp_phil:
                    if 'model ' in line:
                        found_model_field = True

                    elif found_model_field and 'file_name = ' in line:
                        line = line.replace('None', pdb)
                        found_model_field = False

                    elif 'alignment ' in line:
                        found_align_field = True

                    elif found_align_field and 'file_name = ' in line:
                        line = line.replace('None', alignment)
                        found_align_field = False

                    elif 'output ' in line:
                        found_outpt_foldr = True

                    elif found_outpt_foldr and 'folder = ' in line:
                        line = line.replace('.', sculptor_out)
                        found_outpt_foldr = False

                    phil.write(line)

    return out_phil, sculptor_out

def run_template_gdt(target, model):

    infile, wrkdir = generate_lga_input(target, model)

    outfile = '{}/{}_aligned.gdt'.format(wrkdir, infile)

    if not os.path.isfile(outfile):
    
        run_gdt = sp.Popen([templates_gdt, wrkdir, infile], stdout=sp.PIPE, stderr=sp.PIPE, stdin=sp.PIPE)
        stdout, stderr = run_gdt.communicate()

        with open(outfile, 'w') as outp:
            outp.write(stdout.decode("utf-8"))
    
    gdt_out = parse_predictions_gdt(outfile)

    return gdt_out

# 10. For CaBLAM analysis

def parse_molprobity_cablam(molpout):

    results = {}
    
    found_table = False

    with open(molpout, 'r') as molp:
        for line in molp:
            if 'residue : outlier_type' in line:
                found_table = True
            elif 'SUMMARY:' in line:
                found_table = False
            
            elif found_table:
                res = line.split(':')[0]
                res = '{}{}'.format(res.split()[-1].strip(), res.split()[-2].strip())
                contour = line.split(':')[2]
                try:
                    results[res] = float(contour) 
                except:
                    pass
    
    return results

def run_cablam(pdb):

    prepared_target = prepare_pdb_for_molprobity(pdb)

    cablam_out = '{}_molprobity_cablam.txt'.format(prepared_target)

    if not os.path.isfile(cablam_out):
        run_molprobity = sp.Popen([molprobity_cablam, prepared_target], stdout=sp.PIPE, stderr=sp.PIPE, stdin=sp.PIPE)
        stdout, stderr = run_molprobity.communicate()

        with open(cablam_out, 'w') as ramaout:
            ramaout.write(stdout.decode("utf-8"))

    cablam = parse_molprobity_cablam(cablam_out)

    return cablam


def run_CaBlamDiff(target, model, target_interval = None, mode = None, make_plot = True):

    if model != None:
    
        target_cablam_out = run_cablam(target)
        model_cablam_out = run_cablam(model)
        
        if target_interval == None:
            residues_in_target = [key for key in target_cablam_out.keys()]
            residues_in_model = [key for key in model_cablam_out.keys()]

        else:
            residues_in_target = [key for key in target_cablam_out.keys() if key in target_interval]
            residues_in_model = [key for key in model_cablam_out.keys() if key in target_interval]

            chidiff = np.nan

        if mode == None:
            modelled_residues = [key for key in residues_in_target if key in residues_in_model]
        elif mode == 'template':
            modelled_residues = align_structures(target, model)

        dipdifferences = []
        accepted_residues = []
        target_scores = []
        model_scores = []
        for residue in modelled_residues:
            if type(modelled_residues) is dict:
                model_residue = modelled_residues[residue]
            else:
                model_residue = residue

            if residue in target_cablam_out and model_residue in model_cablam_out:
                dipdifferences.append(model_cablam_out[model_residue] - target_cablam_out[residue])
                accepted_residues.append(residue)
                target_scores.append(target_cablam_out[residue])
                model_scores.append(model_cablam_out[model_residue])        

        if make_plot and len(accepted_residues) > 0:

            out_label = '{}_{}_CaBLAM'.format(target.replace('.pdb', ''), model.split('/')[-1].replace('.pdb', ''))

            # plot dipscores
            plt.clf()
            fig = plt.figure(figsize = (10, 2))
            plt.plot(accepted_residues, target_scores, c = 'lightblue', label = '{}'.format(target.split('/')[-1]))
            plt.plot(accepted_residues, model_scores, c = 'steelblue', label = '{}'.format(model.split('/')[-1]))
            plt.xticks([], [])
            plt.ylabel('CaBLAM contour/score')
            plt.xlabel('Residues')
            plt.savefig('{}_CaBLAM_scores.png'.format(out_label))
            plt.close()
            
            # plot dipscore difference
            plt.clf()
            fig = plt.figure(figsize = (10, 2))
            plt.plot(accepted_residues, dipdifferences, c = 'steelblue', label = '{}'.format(model.split('/')[-1]))
            plt.hlines(y=0, xmin=accepted_residues[0], xmax=accepted_residues[-1], linestyle = ':')
            plt.xticks([], [])
            plt.ylim(-1.1,1.1)
            plt.ylabel('CaBLAM difference')
            plt.xlabel('Residues')
            plt.savefig('{}_CaBLAM_diff.png'.format(out_label))
            plt.close()
            
            # plot histogram of dipscore differences
            plt.clf()
            fig = plt.figure(figsize = (4, 3))
            plt.hist(dipdifferences, color = 'steelblue', bins = 100 )
            plt.xlabel('CaBLAM difference')
            plt.ylabel('Count')
            plt.vlines(x=np.mean([i for i in dipdifferences if not np.isnan(i)]), ymin = 0, ymax = 50, color = 'red')
            plt.title('Mean: {}'.format(np.mean([i for i in dipdifferences if not np.isnan(i)])))
            plt.xlim(-1.1, 1.1)
            plt.savefig('{}_CaBLAM_diff_hist.png'.format(out_label))
            plt.close()

        dipdiffs = {'CaBLAMDiff': round(np.mean([i for i in dipdifferences if not np.isnan(i)]), 5)}

    else:
        dipdiffs = {'CaBLAMDiff': np.nan}

    return dipdiffs

def compare_CaBlam_to_DipScore(target, model, target_interval = None, mode = None, make_plot = True):

    if model != None:
    
        target_cablam_out = run_cablam(target)
        model_cablam_out = run_cablam(model)

        target_dipcheck_out, target_chiscore = run_dipcheck(target)
        model_dipcheck_out, model_chiscore = run_dipcheck(model)

        if target_interval == None:
            residues_in_target = [key for key in target_cablam_out.keys()]
            residues_in_model = [key for key in model_cablam_out.keys()]

        else:
            residues_in_target = [key for key in target_cablam_out.keys() if key in target_interval]
            residues_in_model = [key for key in model_cablam_out.keys() if key in target_interval]

            chidiff = np.nan

        if mode == None:
            modelled_residues = [key for key in residues_in_target if key in residues_in_model]
        elif mode == 'template':
            modelled_residues = align_structures(target, model)

        accepted_residues = []
        dipdifferences = []
        cablam_differences = []

        for residue in modelled_residues:
            if type(modelled_residues) is dict:
                model_residue = modelled_residues[residue]
            else:
                model_residue = residue

            if residue in target_cablam_out and model_residue in model_cablam_out and residue in target_dipcheck_out and model_residue in model_dipcheck_out:
                cablam_differences.append(model_cablam_out[model_residue] - target_cablam_out[residue])
                dipdifferences.append(model_dipcheck_out[model_residue] - target_dipcheck_out[residue])

                accepted_residues.append(residue)

        if len(cablam_differences) >= 2:
            correlation = pearsonr(cablam_differences,dipdifferences)[0]
        else:
            correlation = np.nan

        if make_plot and len(accepted_residues) > 0:

            out_label = '{}_{}_CaBLAM_vs_DipScore'.format(target.replace('.pdb', ''), model.split('/')[-1].replace('.pdb', ''))

            # plot dipscores vs CaBLAM            
            plt.clf()
            fig = plt.figure(figsize = (4, 4))
            plt.scatter(cablam_differences, dipdifferences, c = 'steelblue', label = '{}'.format(model.split('/')[-1]))
            plt.hlines(y=0, xmin=-1.1, xmax=1.1, linestyle = ':')
            plt.vlines(x=0, ymin=-1.1, ymax=1.1, linestyle = ':')
            plt.ylim(-1.1,1.1)
            plt.xlim(-1.1,1.1)
            plt.ylabel('DipScore difference')
            plt.xlabel('CaBLAM difference')
            plt.title('Pearson corr: {}'.format(correlation))
            plt.tight_layout()
            plt.savefig('{}_correlation.png'.format(out_label))
            plt.close()

        dipdiffs = {'DipCaBlam_corr': correlation}

    else:
        dipdiffs = {'DipCaBlam_corr': np.nan}

    return dipdiffs

# MAIN FUNCTION FOR MULTITHREADING

def analyse_tables(arguments):

    job                    = arguments[0]
    tables                 = arguments[1]
    templates_folder       = arguments[2]
    targets_classification = arguments[3]
    targets_deposited_pdbs = arguments[4]
    targets_neff           = arguments[5]
    out_folder             = arguments[6]

    for j, table in enumerate(tables):
        target = table.split('/')[-1].split('.txt')[0]
        target_structure = '{}/{}.pdb'.format(targets_folder, target)

        jobID = 'Job.{}: {} ({} of {})'.format(job, target, j+1, len(tables))

        print('\n{}'.format(jobID))

        #if target == 'T1092-D1':

        if not os.path.isfile(target_structure) and os.path.isdir('{}/Targets_submitted'.format(targets_folder)):
            target_submitted = ['{}/Targets_submitted/{}'.format(targets_folder, f) for f in os.listdir('{}/Targets_submitted'.format(targets_folder)) if target in f and f.endswith('.pdb')]
            if len(target_submitted) > 0:
                target_structure = target_submitted[0]

        if os.path.isfile(target_structure):

            target_start = time.time()

            try:
                target_deposited_pdb = targets_deposited_pdbs[target.split('-')[0]]['PDB_code']
                target_title = targets_deposited_pdbs[target.split('-')[0]]['Title']
                target_status = targets_deposited_pdbs[target.split('-')[0]]['Status']
                target_exp_date = targets_deposited_pdbs[target.split('-')[0]]['Exp.Date']
            except:
                target_deposited_pdb = targets_deposited_pdbs[target]['PDB_code']
                target_title = targets_deposited_pdbs[target]['Title']
                target_status = targets_deposited_pdbs[target]['Status']
                target_exp_date = targets_deposited_pdbs[target]['Exp.Date']                

            if target_status == 'Active':

                # parse casp table
                target_table = parse_table(table)
                # download all its models
                models_folder, full_models_folder = download_all_models(target, out_folder = '{}/models'.format(tmp_folder))
                # compute secondary structure and get residue numbering
                target_dssp = run_dssp(target_structure, out_label = target)
                # calculate target difficulty
                target_difficulty, target_classification = compute_difficulty(target, templates_folder, target_deposited_pdb, target_exp_date, res_interval = [int(target_dssp['ResNum'][0]), int(target_dssp['ResNum'][-1])], targets_classification = targets_classification)
                # get target Neff_max
                if target in targets_neff:
                    target_neff = targets_neff[target]
                else:
                    target_neff = np.nan

                # if is not a refinement case, get the target top N templates
                if not target.startswith('R'):
                    templates_table = get_target_templates(target, templates_folder, target_deposited_pdb, target_exp_date, n = top_n_templates, min_cov = min_cov, min_prob = 70, copy_to = mr_folder)
                else:
                    templates_table = None

                if templates_table is not None and len(templates_table) > 0:

                    if not os.path.isfile('{}/{}_templates_updated_data.csv'.format(out_folder, target)):
                    
                        # Analyse the templates
                        print(' ... Analysing templates')
                        print(' ... ... There are {} templates'.format(len(templates_table)))

                        # add new columns to the table where new scores will be stored
                        template_count = 0
                        templates_table['Target_difficulty']     = [target_difficulty for i in range(len(templates_table))]
                        templates_table['Target_Neff']           = [target_neff for i in range(len(templates_table))]
                        templates_table['Target_classification'] = [target_classification for i in range(len(templates_table))]
                        templates_table['Target_title']       = [target_title for i in range(len(templates_table))]
                        templates_table['Target_PDB_code']    = [target_deposited_pdb for i in range(len(templates_table))]
                        templates_table['Template_DipDiff']      = [np.nan for i in range(len(templates_table))]
                        templates_table['Template_ChiDiff']      = [np.nan for i in range(len(templates_table))]
                        templates_table['Template_BBscore']      = [np.nan for i in range(len(templates_table))]
                        templates_table['Template_SCscore']      = [np.nan for i in range(len(templates_table))]
                        templates_table['Template_CaBLAMdiff'] = [np.nan for i in range(len(templates_table))]
                        templates_table['Template_DipCaBLAM_corr'] = [np.nan for i in range(len(templates_table))]
                        templates_table['Template_GDT_TS']       = [np.nan for i in range(len(templates_table))]
                        templates_table['Template_GDT_HA']       = [np.nan for i in range(len(templates_table))]
                        templates_table['Template_GDC_SC']       = [np.nan for i in range(len(templates_table))]

                        for index, row in templates_table.iterrows():
                            template_count += 1

                            curr_template = row.Truncated_pdb
                            #curr_template = '{}/structures/{}_truncated_to_{}.pdb'.format(templates_folder, row.pdb, target)

                            if os.path.isfile(curr_template):

                                print(' ... ... {}.{} | {}.{} Checking template {}'.format(job, len(tables), template_count, len(templates_table), row.pdb))

                                # compute dipdiff of model
                                template_dipdiff = run_dipdiff(target_structure, curr_template,  mode = 'template')
                                # now check Tristan's and Randy's geometry scores
                                template_geometry = run_tristan_geometry_scores(target_structure, curr_template,  mode = 'template')
                                # now check the CaBLAM differences
                                template_cablam = run_CaBlamDiff(target_structure, curr_template,  mode = 'template', make_plot = True)
                                # now compare the CaBLAM and the DipScores
                                template_dipcablam_corr = compare_CaBlam_to_DipScore(target_structure, curr_template,  mode = 'template', make_plot = True)
                                # compute GDT scores
                                template_gdt = run_template_gdt(target_structure, curr_template)

                                # update table
                                templates_table.at[index, 'Template_DipDiff']    = template_dipdiff['DipDiff']
                                templates_table.at[index, 'Template_ChiDiff']    = template_dipdiff['ChiDiff']
                                templates_table.at[index, 'Template_BBscore']    = template_geometry['BBscore']
                                templates_table.at[index, 'Template_SCscore']    = template_geometry['SCscore']
                                templates_table.at[index, 'Template_CaBLAMdiff'] = template_cablam['CaBLAMDiff']
                                templates_table.at[index, 'Template_DipCaBLAM_corr'] = template_dipcablam_corr['DipCaBlam_corr']

                                for score in template_gdt:
                                    templates_table.at[index, 'Template_{}'.format(score)] = template_gdt[score]

                            else:
                                print(" ... ... {}.{} | {}.{} {} ATTENTION: not found!".format(job, len(tables), template_count, len(templates_table), curr_template))

                            #save tables
                            templates_table.to_csv('{}/{}_templates_updated_data.csv'.format(out_folder, target), sep = '\t')
                
                    else:
                        print(' ... Already analysed templates for {}'.format(target))
                else:
                    print(' ... Found no templates for {}'.format(target))


                # Now analyse the models
                print(' ... Analysing models')
                print(' ... ... There are {} models'.format(len(target_table)))

                # add new columns to the table where new scores will be stored
                target_table['Target_difficulty']     = [target_difficulty for i in range(len(target_table))]
                target_table['Target_Neff']           = [target_neff for i in range(len(target_table))]
                target_table['Target_classification'] = [target_classification for i in range(len(target_table))]
                target_table['Target_title']       = [target_title for i in range(len(target_table))]
                target_table['Target_PDB_code']    = [target_deposited_pdb for i in range(len(target_table))]
                target_table['Model_DipDiff']         = [np.nan for i in range(len(target_table))]
                target_table['Model_ChiDiff']         = [np.nan for i in range(len(target_table))]
                target_table['Model_CaBLAMdiff']     = [np.nan for i in range(len(target_table))]
                target_table['Model_DipCaBLAM_corr'] = [np.nan for i in range(len(target_table))]
                target_table['Model_BBscore']         = [np.nan for i in range(len(target_table))]
                target_table['Model_SCscore']         = [np.nan for i in range(len(target_table))]

                model_count = 0
                loop_count = 0
                for index, row in target_table.iterrows():

                    model_start = time.time()

                    model_count += 1

                    curr_model = '{}/{}'.format(models_folder, row.Model)

                    if not os.path.isfile(curr_model) and '-D' in target:
                        curr_model = '{}/{}'.format(full_models_folder,row.Model.split('-D')[0])

                    if os.path.isfile(curr_model):

                        print(' ... ... {}.{} | {}.{} Checking model'.format(job, len(tables), model_count, len(target_table)))
                        # compute dipdiff of model
                        model_dipdiff = run_dipdiff(target_structure, curr_model, make_plot = True)
                        # now check Tristan's and Randy's geometry scores
                        model_geometry = run_tristan_geometry_scores(target_structure, curr_model)
                        # now check the CaBLAM differences
                        model_cablam = run_CaBlamDiff(target_structure, curr_model, make_plot = True)
                        # now compare the CaBLAM and the DipScores
                        model_dipcablam_corr = compare_CaBlam_to_DipScore(target_structure, curr_model, make_plot = True)

                        # update table
                        target_table.at[index, 'Model_DipDiff']    = model_dipdiff['DipDiff']
                        target_table.at[index, 'Model_ChiDiff']    = model_dipdiff['ChiDiff']
                        target_table.at[index, 'Model_BBscore']    = model_geometry['BBscore']
                        target_table.at[index, 'Model_SCscore']    = model_geometry['SCscore']
                        target_table.at[index, 'Model_CaBLAMdiff']    = model_cablam['CaBLAMDiff']
                        target_table.at[index, 'Model_DipCaBLAM_corr'] = model_dipcablam_corr['DipCaBlam_corr']

                        # check how long it took
                        model_end = time.time()
                        numb_seconds = model_end - model_start
                        print(" ... ... {}.{} | {}.{} {} Took {} days {}".format(job, len(tables), model_count, len(target_table), row.Model, int(time.strftime('%d', time.gmtime(numb_seconds)))-1, time.strftime('%H hours %M min %S sec', time.gmtime(numb_seconds))))

                        # save tables
                        target_table.to_csv('{}/{}_predictions_updated_data.csv'.format(out_folder, target), sep = '\t')
                    
                    else:
                        print(" ... ... {}.{} | {}.{} {} ATTENTION: not found!".format(job, len(tables), model_count, len(target_table), row.Model))


                 # check how long it took
                target_end = time.time()
                numb_seconds = target_end - target_start
                print(" ... {} Took {} days {}".format(target, int(time.strftime('%d', time.gmtime(numb_seconds)))-1, time.strftime('%H hours %M min %S sec', time.gmtime(numb_seconds))))

            else:
                print(' -- ATTENTION: target {} was cancelled --> skipping it!'.format(target))
        else:
            print(' -- ATTENTION: no structure was found for {} --> skipping it!'.format(target))

    return out_folder

# MAIN CODE

tables                 = download_tables(output_folder = '{}/tables'.format(tmp_folder))
targets_folder         = download_targets(output_folder = '{}/targets'.format(tmp_folder))
templates_folder       = download_templates_information(output_folder = '{}/templates'.format(tmp_folder))
targets_classification = parse_targets_classification() 
targets_deposited_pdbs = parse_targets_pdbs()
targets_neff           = parse_templates_Neff()

output_folder          = '{}/out_tables'.format(tmp_folder)

if not os.path.isdir(output_folder):
    os.mkdir(output_folder)

unique_targets = sorted(list(set([table.split('-')[0].split('/')[-1].split('.txt')[0] for table in tables])))

# setup and run threads
separated_jobs = chunk_list(tables, n_threads, unique_targets)
list_arguments = [i for i in zip(range(n_threads), separated_jobs, [templates_folder for job in separated_jobs], [targets_classification for job in separated_jobs], [targets_deposited_pdbs for job in separated_jobs], [targets_neff for job in separated_jobs], [output_folder for job in separated_jobs])]

pool = mp.Pool(n_threads)
collected_data = pool.map(analyse_tables, list_arguments)
