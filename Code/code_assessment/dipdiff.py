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
import matplotlib.colors as colors

from matplotlib import rc

rc('text', usetex=False)
rc('font', family='Arial')

# Get inputs
parser = argparse.ArgumentParser('DipDiff calculator. If running with a template, you have to set the memory limit to unlimited due to LGA.\nEXAMPLE USAGE: ulimit -s unlimited; python3 dipdiff.py -target xxx -model xxx')
requiredNamed = parser.add_argument_group('required arguments')
optionalNamed = parser.add_argument_group('optional arguments with defaults')

# required inputs
requiredNamed.add_argument('-target', dest='target', type=str, required=True, help='The target structure the model is to be compared to')
requiredNamed.add_argument('-model', dest='model', type=str, required=True, help='The mode structure to be evaluated')
# optional inputs
optionalNamed.add_argument('-mode', dest='mode', type=str, default = None, help='The mode to be run. If the model is a template, use mode "template" (default: None)')

# Define inputs
args = parser.parse_args()

# tools needed 
dipcheck = 'dipcheck'
templates_lga = 'templates_lga.sh'


# 1. Helping routines

def set_occ(pdb, occ = 1.0, out_folder = None):

    run_pdbtools = sp.Popen(['pdb_occ', '-{}'.format(occ), pdb], stdout=sp.PIPE, stderr=sp.PIPE, stdin=sp.PIPE)
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

# 2. For structure alignment 

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
    
    if len(target.split('/')) == 1:
    	out_folder = 'dipcheck'
    else:
    	out_folder = '{}/dipcheck'.format('/'.join(target.split('/')[:-1]))
    #target = set_bfactor(target, out_folder = out_folder)
    target = set_occ(target, out_folder = out_folder)

    out_dipcheck = '{}.dipout'.format(target)
    out_pdb = '{}.dipout.pdb'.format(target)

    if not os.path.isfile(out_dipcheck):

        run_dipcheck = sp.Popen([dipcheck, 'xyzin', target, 'xyzout', out_pdb], stdout=sp.PIPE, stderr=sp.PIPE, stdin=sp.PIPE)
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


        if len(accepted_residues) > 0:

	        if make_plot:

	            out_label = '{}_{}_DipDiff'.format(target.replace('.pdb', ''), model.split('/')[-1].replace('.pdb', ''))

	            # plot dipscores
	            plt.clf()
	            fig = plt.figure(figsize = (10, 2))
	            # plt.plot(accepted_residues, target_scores, c = 'lightblue', label = '{}'.format(target.split('/')[-1].strip('.pdb')))
	            # plt.plot(accepted_residues, model_scores, c = 'steelblue', label = '{}'.format(model.split('/')[-1].strip('.pdb')))
	            gradient = np.linspace(0, 1, 100).reshape(-1, 1)
	            plt.imshow(gradient , extent=[-1, len(accepted_residues), -0.1,1.1], aspect='auto', cmap=colors.LinearSegmentedColormap.from_list('wr',["w", "r"], N=256), alpha=0.2)
	            plt.plot(accepted_residues, target_scores, c = 'grey', label = 'Target')
	            plt.plot(accepted_residues, model_scores, c = 'black', label = 'Model')
	            plt.xticks([], [])
	            plt.ylim(-0.1,1.1)
	            plt.ylabel('DipScore')
	            plt.xlabel('Residues')
	            plt.tight_layout()
	            plt.legend()
	            plt.savefig('{}_DipCheck_scores.pdf'.format(out_label))
	            plt.close()
	            
	            # plot dipscore difference
	            plt.clf()
	            fig = plt.figure(figsize = (10, 2))
	            gradient = np.linspace(0, 1, 100).reshape(-1, 1)
	            plt.imshow(gradient , extent=[-1, len(accepted_residues), -1.1,1.1], aspect='auto', cmap=colors.LinearSegmentedColormap.from_list('bwr',["b","w", "r"], N=256), alpha=0.2)
	            plt.plot(accepted_residues, dipdifferences, c = 'dimgrey', label = '{}'.format(model.split('/')[-1]))
	            plt.hlines(y=0, xmin=accepted_residues[0], xmax=accepted_residues[-1], linestyle = ':')
	            plt.xticks([], [])
	            plt.ylim(-1.1,1.1)
	            plt.ylabel('DipScore difference')
	            plt.xlabel('Residues')
	            plt.tight_layout()
	            plt.savefig('{}_DipScore_diff.pdf'.format(out_label))
	            plt.close()
	            
	            # plot histogram of dipscore differences
	            plt.clf()
	            fig = plt.figure(figsize = (3, 4))
	            h = plt.hist(dipdifferences, color = 'dimgrey', bins = 100)
	            plt.xlabel('DipScore difference')
	            plt.ylabel('Count')
	            plt.vlines(x=np.mean([i for i in dipdifferences if not np.isnan(i)]), ymin = 0, ymax = 50, color = 'k', linestyle=':')
	            plt.title('Mean/DipDiff: {}'.format(round(np.mean([i for i in dipdifferences if not np.isnan(i)]), 3)))
	            plt.ylim(0, max(h[0])+0.1*max(h[0]))
	            plt.xlim(-1.1, 1.1)
	            plt.tight_layout()
	            plt.savefig('{}_DipScore_diff_hist.pdf'.format(out_label))
	            plt.close()


        dipdiffs = {'DipDiff': round(np.mean([i for i in dipdifferences if not np.isnan(i)]), 3), 'ChiDiff': chidiff}

    else:
        dipdiffs = {'DipDiff': np.nan, 'ChiDiff': np.nan}

    return dipdiffs


if __name__ == '__main__':

	main_start = time.time()

	print('INPUTS:\n- Target: {}\n- Model : {}\n'.format(args.target, args.model))

	dipdiffs = run_dipdiff(args.target, args.model, mode = args.mode, make_plot = True)
	print(' ... Model DipDiff:', dipdiffs['DipDiff'])

	main_end = time.time()
	main_numb_seconds = main_end - main_start
	print("\n#### Finished job after: {} \n".format(time.strftime('%H hours %M min %S sec', time.gmtime(main_numb_seconds))))

