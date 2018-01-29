#! /usr/bin/env python
# -*- coding: utf-8 -*

##########################################################################
# @author: yann.leguen@cea.fr
# Copyright (C) CEA, 2017
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html for details.
##########################################################################

import os
import json
import pandas as pd
import numpy as np
from scipy.stats import pearsonr
import nibabel.gifti.giftiio as gio

import shutil
from heritability_analysis import create_solar_pedigree
from heritability_analysis import convert_pheno_to_solar


def compute_correlation_with_scipy(indir_cog,
                                   columns, areals_org, phen_dir,
                                   areals_name, sd, fil, outdir):
    dict_rhop = {}
    dict_rhop_pval = {}
    outdir2 = os.path.join(outdir, 'tmp_phenotypes_csv')
    if not os.path.isdir(outdir2):
        os.makedirs(outdir2)
    output = os.path.join(outdir, sd)

    # Loop through all selected colums of cognitive traits
    for i, col in enumerate(columns): 
        fname = os.path.join(indir_cog, col+'.csv')
        df_cog = pd.read_csv(fname)
        df_cog.index = df_cog['IID'] 
        # Loop through all areals
        for k, areal in enumerate(areals_org):
            fname2 = os.path.join(phen_dir, areals_name[k]+'_'+fil+'.csv')
            df_beta = pd.read_csv(fname2)
            df_beta.index = df_beta['IID']
            df_beta = df_beta.sort_index()
            if k == 0:
                # Align the index of cognitive trait with the ones of the beta
                # only need to do it once because all areals have same subjects
                df_cog = df_cog.loc[df_beta.index]
                df_cog = df_cog.sort_index()
                df_cog = df_cog.dropna()
            df_beta = df_beta.loc[df_cog.index]
            # Compute correlation between activation and cognitive trait
            corr, p = pearsonr(df_cog[col], df_beta[fil])
            dict_rhop[str(areal)+'_'+col] = corr
            dict_rhop_pval[str(areal)+'_'+col] = p
            # Create a csv for the bivariate genetic analysis
            df = pd.DataFrame()
            df = df_cog
            df[fil] = df_beta[fil]
            output2 = os.path.join(outdir2, str(areal)+'_'+col+'.csv')
            df.to_csv(output2, header=True, index=False)

    encoded = json.dumps(dict_rhop)
    with open(output+'scipy_rhop_dict.json', 'w') as f:
        json.dump(encoded, f)

    encoded = json.dumps(dict_rhop_pval)
    with open(output+'scipy_rhop_pval_dict.json', 'w') as f:
        json.dump(encoded, f)

def run_solar_bivariate(work_dir, pheno, trait1, trait1_, trait2, outdir):
    """
    Parameters
    work_dir: SOLAR working directory containing pedigree
    pheno: phenotype file containing the two traits
    trait1: areal number
    trait1_ : name of first trait in the dataframe
    
    trait2: name second trait
    outdir: output directory for SOLAR analysis
    """

    outdir = os.path.join(outdir, trait1, trait2)
    if not os.path.isdir(outdir):
        os.makedirs(outdir)
    df_phen = convert_pheno_to_solar(pheno, False, restrict_to_white=True)
    pheno = os.path.join(os.path.dirname(pheno), "solar_ready",
                         os.path.basename(pheno))
    cmd = " ".join(["solar bivariate_analysis", work_dir,
                    outdir, trait1_, trait2, pheno])
                    
    print cmd
    # Note: The code for Group covariates is not self adaptive yet
    # i.e. need to check pheno_analysis.tcl file to adjust accordingly
    os.system(cmd)

def preprocess_solar_bivariate(columns, areals_org, phen_dir,
                               fil, work_dir, outdir):

    if not os.path.isdir(outdir):
        os.makedirs(outdir)

    for i, col in enumerate(columns): 
        for k, areal in enumerate(areals_org):
            pheno = os.path.join(phen_dir, str(areal)+'_'+col+'.csv')
            run_solar_bivariate(work_dir, pheno, str(areal), fil, col, outdir)
            
def parse_bivariate_analysis(indir, columns, areals_org, sd, outdir):
    SANITY_CHECK = False
    verbose = False
    nb = ['.', '-', 'e', '0', '1', '2', '3', '4',
          '5', '6', '7', '8', '9', '10']
    output = os.path.join(outdir, sd)


    dict_rhoe = {}
    dict_rhoe_pval = {}
    dict_rhog = {}
    dict_rhog_std = {}
    dict_rhog_pval = {}
    dict_rhop = {}
    dict_rhop_pval = {}
    # Loop through all selected colums of cognitive traits
    for i, col in enumerate(columns):
        # Loop through all areals
        for k, areal in enumerate(areals_org):          
            file_path = os.path.join(indir, str(areal), col, 'polygenic.out')
            if os.path.isfile(file_path):
                for line in open(file_path, 'r'):
                    if 'RhoE is' in line:
                        if verbose:
                            print line[4:len(line)-1]
                        if len(line) > 30:
                            h2 = line[12:len(line)-16] 
                            p = line[28:len(line)-1] 
                        else:
                            h2 = line[12:len(line)-8] 
                            p = '1'
                        for z,l in enumerate(h2):
                            if not (l  in nb):
                                break
                        h2 = float(h2[:z])
                        p = float(p)
                        if SANITY_CHECK:
                            if p > 1 or p < 0 or abs(h2) > 1:
                                print line[4:len(line)-1]
                                print ("We extracted rhoe: "+str(h2)+
                                       " pval: "+str(p))
                        if verbose:
                            print "We extracted rhoe: "+str(h2)+" pval: "+str(p)
                        dict_rhoe[str(areal)+'_'+col] = h2
                        dict_rhoe_pval[str(areal)+'_'+col] = p
                    if 'RhoG is' in line:
                        if verbose:
                            print line[4:len(line)-1]
                        h2 = line[12:len(line)]                
                        for z,l in enumerate(h2):
                            if not (l  in nb):
                                break
                        h2 = float(h2[:z])
                        if verbose:
                            print "We extracted rhog: "+str(h2)
                        if SANITY_CHECK:
                            if  abs(h2) > 1:
                                print line[4:len(line)-1]
                                print "We extracted rhog: "+str(h2)
                        dict_rhog[str(areal)+'_'+col] = h2

                    if 'RhoG Std. Error:' in line:
                        if verbose:
                            print line[4:len(line)-1]
                        h2 = line[26:len(line)-1]
                        if verbose:
                            print "RHOG STD :"+str(h2)+"\n"
                        if h2 != 'Not Computable':
                            for z,l in enumerate(h2):
                                if not (l  in nb):
                                    break
                            h2 = float(h2[:z])
                            if verbose:
                                print "We extracted rhog std: "+str(h2)
                            if SANITY_CHECK:
                                if  abs(h2) > 1:
                                    print line[4:len(line)-1]
                                    print "We extracted rhog std: "+str(h2)
                            dict_rhog_std[str(areal)+'_'+col] = h2
                        else:
                            dict_rhog_std[str(areal)+'_'+col] = 1
                    if 'RhoG different from zero' in line:
                        if verbose:
                            print line[4:len(line)-1]
                        p = line[38:len(line)-1]                
                        p = float(p)
                        if verbose:
                            print "We extracted p-rhog: "+str(p)
                        if SANITY_CHECK:
                            if p > 1 or p < 0 :
                                print line[4:len(line)-1]
                                print "We extracted rhog pval: "+str(p)
                        dict_rhog_pval[str(areal)+'_'+col] = p 
                    if 'Derived Estimate of RhoP is' in line:
                        h2 = line[36:len(line)]                
                        for z,l in enumerate(h2):
                            if not (l  in nb):
                                break
                        h2 = float(h2[:z])
                        if verbose:
                            print line[4:len(line)-1]
                            print "We extracted rhop: "+str(h2)
                        if SANITY_CHECK:
                            if  abs(h2) > 1:
                                print line[4:len(line)-1]
                                print "We extracted rhog: "+str(h2)
                        dict_rhop[str(areal)+'_'+col] = h2
                    if 'RhoP different from zero' in line:
                        p = line[38:len(line)-1]                
                        p = float(p)
                        if verbose:
                            print line[4:len(line)-1]
                            print "We extracted p-rhop: "+str(p)
                        if SANITY_CHECK:
                            if p > 1 or p < 0 :
                                print line[4:len(line)-1]
                                print "We extracted rhog pval: "+str(p)
                        dict_rhop_pval[str(areal)+'_'+col] = p  
            else:
                if True:
                    print "Check why "+file_path+" doesn't exist !"
    encoded = json.dumps(dict_rhoe)
    with open(output+'rhoe_dict.json', 'w') as f:
        json.dump(encoded, f)

    encoded = json.dumps(dict_rhoe_pval)
    with open(output+'rhoe_pval_dict.json', 'w') as f:
        json.dump(encoded, f)
    encoded = json.dumps(dict_rhog)
    with open(output+'rhog_dict.json', 'w') as f:
        json.dump(encoded, f)

    encoded = json.dumps(dict_rhog_std)
    with open(output+'rhog_std_dict.json', 'w') as f:
        json.dump(encoded, f)
        
    encoded = json.dumps(dict_rhog_pval)
    with open(output+'rhog_pval_dict.json', 'w') as f:
        json.dump(encoded, f)
    encoded = json.dumps(dict_rhop)
    with open(output+'rhop_dict.json', 'w') as f:
        json.dump(encoded, f)

    encoded = json.dumps(dict_rhop_pval)
    with open(output+'rhop_pval_dict.json', 'w') as f:
        json.dump(encoded, f)

        
def table_bivariate_cog_traits(indir, areals_org, col,
                               areals_name, THR, output):
    verbose = False
    nb = ['.', '-', 'e', '0', '1', '2', '3',
          '4', '5', '6', '7', '8', '9', '10']

    elements = ['trait', 'rhop', 'rhop_pval',
                'rhog', 'rhog_std', 'rhog_pval',
                'rhoe', 'rhoe_pval']

    df = pd.DataFrame()
    dict_df = {}
    for elem in elements:
        dict_df[elem] = []
    
    with open(indir+'rhog_pval_dict.json', 'r') as f:
        data = json.load(f)
    dict_rhog_pval = json.loads(data)
    with open(indir+'rhog_dict.json', 'r') as f:
        data = json.load(f)
    dict_rhog = json.loads(data)
    with open(indir+'rhog_std_dict.json', 'r') as f:
        data = json.load(f)
    dict_rhog_std = json.loads(data)
    with open(indir+'rhoe_pval_dict.json', 'r') as f:
        data = json.load(f)
    dict_rhoe_pval = json.loads(data)
    with open(indir+'rhoe_dict.json', 'r') as f:
        data = json.load(f)
    dict_rhoe = json.loads(data)
    with open(indir+'scipy_rhop_dict.json', 'r') as f:
        data = json.load(f)
    dict_scipy_rhop = json.loads(data)
    with open(indir+'scipy_rhop_pval_dict.json', 'r') as f:
        data = json.load(f)
    dict_scipy_pval = json.loads(data)
    
    for k, areal in enumerate(areals_org):
        if dict_rhog.has_key(str(areal)+'_'+col):
            # and areal in list_common_areals[side]:
            if dict_scipy_pval[str(areal)+'_'+col] < THR:              
                dict_df['trait'].append(areals_name[k])
                dict_df['rhog'].append(round(dict_rhog[str(areal)+'_'+col],2))
                dict_df['rhog_std'].append(round(dict_rhog_std[str(areal)+
                                                               '_'+col],2))
                p = dict_rhog_pval[str(areal)+'_'+col]
                if p<0.01:
                    p_temp = '%.1e'% p
                    p_temp= p_temp.replace('e-0', '·10-')
                    p_temp = p_temp.replace('e-', '·10-')
                    dict_df['rhog_pval'].append(p_temp)
                else:
                    dict_df['rhog_pval'].append(str(round(p,2)))

                dict_df['rhop'].append(round(dict_scipy_rhop[str(areal)+
                                                             '_'+col],2))
                p = dict_scipy_pval[str(areal)+'_'+col]
                if p<0.01:
                    p_temp = '%.1e'% p
                    p_temp= p_temp.replace('e-0', '·10-')
                    p_temp = p_temp.replace('e-', '·10-')
                    dict_df['rhop_pval'].append(p_temp)
                else:
                    dict_df['rhop_pval'].append(str(round(p,2)))

                dict_df['rhoe'].append(round(dict_rhoe[str(areal)+'_'+col],2))
                p = dict_rhoe_pval[str(areal)+'_'+col]
                if p<0.01:
                    p_temp = '%.1e'% p
                    p_temp= p_temp.replace('e-0', '·10-')
                    p_temp = p_temp.replace('e-', '·10-')
                    dict_df['rhoe_pval'].append(p_temp)
                else:
                    dict_df['rhoe_pval'].append(str(round(p,2)))


    for key in dict_df.keys():
        if verbose:
            print key
            print len(dict_df[key])
        df[key] = dict_df[key]

    columns_order = ['trait', 'rhop', 'rhop_pval', 'rhog',
                     'rhog_std', 'rhog_pval', 'rhoe', 'rhoe_pval']
    df = df[columns_order]
    df.columns = ['Trait', 'rhop', 'rhop_p', 'rhog',
                  'rhog_std', 'rhog_p', 'rhoe', 'rhoe_p']

    array_result_rhop = []
    array_result_rhog = []
    array_result_rhoe = []
    for k in df.index:
         array_result_rhop.append(str(df.loc[k]['rhop'])+
                                  ' ('+str(df.loc[k]['rhop_p'])+')')
         array_result_rhog.append(str(df.loc[k]['rhog'])+
                                  '±'+str(df.loc[k]['rhog_std'])+
                                  ' ('+str(df.loc[k]['rhog_p'])+')')
         array_result_rhoe.append(str(df.loc[k]['rhoe'])+
                                  ' ('+str(df.loc[k]['rhoe_p'])+')')

    df['rhog±SE (p)'] = array_result_rhog
    df['rhop (p)'] = array_result_rhop
    df['rhoe (p)'] = array_result_rhoe
    df =df[['Trait', 'rhop (p)', 'rhog±SE (p)', 'rhoe (p)']]
    df.index = df['Trait']
    df= df.sort_index()
    df['Trait'] = [trait.replace('_', ' ') for trait in df['Trait']]
    df.to_csv(output, sep=',', header=True, index=False)
    
if __name__== '__main__':

    ROOT_DIR= ""
    indir_cog = os.path.join(ROOT_DIR
                             '2017_HCP_tfMRI_LANGUAGE',
                             'cognitive_traits')


    columns_lang = ['Language_Task_Acc', 'Language_Task_Median_RT',
	            'Language_Task_Story_Acc', 'Language_Task_Story_Median_RT',
                    'Language_Task_Story_Avg_Difficulty_Level',
                    'Language_Task_Math_Acc', 'Language_Task_Math_Median_RT',
                    'Language_Task_Math_Avg_Difficulty_Level']
    
    columns_NIH = ['PMAT24_A_CR', 'PMAT24_A_SI', 'PMAT24_A_RTCR',
                   'PicVocab_Unadj', 'PicVocab_AgeAdj', 'ProcSpeed_Unadj',
                   'ProcSpeed_AgeAdj', 'IWRD_TOT', 'IWRD_RTC',
                   'ListSort_Unadj', 'ListSort_AgeAdj', 'PicSeq_AgeAdj',
                   'PicSeq_Unadj', 'CardSort_AgeAdj', 'CardSort_Unadj',
                   'Flanker_AgeAdj', 'Flanker_Unadj', 'ReadEng_AgeAdj',
                   'ReadEng_Unadj']
    columns_NIH_AgeAdj = [p for p in columns_NIH if ('Unadj' not in p)]

    columns = columns_lang+columns_NIH_AgeAdj

    outdir = os.path.join(ROOT_DIR,
                          '2017_HCP_tfMRI_LANGUAGE',
                          'bivariate_genetic_analysis')

    
    # Original HCP files restricted and unrestricted data
    restrict = os.path.join(ROOT_DIR, '2017_HCP_tfMRI_LANGUAGE',
                            'pedigree', 'RESTRICTED_HCP_file.csv')
    unrestrict = os.path.join(ROOT_DIR, '2017_HCP_tfMRI_LANGUAGE',
                              'pedigree', 'unrestricted_HCP_file.csv')
    # Pedigree output formatted for SOLAR
    pedigree = os.path.join(ROOT_DIR, '2017_HCP_tfMRI_LANGUAGE',
                            'pedigree', 'HCP_S1200_pedigree.csv')
    
    # Create the pedigree
    df_ped = create_solar_pedigree(restrict, unrestrict, pedigree)
    work_dir = os.path.join(outdir, 'SOLAR')
    if not os.path.isdir(work_dir):
        os.makedirs(work_dir)
    os.system("solar makeped "+work_dir+" "+pedigree)

    # Path to the areals on the 41k template
    group_path = os.path.join(ROOT_DIR, '2017_HCP_tfMRI_LANGUAGE',
                              'HCP_MMP1.0')
    # Filename containing the areals labels
    filename = os.path.join(group_path, 'areal_names.csv')
    df_labels = pd.read_csv(filename)
    df_labels.index = df_labels['Index']
    areals_name = [name.replace('\n', '').replace('/', ' ').replace(' ','_')
                    for name in df_labels['Area Description']]
    
    reg = '_MSMAll'
    tasks = ['LANGUAGE']
    COPE_NUMS = [4]
    indir = os.path.join(ROOT_DIR, '2017_HCP_tfMRI_LANGUAGE',
                         'PHENOTYPES')
    fil = 'pe1'
    median = True

    columns_add = sorted(list(set(columns_NIH)-set(columns)))
    #columns = columns+columns_add
    selected = ['Language_Task_Acc',
	        'Language_Task_Story_Acc',
                'Language_Task_Story_Avg_Difficulty_Level',
                'Language_Task_Math_Acc',
                'Language_Task_Math_Median_RT',
                'Language_Task_Math_Avg_Difficulty_Level',
                'ListSort_AgeAdj',
                'ReadEng_AgeAdj',
                'PMAT24_A_CR',
                'PicVocab_AgeAdj']
    columns = selected
    
    for i, task in enumerate(tasks):
        for COP in COPE_NUMS:
            for sd in ['L', 'R']:
                file_areals = os.path.join(group_path,
                                            sd+'.fsaverage41k.label.gii')
                array_areals = gio.read(file_areals).darrays[0].data
                areals_org =  np.unique(array_areals)
                areals_org = areals_org[1:]
                file_areals = os.path.join(group_path,
                                            sd+'.fsaverage41k.label.gii')


                if median:
                    phen_dir = os.path.join(indir, task+'_'+str(COP),
                                            fil, sd,
                                            'pheno_median_value')
                    outdir2 = os.path.join(outdir, task+'_'+str(COP),
                                            fil, sd,
                                            'pheno_median_value')
                    outdir3  = os.path.join(work_dir, task+'_'+str(COP),
                                            fil, sd,
                                            'pheno_median_value')
                else:
                    phen_dir = os.path.join(indir, task+'_'+str(COP),
                                            fil, sd,
                                            'pheno_mean_value')
                    outdir2 = os.path.join(outdir, task+'_'+str(COP),
                                            fil, sd,
                                            'pheno_mean_value')
                    outdir3  = os.path.join(work_dir, task+'_'+str(COP),
                                            fil, sd,
                                            'pheno_mean_value')
                    
                compute_correlation_with_scipy(indir_cog,
                                               columns, areals_org,
                                               phen_dir, areals_name,
                                               sd, fil, outdir2)

                phen_dir = os.path.join(outdir2, 'tmp_phenotypes_csv')
                preprocess_solar_bivariate(columns_add, areals_org,
                                           phen_dir, fil, work_dir, outdir3)

                shutil.rmtree(phen_dir)

                # Input directory for SOLAR bivariate files
                indir = outdir3
                parse_bivariate_analysis(indir, columns, areals_org, sd, outdir2)

                THR = 5e-2/360.
                indir = outdir2+'/'+sd
                out_tab = os.path.join(ROOT_DIR, '2017_HCP_tfMRI_LANGUAGE'
                                       'table_results_selected')
                if not os.path.isdir(out_tab):
                    os.makedirs(out_tab)
                
                for col in columns:

                    
                    output = os.path.join(out_tab, 'bivariate_analysis_'+
                                          col+'_hem'+sd+'.csv')
                    table_bivariate_cog_traits(indir, areals_org, col,
                                               areals_name, THR, output)
