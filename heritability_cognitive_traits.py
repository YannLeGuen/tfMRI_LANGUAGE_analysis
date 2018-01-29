#! /usr/bin/env python
# -*- coding: utf-8 -*

##########################################################################
# Copyright (C) CEA, 2017
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html for details.
##########################################################################

import os
import json
import pandas as pd
from scipy.stats import pearsonr

from heritability_analysis import create_solar_pedigree
from heritability_analysis import convert_pheno_to_solar

def select_column(unrestrict, columns, outdir):
    """
    Create phenotype with a single phenotype
    """
    if not os.path.isdir(outdir):
        os.makedirs(outdir)
    df = pd.read_csv(unrestrict)
    df['IID'] = df['Subject']
    for col in columns:
        df_sub = df[['IID', col]].dropna()
        output = os.path.join(outdir, col)
        df_sub.to_csv(output+'.csv',  header=True, index=False)

        
def run_solar(trait, pheno, work_dir, cc):
    """
    Parameters
    trait: considered trait in phenotype file
    task: task and contrast considered
    areal: areal name considered
    median: median or mean value directory
    pheno: file path containing the phenotype
    work_dir: root folder in which solar pedigree and outputs are written
    cc: boolean stating if it's a case/control phenotype
    """

    outdir = os.path.join(work_dir, 'cognitive_trait', trait)
    print "\n\n\n"
    print outdir
    if not os.path.isdir(outdir):
        os.makedirs(outdir)

    # Convert phenotype to solar
    df_phen = convert_pheno_to_solar(pheno, cc, restrict_to_white=True)
    pheno = os.path.join(os.path.dirname(pheno), "solar_ready",
                         os.path.basename(pheno))
    # Run solar
    if cc:
        cmd = ("solar case_control_analysis "+work_dir+" "+
               outdir+" "+trait+" "+pheno)
    else:
        cmd = "solar pheno_analysis "+work_dir+" "+outdir+" "+trait+" "+pheno

    print cmd
    # Note: The code for Group covariates is not self adaptive yet
    # i.e. need to check pheno_analysis.tcl file to adjust accordingly
    os.system(cmd)
    return df_phen

def table_result(work_dir, columns, prefix, outdir='/tmp/'):
    """
    Parameters
    work_dir: root folder in which solar pedigree and outputs are written
    task: task and contrast considered
    fil: feature used to approximate activation (pe1, cope1, zstat1)
    median: boolean stating median or mean value
    areals: list of areals
    outdir: output directory for h2 and pval dictionaries
    """
    if not os.path.isdir(outdir):
        os.makedirs(outdir)
    # Reference list of numbers for parsing
    verbose = False
    cnt = 0
    nb = ['.', '-', 'e', '0', '1', '2', '3', '4',
          '5', '6', '7', '8', '9', '10']
    covariates = ['age*sex', 'age^2*sex', 'age^2',
                  'age',  'sex', 'group', 'educ']
    elements = ['trait', 'nb_subj', 'h2', 'h2std', 'h2cov', 'pval']
    df = pd.DataFrame()
    dict_df = {}
    for elem in elements:
        dict_df[elem] = []
    for cov in covariates:
        dict_df[cov] = []

    sub_dir = os.path.join(work_dir, 'cognitive_trait')
    for col in columns:
        solar_output = os.path.join(sub_dir, col,
                                    'polygenic.out')
        for line in open(solar_output, 'r'):
            if 'Trait' in line:
                trait = col #line[14:24]
                dict_df['trait'].append(col)
                ref = len('Trait:       analyze_'+col+'  Individuals:  ')
                if  ref < 50:
                    nb_subj = int(line[50:len(line)-1])
                else:
                    nb_subj = int(line[ref:len(line)])
                dict_df['nb_subj'].append(nb_subj)
                if verbose:
                    print line
                    print 'Trait extracted: '+str(trait)
                    print 'Subjects extracted: '+ str(nb_subj)
            # Find the heritability estimate and pval line
            elif 'H2r is' in line:#and '(Significant)' in line:
                if verbose:
                    print line[4:len(line)-15]
                h2 = line[11:len(line)-30] 
                if 'Not Significant' in line:
                    p = float(line[26:len(line)-19])
                else:
                    p = float(line[26:len(line)-15])
                for k,l in enumerate(h2):
                    if not (l  in nb):
                        break
                h2 = float(h2[:k])
                p = float(p)
                if verbose:
                    print "We extracted h2: "+str(h2)+" pval: "+str(p)
                dict_df['h2'].append(round(h2, 2))
                if p < 0.01:
                    p_temp = '%.1e'% p
                    p_temp= p_temp.replace('e-0', '·10-')
                    p_temp = p_temp.replace('e-', '·10-')
                    dict_df['pval'].append(p_temp)
                else:
                    dict_df['pval'].append(str(round(p,2)))

                if h2 == 0:
                    dict_df['h2std'].append(0)
            # Find the standard deviation of heritability estimate
            elif 'H2r Std.' in line:
                h2std= float(line[25:len(line)])                        
                if verbose:
                    print line
                    print 'Std extracted: '+str(h2std)
                dict_df['h2std'].append(round(h2std,2))
            elif 'Proportion of Variance' in line:
                if verbose:
                    print line
                cnt = 1
            elif cnt == 1:
                h2cov = float(line[6:15])
                if verbose:
                    print line
                    print "we found h2cov: "+str(h2cov)
                cnt = 0
                dict_df['h2cov'].append(round(h2cov*100,1))
            # Find the p-val of covariates
            else:
                for cov in covariates:
                    if cov in line and 'Significant' in line:
                        if 'Not Significant' in line:
                            p = float(line[47:len(line)-20])
                        else:
                            p = float(line[47:len(line)-16])
                        if verbose:
                            print line
                            print 'pval found for '+cov+ ' '+str(p)
                        if p < 0.01:
                            p_temp = '%.1e'% p
                            """p_temp0 = p_temp[:5]
                            for k,c in enumerate(p_temp[5:]):
                            if k==0 and int(c)==0:
                            pass
                            else:
                            p_temp0+=list_exponent[int(c)]"""
                            p_temp = p_temp.replace('e-0', '·10-')
                            p_temp = p_temp.replace('e-', '·10-')
                            dict_df[cov].append(p_temp)
                        else:
                            dict_df[cov].append(str(round(p,2)))
                        break
        # case when covariates explain 0% of variance
        if len(dict_df['h2std']) > len(dict_df['h2cov']):
            dict_df['h2cov'].append(0)
    for key in dict_df.keys():
        if verbose:
            print key
            print len(dict_df[key])
        df[key] = dict_df[key]


    columns_order = ['trait', 'h2', 'h2std', 'pval', 'age', 'age^2', 'sex',
                     'age*sex', 'age^2*sex', 'group', 'educ', 'h2cov', 'nb_subj']
    df = df[columns_order]
    df.columns = ['Trait', 'h²', 'std', 'p', 'Age(p)', 'Age²(p)', 'Sex(p)',
                  'Age*Sex(p)', 'Age²*Sex(p)', 'Group(p)',
                  'Educ(p)', 'h²cov(%)', 'Subjects']
    array_result = []
    for k in df.index:
         array_result.append((str(df.loc[k]['h²'])+'±'+str(df.loc[k]['std'])+
                              ' ('+str(df.loc[k]['p'])+')'))
    df['h²±SE(p)'] = array_result
    df =df[['Trait', 'h²±SE(p)', 'Age(p)', 'Age²(p)', 'Sex(p)', 'Age*Sex(p)',
            'Age²*Sex(p)', 'Group(p)', 'Educ(p)', 'h²cov(%)', 'Subjects']]

    n_traits = []
    for key in df['Trait']:
        nkey = key.replace('_', ' ')
        n_traits.append(nkey)
    df['Trait'] = n_traits
    df.index = df['Trait']
    df= df.sort_index()
    output = os.path.join(outdir, prefix+'_table_solar_result.csv')
    df.to_csv(output, sep=',', header=True, index=False)

    print df

    
def compute_correlation_with_scipy(indir, traits, outdir):

    dict_rhop = {}
    dict_rhop_pval = {}
    outdir2 = os.path.join(outdir, 'tmp_phenotypes_csv')
    if not os.path.isdir(outdir2):
        os.makedirs(outdir2)
    output = outdir+'/'
    traits = sorted(traits)
    for k, col1 in enumerate(traits):
        fname = os.path.join(indir, col1+'.csv')
        df1 = pd.read_csv(fname, sep=',')
        df1.index = df1['IID']
        for col2 in traits[k+1:]:
            fname2 = os.path.join(indir, col2+'.csv')
            df2 = pd.read_csv(fname2, sep=',')
            df2.index = df2['IID']
            df = pd.DataFrame()
            df[col1] = df1[col1]
            df[col2] = df2[col2]
            df = df.dropna()
            df['IID'] = df.index
            output2 = os.path.join(outdir2, col1+'_'+col2+'.csv')
            df.to_csv(output2, header=True, index=False)
            # Compute correlation between activation and cognitive trait
            corr, p = pearsonr(df[col1], df[col2])
            dict_rhop[col1+'_'+col2] = corr
            dict_rhop_pval[col1+'_'+col2] = p
            
    encoded = json.dumps(dict_rhop)
    with open(output+'scipy_rhop_dict.json', 'w') as f:
        json.dump(encoded, f)

    encoded = json.dumps(dict_rhop_pval)
    with open(output+'scipy_rhop_pval_dict.json', 'w') as f:
        json.dump(encoded, f)

def preprocess_solar_bivariate(traits, phen_dir, work_dir, outdir):
    traits = sorted(traits)
    
    if not os.path.isdir(outdir):
        os.makedirs(outdir)
        
    for k, col1 in enumerate(traits):
        for col2 in traits[k+1:]:
            pheno = os.path.join(phen_dir, col1+'_'+col2+'.csv')
            run_solar_bivariate(work_dir, pheno, col1, col2 , outdir)

            
def run_solar_bivariate(work_dir, pheno, trait1, trait2, outdir):
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
                    outdir, trait1, trait2, pheno])
                    
    print cmd
    # Note: The code for Group covariates is not self adaptive yet
    # i.e. need to check pheno_analysis.tcl file to adjust accordingly
    os.system(cmd)


def parse_bivariate_analysis(indir, traits, outdir):
    """
    Parse SOLAR bivariate analysis outputs
    indir: SOLAR working directory where each bivariate folder is contained
    traits: all the cognitive traits considered
    outdir: output directory for the json dictionnaries
    """
    
    SANITY_CHECK = False
    verbose = False
    nb = ['.', '-', 'e', '0', '1', '2', '3', '4',
          '5', '6', '7', '8', '9', '10']
    output = outdir +'/'

    traits = sorted(traits)

    dict_rhoe = {}
    dict_rhoe_pval = {}
    dict_rhog = {}
    dict_rhog_std = {}
    dict_rhog_pval = {}
    dict_rhop = {}
    dict_rhop_pval = {}
    # Loop through all selected colums of cognitive traits
    for k, col1 in enumerate(traits):
        print k, col1
        for col2 in traits[k+1:]:
            file_path = os.path.join(indir, col1, col2, 'polygenic.out')
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
                                print "We extracted rhoe: "+str(h2)+" pval: "+str(p)
                        if verbose:
                            print "We extracted rhoe: "+str(h2)+" pval: "+str(p)
                        dict_rhoe[col1+'_'+col2] = h2
                        dict_rhoe_pval[col1+'_'+col2] = p
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
                        dict_rhog[col1+'_'+col2] = h2

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
                            dict_rhog_std[col1+'_'+col2] = h2
                        else:
                            dict_rhog_std[col1+'_'+col2] = 1
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
                        dict_rhog_pval[col1+'_'+col2] = p 
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
                        dict_rhop[col1+'_'+col2] = h2
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
                        dict_rhop_pval[col1+'_'+col2] = p  
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

def table_bivariate_cog_traits(indir, rows, columns, output):

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

    df = pd.DataFrame(columns=columns, index=rows)
    for k1, col1 in enumerate(rows):
        line = []
        for k2, col2 in enumerate(columns):
            cols = [col1, col2]
            # neeed the columns to be in alphabetical order
            c1, c2 = sorted(cols)

            if col2 == col1:
                line.append('*')
            elif k1 > k2:
                text = str(round(dict_scipy_rhop[c1+'_'+c2],2))
                p = dict_scipy_pval[c1+'_'+c2]
                if p<0.01:
                    p_temp = '%.0e'% p
                    p_temp = p_temp[1:]
                    p_temp= p_temp.replace('e-0', '10-')
                    p_temp = p_temp.replace('e-', '10-')
                else:
                    p_temp = str(round(p,2))
                text = text+' ('+p_temp+')'
                line.append(text)
            elif k2> k1:
                if dict_rhog.has_key(c1+'_'+c2):
                    text = str(round(dict_rhog[c1+'_'+c2],2))
                    text = text+'±'+str(round(dict_rhog_std[c1+'_'+c2],2))
                    p = dict_rhog_pval[c1+'_'+c2]
                    if p<0.01:
                        p_temp = '%.0e'% p
                        p_temp = p_temp[1:]
                        p_temp= p_temp.replace('e-0', '10-')
                        p_temp = p_temp.replace('e-', '10-')
                    else:
                        p_temp = str(round(p,2))
                    text = text+' ('+p_temp+')'
                else:
                    text = 'Rhog ERROR'
                    print "Rhog ERROR: "+c1+'_'+c2
                line.append(text)
        #print df.head
        #print line
        df.loc[col1] = line
        df.loc[col1] = df.loc[col1].astype('str')

    #print df
    df.to_csv(output, header=True, index=True)
if __name__== '__main__':

    ROOT_DIR = ""
    
        # Original HCP files restricted and unrestricted data
    restrict = os.path.join(ROOT_DIR, '2017_HCP_tfMRI_LANGUAGE',
                            'pedigree', 'RESTRICTED_HCP_file.csv')
    unrestrict = os.path.join(ROOT_DIR, '2017_HCP_tfMRI_LANGUAGE',
                              'pedigree', 'unrestricted_HCP_file.csv')
    # Pedigree output formatted for SOLAR
    pedigree = os.path.join(ROOT_DIR, '2017_HCP_tfMRI_LANGUAGE',
                            'pedigree', 'HCP_S1200_pedigree.csv')

    outdir = os.path.join(ROOT_DIR,
                          '2017_HCP_tfMRI_LANGUAGE',
                          'cognitive_traits')

    columns_lang = ['Language_Task_Acc', 'Language_Task_Median_RT',
	            'Language_Task_Story_Acc', 'Language_Task_Story_Median_RT',
                    'Language_Task_Story_Avg_Difficulty_Level',
                    'Language_Task_Math_Acc', 'Language_Task_Math_Median_RT',
                    'Language_Task_Math_Avg_Difficulty_Level']
    select_column(unrestrict, columns_lang, outdir)
    
    columns_NIH = ['PMAT24_A_CR', 'PMAT24_A_SI', 'PMAT24_A_RTCR',
                   'PicVocab_Unadj', 'PicVocab_AgeAdj', 'ProcSpeed_Unadj',
                   'ProcSpeed_AgeAdj', 'IWRD_TOT', 'IWRD_RTC',
                   'ListSort_Unadj', 'ListSort_AgeAdj', 'PicSeq_AgeAdj',
                   'PicSeq_Unadj', 'CardSort_AgeAdj', 'CardSort_Unadj',
                   'Flanker_AgeAdj', 'Flanker_Unadj', 'ReadEng_AgeAdj',
                   'ReadEng_Unadj']
    select_column(unrestrict, columns_NIH, outdir)
    columns_NIH_AgeAdj = [p for p in columns_NIH if ('Unadj' not in p)] 
    columns = columns_lang+columns_NIH_AgeAdj

    ########## Univariate genetic analysis of the cognitive traits ############
    work_dir = os.path.join(ROOT_DIR,
                            '2017_HCP_tfMRI_LANGUAGE',
                            'SOLAR_PHEN', 'solar_work_dir')
    if not os.path.isdir(work_dir):
        os.makedirs(work_dir)

    os.system("solar makeped "+work_dir+" "+pedigree)

    for col in columns_lang+columns_NIH:
        pheno = os.path.join(outdir, col+'.csv')
        run_solar(col, pheno, work_dir, False)


    outdir = os.path.join(ROOT_DIR,
                          '2017_HCP_tfMRI_LANGUAGE',
                          'dictionaries_heritability',
                          'cognitive_traits')

    for cols, prefix in zip([columns_lang, columns_NIH_AgeAdj],
                            ['Language_tfMRI_cognitive', 'NIH_cognitive']):
        table_result(work_dir, cols, prefix, outdir=outdir)



    ########### Bivariate genetic analysis of the cognitive traits ############
    work_dir = os.path.join(ROOT_DIR,
                            '2017_HCP_tfMRI_LANGUAGE',
                            'SOLAR_PHEN', 'solar_work_dir')
    if not os.path.isdir(work_dir):
        os.makedirs(work_dir)

    os.system("solar makeped "+work_dir+" "+pedigree)

    
    indir = os.path.join(ROOT_DIR '2017_HCP_tfMRI_LANGUAGE'
                         'cognitive_traits')
    outdir =  os.path.join(indir, 'bivariate_analysis')

    compute_correlation_with_scipy(indir, columns, outdir)

    phen_dir = os.path.join(outdir, 'tmp_phenotypes_csv')

    work_dir = os.path.join(outdir, 'solar_work_dir')
    if not os.path.isdir(work_dir):
        os.makedirs(work_dir)

    outdir =  os.path.join(indir, 'bivariate_analysis', 'SOLAR')
    preprocess_solar_bivariate(columns, phen_dir, work_dir, outdir)

    solar_dir = outdir 
    outdir =  os.path.join(indir, 'bivariate_analysis')
    parse_bivariate_analysis(solar_dir, columns, outdir)
    output = os.path.join(outdir, 'table_cognitive_traits_rho.csv')
    table_bivariate_cog_traits(outdir+'/', columns, columns, output)
