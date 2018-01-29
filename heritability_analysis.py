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
import glob
import re
import numpy as np
import pandas as pd
import nibabel.gifti.giftiio as gio

# verbose for details of race and ethnicity
verbose = False

def create_solar_pedigree(org_ped_restrict, org_ped_unrestrict, output):
    """
    Parameters
    org_ped_restricted: file containing the restricted pedigree from HCP
    org_ped_unrestricted: file containing the unrestricted pedigree from HCP
    output: output file containing SOLAR formatted pedigree
    """
    outdir = os.path.dirname(output)
    if not os.path.isdir(outdir):
        os.makedirs(outdir)
    # SOLAR pedigree required columns
    columns = ['ID', 'FA', 'MO', 'Age', 'SEX']
    # Load the full restricted HCP data
    df_src = pd.read_csv(org_ped_restrict)
    # Select columns of interest
    df = df_src[['Subject', 'Father_ID', 'Mother_ID', 'Age_in_Yrs']]
    df.index = df['Subject']
    # Gender is in unrestricted data
    df_src2 = pd.read_csv(org_ped_unrestrict)
    df2 = df_src2[['Subject', 'Gender']]
    df2.index = df2['Subject']
    # Merge the gender in first dataframe
    df['Gender']  = df2['Gender']
    # Assign the SOLAR expected column names
    df.columns = columns

    # Merge zygosity genetically tested (GT) and self reported (SR)
    df_src.index = df_src['Subject']

    SR = list(df_src['ZygositySR'])
    Zyg_corrected = [p if p != ' '  else SR[k]
                     for k, p in enumerate(df_src['ZygosityGT'])]
    Zyg_corrected = [p if p != 'DZ' else 'NotMZ' for p in Zyg_corrected]
    df_src['Zygosity_corr'] = Zyg_corrected
    
    # SOLAR expect a column with the twin pair number for MZ twins
    # The code below prepares this column (not optimized)
    # We also prepare household ID column for the pair of twin only if needed
    df_short = df_src[['Subject', 'Zygosity_corr', 'Mother_ID', 'Father_ID']]
    df_short.index = df_short['Subject']
    mz_zyg = np.zeros(len(df_short.index))
    twin_zyg_hid = np.zeros(len(df_short.index))
    mztwins_done = []
    dztwins_done = []
    count = 1
    count_hid = 1
    for i,zyg in enumerate(df_short['Zygosity_corr']):
        if zyg == 'MZ' and df_short.index[i] not in mztwins_done:
            MO = df_short.loc[df_short.index[i]]['Mother_ID']
            FA = df_short.loc[df_short.index[i]]['Father_ID']
            df_sub = df_short.loc[(df_short['Mother_ID'] == MO) &
                                  (df_short['Father_ID'] == FA) ]

            df_sub = df_sub.loc[df_sub['Zygosity_corr'] == 'MZ']
            # Check if at least 2 twins in the pair,
            # because sometimes HCP has singleton pair (at least in S900)..
            if len(df_sub.index) >= 2:
                for subj in df_sub.index:
                    mztwins_done.append(subj)
                    ind = list(df_short.index).index(subj)
                    mz_zyg[ind] = count
                    twin_zyg_hid[ind] = count_hid
                count += 1
                count_hid += 1
        elif zyg == 'NotMZ' and df_short.index[i] not in dztwins_done:
            MO = df_short.loc[df_short.index[i]]['Mother_ID']
            df_sub = df_short.loc[df_short['Mother_ID'] == MO ]
            df_sub = df_sub.loc[df_sub['Zygosity_corr'] == 'NotMZ']
            if len(df_sub.index) >= 2:
                for subj in df_sub.index:
                    mztwins_done.append(subj)
                    ind = list(df_short.index).index(subj)
                    twin_zyg_hid[ind] = count_hid
                count_hid += 1
        else:
            pass
    df['MZTWIN'] = mz_zyg
    df['MZTWIN'] = [ int(i) for i in df['MZTWIN'] ]
    # Format the columns ID, FA and MO
    for column in ['ID', 'FA', 'MO']:
        df[column] = ['%06d' % int(i) if not np.isnan(i)
                      else 0  for i in df[column]]

    # Create required family ID column by concatenating MO and FA
    l  = list(df['MO'].astype(str)+df['FA'].astype(str))
    df['FAMID'] = l
    for column in ['FAMID']:
        df[column] = ['%012d' % i for i in df[column].astype(int)]

    # Save SOLAR formatted pedigree (we did not include HID yet)
    df.to_csv(output, header=True, index=False)
    return df

def convert_pheno_to_solar(pheno, case_control, org_ped_restricted,
                           restrict_to_white=True):
    """
    Parameters
    pheno: file path containing the phenotype
    case_control: need to convert float phenotype to int
    org_ped_restricted: file containing the restricted pedigree from HCP
    restrict_to_white: boolean to restrict the analysis to white
    """
    cov = org_ped_restricted
   
    df_cov = pd.read_csv(cov)
    df_cov  = df_cov[['Subject', 'Age_in_Yrs', 'Race',
                      'Ethnicity', 'SSAGA_Educ']]
    df_cov.index = df_cov['Subject']
    # Identify all the available Race and Ethinicity labels
    race_types = list(set(df_cov['Race']))
    ethnicity_types = list(set(df_cov['Ethnicity']))
    
    # Print some statistics for the whole HCP dataset
    if verbose:
        print "\nTotal individual "+str(df_cov.shape[0])
        print " RACE"
        for race in race_types:
            print race +" "+str(df_cov.loc[df_cov['Race']== race].shape[0])
        print "\n ETHNI"
        for eth in ethnicity_types:
            print eth+" "+str(df_cov.loc[df_cov['Ethnicity']== eth].shape[0])


    df_phen = pd.read_csv(pheno)
    if 'IID' not in df_phen.columns:
        df_phen = df_phen.rename(columns={'ID': 'IID'})
    columns_phen = list(set(df_phen.columns)-set(['IID']))
    df_phen = df_phen.dropna()
    df_phen.index = df_phen['IID'].astype('int')
    # Add the Age column as covariate
    df_phen['Age'] = df_cov['Age_in_Yrs']
    # Add eTIV column as covariate
    df_phen['Educ'] = df_cov['SSAGA_Educ']
    # Find the phenotype column name
    df_phen = df_phen[[u'IID', 'Age', 'Educ']+columns_phen]
                         
    # Exclude individuals with race and ethnicity not included
    df_cov = df_cov.drop(df_cov.loc[df_cov['Race']==
                                    'Unknown or Not Reported'].index)
    df_cov = df_cov.drop(df_cov.loc[df_cov['Ethnicity']==
                                    'Unknown or Not Reported'].index)
    df_cov = df_cov.drop(df_cov.loc[df_cov['Race']==
                                    'More than one'].index)
    df_cov = df_cov.drop(df_cov.loc[df_cov['Race']==
                                    'Am. Indian/Alaskan Nat.'].index)
    df_phen = df_phen.loc[df_cov.index]

    # Bit of code to restrict the analysis to white not hispanic
    if restrict_to_white:
        df_phen = df_phen.loc[df_cov.loc[df_cov['Race']=='White'].index]
        df_cov = df_cov.loc[df_phen.index]


    # Update the list of races and ethinicities
    race_types = list(set(df_cov['Race']))
    ethnicity_types = list(set(df_cov['Ethnicity']))

    if verbose:
        # Display statistic for this phenotype
        print "\nTotal individual "+str(df_cov.shape[0])
        print " RACE"
        for race in race_types:
            print race+" "+str(df_cov.loc[df_cov['Race']== race].shape[0])
        print " ETHNI"
        for eth in ethnicity_types:
            print eth+" "+str(df_cov.loc[df_cov['Ethnicity']== eth].shape[0])
        
    count = 0
    count_max = None
    max_val = 0
    #print "\n"
    for race in race_types:
        for ethni in ethnicity_types:
            if df_cov.loc[ (df_cov['Ethnicity']== ethni) &
                           (df_cov['Race']== race) ].shape[0] != 0:
                df_phen['Group'+str(count)] = np.asarray(
                    ((df_cov['Ethnicity']== ethni) &
                     (df_cov['Race']== race))).astype('int')

                if max_val < sum(df_phen['Group'+str(count)]):
                    max_val = sum(df_phen['Group'+str(count)])
                    count_max = count
                count += 1
    # We withdraw the main group
    df_phen = df_phen.drop('Group'+str(count_max),1)
    # If any drop Nan values
    df_phen = df_phen.dropna()
    df_phen['IID'] = df_phen['IID'].astype('int')
    # Do not overwrite the previous pheno file
    dirout = os.path.join(os.path.dirname(pheno), "solar_ready")
    if not os.path.isdir(dirout):
        os.makedirs(dirout)    
    pheno = os.path.join(dirout, os.path.basename(pheno))
    # For plis de passage phenotype
    if case_control:
        for col in columns_phen:
            df_phen[col] = df_phen[col].astype('int')
    df_phen.to_csv(pheno, header=True, index=False)
    print df_phen.head()
    return df_phen, pheno
    
def run_solar(trait, task, sd, median, areal, pheno, work_dir, cc,
              org_ped_restricted):
    """
    Parameters
    trait: considered trait in phenotype file
    task: task and contrast considered
    areal: areal name considered
    median: median or mean value directory
    pheno: file path containing the phenotype
    work_dir: root folder in which solar pedigree and outputs are written
    cc: boolean stating if it's a case/control phenotype
    org_ped_restricted: file containing the restricted pedigree from HCP
    """

    outdir = os.path.join(work_dir, task, sd, median, areal, trait)
    print "\n\n\n"
    print outdir
    if not os.path.isdir(outdir):
        os.makedirs(outdir)

    # Convert phenotype to solar
    df_phen, file_phen = convert_pheno_to_solar(pheno, cc, org_ped_restricted,
                                                restrict_to_white=True)
    pheno = os.path.join(os.path.dirname(pheno), "solar_ready",
                         os.path.basename(pheno))
    # Run solar
    if cc:
        cmd = "solar case_control_analysis "+work_dir+" "+outdir+" "+trait+" "+pheno
    else:
        cmd = "solar pheno_analysis "+work_dir+" "+outdir+" "+trait+" "+pheno

    print cmd
    # Note: The code for Group covariates is not self adaptive yet
    # i.e. need to check pheno_analysis.tcl file to adjust accordingly
    os.system(cmd)
    return df_phen, file_phen

def parse_solar_out(work_dir, task, fil, df_labels,
                    median, areals, outdir='/tmp/'):
    """
    Parameters
    work_dir: root folder in which solar pedigree and outputs are written
    task: task and contrast considered
    fil: feature used to approximate activation (pe1, cope1, zstat1)
    median: boolean stating median or mean value
    areals: list of areals
    outdir: output directory for h2 and pval dictionaries
    """
    if median:
        outdir = os.path.join(outdir, fil+'_median_value', task)
    else:
        outdir = os.path.join(outdir, fil+'_mean_value', task)
    if not os.path.isdir(outdir):
        os.makedirs(outdir)
    
    # Dictionaries containing h2 and associated p-values estimates
    dict_h2 = {}
    dict_pval = {}
    # Reference list of numbers for parsing
    nb = ['.', '-', 'e', '0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10']
    for sd in ['R', 'L']:
        dict_h2[sd] = {}
        dict_pval[sd] = {}
        if median:
            sub_dir = os.path.join(work_dir, task, sd, 'median_value')
        else:
            sub_dir = os.path.join(work_dir, task, sd, 'mean_value')
        for areal in areals:
            solar_output = os.path.join(sub_dir, areal, fil,
                                        'polygenic.out')
            for line in open(solar_output, 'r'):
                # Find the heritability estimate and pval line
                if 'H2r is' in line and '(Significant)' in line:
                    if verbose:
                        print line[4:len(line)-15]
                    h2 = line[11:len(line)-30] 
                    p = line[26:len(line)-15]                
                    for k,l in enumerate(h2):
                        if not (l  in nb):
                            break
                    h2 = float(h2[:k])
                    p = float(p)
                    if verbose:
                        print "We extracted h2: "+str(h2)+" pval: "+str(p)
                    dict_h2[sd][df_labels.loc[areal]['Index']] = h2
                    dict_pval[sd][df_labels.loc[areal]['Index']] = p
                    """
                    Code below parse output in case of common environment model
                    It needs to be adapted to new variables

                    if 'C2 is' in line and '(Significant)' in line:
                        #print line[4:len(line)-15]
                        c2 = line[12:len(line)-30] 
                        p = line[26:len(line)-15]                
                        for k,l in enumerate(c2):
                            if not (l  in nb):
                                break
                        c2 = float(c2[:k])
                        p = float(p)
                        #if p<5e-2/10.0:
                        #print "We extracted c2: "+str(c2)+" pval: "+str(p)
                        #dict_h2[pheno][areal][i] = h2
                        #dict_pval[pheno][areal][i] = p
                        dict_h2[pheno][trait] = c2
                        dict_pval[pheno][trait] = p

                    elif 'C2 is' in line:# and '(Significant)' in line:
                        #print line
                        if len(line) < 22:
                            c2 = line[11:len(line)-1]
                            p = 1
                            #print "We extracted c2: "+str(c2)+" pval: "+str(p)
                        else:
                            c2 = line[12:len(line)-30] 
                            if 'Not Significant' in line:
                                p = float(line[26:len(line)-19])
                            else:
                                p = float(line[26:len(line)-15])
                            for k,l in enumerate(c2):
                                if not (l  in nb):
                                    break
                            c2 = float(c2[:k])
                            p = float(p)
                            #print "We extracted c2: "+str(c2)+" pval: "+str(p)
                        dict_h2[pheno][trait] = c2
                        dict_pval[pheno][trait] = p
                    """
    encoded = json.dumps(dict_h2)
    output_h2 = os.path.join(outdir, 'h2_dict.json')
    with open(output_h2, 'w') as f:
        json.dump(encoded, f)

    encoded = json.dumps(dict_pval)
    output_pval = os.path.join(outdir, 'pval_dict.json')
    with open(output_pval, 'w') as f:
        json.dump(encoded, f)



def table_result(work_dir, task, fil, median, areals, thr, outdir='/tmp/'):
    """
    Parameters
    work_dir: root folder in which solar pedigree and outputs are written
    task: task and contrast considered
    fil: feature used to approximate activation (pe1, cope1, zstat1)
    median: boolean stating median or mean value
    areals: list of areals
    outdir: output directory for h2 and pval dictionaries
    """
    if median:
        outdir = os.path.join(outdir, fil+'_median_value', task)
    else:
        outdir = os.path.join(outdir, fil+'_mean_value', task)
    if not os.path.isdir(outdir):
        os.makedirs(outdir)
    # Reference list of numbers for parsing
    verbose = True
    nb = ['.', '-', 'e', '0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10']
    covariates = ['age*sex', 'age^2*sex', 'age^2', 'age',  'sex', 'group', 'educ']
    elements = ['trait', 'nb_subj', 'h2', 'h2std', 'pval']
    for sd in ['R', 'L']:
        df = pd.DataFrame()
        dict_df = {}
        for elem in elements:
            dict_df[elem] = []
        for cov in covariates:
            dict_df[cov] = []
        if median:
            sub_dir = os.path.join(work_dir, task, sd, 'median_value')
        else:
            sub_dir = os.path.join(work_dir, task, sd, 'mean_value')
        for areal in areals:
            solar_output = os.path.join(sub_dir, areal, fil,
                                        'polygenic.out')
            for line in open(solar_output, 'r'):
                if 'H2r is' in line and '(Significant)' in line: 
                    p = float(line[26:len(line)-15])
                    if p < thr:
                        significant = True
                    else:
                        significant = False
                    break
                else:
                    significant = False
            if significant:
                for line in open(solar_output, 'r'):
                    if 'Trait' in line:                            
                        trait = areal #line[14:24]
                        dict_df['trait'].append(areal)                            
                        nb_subj = int(line[50:len(line)])
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
                    # Find the p-val of covariates
                    else:
                        for cov in covariates:
                            if cov in line and 'Significant' in line:
                                if 'Not Significant' in line:
                                    if verbose:
                                        print sd
                                        print line
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
        for key in dict_df.keys():
            if verbose:
                print key
                print len(dict_df[key])
            df[key] = dict_df[key]


        columns_order = ['trait', 'h2', 'h2std', 'pval', 'age', 'age^2',
                         'sex', 'age*sex', 'age^2*sex', 'group',
                         'educ', 'nb_subj']
        df = df[columns_order]
        df.columns = ['Trait', 'h²', 'std', 'p', 'Age(p)', 'Age²(p)',
                      'Sex(p)', 'Age*Sex(p)', 'Age²*Sex(p)', 'Group(p)',
                      'Educ(p)', 'Subjects']
        array_result = []
        for k in df.index:
             array_result.append((str(df.loc[k]['h²'])+'±'+
                                  str(df.loc[k]['std'])+' ('+
                                  str(df.loc[k]['p'])+')'))
        df['h²±SE(p)'] = array_result
        df =df[['Trait', 'h²±SE(p)', 'Age(p)', 'Age²(p)',
                'Sex(p)', 'Age*Sex(p)', 'Age²*Sex(p)', 'Group(p)',
                'Educ(p)', 'Subjects']]

        n_traits = []
        for key in df['Trait']:
            nkey = key.replace('_', ' ')
            n_traits.append(nkey)
        df['Trait'] = n_traits
        df.index = df['Trait']
        df= df.sort_index()
        output = os.path.join(outdir, 'table_solar_signif_result'+sd+'.csv')
        df.to_csv(output, sep=',', header=True, index=False)

        print df
        
def print_dataset_stats(org_ped_restrict, df_phen, df_ped):
    """
    Print the dataset characteristics including
    the number of twins MZ/DZ, siblings, men/women ...
    """
    # Load the full restricted HCP data
    df_src = pd.read_csv(org_ped_restrict)
    # Merge zygosity genetically tested (GT) and self reported (SR)
    df_src.index = df_src['Subject']

    SR = list(df_src['ZygositySR'])
    Zyg_corrected = [p if p != ' '  else SR[k]
                     for k, p in enumerate(df_src['ZygosityGT'])]
    Zyg_corrected = [p if p != 'DZ' else 'NotMZ' for p in Zyg_corrected]
    df_src['Zygosity_corr'] = Zyg_corrected
    
    df_ped = df_ped.loc[df_phen.index]
    print "\n\n Number of subjects: "+str(df_phen.shape[0])
    print "Hispanic "+str(sum(df_phen['Group1']))
    count_mz_twin = 0
    count_mz_twin_siblings = 0
    count_mz_twin_half_siblings = 0
    mz_done = []
    list_mz = []
    s_ids_done = []
    for i in range(1, max(df_ped['MZTWIN'])):
        if df_ped.loc[df_ped['MZTWIN'] == i].shape[0] >= 2:
            list_mz.append(list(df_ped.loc[df_ped['MZTWIN'] == i].index))
            count_mz_twin+=1
            mo = list(df_ped.loc[df_ped['MZTWIN'] == i]['MO'])[0]
            fa = list(df_ped.loc[df_ped['MZTWIN'] == i]['FA'])[0]
            mz_done.append(mo)
            count_siblings = df_ped.loc[((df_ped['MO'] ==  mo) &
                                         (df_ped['FA'] ==  fa))].shape[0]-2
            count_half_siblings = df_ped.loc[(df_ped['MO'] ==
                                              mo)].shape[0]-2-count_siblings
            count_mz_twin_siblings += count_siblings
            count_mz_twin_half_siblings += count_half_siblings                          
            s_ids_done+list(df_ped.loc[df_ped['MO'] ==  mo].index)
        elif df_ped.loc[df_ped['MZTWIN'] == i].shape[0] == 1:
            pass

    print ("Included in our analysis MZ twin pair "+ str(count_mz_twin)+
           ", their siblings "+str(count_mz_twin_siblings)+" and half sibs "+
           str(count_mz_twin_half_siblings))

    df_short = df_src[['Subject', 'Zygosity_corr', 'Mother_ID', 'Father_ID']]
    df_short.index = df_short['Subject']
    df_ped_dz = df_short.loc[df_short['Zygosity_corr'] == 'NotMZ']
    dz_id = sorted(set(df_ped_dz.index)&set(df_phen.index))
    df_ped_dz = df_ped_dz.loc[dz_id]
    dztwins_done = []
    count_dz = 0
    count_dz_twin_siblings = 0
    count_dz_twin_half_siblings = 0
    dz_done = []
    list_dz = []
    for s_id in df_ped_dz.index:
        mo = df_ped_dz.loc[s_id]['Mother_ID']
        fa = df_ped_dz.loc[s_id]['Father_ID']
        if mo not in dztwins_done:
            if df_ped_dz.loc[df_ped_dz['Mother_ID'] == mo].shape[0] >= 2:
                list_dz.append(list(df_ped_dz.loc[df_ped_dz['Mother_ID'] ==
                                                  mo].index))
                count_dz +=1
                dz_done.append(str(0)+str(mo))            

                count_siblings = df_ped.loc[((df_ped['MO'] ==
                                              '0'+str(mo)) &
                                             (df_ped['FA'] ==
                                              str(0)+str(fa)))].shape[0]-2
                count_half_siblings = df_ped.loc[(df_ped['MO'] ==
                                                  '0'+str(mo))].shape[0]-2-count_siblings
                count_dz_twin_siblings += count_siblings
                count_dz_twin_half_siblings += count_half_siblings  
                s_ids_done+list(df_ped.loc[df_ped['MO'] ==  mo].index)
                dztwins_done.append(mo)
            else:
                pass
    print ("We included "+str(len(dztwins_done))+" DZ pair and "+
           str(count_dz_twin_siblings)+" of their siblings and "+
           str(count_dz_twin_half_siblings)+" half siblings")
    mo_id_done= dz_done+mz_done
    
    count_only_siblings = 0
    count_only_half_siblings = 0
    count_alone = 0
    for s_id in df_phen.index:
        if not df_ped.loc[s_id]['MO'] in mo_id_done:
            mo = df_ped.loc[s_id]['MO']
            fa= df_ped.loc[s_id]['FA']
            mo_id_done.append(mo)
            if df_ped.loc[df_ped['MO'] == mo].shape[0] >= 2:
                count_siblings = df_ped.loc[((df_ped['MO'] ==  mo) &
                                             (df_ped['FA'] ==  fa))].shape[0]
                count_half_siblings = df_ped.loc[(df_ped['MO'] ==
                                                  mo)].shape[0]-count_siblings
                count_only_siblings += count_siblings
                count_only_half_siblings += count_half_siblings
                s_ids_done+list(df_ped.loc[df_ped['MO'] ==  mo].index)
            else:
               count_alone +=1 
               s_ids_done+list(df_ped.loc[df_ped['MO'] ==  mo].index)
    print ("We found "+str(count_only_siblings)+" siblings "+
           str(count_only_half_siblings)+" half siblings without twins included "+
           str(count_alone)+" unrelated subjects")
    print df_ped.head()
    print df_ped.shape
    print ("Min age "+str(np.min(df_ped['Age']))+
           ", Max age"+str(np.max(df_ped['Age'])))
    print ("Mean age "+str(np.mean(df_ped['Age']))+
           " std "+str(np.std(df_ped['Age'])))
    print ("Number men "+str(sum(df_ped['SEX'] == 'M'))+
           " women "+str(sum(df_ped['SEX'] == 'F')))


def kinship_of_individuals(pedigree, phen_file):

    """
    Methods to create usable kinship matrix for R package RSP

    Parameters:
    phen_file: take one phenotype for the task considered,
       no matter each areal they all contain the same individuals
    pedigree: HCP pedigree
    """
    
    df_phen = pd.read_csv(phen_file)
    df_phen.index = df_phen['IID']
    df_ped = pd.read_csv(pedigree)
    df_ped.index = df_ped['ID']
    df = df_ped.loc[df_phen.index]
    for dad in list(set(df_ped['FA'])):
        df.loc[dad] = [dad, 0, 0, 0, 'M', 0, 0]
    for mum in list(set(df_ped['MO'])):
        df.loc[mum] = [mum, 0, 0, 0, 'F', 0, 0]
    df = df.replace('M', 1)
    df = df.replace('F', 2)
    DIR_ped = os.path.join(ROOT_DIR, '2017_HCP_tfMRI_LANGUAGE'
                           'pedigree')
    output = os.path.join(DIR_ped, 'RSP_kinship_format.csv')
    df.to_csv(output, header=True, index=False)

    cmd = 'Rscript equivalent_power_calculation.R '+output
    os.system(cmd)
    
    
if __name__== '__main__':

    ROOT_DIR= ""
    
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
    
    # Folder in which solar pedigree and outputs are written.
    # cf scripts makeped.tcl & pheno_analysis.tcl
    work_dir = os.path.join(ROOT_DIR, '2017_HCP_tfMRI_LANGUAGE',
                            'SOLAR_PHEN', 'solar_work_dir')
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
    df_labels['Area Description'] = areals_name
    df_labels.index = df_labels['Area Description']    

    reg = '_MSMAll'
    tasks = ['LANGUAGE']
    COPE_NUMS = [1, 2]
    indir = os.path.join(ROOT_DIR,
                         '2017_HCP_tfMRI_LANGUAGE',
                         'PHENOTYPES')
    fil = 'pe1'
    median = True
    
    for i, task in enumerate(tasks):
        for COP in COPE_NUMS:
            for sd in ['R', 'L']:
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
                else:
                    phen_dir = os.path.join(indir, task+'_'+str(COP),
                                            fil, sd,
                                            'pheno_mean_value')

                for k, areal in enumerate(areals_org):
                    if median:
                        fname = os.path.join(phen_dir,
                                             areals_name[k]+'_'+fil+'.csv')
                        phen_file = os.path.join(phen_dir, 'solar_ready',
                                                 areals_name[k]+'_'+fil+'.csv')
                        df_phen = run_solar(fil, task+'_'+str(COP), sd,
                                            'median_value', areals_name[k],
                                            fname, work_dir, False, restrict)
                    else:
                        fname = os.path.join(phen_dir,
                                             areals_name[k]+'_'+fil+'.csv')
                        df_phen = run_solar(fil, task+'_'+str(COP), sd,
                                            'mean_value', areals_name[k],
                                            fname, work_dir, False, restrict)
    
    print_dataset_stats(restrict, df_phen, df_ped)

    kinship_of_individuals(pedigree, phen_file)
    
    df_labels.index = df_labels['Area Description']
    outdir = os.path.join(ROOT_DIR,
                          '2017_HCP_tfMRI_LANGUAGE',
                          'dictionaries_heritability')
    for i, task in enumerate(tasks):
        for COP in COPE_NUMS:
            parse_solar_out(work_dir, task+'_'+str(COP), fil, df_labels,
                            median, areals_name, outdir=outdir)

            thr = 5e-2/(360)
            table_result(work_dir, task+'_'+str(COP), fil,
                         median, areals_name, thr, outdir=outdir)

    
