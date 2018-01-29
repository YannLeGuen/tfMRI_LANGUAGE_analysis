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
import numpy as np
import pandas as pd
import nibabel.gifti.giftiio as gio
from multiprocessing import cpu_count
from multiprocessing import Pool
import time


def median_value_phenotyping(parameters):
    """
    Parameters
    task: HCP task name
    reg: type of registration either 'MSMAll' or '' for MSMSulc
    COPE: COPE number cf https://wiki.humanconnectome.org/display/..
                                ..PublicData/Task+fMRI+Contrasts
    s_ids: list of subject ids
    sd: hemisphere side
    areal: areal number on the texture of the HCP multimodal parcellation
    k: numero of the areal in areals_name table
    outdir: output directory
    """
    task, reg, COPE, s_ids, file_areals, sd, areal, k, outdir = parameters
    array_areals = gio.read(file_areals).darrays[0].data
    index = np.where(array_areals == areal)[0]
    array_s_ids = []
    array_mu = []
    for s_id in s_ids:
        f_path = os.path.join(path_org, '3T_tfMRI_'+task, s_id,
                              'MNINonLinear/Results',
                              'tfMRI_'+task,
                              'tfMRI_'+task+'_hp200_s2_level2'+reg+'.feat',
                              'GrayordinatesStats/cope'+str(COPE)+'.feat')
        # For some subjects the MSMAll reg is not available in HCP database
        # this does not apply to the LANGUAGE tasks
        if not os.path.isdir(f_path):
            reg_bis = ''
            f_path = os.path.join(path_org, '3T_tfMRI_'+task, s_id,
                                  'MNINonLinear/Results',
                                  'tfMRI_'+task,
                                  'tfMRI_'+task+'_hp200_s2_level2'+
                                  reg_bis+'.feat',
                                  'GrayordinatesStats/cope'+str(COPE)+
                                  '.feat')
            
        metric_out = os.path.join(f_path, sd+"_"+fil+".41k_fsavg_"+
                                  sd+".func.gii")
        if os.path.isfile(metric_out):
            array_s_ids.append(s_id)
            array_fil = gio.read(metric_out).darrays[0].data
            if median:
                #array_mu.append(np.amax(array_fil[index]))
                array_mu.append(np.median(array_fil[index]))
            else:
                array_mu.append(np.mean(array_fil[index]))
        else:
            pass
            #print metric_out+" doesn't exist for "+s_id
    df = pd.DataFrame()
    df['IID'] = array_s_ids
    df[fil] = array_mu
    output = os.path.join(outdir, areals_name[k]+'_'+fil)
    print "saving "+output
    df.to_csv(output+'.csv',  header=True, index=False)






if __name__ == '__main__': 

    ROOT_DIR = ""
    # Filename containing the areals labels
    filename = os.path.join(ROOT_DIR, '2017_HCP_tfMRI_LANGUAGE',
                            'HCP_MMP1.0/areal_names.csv')
    df_labels = pd.read_csv(filename)
    df_labels.index = df_labels['Index']
    areals_name = [name.replace('\n', '').replace('/', ' ').replace(' ','_')
                    for name in df_labels['Area Description']]
    df_labels['Area Description'] = areals_name

    # Path to the database containing the functional activation in gii format

    path_org = os.path.join(ROOT_DIR, 'HCP_tfMRI')
    # Path where the phenotypes files are saved
    p_dest = os.path.join(ROOT_DIR,
                          '2017_HCP_tfMRI_LANGUAGE',
                          'PHENOTYPES')
   
    
    tasks = ['LANGUAGE']
    COPE_NUMS = [[1,4]]
    registrations = ['_MSMAll']#['', '_MSMAll']
    # We consider the parameter estimates (beta) as proxy to our activations
    #fils = ['pe1', 'cope1', 'zstat1'] # cope1 and pe1 are identical
    fil = 'pe1'
    # Set median to True to use median pe value in the areal as phenotype,
    # else set to False to use mean pe value.
    median = True

    # Path to the areals on the 41k template
    group_path = os.path.join(ROOT_DIR, '2017_HCP_tfMRI_LANGUAGE'
                              'functional_analysis/HCP_MMP1.0')
    file_areals = os.path.join(group_path, 'L'+'.fsaverage41k.label.gii')
    array_areals = gio.read(file_areals).darrays[0].data
    areals_org =  np.unique(array_areals)
    areals_org = areals_org[1:]

    parameters = []
    # Select the task
    for i, task in enumerate(tasks):
        path_org_task = os.path.join(path_org, '3T_tfMRI_'+task)
        s_ids = [s_id for s_id in os.listdir(path_org_task)
                 if os.path.isdir(os.path.join(path_org_task, s_id))]
        print task+' '+str(len(s_ids))
        # Select registration type MSMAll or MSMSulc = ''
        for reg in registrations:
            for COP in range(COPE_NUMS[i][0],COPE_NUMS[i][1]+1):
                for sd in ['R', 'L']:
                    file_areals = os.path.join(group_path,
                                                sd+'.fsaverage41k.label.gii')

                    for k, areal in enumerate(areals_org):
                        if median:
                            outdir = os.path.join(p_dest, task+'_'+str(COP),
                                                  fil, sd, 'pheno_median_value')
                        else:
                            outdir = os.path.join(p_dest, task+'_'+str(COP),
                                                  fil, sd, 'pheno_mean_value')
                        if not os.path.isdir(outdir):
                            os.makedirs(outdir)
                        parameters.append([task, reg, COP, s_ids, file_areals,
                                           sd, areal, k, outdir])
    t0 = time.time()
    number_CPU = cpu_count()-4
    pool = Pool(processes = number_CPU)
    pool.map(median_value_phenotyping, parameters)
    pool.close()
    pool.join()
    print "Elapsed time for all areals phenotyping: "+str(time.time()-t0)
