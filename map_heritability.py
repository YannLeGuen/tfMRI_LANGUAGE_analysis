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
import nibabel.gifti.giftiio as gio
# Import to create a poster figure with the snapshots
import PIL.Image
# Import related to Anatomist object manipulation
from soma import aims
import anatomist.api 
ana = anatomist.api.Anatomist()
import paletteViewer
# Define the referential and transformation necessary,
# because Freesurfer meshes are indexed in the reverse order.
ref = ana.createReferential()
tr = ana.createTransformation([0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, -1],
                              ref, ana.centralReferential())

# Load seaborn colorbars in Anatomist
DIR_CB = ""
cbar = ['YlOrRd', 'RdYlGn']
for cb in cbar:
    directory = os.path.join(DIR_CB, 'COLORBAR')
    filename = os.path.join(directory, cb+'.txt')
    c_sns = np.loadtxt(filename)
    c_sns = c_sns.astype('int')

    if cb == 'RdYlGn':
        c_sns = np.concatenate((c_sns,[[200,200,200]]))
        c_sns = np.flipud(c_sns)
    elif cb == 'YlOrRd':
        c_sns = np.concatenate(([[200,200,200]], c_sns))
        
    pal = ana.createPalette(cb)
    pal.setColors(c_sns.ravel(), 'RGB')


# Method to notify the Anatomist window of changes and remove cursor
def updateWindow(window, obj):
    window.addObjects(obj)
    obj.setChanged()
    obj.notifyObservers()
    ana.execute('WindowConfig', windows=[window], cursor_visibility=0)



def map_heritability(template, sd, areals, indir, cb, interval_cb,
                          pval_thresh, outdir):
    """
    Parameters
    template: path to template mesh file
    areals: path to gii file containing the texture of areals
    sd: hemisphere side from which SOLAR estimates are displayed
    indir: directory containing h2 and pval dictionaries
    cb: colorbar names must be an array with [cb_h2, cb_pval]
    interval_cb: number boundaries for the colobar display
                 must be a double array [inter_h2[0,1], inter_pval[0,1]]
    pval_thresh: p-value threshold
    outdir: output directory for the snapshots
    """
    if not os.path.isdir(outdir):
        os.makedirs(outdir)
    # Define the quaternions that are gonna be used for the snapshots
    if (sd == "L"):
        view_quaternions = {'intern' : [0.5, 0.5, 0.5, 0.5],
                            'extern' : [0.5, -0.5, -0.5, 0.5],
                            'bottom' : [1, 0, 0, 0],
                            'top': [0, 0, 1, 0]}
    else:
        view_quaternions = {'intern' : [0.5, -0.5, -0.5, 0.5],
                            'extern' : [0.5, 0.5, 0.5, 0.5],
                            'bottom' : [1, 0, 0, 0],
                            'top': [0, 0, 1, 0]}

    # Load h2 and pval dictionaries
    path_h2 = os.path.join(indir, 'h2_dict.json')
    with open(path_h2, 'r') as f:
        data = json.load(f)
    dict_h2 = json.loads(data)
    path_pval = os.path.join(indir, 'pval_dict.json')
    with open(path_pval, 'r') as f:
        data = json.load(f)
    dict_pval = json.loads(data)

    # Areal value in each vertex
    array_areals = gio.read(areals).darrays[0].data
    # Create arrays that will contain h2 and pval estimates on each vertex
    array_h2 = np.zeros(len(array_areals))
    array_pval = np.zeros(len(array_areals))
    list_areals =  np.unique(array_areals)
    
    # Load the mesh in anatomist
    mesh = ana.loadObject(template)

    for areal in list_areals:
        # make sure areal has the right format
        areal = str(int(areal))
        if dict_h2[sd].has_key(areal):
            # Finally, checl if the p-value pass the pval threshold
            if dict_pval[sd][areal] < pval_thresh :
                # Select the index of vertices belonging to the areal
                ind = np.where(array_areals == float(areal))[0]
                # Give them the h2 estimate of this areal
                array_h2[ind] = dict_h2[sd][areal]
                # And the log10 pvalue
                array_pval[ind] = -np.log10(dict_pval[sd][areal])

    ### Block performing the display of h2 estimates ###
    # Create anatomist window
    window = ana.createWindow('3D', geometry=[0, 0, 584, 584])
    tex = aims.TimeTexture(dtype='FLOAT')
    tex[0].assign(array_h2)                  
    pits_tex = ana.toAObject(tex)
    tex_mesh = ana.fusionObjects([mesh, pits_tex], method='FusionTexSurfMethod')
    tex_mesh.assignReferential(ref)
    tex_mesh.setMaterial(front_face='counterclockwise')
    pits_tex.setPalette(cb[0], minVal=interval_cb[0][0],
                        maxVal=interval_cb[0][1], absoluteMode=True)
    updateWindow(window, tex_mesh)

    # Loop through the quaternions and do the snapshot
    for vw in view_quaternions.keys():
            q = aims.Quaternion(view_quaternions[vw])
            window.camera(view_quaternion=view_quaternions[vw],
                                 zoom=0.65)
            # Snapshot file output
            output = os.path.join(outdir, 'snapshot_h2_'+sd+'_'+vw+'.png')
            ana.execute('WindowConfig', windows=[window], snapshot=output)


    #### Block performing the display of pvalues estimates ###
    # Create anatomist window, second time (else zoom problem for snapshots)
    window = ana.createWindow('3D', geometry=[0, 0, 584, 584])
    tex = aims.TimeTexture(dtype='FLOAT')
    tex[0].assign(array_pval)                  
    pits_tex = ana.toAObject(tex)
    tex_mesh = ana.fusionObjects([mesh, pits_tex],
                                 method='FusionTexSurfMethod')
    tex_mesh.assignReferential(ref)
    tex_mesh.setMaterial(front_face='counterclockwise')
    pits_tex.setPalette(cb[1], minVal=interval_cb[1][0],
                        maxVal=interval_cb[1][1], absoluteMode=True)
    updateWindow(window, tex_mesh)
    
    # Loop through the quaternions and do the snapshot
    for vw in view_quaternions.keys():
            q = aims.Quaternion(view_quaternions[vw])
            window.camera(view_quaternion=view_quaternions[vw],
                                 zoom=0.65)
            # Snapshot file output
            output = os.path.join(outdir, 'snapshot_pval_'+sd+'_'+vw+'.png')
            ana.execute('WindowConfig', windows=[window], snapshot=output)

    ana.releaseObject(pits_tex)
    pits_tex = None
    ana.releaseObject(tex_mesh)
    tex_mesh = None


    
def create_poster(images, output_image, tile=4, padding=4, bgcolor=(255,255,255),
                  align_grid=False):
    # print 'calculating poster size...'
    rowsizes = []
    maxsize = [0, 0]
    n = 0
    for imagename in images:
        image = PIL.Image.open( imagename )
        if n % tile == 0:
            rowsizes.append( ( 0, 0 ) )
        rowsizes[-1] = ( rowsizes[-1][0] + image.size[0] + padding,
            max( rowsizes[-1][1], image.size[1] + padding ) )
        maxsize = (max((maxsize[0], image.size[0])),
                   max((maxsize[1], image.size[1])))
        n += 1
    print 'sizes:', rowsizes
    if align_grid:
        size = ((maxsize[0] + padding) * min(tile, len(images)) - padding,
                (maxsize[1] + padding) * len(rowsizes) - padding)
    else:
        size = ( max( [ rs[0] for rs in rowsizes ] ) +padding,
                 sum( [ rs[1] for rs in rowsizes ] ) +padding)

    print 'size:', size
    outimage = PIL.Image.new( 'RGB', size, bgcolor )
    print 'creating image...'
    n = 0
    line = 0
    xpos = 0
    ypos = 0
    #rowsizes.insert( 0, ( 0, 0 ) )
    for imagename in images:
        image = PIL.Image.open( imagename )
        if align_grid:
            bbox = ((maxsize[0] + padding) * n, (maxsize[1] + padding) * line,
                    0, 0)
            bbox = (bbox[0], bbox[1],
                    bbox[0] + image.size[0], bbox[1] + image.size[1])
        else:
            y = ypos + ( rowsizes[line][1] - image.size[1] ) / 2
            bbox = ( xpos, y, xpos+image.size[0], y + image.size[1] )
        outimage.paste( image, bbox )
        xpos = bbox[2] + padding
        n += 1
        if n == tile:
            n = 0
            line += 1
            xpos = 0
            ypos = bbox[3] + padding
    path = os.path.dirname(output_image)
    outimage.save( open( output_image, 'w' ) )
    # print 'done.'

def poster_h2_pval(indir, output):
    """
    Parameters
    indir: directory containing h2 and pval snapshots
    output: filename of the poster output
    """
    outdir = os.path.dirname(output)
    if not os.path.isdir(outdir):
        os.makedirs(outdir)
    padding = 6
    # List containing the snapshots filename ordered for the poster display
    images = []
    view = ['extern', 'intern']
    for ft in ['h2', 'pval']:
        for sd in ['L', 'R']:
            for v in view:
                image = os.path.join(indir, 'snapshot_'+ft+'_'+sd+'_'+v+'.png')
                images.append(image)
        v = 'bottom'
        for sd in ['L', 'R']:
            image = os.path.join(indir, 'snapshot_'+ft+'_'+sd+'_'+v+'.png')
            images.append(image)
        v = 'top'
        for sd in ['R', 'L']:
            image = os.path.join(indir, 'snapshot_'+ft+'_'+sd+'_'+v+'.png')
            images.append(image)
    # Create a poster with the snapshots
    create_poster(images,  output, padding=padding, tile=8, align_grid=False)

    # Code section to add the colorbars to the image
    # The colorbars need to have previously been saved in a file
    # Unfortunately Anatomist does not enable to do it through line of code
    # So it needs to be done manually
    image = PIL.Image.open(output)
    fts = ['herit45', 'pval4']
    fts = ['cbar_herit01_05', 'cbar_pval3_10']
    colorbar = {}
    max_large = 0
    for ft in fts:
        colorbar[ft] = PIL.Image.open(DIR_CB+'/PALETTES/'+
                                      ft+'.png')
        if int(colorbar[ft].size[0]) > max_large:
            max_large = int(colorbar[ft].size[0])
        colorbar[ft] = colorbar[ft].resize((int(colorbar[ft].size[0]),
                                            int((image.size[1]-
                                                 (len(fts)-1)*padding)
                                                /len(fts))))
    bgcolor = (255,255,255)
    size =( image.size[0]+padding+max_large, 
            image.size[1])
    outimage = PIL.Image.new( 'RGB', size, bgcolor)
    bbox = (0, 0, image.size[0], image.size[1])
    outimage.paste( image, bbox )
    for k,ft in enumerate(fts):
        if ft != None:
            bbox = (image.size[0]+padding,
                    int(image.size[1]*k/len(fts))+k*padding, size[0],
                    int(image.size[1]*k/len(fts))+k*padding+
                    int((image.size[1]-(len(fts)-1)*padding)/len(fts)))
            outimage.paste(colorbar[ft], bbox)

    path = os.path.dirname(output)
    outimage.save(open(output, 'w'))

if __name__ == '__main__':
    # Careful when running this script with colorbar interval between 0 and 1
    # You must use brainvisa bug fix and not the current release as of 06/2017
    tasks = ['LANGUAGE']
    COPE_NUMS = [1, 2]
    fil = 'pe1'
    median = True

    ROOT_DIR = ""
    indir = os.path.join(ROOT_DIR,
                         '2017_HCP_tfMRI_LANGUAGE',
                         'dictionaries_heritability')

    outdir = os.path.join(ROOT_DIR,
                          '2017_HCP_tfMRI_LANGUAGE',
                          'snapshots')
    outdir_poster = os.path.join(ROOT_DIR,
                                 '2017_HCP_tfMRI_LANGUAGE',
                                 'poster_results')

    group_path = os.path.join(ROOT_DIR,
                              '2017_HCP_tfMRI_LANGUAGE',
                              'HCP_MMP1.0')

    # Define colorbar names for h2 and pval respectively
    cbar = ['YlOrRd', 'RdYlGn']
    interval_cb = [[0.1,0.5], [3, 10]]
    # Divide by nb areals (58)* nb hemispheres (2) to get Bonferroni threshold
    PVAL_THRESH = 5e-2/(180*2)


    if median:
        indir = os.path.join(indir, fil+'_median_value')
        outdir = os.path.join(outdir, fil+'_median_value')
        outdir_poster = os.path.join(outdir_poster, fil+'_median_value')
    else:
        indir = os.path.join(indir, fil+'_mean_value')
        outdir = os.path.join(outdir, fil+'_mean_value')
        outdir_poster = os.path.join(outdir_poster, fil+'_mean_value')


    for i, task in enumerate(tasks):
        for COP in COPE_NUMS:
            indir2 = os.path.join(indir, task+'_'+str(COP))
            outdir2 = os.path.join(outdir, task+'_'+str(COP))
            indir3 = outdir2
            for sd in ['L', 'R']:
                areals = os.path.join(group_path,
                                      sd+'.fsaverage164k.label.gii')
                template  = os.path.join(ROOT_DIR, 'folder_gii',
                                         sd.lower()+'h.inflated.white.gii')
                map_heritability(template, sd, areals, indir2, cbar,
                                      interval_cb, PVAL_THRESH, outdir2)



            output = os.path.join(outdir_poster,
                                  task+'_'+str(COP)+'_poster.png')
            # Nicely crop the previous png snapshots
            os.system('ls '+indir3+'/snapshot*.png|while read input;'
                      'do convert  $input -trim /tmp/toto.png;'
                      'convert  /tmp/toto.png -transparent white $input;'
                      'done')
            # Remove intermediate file
            os.system('rm /tmp/toto.png')
            # Create the poster with snapshots and colorbars
            # Need to have saved the colorbar snapshots
            # manually before see code above
            poster_h2_pval(indir3, output)
