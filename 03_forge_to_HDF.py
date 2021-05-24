# -*- coding: utf-8 -*-
"""
Created on Fri Mar 19 03:25:16 2021

@author: arger

Once you run forge up on the OSC, this takes the downloaded file folders and
consolidates the binaries into human readable stuff.

NOTE TO OTHERS USING THIS FOR THEIR OWN PURPOSES OUTSIDE OF THIS PITTING MODEL:
    Writing a totally generalized version of this problem is difficult, so 
    instead i am just going to write the important decoders/translators, then
    smash them together in the way i need them. feel free to mix and match 
    whatever for your own uses, but this will not work for everyone as written.
    
Sort of curious if anyone but me reads these info files. Is it just a message to
myself? am i just looking forward in time to speak through text to a future 
version of myself? If so, remember to edit the timestamps so its not so obvious
you wote this at 3am the day before presenting it. 

NOTE TO SELF:conda.black before commit.

Also add asserts
"""

import numpy as np
import glob, os, re,h5py
import functions as fn
from pyevtk.hl import gridToVTK
from zipfile import ZipFile



def folder_to_vtr(location, overwrite=False):
    # NOTE: it makes a lot of sense to do generalized h5 files at some point,
    # but this is a binary VTK, which is essentially equivilant for this basic
    # of a file structure.
    foldername = location.split("/")[-1].split("\\")[-1]
    vtk_fname = "{}/{}.vtr".format(location, foldername)
    if np.all([len(glob.glob(vtk_fname)) != 0, overwrite == False]) == True:
        print( """
vtr file already exists in folder,skipping. Pass 'overwrite = True' into this
function to ignore this warning and overwrite the previous vtk files""")
    else:
        sig = fn.sig6_binary_to_array(location)
#        sig_xyz = fn.sig6_binary_to_array_zxy(location)
#        sig_xzy = fn.sig6_binary_to_array_xzy(location)
        sig_zyx = fn.sig6_binary_to_array_zyx(location)
        els = fn.els_binary_to_array(location)
        gID = fn.gID_binary_to_array(location)
        dct = els + sig+sig_zyx
        dct.append(gID)
        to_VTK = dict(
            zip(
            [
                "e_xx",
                "e_yy",
                "e_zz",
                "e_yz",
                "e_xz",
                "e_zx",
                "sig_xx",
                "sig_yy",
                "sig_zz",
                "sig_yz",
                "sig_xz",
                "sig_zx",
                "Grain_ID",
            ],
            dct,
        )
    )
        x1 = gID.shape[0]+1
        y1 = gID.shape[1]+1
        z1 = gID.shape[2]+1
        gridToVTK(vtk_fname[:-4], np.arange(x1), np.arange(y1), np.arange(z1), cellData=to_VTK)
        
    return ()

#line to do one vtr this way
folder_to_vtr("Simulations/Experiment_01")
# #Line to do all the folders in one loops
# for folder in glob.glob("Simulations/*")[1:]:
#     print("writing vtr file for {}...".format(folder.split("/")[-1]))
#     folder_to_vtr(folder, overwrite= False)
#     print("zipping up iout files...")
#     os.chdir(folder)
#     zf =ZipFile( "iout.zip",'w')
#     for file in glob.glob("*.iout*"):
#         zf.write(file)
#         os.remove(file)
#     zf.close()
#     os.chdir("../../../")
# print("Done\n\n")
