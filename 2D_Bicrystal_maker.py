# -*- coding: utf-8 -*-
"""
Created on Mon Apr 19 12:47:55 2021

@author: arger
===============================================================================
Script for making 2D bicrystals. Borrows heavily from Pitting model made for
Alana.
===============================================================================

"""

import numpy as np
import functions as fn

# =========================================================================== #
#    --------------   BEGINNING OF USER-DEFINED VARIABLES   --------------    #
# =========================================================================== #
# Simulation Dimensions
xdim,ydim =256,256
Bicrystal_1_length = 120
Bicrystal_1_ori = np.array([1,1,1]) # in degrees
Bicrystal_2_ori = np.array([30,60,90])# in degrees
# =========================================================================== #
#    -----------------   END OF USER-DEFINED VARIABLES   -----------------    #
# =========================================================================== #
    
    
#Step 1: create initial simulation
xx,yy = np.meshgrid(np.arange(xdim),np.arange(ydim))
xtal_mask = yy<Bicrystal_1_length

Eul_1 = (xtal_mask*Bicrystal_1_ori[0])+((xtal_mask == 0)*Bicrystal_2_ori[0])
Eul_2 = (xtal_mask*Bicrystal_1_ori[1])+((xtal_mask == 0)*Bicrystal_2_ori[1])
Eul_3 = (xtal_mask*Bicrystal_1_ori[2])+((xtal_mask == 0)*Bicrystal_2_ori[2])

i_loc = yy +1
j_loc = xx +1
k_loc = (xx*0)+1

Grain_ID = xtal_mask+1
Phase_ID = Grain_ID*1

# In the steps above, I made each of the colums we need as 2D arrays. this next
# line takes all those variables and unravels them into columns, then
# stacks them next to each other into the 2D microstructure file you need.
dat = np.vstack([Eul_1.ravel(),
                 Eul_2.ravel(),
                 Eul_3.ravel(),
                 i_loc.ravel(),
                 j_loc.ravel(),
                 k_loc.ravel(),
                 Grain_ID.ravel(),
                 Phase_ID.ravel()]).T

# comment out these next two lines if you want a truly 2D simulation with no 
# air cushion to prevent wacky stress/strain stuff
#air_cushion = dat *np.array([0,0,0,0,1,1,0,0])+np.array([0,0,0,2,0,0,3,2])
#dat = np.vstack([dat,air_cushion])

init_savename = "input.init"


# ========================================================================== #
# ------------ Function that actually writes out the experiment ------------ #
# ========================================================================== #

fn.Create_Simulation('Example','initial_MS.init',dat,xdim,ydim)
