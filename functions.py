# -*- coding: utf-8 -*-
"""
Created on Tue Mar 16 10:05:39 2021

@author: arger
Just a file for holding extra functions. Makes the other scripts more compact

"""
import glob,os,time,h5py,scipy,re
import numpy as np




def strip_metadata_from_input(input_file_location):
    with open(input_file_location, "r") as infile:
        text = infile.read()
        infile.close()
    # make a copy of the text sans comments for easy parsing
    txt = "\n".join([x.split("#")[0] for x in text.split("\n")])
    # find the celldim line and extract the values from it (this is a really dumb way
    # to do this btw, but it works. google Regex for the 'right' method)
    x, y, z = np.array(
        re.sub(" +", " ", txt.split("CellDim")[1].split("\n")[0]).split(" ")[1:]
    ).astype(np.int)
    n = x * y * z
    # add more data grabbing here if desired

    # now make dictionary
    Meta_Dictionary = {"x": x, "y": x, "z": x, "xyz": np.array([x, y, z]), "n": n}
    return Meta_Dictionary


def get_xyz(MD_dict, dimensions, location):
    if type(MD_dict) == dict:
        try:
            dimensions = MD_dict["xyz"]
        except:
            print(
                """Import from provided dict failed, skipping to next
method for determining dimensions"""
            )
    if dimensions == "infer":
        input_file = glob.glob(location + "/*.in")
        assert (
            len(input_file) == 1
        ), """
error: more than one input file in folder. Either remove duplicate input files,
or explicitly pass in xyz dimensions"""
        MD_dict = strip_metadata_from_input(input_file[0])
        dimensions = MD_dict["xyz"]
    x, y, z = dimensions
    return (x, y, z)


def sig6_binary_to_array(location, prefix="sig_", dimensions="infer", MD_dict="N/A"):
    assert os.path.isdir(location), " Chosen Folder path does not exist"
    filelist = np.array(glob.glob("{}/{}*.iout*".format(location, prefix)))
    assert (
        len(filelist) == 6
    ), """
Whoops, looks like you dont have exactly 6 stress file binaries! Double check
to make sure you have correctly downloaded all files, that your location variable
is correct, that that you passed in a binary prefix unique to the sigma files
in question."""
    filelist.sort()
    # Now a little helper function for getting dimensions of simulation.
    x, y, z = get_xyz(MD_dict, dimensions, location)
    # this is the actual sigma file reader. reads binaries into 6 numpy arrrays
    sig = [np.frombuffer(open(file, "rb").read(), np.double) for file in filelist]
    # This reshapes sig according to x,y,z
    sigmas = [s.reshape(z, y, x) for s in sig]
    return sigmas

def sig6_binary_to_array_zxy(location, prefix="sig_", dimensions="infer", MD_dict="N/A"):
    assert os.path.isdir(location), " Chosen Folder path does not exist"
    filelist = np.array(glob.glob("{}/{}*.iout*".format(location, prefix)))
    assert (
        len(filelist) == 6
    ), """
Whoops, looks like you dont have exactly 6 stress file binaries! Double check
to make sure you have correctly downloaded all files, that your location variable
is correct, that that you passed in a binary prefix unique to the sigma files
in question."""
    filelist.sort()
    # Now a little helper function for getting dimensions of simulation.
    x, y, z = get_xyz(MD_dict, dimensions, location)
    # this is the actual sigma file reader. reads binaries into 6 numpy arrrays
    sig = [np.frombuffer(open(file, "rb").read(), np.double) for file in filelist]
    # This reshapes sig according to x,y,z
    sigmas = [s.reshape(z, x, y) for s in sig]
    return sigmas
def sig6_binary_to_array_xzy(location, prefix="sig_", dimensions="infer", MD_dict="N/A"):
    assert os.path.isdir(location), " Chosen Folder path does not exist"
    filelist = np.array(glob.glob("{}/{}*.iout*".format(location, prefix)))
    assert (
        len(filelist) == 6
    ), """
Whoops, looks like you dont have exactly 6 stress file binaries! Double check
to make sure you have correctly downloaded all files, that your location variable
is correct, that that you passed in a binary prefix unique to the sigma files
in question."""
    filelist.sort()
    # Now a little helper function for getting dimensions of simulation.
    x, y, z = get_xyz(MD_dict, dimensions, location)
    # this is the actual sigma file reader. reads binaries into 6 numpy arrrays
    sig = [np.frombuffer(open(file, "rb").read(), np.double) for file in filelist]
    # This reshapes sig according to x,y,z
    sigmas = [s.reshape(x, z, y) for s in sig]
    return sigmas
def sig6_binary_to_array_zyx(location, prefix="sig_", dimensions="infer", MD_dict="N/A"):
    assert os.path.isdir(location), " Chosen Folder path does not exist"
    filelist = np.array(glob.glob("{}/{}*.iout*".format(location, prefix)))
    assert (
        len(filelist) == 6
    ), """
Whoops, looks like you dont have exactly 6 stress file binaries! Double check
to make sure you have correctly downloaded all files, that your location variable
is correct, that that you passed in a binary prefix unique to the sigma files
in question."""
    filelist.sort()
    # Now a little helper function for getting dimensions of simulation.
    x, y, z = get_xyz(MD_dict, dimensions, location)
    # this is the actual sigma file reader. reads binaries into 6 numpy arrrays
    sig = [np.frombuffer(open(file, "rb").read(), np.double) for file in filelist]
    # This reshapes sig according to x,y,z
    sigmas = [s.reshape(z, y,x) for s in sig]
    return sigmas

def els_binary_to_array(location, prefix="els_", dimensions="infer", MD_dict="N/A"):
    assert os.path.isdir(location), " Chosen Folder path does not exist"
    filelist = np.array(glob.glob("{}/{}*.iout*".format(location, prefix)))
    assert (
        len(filelist) == 6
    ), """
Whoops, looks like you dont have exactly 6 els binaries! Double check
to make sure you have correctly downloaded all files, that your location variable
is correct, that that you passed in a binary prefix unique to the els files
in question."""
    filelist.sort()
    # Now a little helper function for getting dimensions of simulation.
    x, y, z = get_xyz(MD_dict, dimensions, location)
    # this is the actual els file reader. reads binaries into 6 numpy arrrays
    els = [np.frombuffer(open(file, "rb").read(), np.double) for file in filelist]
    # This reshapes sig according to x,y,z
    elss = [s.reshape(z, y,x) for s in els]
    return elss

def gID_binary_to_array(
    location, prefix="Aus_Grain", dimensions="infer", MD_dict="N/A"
):
    assert os.path.isdir(location), " Chosen Folder path does not exist"
    filelist = np.array(glob.glob("{}/{}*.iout*".format(location, prefix)))
    assert (
        len(filelist) == 1
    ), """
Whoops, there should only be a single gID file, but you have chosen {}""".format(
        len(filelist)
    )
    # Now a little helper function for getting dimensions of simulation.
    x, y, z = get_xyz(MD_dict, dimensions, location)
    # this is the actual els file reader. reads binaries into 6 numpy arrrays
    gID = np.frombuffer(open(filelist[0], "rb").read(), np.int).reshape(z, y, x)
    return gID

def orientation_tex_to_array(
    location, prefix="tex_", dimensions="infer", MD_dict="N/A"
):
    assert os.path.isdir(location), " Chosen Folder path does not exist"
    filelist = np.array(glob.glob("{}/{}*.iout*".format(location, prefix)))
    assert (
        len(filelist) == 1
    )
    file = filelist[0]
    # Now a little helper function for getting dimensions of simulation.
    x, y, z = get_xyz(MD_dict, dimensions, location)
    # this is the actual els file reader. reads binaries into 6 numpy arrrays
    tex = np.frombuffer(open(file, "rb").read(), np.double)
    # This reshapes sig according to x,y,z
    return tex

def input_text_gen(x,y,z,init_name):
    date = "_".join(np.array([i for i in time.localtime()]).astype(str))
    return( """
# ========================================================================== #
# Dummy Script Auto-generated on {}
# ========================================================================== #
#	Materials types, properties, phases
# -------------------------------------------------------------------------- #
N_phases	2

# Set initial cell dimensions and choose init file to load
CellDim    {} {} {}
initial_ms	  {}

# Type of phase-1 and phase-2 (in accordance with the following order)
# 0(xtal) 1(gas) ie, the 1st phase is PhaseID 0 and the 2nd phase is PhaseID 1
Type_phases		0	0

# Elastic constants of phase 1: c11, c12, c44, nu
ElastConst_PhaseI	162700.	119000.	71150.	0.35

# Data file containing info of phase 1
Slip_PhaseI		Cu_p_hard.sx

# Elastic constants of phase 1: c11, c12, c44, nu
# Gas, ignored!
ElastConst_PhaseII	162700.	119000.	71150.	0.35

# Data file containing info of phase 2
Slip_PhaseII	Cu_p_soft.sx


#====================================================
#	Boundary conditions
#----------------------------------------------------ss
# Flags for velocity gradient (0:unknown, 1:known) Voigt notation
# 1 is where we are SETTING the velocity of the deformation
# 0 is where we are letting the code solve for what the deformation should be
VelGrad_BC_Flag		0 0 0 0 0 1


# Velocity gradient
VelGrad_BC		0.0 0.0 0.0 0.0 0.0 1E-3

# Flags for Cauchy stress
# 1 is where we are SETTING the stress boundary condition
# 0 is where we are letting the code solve for what the stress should be
Stress_BC_Flag		1 1 1 1 1 0

# Applied Cauchy stress
# Setting these to zero says no applied stress other than what 
Stress_BC	 1.0 1.0 1.0 1.0 1.0 1.0 

# elasticity BC: 0 -(e_homo=0), 1 -(relaxed bc)
ElastBC		1

#====================================================
#	Simulation setup
#----------------------------------------------------
# Time increment
TimeStep 1

# Total simulation steps
N_steps		10

# Error tolerance
Err		5E-5

# Maximum iteration steps allowed to solve micromechanical problem
IterMax		20

# write fields?(0:No, 1:Yes)	write step interval
PrintControl	1     1

# Update texture flag
Update_Flag		1

# Hardening flag
Hard_Flag		1

# Texture output flag
Tex_Flag		1


#====================================================
#	Dislocation density based constitutive model
#----------------------------------------------------
# A. Ma et al. Acta Mater. 54(2006) 2169
# lattice parameters for each phase [nm]
a0	0.362	0.0
# voxel length used in GND [mm]
#L0	0.03
L0     0.003
# Used to scale L0 to saturate GND
#GND_ScaleRatio	1000.0
GND_ScaleRatio	6000.0
# Activation energy for slip [1E-19J]
Q_slip		20000000000.0	0.0
# Activation energy for climb [1E-19J]
Q_bulk		300000000000.51 0.0
# Initial SSD density [1E12/m^2]
rho_SSD_initial	   0.05
# Temperature [K]
T__K  30
# Constant for passing stress
C_1		0.48	0.0
# Constant for jump width
C_2		3.0	0.0
# Constant for obstacle width
C_3		3.0	0.0
# Constant for lock forming rate [1E6/m]
C_4 	  800.0 	0.0
## Constant for 6thermal annihilation rate 
C_5	16.0	0.0
## Constant for thermal annihilation rate [1/nm]
C_6	1100E0	0.0

## Constant for dipole forming rate [1E-37m^5*s^C_8/s]
C_7		1.25E-8	0.0
#C_7		0.0	0.0
# Constant for nonlinear climb of edge dislocations
C_8		0.24	0.0
# Self interaction
Selfinter	0.125	0.0
# coplanar interaction
Coplanar	0.125	0.0
# cross slip
CrossSlip	0.725	0.0
#CrossSlip	0.8	0.0
# glissile junction
GlissileJunction	0.125 0.0
#GlissileJunction	0.8 0.0
# Hirth loc
HirthLock	0.05	0.0
#HirthLock	0.8	0.0
# Lomer-Cottrell Lock
LomerCottrellLock 0.18 0.0
#LomerCottrellLock	0.8 0.0
#Diffusion coefficient for climb equation
D_0     4E-05
#nucleus_radius
nucleus_radius  4
#====================================================
#	Phase field of recrystallization
#----------------------------------------------------
# System size for Phase-field (not need to be equal to FFDynamic recrystallization of copper polycrystals with different puritiesT size)
CellDim_pf	64 64 64
#CellDim_pf 512 512 1
# Seed f32 GSL random number generator
RandSeed 123
# GB energy [J/m^2]
E_gb	0.625
# GB mob2lity [m^4/(MJ*s)]
M_gb  7.4E-11
# Burgers vector [nm]
bb_len		0.256
# Scaler of stored energy [1]
Scale_Fdeform	0.5
# Characteristic nucleation strength (disl. density difference) for DRX [1E12/m^2]
k_c 320000.0
# Expon0nt in the statistical model for DRX nucleation [1]
alpha 8.5
# Norm3lized GB mobility in phase-field [in grid unit]
# inversely proportional to the # of phase-field steps
#M_bar   12.0
# Coefficient to be multiplied to k_c for re-nucleation
zeta_kc 1.0
# Aging time for DRX re-nucleation
age_drx 1E0
# Scaler used to be multiplied to SSDs for DRX grains after pahse-field
ssdScaler_DRX 0.74
gndScaler_DRX  0.66
sigScaler_DRX	0.66
# Static or dynamic implementation of nucleation? (1: static, 0: dynamic)
Nucl_Static_Flag	1
# Characteristic aging time for newly-formed DRX grains, normalized by TimeStep
tau_DRX	1000.0
# Additional barrier for slip for new DRX grains [J/m^2]
DeltaQ_DRX	0.0
# Additional passing stress for new DRX grains 
DeltaTau_DRX	0.0
#pf_length_scale
pf_length_scale     2.5E-06
#pf_time_step
pf_time_step        0.004

""".format(date,x,y,z,init_name))

def SX_hard_text_gen():
    return("""# Crystal system (XtalSys):
CUBIC

# Crystal axex (XtalAxis):
1.0 1.0 1.0

# Total number of modes listed in the file (N_modes_max):
3

# Number of modes to be used in the calculation (N_modes):
1

# Labels of the modes to be used (iMode):
1

# Mode# 7 -- <110>(111) SLIP (same as Mode#1, used to check PlasticInit())
	nsmx	nrsx	gamd0x	twshx	isectwx
	6		10.		1.0		0.0		0
	tau0xf	tau0xb	tau1x	thet0	thet1
	150.0		150.0		5.		400.	250.
	hselfx	hlatex
	1.0		1.0
	Slip (n-b)
	1 1 -1		0 1 1
	1 1 -1		1 0 1
	1 1 -1		1 -1 0
	1 -1 -1		0 1 -1
	1 -1 -1		1 0 1
	1 -1 -1		1 1 0
# Mode# 1 -- <110>(111) SLIP
	nsmx	nrsx	gamd0x	twshx	isectwx
	12		10.		1.0		0.0		0
	tau0xf	tau0xb	tau1x	thet0	thet1
	150.0		150.0		5.		400.	250.
	hselfx	hlatex
	1.0		1.0
	Slip (n-b)
	1 1 -1		0 1 1
	1 1 -1		1 0 1
	1 1 -1		1 -1 0
	1 -1 -1		0 1 -1
	1 -1 -1		1 0 1
	1 -1 -1		1 1 0
	1 -1 1		0 1 1
	1 -1 1		1 0 -1
	1 -1 1		1 1 0
	1 1 1		0 1 -1
	1 1 1		1 0 -1
	1 1 1		1 -1 0
# Mode# 100 <110>(111) SLIP
	nsmx	nrsx	gamd0x	twshx	isectwx
	12		10.		1.0		0.0		0
	tau0xf	tau0xb	tau1x	thet0	thet1
	9.0		9.0		5.		400.	250.
	hselfx	hlatex
	1.0		1.0
	Slip (n-b)
	1 1 -1		0 1 1
	1 1 -1		1 0 1
	1 1 -1		1 -1 0
	1 -1 -1		0 1 -1
	1 -1 -1		1 0 1
	1 -1 -1		1 1 0
	1 -1 1		0 1 1
	1 -1 1		1 0 -1
	1 -1 1		1 1 0
	1 1 1		0 1 -1
	1 1 1		1 0 -1
	1 1 1		1 -1 0
""")

def SX_soft_text_gen():
    return("""# Crystal system (XtalSys):
CUBIC

# Crystal axex (XtalAxis):
1.0 1.0 1.0

# Total number of modes listed in the file (N_modes_max):
3

# Number of modes to be used in the calculation (N_modes):
1 

# Labels of the modes to be used (iMode):
1 

# Mode# 1 -- <110>(111) SLIP (same as Mode#1, used to check PlasticInit())
	nsmx	nrsx	gamd0x	twshx	isectwx
	12		20.		1.0		0.0		0
	tau0xf	tau0xb	tau1x	thet0	thet1
	14.25		14.25		35.0		90.0	14.0
	hselfx	hlatex2
	1.0	1.5
	Slip (n-b)
            1  1  1        0  1 -1           
	    1  1  1        1  0 -1
	    1  1  1        1 -1  0
	   -1  1  1        0  1 -1
	   -1  1  1        1  0  1
	   -1  1  1        1  1  0
	   -1 -1  1        0  1  1
	   -1 -1  1        1  0  1
	   -1 -1  1        1 -1  0
	    1 -1  1        0  1  1
	    1 -1  1        1  0 -1
	    1 -1  1        1  1  0     
# Mode# 2 -- <112>(111) SLIP
	nsmx	nrsx	gamd0x	twshx	isectwx
	12		10.		1.0		0.0		0
	tau0xf	tau0xb	tau1x	thet0	thet1
	9.0		9.0		5.		400.	250.
	hselfx	hlatex
	1.0		1.0
	Slip (n-b)
	1  1  1       -2  1  1
   1  1  1        1 -2  1
   1  1  1        1  1 -2
  -1  1  1        2  1  1
  -1  1  1       -1 -2  1
  -1  1  1       -1  1 -2
  -1 -1  1        2 -1  1
  -1 -1  1       -1  2  1
  -1 -1  1       -1 -1 -2
   1 -1  1       -2 -1  1
   1 -1  1        1  2  1
   1 -1  1        1 -1 -2
# Mode# 100 <110>(111) SLIP
	nsmx	nrsx	gamd0x	twshx	isectwx
	12		10.		1.0		0.0		0
	tau0xf	tau0xb	tau1x	thet0	thet1
	9.0		9.0		5.		400.	250.
	hselfx	hlatex
	1.0		1.0
	Slip (n-b)
	1 1 -1		0 1 1
	1 1 -1		1 0 1
	1 1 -1		1 -1 0
	1 -1 -1		0 1 -1
	1 -1 -1		1 0 1
	1 -1 -1		1 1 0
	1 -1 1		0 1 1
	1 -1 1		1 0 -1
	1 -1 1		1 1 0
	1 1 1		0 1 -1
	1 1 1		1 0 -1
	1 1 1		1 -1 0

""")

def Create_Simulation(folder_name,init_name,dat,xlen,ylen):
    folder_location ="Simulations/"+folder_name 
    
    #Initial assertion checks to prevent overwriting data
    if len(glob.glob(folder_location)) == 0:
        os.mkdir(folder_location)
    else:
        errortxt = """
ERROR: The chosen folder has completed results in it. This script will not
write to a folder with completed results. To write out folder, either change
the foldername for this function, or move the existing data to a different
location"""
        assert len(glob.glob(folder_location+"/*.iout")) == 0, errortxt
        assert len(glob.glob(folder_location+"/*.iout*")) == 0, errortxt
        assert len(glob.glob(folder_location+"/*.out")) == 0, errortxt
    
    #Write Text files
    with open(folder_location+"/input.in",'w', newline ='') as infile:
        text = input_text_gen(xlen, ylen, 1,init_name)
        infile.write(text)
        infile.close()
    
    with open(folder_location+"/PBS_instructions.txt",'w', newline ='') as infile:
        infile.write("""Code for interactive Batch:
sinteractive -N 1 -n 48 -p largemem -t 02:00:00 -A PAA0023 bigmem_Forge_{}
Alternate for if the largemem nodes are taken/slow:
sinteractive -N 2 -n 48 -p parallel -t 8:00:00 -A PAA0023 large_boi_{}

Command for running forge (works for Austin anyway, may not for you):
mpiexec ~/bin/forge/ ./input.in
""".format(folder_name,folder_name))
        infile.close()        
    
    with open(folder_location+"/Cu_p_soft.sx",'w', newline ='') as infile:
        text = SX_soft_text_gen()
        infile.write(text)
        infile.close()
    
    with open(folder_location+"/Cu_p_hard.sx",'w', newline ='') as infile:
        text = SX_hard_text_gen()
        infile.write(text)
        infile.close()
    
    #Create initial Simulation
    MS_savename = folder_location+"/"+init_name
    print(MS_savename)
    np.savetxt(MS_savename,dat,
           fmt='%0.1f %0.1f %0.1f  %i %i %i  %i %i',
           delimiter= ' ', newline ='\n')
    with open(MS_savename,'rb') as file:
        text = file.read()
        file.close()
    unix_text = b'\n'.join(text.split(b'\r\n'))    
    with open(MS_savename,'wb') as file:
        file.write(unix_text)
    return()
