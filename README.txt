Basic How-To:

2D_Bicrtystal_maker.py will create a 256x256 simulation for you. you can open 
the file up in either spyder or just a text editor, and change several of the 
parameters (orientation, size, dividing line between grain, etc). If you install
anaconda in your OSC folder, you can also run it there as well. 

The generated experiment is in Simulations/(The-name-you-gave-it-in-2D_Bicrystal_maker.py)
put that file in an OSC folder, run forge, then download the result back to your computer.

once you have the results with all the iout files included, run forge_to_hdf.py, which
will generate a vtk file. THAT file is readable by paraview natively. also, since your data
is 2D, you can also just read that vtk (which is just an hdf5 under a different name) into 
any software (Matlab, python, jupyter, etc) and look at it using image viewing software


-Austin 