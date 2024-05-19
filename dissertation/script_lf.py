#Final Year Project - Making Galaxies - Peter Thomas
#Exercise: produce a plot of the stellar-to-halo mass ratio for the galaxies 
#(z=0)

# Python 3 script to plot stellar mass function

import matplotlib
from matplotlib.pyplot import *
import numpy as np

# Magnitude that is desired: 0 - u, 1 - g, 2 - r, 3 - i, 4 - v
mag_type=2

# Parameters for Hen14
hubble=0.673

# MR
boxside_MR=480.28  # Units Mpc/h
datafile_MR='data/Hen14_sfh2/snap_58.pkl'

# Define limits of plot
# Because I want to use these as plotting data, they need to be arrays,
# not lists.  The following seems to work
xrange=np.array([-27,-17])
binperdex=5
nbin=binperdex*(xrange[1]-xrange[0])

#--------------------------------------------------------------------

if mag_type==0: Xlabel='u'
if mag_type==1: Xlabel='g'
if mag_type==2: Xlabel='r'
if mag_type==3: Xlabel='i'
if mag_type==4: Xlabel='z'
Ylabel=r'$N/$dex$\,Mpc^3$'
Title='Luminosity function'
pngfile='lf.png'
