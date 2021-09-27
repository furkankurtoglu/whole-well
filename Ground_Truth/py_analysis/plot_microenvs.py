# -*- coding: utf-8 -*-
"""
Created on Fri Aug 27 14:41:55 2021

@author: Furkan
"""


import importlib.machinery

pyMCDS = importlib.machinery.SourceFileLoader('pyMCDS','./analysis/pyMCDS.py').load_module()

import os.path
from os import path

from pathlib import Path
import matplotlib.pyplot as plt
import matplotlib.animation

import numpy as np
import pandas as pd
from fury import window, actor, utils, primitive, io, ui
from fury.data import read_viz_textures, fetch_viz_textures
import itertools
import vtk
import glob
import time
import random
import scipy.io as sio
import xml.etree.ElementTree as ET
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

saving_times = np.array([0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 20.0])

def data_parser (time_point):
    # Fine MicroEnv Data Parsing
    fine_tuple = []
    coarse_tuple = []
    transfer_tuple = []
    if path.exists(time_point + "_microenvironment0.mat"):
        fine_data = sio.loadmat(time_point + "_microenvironment0.mat")['multiscale_microenvironment']
        dx = fine_data[0,1]-fine_data[0,0]
        fine_x = np.unique(fine_data[0,:])
        fine_y = np.unique(fine_data[1,:])
        fine_X, fine_Y = np.meshgrid(fine_x, fine_y)
        fine_oxy = fine_data[4,np.where(fine_data[2,:] == 16)]
        fine_oxy = fine_oxy.reshape((len(fine_y),len(fine_x)))
        fine_glu = fine_data[5,np.where(fine_data[2,:] == 16)]
        fine_glu = fine_glu.reshape((len(fine_y),len(fine_x)))
        fine_chem = fine_data[6,np.where(fine_data[2,:] == 16)]
        fine_chem = fine_chem.reshape((len(fine_y),len(fine_x)))
        fine_oxy_tuple = (fine_X, fine_Y, fine_oxy)
        fine_glu_tuple = (fine_X, fine_Y, fine_glu)
        fine_chem_tuple = (fine_X, fine_Y, fine_chem)
        fine_tuple = (fine_oxy_tuple, fine_glu_tuple, fine_chem_tuple)
    
    # Coarse MicroEnv Data Parsing
    if path.exists(time_point + "_microenvironment1.mat"):
        coarse_data = sio.loadmat(time_point + "_microenvironment1.mat")['multiscale_microenvironment']


    if path.exists(time_point + "_microenvironment2.mat"):
        coarse_data = sio.loadmat(time_point + "_microenvironment2.mat")['multiscale_microenvironment']

    return fine_tuple, coarse_tuple, transfer_tuple
    



def get_subs_name ():
    tree = ET.parse("initial.xml")
    root = tree.getroot()
    
    subs_names = []
    for substrate in root.iter('variable'):
        subs_names.append(substrate.attrib['name'])
    
    return subs_names
        

    

main_path = Path(os.getcwd()).parent
out_path = os.path.join(main_path, "output")

os.chdir(out_path)

time_point = "output000000"
number_of_frames = len(saving_times)

subs_list = get_subs_name()

#%% Temporospatial Plot
fig, axs = plt.subplots()

# color bar
tp = "final"
ft, ct, tt = data_parser(tp)
fine_X, fine_Y, fine_oxy = ft[0]
zmin = min([min(zl) for zl in fine_oxy])
zmax = max([max(zl) for zl in fine_oxy])
levels = np.linspace(zmin, 38,41)
kw = dict(levels=levels, vmin=zmin, vmax=38.000, origin='lower')
cp = axs.contourf(fine_Y,fine_X,fine_oxy, **kw)
cbar = plt.colorbar(cp,format='%0.2f')
axs.clear()



def animate(i):
    time_p= time_point + '%02d'%(i)
    ft, ct, tt = data_parser(time_p)
    fine_X, fine_Y, fine_oxy = ft[0]
    axs.clear()
    axs.contourf(fine_Y,fine_X,fine_oxy, **kw)
    axs.set_title('Oxygen, Z=16 um, time = ' +str(saving_times[i])+ ' minutes') 
    axs.invert_xaxis()
    axs.axis('scaled')
    
    
    
number_of_frames = len(saving_times)

ani = matplotlib.animation.FuncAnimation(fig,animate,blit=False, frames=number_of_frames,repeat=False)

plt.show()

ani.save('./oxygen.gif', writer='imagemagick', fps=4)



fig2, ax = plt.subplots()

# color bar
tp = "final"
ft, ct, tt = data_parser(tp)
fine_X, fine_Y, fine_glu = ft[1]
zmin2 = min([min(zl) for zl in fine_glu])
zmax2 = max([max(zl) for zl in fine_glu])
levels2 = np.linspace(zmin2, 16.897255)
kw2 = dict(levels=levels2, vmin=zmin2, vmax=16.897255, origin='lower')
cp2 = ax.contourf(fine_Y,fine_X,fine_glu, **kw2)
cbar2 = plt.colorbar(cp2,format='%0.2f')
ax.clear()

def animate2(i):
    time_p= time_point + '%02d'%(i)
    ft, ct, tt = data_parser(time_p)
    fine_X, fine_Y, fine_glu = ft[1]
    ax.clear()
    ax.contourf(fine_Y,fine_X,fine_glu, **kw2)
    ax.set_title('Glucose, Z=16 um, time = ' +str(saving_times[i])+ ' minutes') 
    ax.invert_xaxis()
    ax.axis('scaled')
    


ani2 = matplotlib.animation.FuncAnimation(fig2,animate2,blit=False, frames=number_of_frames,repeat=False)

plt.show()

ani2.save('./glucose.gif', writer='imagemagick', fps=4)






fig3, ax3 = plt.subplots()

# color bar
tp = "final"
ft, ct, tt = data_parser(tp)
fine_X, fine_Y, fine_chem = ft[2]
zmin3 = min([min(zl) for zl in fine_chem])
zmax3 = max([max(zl) for zl in fine_chem])
levels3 = np.linspace(0, zmax3)
kw3 = dict(levels=levels3, vmin=0, vmax=zmax3, origin='lower')
cp3 = ax3.contourf(fine_Y,fine_X,fine_chem, **kw3)
cbar3 = plt.colorbar(cp3,format='%0.5f')
ax3.clear()

def animate3(i):
    time_p= time_point + '%02d'%(i)
    ft, ct, tt = data_parser(time_p)
    fine_X, fine_Y, fine_chem = ft[2]
    ax3.clear()
    ax3.contourf(fine_Y,fine_X,fine_chem, **kw3)
    ax3.set_title('Chemokine, Z=16 um, time = ' +str(saving_times[i])+ ' minutes') 
    ax3.invert_xaxis()
    ax3.axis('scaled')
    


ani3 = matplotlib.animation.FuncAnimation(fig3,animate3,blit=False, frames=number_of_frames,repeat=False)

plt.show()

ani3.save('./chemokine.gif', writer='imagemagick', fps=4)
#%% Temporal Plot

# fig, axs = plt.subplots()
# time_point = "output00000001"

# fine_data = sio.loadmat(time_point + "_microenvironment0.mat")['multiscale_microenvironment']
# dx = fine_data[0,1]-fine_data[0,0]
# fine_x = np.unique(fine_data[0,:])
# fine_y = np.unique(fine_data[1,:])
# fine_X, fine_Y = np.meshgrid(fine_x, fine_y)
# fine_oxy = fine_data[4,np.where(fine_data[2,:] == 16)]
# fine_oxy = fine_oxy.reshape((len(fine_y),len(fine_x)))
# fine_tuple = (fine_X, fine_Y, fine_oxy)

# ft, ct, tt = data_parser(time_point)
# fine_X, fine_Y, fine_oxy = ft
# zmin = min([min(zl) for zl in fine_oxy])
# zmax = max([max(zl) for zl in fine_oxy])
# levels = np.linspace(zmin, 38,41)
# kw = dict(levels=levels, vmin=zmin, vmax=38.000, origin='lower')
# cp = axs.contourf(fine_Y,fine_X,fine_oxy, **kw)
# cbar = plt.colorbar(cp)
# plt.show()

    


