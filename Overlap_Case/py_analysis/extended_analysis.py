# -*- coding: utf-8 -*-
"""
Created on Wed Oct 13 10:26:36 2021

@author: fkurtog
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
#from fury import window, actor, utils, primitive, io, ui
#from fury.data import read_viz_textures, fetch_viz_textures
import itertools
# import vtk
import glob
import time
import random
import scipy.io as sio
import xml.etree.ElementTree as ET
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

saving_times = np.array([0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 20.0])


main_path = Path(os.getcwd()).parent
out_path = os.path.join(main_path, "output")

os.chdir(out_path)

time_point = "output000000"
number_of_frames = len(saving_times)


fine_data_before_diffusion = sio.loadmat(time_point + "_microenvironment0_before_diffusion.mat")['multiscale_microenvironment']
fine_data_before_diffusion_oxy = fine_data_before_diffusion[4,np.where(fine_data_before_diffusion[2,:] == 16)]


coarse_before_coarsening = sio.loadmat(time_point + "_microenvironment1_before_coarsening.mat")['multiscale_microenvironment']
coarse_before_coarsening_oxy = coarse_before_coarsening[4,:]


coarse_after_coarsening = sio.loadmat(time_point + "_microenvironment1_after_coarsening.mat")['multiscale_microenvironment']
coarse_after_coarsening_oxy = coarse_after_coarsening[4,:]

coarse_after_diffusion = sio.loadmat(time_point + "_microenvironment1_after_diffusion.mat")['multiscale_microenvironment']
coarse_after_diffusion_oxy = coarse_after_diffusion[4,:]

coarse_after_diffusion_with_DC = sio.loadmat(time_point + "_microenvironment1_after_diffusion_with_DC.mat")['multiscale_microenvironment']
coarse_after_diffusion_with_DC_oxy = coarse_after_diffusion_with_DC[4,:]


fine_after_projection = sio.loadmat(time_point + "_microenvironment1_after_coarsening.mat")['multiscale_microenvironment']
fine_after_projection_oxy = fine_after_projection[4,:]


#%%
# Coarse Diffusion Total
tot_coarse_bef_dif = np.sum(coarse_after_coarsening_oxy)
tot_coarse_aft_dif = np.sum(coarse_after_diffusion_oxy)


#%% 
# Dirichlet BC update
tot_coarse_aft_dif_DC = np.sum(coarse_after_diffusion_with_DC_oxy)



#%%
# 