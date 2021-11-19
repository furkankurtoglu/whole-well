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

time_point = "output00000016"
number_of_frames = len(saving_times)

#%% Fine Microenvironment Diffusion
fine_data_before_diffusion = sio.loadmat(time_point + "_microenvironment0_before_diffusion.mat")['multiscale_microenvironment']
aa_fine_data_before_diffusion_glu = np.sum(fine_data_before_diffusion[5,:])

fine_data_after_diffusion = sio.loadmat(time_point + "_microenvironment0_after_diffusion.mat")['multiscale_microenvironment']
ab_fine_data_after_diffusion_glu = np.sum(fine_data_before_diffusion[5,:])


#%% 
coarse_before_coarsening = sio.loadmat(time_point + "_microenvironment1_before_coarsening.mat")['multiscale_microenvironment']
ba_coarse_before_coarsening_glu = np.sum(coarse_before_coarsening[5,:])


coarse_after_coarsening = sio.loadmat(time_point + "_microenvironment1_after_coarsening.mat")['multiscale_microenvironment']
bb_coarse_after_coarsening_glu = np.sum(coarse_after_coarsening[5,:])

coarse_after_diffusion = sio.loadmat(time_point + "_microenvironment1_after_diffusion.mat")['multiscale_microenvironment']
bc_coarse_after_diffusion_glu = np.sum(coarse_after_diffusion[5,:])

coarse_after_diffusion_with_DC = sio.loadmat(time_point + "_microenvironment1_after_diffusion_with_DC.mat")['multiscale_microenvironment']
bd_coarse_after_diffusion_with_DC_glu = np.sum(coarse_after_diffusion_with_DC[5,:])


fine_before_projection = sio.loadmat(time_point + "_microenvironment0_before_projection.mat")['multiscale_microenvironment']
ac_fine_before_projection_glu = np.sum(fine_before_projection[5,:])

fine_after_projection = sio.loadmat(time_point + "_microenvironment0_after_projection.mat")['multiscale_microenvironment']
ad_fine_after_projection_glu = np.sum(fine_after_projection[5,:])


#%%
fine_gain = ad_fine_after_projection_glu - ac_fine_before_projection_glu
coarse_lost = ba_coarse_before_coarsening_glu - bb_coarse_after_coarsening_glu
Net_Difference = (fine_gain - coarse_lost) * 2880 *2880 * 512