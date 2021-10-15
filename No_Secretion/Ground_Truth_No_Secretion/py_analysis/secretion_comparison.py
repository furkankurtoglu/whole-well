# -*- coding: utf-8 -*-
"""
Created on Wed Oct 13 10:19 2021

@author: Kali
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
import itertools
import vtk
import glob
import time
import random
import scipy.io as sio
import xml.etree.ElementTree as ET
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

saving_times = np.array([0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0])


main_path = Path(os.getcwd()).parent
out_path = os.path.join(main_path, "output")

os.chdir(out_path)

time_point = "output000000"
number_of_frames = len(saving_times)

total_o2_before_secretion = []
total_glu_before_secretion = []
total_chem_before_secretion = []

for i in range(len(saving_times)):
    time_p = time_point + '%02d'%(i)
    if path.exists(time_p + "_microenvironment0_before_secretion.mat"):
        fine_data = sio.loadmat(time_p + "_microenvironment0_before_secretion.mat")['multiscale_microenvironment']
        micEnv_O2_before =sum(fine_data[4,:])
        micEnv_glu_before = sum(fine_data[5,:])
        micEnv_chem_before = sum(fine_data[6,:])
        total_o2_before_secretion.append(micEnv_O2_before)
        total_glu_before_secretion.append(micEnv_glu_before)
        total_chem_before_secretion.append(micEnv_chem_before)

total_o2_after_secretion = []
total_glu_after_secretion = []
total_chem_after_secretion = []

for i in range(len(saving_times)):
    time_p = time_point + '%02d'%(i)
    if path.exists(time_p + "_microenvironment0_after_secretion.mat"):
        fine_data = sio.loadmat(time_p + "_microenvironment0_after_secretion.mat")['multiscale_microenvironment']
        micEnv_O2_after =sum(fine_data[4,:])
        micEnv_glu_after = sum(fine_data[5,:])
        micEnv_chem_after = sum(fine_data[6,:])
        total_o2_after_secretion.append(micEnv_O2_after)
        total_glu_after_secretion.append(micEnv_glu_after)
        total_chem_after_secretion.append(micEnv_chem_after)

total_o2_before_secretion = np.asarray(total_o2_before_secretion)*(2880*2880*5360)
total_glu_before_secretion = np.asarray(total_glu_before_secretion)*(2880*2880*5360)
total_chem_before_secretion = np.asarray(total_chem_before_secretion)*(2880*2880*5360)

total_o2_after_secretion = np.asarray(total_o2_after_secretion)*(2880*2880*5360)
total_glu_after_secretion = np.asarray(total_glu_after_secretion)*(2880*2880*5360)
total_chem_after_secretion = np.asarray(total_chem_after_secretion)*(2880*2880*5360)

oxy_difference = total_o2_after_secretion - total_o2_before_secretion
glu_difference = total_glu_after_secretion - total_glu_before_secretion
chem_difference = total_chem_after_secretion - total_chem_before_secretion

plt.plot(saving_times[1:], total_o2_after_secretion, c = "blue")
plt.show()
plt.plot(saving_times[1:], total_glu_after_secretion, c = "orange")
plt.show()
plt.plot(saving_times[1:], total_chem_after_secretion, c = "green")
plt.show()