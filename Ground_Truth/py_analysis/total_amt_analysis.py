# -*- coding: utf-8 -*-
"""
Created on Fri Oct 15 0:00 2021

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

saving_times = np.array([0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 20.0])
dt = np.array([0.0, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0])

main_path = Path(os.getcwd()).parent
out_path = os.path.join(main_path, "output")

os.chdir(out_path)

fine_data = sio.loadmat("initial_microenvironment0.mat")['multiscale_microenvironment']
init_total_oxy =sum(fine_data[4,:])*(32*32*32)
init_total_glu = sum(fine_data[5,:])*(32*32*32)
init_total_chem = sum(fine_data[6,:])*(32*32*32)
time_point = "output000000"
number_of_frames = len(saving_times)

total_oxy = []
total_glu = []
total_chem = []

o2 = []
g = []
c = []

total_oxy_diff = []
total_glu_diff = []
total_chem_diff = []

prev_oxy = init_total_oxy
prev_glu = init_total_glu
prev_chem = init_total_chem

o2_uptake_rate_per_cell = 0.005
glu_uptake_rate_per_cell = 0.010
chem_secretion_rate_per_cell_per_min = 0.01
number_of_cells = 1

uptaken_o2_arr = []
uptaken_glu_arr = []

cell_volume = 2494

for i in range(len(saving_times)):
    time_p = time_point + '%02d'%(i)
    if path.exists(time_p + "_microenvironment0.mat"):
        fine_data = sio.loadmat(time_p + "_microenvironment0.mat")['multiscale_microenvironment']
        oxy = sum(fine_data[4,:])*(32*32*32)
        o2.append(oxy)
        glu = sum(fine_data[5,:])*(32*32*32)
        g.append(glu)
        chem = sum(fine_data[6,:])*(32*32*32)
        c.append(chem)
        uptaken_o2 = o2_uptake_rate_per_cell * number_of_cells * cell_volume * saving_times[i] * 0.285
        uptaken_glu = glu_uptake_rate_per_cell * number_of_cells * cell_volume * saving_times[i] * 16.897255
        uptaken_o2_arr.append(uptaken_o2)
        uptaken_glu_arr.append(uptaken_glu)
        total_oxy_diff.append(init_total_oxy - oxy)
        prev_oxy = oxy
        total_glu_diff.append(init_total_glu - glu)
        prev_glu = glu
        total_chem_diff.append(init_total_chem - chem)
        prev_chem = chem
        amt = oxy + uptaken_o2
        total_oxy.append(amt)
        total_glu.append(glu + uptaken_glu)
        total_chem.append(chem)

total_oxy = np.asarray(total_oxy)
total_glu = np.asarray(total_glu)
total_chem = np.asarray(total_chem)
uptaken_o2_arr = np.asarray(uptaken_o2_arr)
print(total_oxy_diff)
print("_______________")
print(uptaken_o2_arr)
print("_______________")
print(uptaken_o2_arr/total_oxy_diff)
print("_______________")
print(total_chem)
#uptaken_o2_arr = np.asarray(uptaken_o2_arr)
#plt.plot(saving_times, uptaken_o2_arr, c = "blue")
#plt.plot(saving_times, total_oxy_diff, c = "red")
#plt.show()
#uptaken_glu_arr = np.asarray(uptaken_glu_arr)*(2880*2880*5360)
#plt.plot(saving_times, uptaken_glu_arr, c = "orange")
#plt.show()
plt.plot(saving_times, total_oxy, c = "blue")
plt.title("Total Oxygen (Mass Conservation)")
plt.xlabel("Time (min)")
plt.ylabel("Total Amount (mmol)")
plt.show()
plt.plot(saving_times, total_glu, c = "orange")
plt.title("Total Glucose (Mass Conservation)")
plt.xlabel("Time (min)")
plt.ylabel("Total Amount (mmol)")
plt.show()
plt.plot(saving_times, total_chem, c = "green")
plt.title("Total Chemokine (Mass Conservation)")
plt.xlabel("Time (min)")
plt.ylabel("Total Amount (mmol)")
plt.show()
plt.plot(saving_times, o2, c = "blue")
plt.title("Total Microenvironment Oxygen")
plt.xlabel("Time (min)")
plt.ylabel("Total Amount (mmol)")
plt.show()
plt.plot(saving_times, g, c = "orange")
plt.title("Total Microenvironment Glucose")
plt.xlabel("Time (min)")
plt.ylabel("Total Amount (mmol)")
plt.show()
plt.plot(saving_times, c, c = "green")
plt.title("Total Microenvironment Chemokine")
plt.xlabel("Time (min)")
plt.ylabel("Total Amount (mmol)")
plt.show()
