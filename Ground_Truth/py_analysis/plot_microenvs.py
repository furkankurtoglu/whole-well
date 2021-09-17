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
        fine_oxy = fine_data[5,np.where(fine_data[2,:] == 16)]
        fine_oxy = fine_oxy.reshape((len(fine_y),len(fine_x)))
        fine_tuple = (fine_X, fine_Y, fine_oxy)
    
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


subs_list = get_subs_name()




fig, axs = plt.subplots()


def animate(i):
    time_p= time_point + '%02d'%(i)
    ft, ct, tt = data_parser(time_p)
    fine_X, fine_Y, fine_oxy = ft
    axs.clear()
    axs.contourf(fine_X,fine_Y,fine_oxy)
    axs.set_title('%02d'%(i)) 
    axs.axis('equal')

ani = matplotlib.animation.FuncAnimation(fig,animate,blit=False, frames=20,repeat=False)

plt.show()







# def plot_micenvs (time_point):
#     fine_data = sio.loadmat(time_point + "_microenvironment0.mat")['multiscale_microenvironment']
#     # coarse_data = sio.loadmat(time_point + "_microenvironment1.mat")['multiscale_microenvironment']
#     dx = fine_data[0,1]-fine_data[0,0]
#     fine_x = np.unique(fine_data[0,:])
#     fine_y = np.unique(fine_data[1,:])
#     fine_X, fine_Y = np.meshgrid(fine_x, fine_y)
#     fine_oxy = fine_data[4,np.where(fine_data[2,:] == 16)]
#     fine_oxy = fine_oxy.reshape((len(fine_y),len(fine_x)))

#     fig,axs = plt.subplots(1,1)
#     cp = axs.contourf(fine_X,fine_Y,fine_oxy)
#     axs.axis('equal')
#     axs.set(xlim=(-3000,3000), ylim=(-500,500))
#     fig.colorbar(cp,format = "%f")
#     fig.tight_layout()
#     return fine_data





# main_path = Path(os.getcwd()).parent
# out_path = os.path.join(main_path, "output")

# os.chdir(out_path)

# time_point = "output0000000"
# # data = plot_micenvs(time_point)


# # data1 = data[4,np.where(data[0,:]) == -2864]
# # print(data1)

# for i in range(0,10):
#     timepoint=time_point+str(i)
#     plot_micenvs(timepoint)
    
