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

Temporospatial_Plotting = 'N'
Total_Amount_Analysis = 'Y'





if Temporospatial_Plotting == 'Y':
    def data_parser (time_point):
        # Fine MicroEnv Data Parsing
        fine_tuple = []
        coarse_tuple = []
        transfer_tuple = []
       
        if path.exists(time_point + "_microenvironment0.mat"):
            fine_data = sio.loadmat(time_point + "_microenvironment0.mat")['multiscale_microenvironment']
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
            coarse_y = coarse_data[0,:]
            coarse_x = np.unique(fine_data[0,:])
            coarse_X, coarse_Y = np.meshgrid(coarse_x, coarse_y)
            coarse_oxy = coarse_data[4,:]
            coarse_oxy = np.transpose(np.tile(coarse_oxy,(90,1)))
            coarse_glu = coarse_data[5,:]
            coarse_glu = np.transpose(np.tile(coarse_glu,(90,1)))
            coarse_chem = coarse_data[6,:]
            coarse_chem = np.transpose(np.tile(coarse_chem,(90,1)))
            coarse_tuple = (coarse_X, coarse_Y, coarse_oxy, coarse_glu, coarse_chem)
            
            
        if path.exists(time_point + "_microenvironment2.mat"):
            transfer_region = sio.loadmat(time_point + "_microenvironment2.mat")['multiscale_microenvironment']
    
        return fine_tuple, coarse_tuple, transfer_tuple
        
    def get_subs_name ():
        tree = ET.parse("initial.xml")
        root = tree.getroot()
        
        subs_names = []
        for substrate in root.iter('variable'):
            subs_names.append(substrate.attrib['name'])
        
        return subs_names
            
    subs_list = get_subs_name()
    
    fig, axs = plt.subplots()
    
    # color bar
    tp = "output00000020"
    ft, ct, tt = data_parser(tp)
    fine_X, fine_Y, fine_oxy = ft[0]
    cX, cY, cOxy, cGlu, cChem = ct
    w_X = np.concatenate((fine_X,cX),axis=0)
    w_Y = np.concatenate((fine_Y,cY),axis=0)
    w_O = np.concatenate((fine_oxy,cOxy),axis=0)
    zmin = min([min(zl) for zl in w_O])
    zmax = max([max(zl) for zl in w_O])
    levels = np.linspace(0, 0.28500001,41)
    kw = dict(levels=levels, vmin=0, vmax=0.28500001, origin='lower')
    cp = axs.contourf(w_Y,w_X,w_O, **kw)
    cbar = plt.colorbar(cp,format='%0.4f')
    axs.clear()
    
    
    
    def animate(i):
        time_p= time_point + '%02d'%(i)
        ft, ct, tt = data_parser(time_p)
        fine_X, fine_Y, fine_oxy = ft[0]
        cX, cY, cOxy, cGlu, cChem = ct
        w_X = np.concatenate((fine_X,cX),axis=0)
        w_Y = np.concatenate((fine_Y,cY),axis=0)
        w_O = np.concatenate((fine_oxy,cOxy),axis=0)
        axs.clear()
        axs.contourf(w_Y,w_X,w_O, **kw)
        axs.set_title('Oxygen, Z=16 um, time = ' +str(saving_times[i])+ ' minutes') 
        axs.invert_xaxis()
        axs.axis('scaled')
        
        
        
        
    number_of_frames = len(saving_times)
    
    ani = matplotlib.animation.FuncAnimation(fig,animate,blit=False, frames=number_of_frames,repeat=False)
    
    plt.show()
    
    ani.save('./oxygen.gif', writer='imagemagick', fps=4)
    
    
    
    fig2, ax = plt.subplots()
    
    # color bar
    tp = "output00000020"
    ft, ct, tt = data_parser(tp)
    fine_X, fine_Y, fine_glu = ft[1]
    cX, cY, cOxy, cGlu, cChem = ct
    w_X = np.concatenate((fine_X,cX),axis=0)
    w_Y = np.concatenate((fine_Y,cY),axis=0)
    w_G = np.concatenate((fine_glu,cGlu),axis=0)
    zmin2 = min([min(zl) for zl in w_G])
    zmax2 = max([max(zl) for zl in w_G])
    levels2 = np.linspace(0, 16.897255)
    kw2 = dict(levels=levels2, vmin=0, vmax=16.897255, origin='lower')
    cp2 = ax.contourf(w_X,w_Y,w_G, **kw2)
    cbar2 = plt.colorbar(cp2,format='%0.2f')
    ax.clear()
    
    def animate2(i):
        time_p= time_point + '%02d'%(i)
        ft, ct, tt = data_parser(time_p)
        fine_X, fine_Y, fine_glu = ft[1]
        cX, cY, cOxy, cGlu, cChem = ct
        w_X = np.concatenate((fine_X,cX),axis=0)
        w_Y = np.concatenate((fine_Y,cY),axis=0)
        w_G = np.concatenate((fine_glu,cGlu),axis=0)
        ax.clear()
        ax.contourf(w_Y,w_X,w_G, **kw2)
        ax.set_title('Glucose, Z=16 um, time = ' +str(saving_times[i])+ ' minutes') 
        ax.invert_xaxis()
        ax.axis('scaled')
        
    
    
    ani2 = matplotlib.animation.FuncAnimation(fig2,animate2,blit=False, frames=number_of_frames,repeat=False)
    
    plt.show()
    
    ani2.save('./glucose.gif', writer='imagemagick', fps=4)
    
    
    
    
    
    
    # fig3, ax3 = plt.subplots()
    
    # # color bar
    # tp = "output00000020"
    # ft, ct, tt = data_parser(tp)
    # fine_X, fine_Y, fine_chem = ft[2]
    # cX, cY, cOxy, cGlu, cChem = ct
    # w_X = np.concatenate((fine_X,cX),axis=0)
    # w_Y = np.concatenate((fine_Y,cY),axis=0)
    # w_C = np.concatenate((fine_chem,cChem),axis=0)
    # zmin3 = min([min(zl) for zl in w_C])
    # zmax3 = max([max(zl) for zl in w_C])
    # levels3 = np.linspace(0, zmax3)
    # kw3 = dict(levels=levels3, vmin=0, vmax=zmax3, origin='lower')
    # cp3 = ax3.contourf(w_X,w_Y,w_C, **kw3)
    # cbar3 = plt.colorbar(cp3,format='%0.5f')
    # ax3.clear()
    
    # def animate3(i):
    #     time_p= time_point + '%02d'%(i)
    #     ft, ct, tt = data_parser(time_p)
    #     fine_X, fine_Y, fine_chem = ft[2]
    #     cX, cY, cOxy, cGlu, cChem = ct
    #     w_X = np.concatenate((fine_X,cX),axis=0)
    #     w_Y = np.concatenate((fine_Y,cY),axis=0)
    #     w_C = np.concatenate((fine_chem,cChem),axis=0)
    #     ax3.clear()
    #     ax3.contourf(w_Y,w_X,w_C, **kw3)
    #     ax3.set_title('Chemokine, Z=16 um, time = ' +str(saving_times[i])+ ' minutes') 
    #     ax3.invert_xaxis()
    #     ax3.axis('scaled')
        
    
    
    # ani3 = matplotlib.animation.FuncAnimation(fig3,animate3,blit=False, frames=number_of_frames,repeat=False)
    
    # plt.show()
    
    # ani3.save('./chemokine.gif', writer='imagemagick', fps=4)
    



#%%



if Total_Amount_Analysis == 'Y':
    o2_uptake_rate_per_cell = 0.005
    glu_uptake_rate_per_cell = 0.01
    chem_secretion_rate_per_cell_per_min = 0.01
    number_of_cells = 170278
    
    
    total_O2 = []
    total_glu = []
    total_chem = []
    
    initial_O2= 0;
    
    previous_data = np.array([0,0,0])
    previous_time = 0;
    
    for i in range(number_of_frames):
        time_p = time_point + '%02d'%(i)
        if path.exists(time_p + "_microenvironment0.mat"):
            fine_data = sio.loadmat(time_p + "_microenvironment0.mat")['multiscale_microenvironment']  
            micEnv_O2 = sum(fine_data[4,:])
            micEnv_glu = sum(fine_data[5,:])
            micEnv_chem = sum(fine_data[6,:])

            coarse_data = sio.loadmat(time_p + "_microenvironment1.mat")['multiscale_microenvironment']
            coarse_oxy = round(sum(coarse_data[4,:]),2)
            coarse_glu = sum(coarse_data[5,:])
            coarse_chem = sum(coarse_data[6,:])
            
            if i == 0:
                initial_O2 = micEnv_O2
                initial_glu = micEnv_glu
                initial_chem = micEnv_chem
    
    
            total_O2.append(micEnv_O2*2880*512*2880 + coarse_oxy*2880*2880*4880)
            total_glu.append(micEnv_glu*2880*512*2880 + coarse_glu*2880*2880*4880)
            total_chem.append(micEnv_chem)
    
            
    total_O2_c = [x / (2880*2880*5392) for x in total_O2]
    total_glu_c = [x / (2880*2880*5392) for x in total_glu]
    
    
    plt.figure()
    plt.plot(saving_times, total_O2_c)
    plt.title('Oxygen')
    plt.xlabel('time(min)')
    plt.ylabel('Concentration(mM)')
    plt.figure()
    plt.plot(saving_times, total_glu_c)
    plt.title('Glucose')
    plt.xlabel('time(min)')
    plt.ylabel('Concentration(mM)')
    
    # plt.figure()
    # plt.plot(saving_times, total_chem)
    # plt.title('Chemokine')
    # plt.xlabel('time(min)')
    # plt.ylabel('Concentration(mM)')
    
    