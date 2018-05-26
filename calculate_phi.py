# -*- coding: utf-8 -*-
"""
Created on Wed Apr  4 11:44:37 2018

@author: ydyoo

Calculate phi from electric field
"""

import functions_2D_periodic as f2d
import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
from mpl_toolkits.mplot3d import Axes3D

exec(open("Parameters.txt").read())

#with open('Ex.txt') as Ex_data:
#    Ex_all = np.loadtxt(Ex_data);
#with open('Ey.txt') as Ey_data:
#    Ey_all = np.loadtxt(Ey_data);
#with open('Ez.txt') as Ez_data:
#    Ez_all = np.loadtxt(Ez_data);
#with open('Bx.txt') as Bx_data:
#    Bx_all = np.loadtxt(Bx_data);
#with open('By.txt') as By_data:
#    By_all = np.loadtxt(By_data);
#with open('Bz.txt') as Bz_data:
#    Bz_all = np.loadtxt(Bz_data);
#with open('uxi.txt') as uxi_data:
#    uxi_all = np.loadtxt(uxi_data);
#with open('uyi.txt') as uyi_data:
#    uyi_all = np.loadtxt(uyi_data);
#with open('uzi.txt') as uzi_data:
#    uzi_all = np.loadtxt(uzi_data);

for t in range(101,102):
    #calculate actual trajectory
    Ex = Ex_all[int(t*Xmax):int((t+1)*Xmax),:];
    Ey = Ey_all[int(t*Xmax):int((t+1)*Xmax),:];
    Ez = Ez_all[int(t*Xmax):int((t+1)*Xmax),:];

    Grad_E = f2d.divergence_periodic_2D(Ex,Ey,Ez,dx,dy);
    
    plt.contourf(Grad_E)
    plt.colorbar
#    phi = f2d.laplacian_scalar_int_2D(-Grad_E,dx,dy,int_error);
    
#    with open('phi.txt','ab') as phi_file:
#        np.savetxt(phi_file,phi);
#    

