# -*- coding: utf-8 -*-
"""
Created on Tue Mar  6 16:56:18 2018

@author: ydyoo

Calculate E from u_i and Q_i by using Generalized Ohm's law
"""
import functions_2D_periodic as f2d
import numpy as np
import matplotlib.pyplot as plt

exec(open("Parameters.txt").read())

with open('Qxe.txt') as Qxe_data:
    Qxe_all = np.loadtxt(Qxe_data);
with open('Qye.txt') as Qye_data:
    Qye_all = np.loadtxt(Qye_data);
with open('Qze.txt') as Qze_data:
    Qze_all = np.loadtxt(Qze_data);
with open('uxe.txt') as uxe_data:
    uxe_all = np.loadtxt(uxe_data);
with open('uye.txt') as uye_data:
    uye_all = np.loadtxt(uye_data);
with open('uze.txt') as uze_data:
    uze_all = np.loadtxt(uze_data);

#using subscript e because I'm lazy and I don't want to change everything from
#electric field calculation from electron motion


#initialize u_i and u_e


for t in range(Tmax-1):
    
    Qxe = Qxe_all[int(t*Xmax):int((t+1)*Xmax),:];
    Qye = Qye_all[int(t*Xmax):int((t+1)*Xmax),:];
    Qze = Qze_all[int(t*Xmax):int((t+1)*Xmax),:];
    uxe = uxe_all[int(t*Xmax):int((t+1)*Xmax),:];
    uye = uye_all[int(t*Xmax):int((t+1)*Xmax),:];
    uze = uze_all[int(t*Xmax):int((t+1)*Xmax),:];
    
    electron_KE = (np.multiply(uxe,uxe)+np.multiply(uye,uye)+np.multiply(uze,uze))/2;
       
    ue_cross_Qe_x,ue_cross_Qe_y,ue_cross_Qe_z=f2d.cross(uxe,uye,uze,Qxe,Qye,Qze);    
    
    uxe_new = uxe_all[int((t+1)*Xmax):int((t+2)*Xmax),:];
    uye_new = uye_all[int((t+1)*Xmax):int((t+2)*Xmax),:];
    uze_new = uze_all[int((t+1)*Xmax):int((t+2)*Xmax),:];
    
    duxedt,duyedt,duzedt = uxe_new-uxe, uye_new-uye, uze_new-uze;
    
    tempx,tempy,tempz = f2d.gradient_2D(electron_KE,dx,dy);
    
    Ex = ue_cross_Qe_x-duxedt-tempx;
    Ey = ue_cross_Qe_y-duyedt-tempy;
    Ez = ue_cross_Qe_z-duzedt-tempz;
    
    print(t);
    
    with open('Ex.txt','ab') as Exi_file:
        np.savetxt(Exi_file,Ex);
    with open('Ey.txt','ab') as Eyi_file:
        np.savetxt(Eyi_file,Ey);
    with open('Ez.txt','ab') as Ezi_file:
        np.savetxt(Ezi_file,Ez);