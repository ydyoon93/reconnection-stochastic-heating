# -*- coding: utf-8 -*-
"""
Created on Wed Apr  4 11:44:37 2018

@author: ydyoo

Calculate ion trajectories and compare with ExB and polarization drift
to see if there is any stochastic heating
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
#with open('phi.txt') as phi_data:
#    phi_all = np.loadtxt(phi_data);
    
mass_ratio = 10.;
for option in range(2):
    t=200.;
    Ex = Ex_all[int(t*Xmax):int((t+1)*Xmax),:];
    Ey = Ey_all[int(t*Xmax):int((t+1)*Xmax),:];
    Ez = Ez_all[int(t*Xmax):int((t+1)*Xmax),:];
    Bx = Bx_all[int(t*Xmax):int((t+1)*Xmax),:];
    By = By_all[int(t*Xmax):int((t+1)*Xmax),:];
    Bz = Bz_all[int(t*Xmax):int((t+1)*Xmax),:];
    x_start = np.arange(100.,130.,0.1);
    y_start = np.arange(300.,302.,0.1);
    z_start = np.zeros(np.size(x_start)*np.size(y_start));
    tempx,tempy = np.meshgrid(x_start,y_start);
    x = np.array([np.float64(np.reshape(tempx,(np.size(x_start)*np.size(y_start))))]);
    y = np.array([np.float64(np.reshape(tempy,(np.size(x_start)*np.size(y_start))))]);
    z = np.array([np.float64(z_start)]);
    ux = np.array([np.zeros(np.size(x_start)*np.size(y_start))]);
    uy = np.array([np.zeros(np.size(x_start)*np.size(y_start))]);
    uz = np.array([np.zeros(np.size(x_start)*np.size(y_start))]);
    
    for t in range(580):
    #    calculate actual trajectory
        if option == 1:
            Ex = Ex_all[int(t*Xmax):int((t+1)*Xmax),:];
            Ey = Ey_all[int(t*Xmax):int((t+1)*Xmax),:];
            Ez = Ez_all[int(t*Xmax):int((t+1)*Xmax),:];
            Bx = Bx_all[int(t*Xmax):int((t+1)*Xmax),:];
            By = By_all[int(t*Xmax):int((t+1)*Xmax),:];
            Bz = Bz_all[int(t*Xmax):int((t+1)*Xmax),:];
        
        Ex_interp = sp.interpolate.RectBivariateSpline(np.arange(Xmax),np.arange(Ymax),Ex).ev(x[-1,:],y[-1,:]);
        Ey_interp = sp.interpolate.RectBivariateSpline(np.arange(Xmax),np.arange(Ymax),Ey).ev(x[-1,:],y[-1,:]);
        Ez_interp = sp.interpolate.RectBivariateSpline(np.arange(Xmax),np.arange(Ymax),Ez).ev(x[-1,:],y[-1,:]);
        Bx_interp = sp.interpolate.RectBivariateSpline(np.arange(Xmax),np.arange(Ymax),Bx).ev(x[-1,:],y[-1,:]);
        By_interp = sp.interpolate.RectBivariateSpline(np.arange(Xmax),np.arange(Ymax),By).ev(x[-1,:],y[-1,:]);
        Bz_interp = sp.interpolate.RectBivariateSpline(np.arange(Xmax),np.arange(Ymax),Bz).ev(x[-1,:],y[-1,:]);
        
        x = np.vstack((x,x[-1,:]+ux[-1,:]/dx));
        y = np.vstack((y,y[-1,:]+uy[-1,:]/dy));
        z = np.vstack((z,z[-1,:]+uz[-1,:]));
        
        ux_old = ux[-1,:];
        uy_old = uy[-1,:];
        uz_old = uz[-1,:];
        Ax,Ay,Az = Bx_interp/2/mass_ratio,By_interp/2/mass_ratio,Bz_interp/2/mass_ratio;
        Cx = (ux_old + (Ex_interp+(uy_old*Bz_interp-uz_old*By_interp)/2)/mass_ratio);
        Cy = (uy_old + (Ey_interp+(uz_old*Bx_interp-ux_old*Bz_interp)/2)/mass_ratio);
        Cz = (uz_old + (Ez_interp+(ux_old*By_interp-uy_old*Bx_interp)/2)/mass_ratio);
        AdotC = Ax*Cx+Ay*Cy+Az*Cz;
        Asquared = Ax**2+Ay**2+Az**2;
        
        ux = np.vstack((ux,(Cx + AdotC*Ax - Ay*Cz + Az*Cy)/(1+Asquared)));
        uy = np.vstack((uy,(Cy + AdotC*Ay - Az*Cx + Ax*Cz)/(1+Asquared)));
        uz = np.vstack((uz,(Cz + AdotC*Az - Ax*Cy + Ay*Cx)/(1+Asquared)));
    
    #    #calculate using ui information
    #    uxi_interp = sp.interpolate.RectBivariateSpline(np.arange(Xmax),np.arange(Ymax),uxi).ev(x_ui[-1],y_ui[-1]);
    #    uyi_interp = sp.interpolate.RectBivariateSpline(np.arange(Xmax),np.arange(Ymax),uyi).ev(x_ui[-1],y_ui[-1]);
    #    uzi_interp = sp.interpolate.RectBivariateSpline(np.arange(Xmax),np.arange(Ymax),uzi).ev(x_ui[-1],y_ui[-1]);
    #    
    #    x_ui = np.vstack((x_ui,x_ui[-1]+uxi_interp/dx));
    #    y_ui = np.vstack((y_ui,y_ui[-1]+uyi_interp/dy));
    #    z_ui = np.vstack((z_ui,z_ui[-1]+uzi_interp));
        if t % 100 == 0:
            print(t)
    if option == 0:
        x_stable = x;
        y_stable = y;
        z_stable = z;
        ux_stable = ux;
        uy_stable = uy;
        uz_stable = uz;
        
    x_real = (x-Xmax/2)*dx;
    y_real = (y-Ymax/2)*dy;
    z_real = z;

    x_real_stable = (x_stable-Xmax/2)*dx;
    y_real_stable = (y_stable-Ymax/2)*dy;
    z_real_stable = z_stable;

        
plt.figure(1);
plot_trajectory(x_real_stable,y_real_stable,z_real_stable,ux_stable,uy_stable,uz_stable);
plt.figure(2);
plot_trajectory(x_real,y_real,z_real,ux,uy,uz);

#sp.io.savemat('trajectory_2.mat',mdict = {'x_stable':x_real_stable,'y_stable':y_real_stable,'z_stable':z_real_stable,'x_stochastic':x_real,'y_stochastic':y_real,'z_stochastic':z_real,'ux_stable':ux_stable,'uy_stable':uy_stable,'uz_stable':uz_stable,'ux_stochastic':ux,'uy_stochastic':uy,'uz_stochastic':uz})

def velocity_sign(array):
    sign_array = array/np.abs(array);
    return sign_array

def plot_trajectory(x,y,z,ux,uy,uz):
    plt.subplot(2,2,1);
    plt.plot(x[:,5700:5710],y[:,5700:5710],linewidth = 0.2,color='r');
    plt.subplot(2,2,2);
    plt.scatter(uy[-1,:],y[-1,:],s=0.1);
    plt.subplot(2,2,3);
    plt.hist(np.abs(uy[-1,:]-np.mean(uy[-1,:]))**2,bins=20,range=(0,0.0001));
#    plt.scatter(ux[-1,:],x[-1,:],s=1);
    plt.subplot(2,2,4)
    plt.scatter(ux[-1,:],x[-1,:],s=0.1);
    

