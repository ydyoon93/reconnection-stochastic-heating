# -*- coding: utf-8 -*-
"""
Created on Wed Mar 14 00:59:24 2018

@author: ydyoo

Simulate Reconnection with periodic boundary conditions
"""

import numpy as np
import matplotlib.pyplot as plt
import time
import functions_2D_periodic as f2d


exec(open("Parameters.txt").read())

#Initialize variables
Bx=np.zeros((Xmax,Ymax),np.float64);
By=np.zeros((Xmax,Ymax),np.float64);
Bz=np.zeros((Xmax,Ymax),np.float64);
Qxi=np.zeros((Xmax,Ymax),np.float64);
Qyi=np.zeros((Xmax,Ymax),np.float64);
Qzi=np.zeros((Xmax,Ymax),np.float64);
Qxe=np.zeros((Xmax,Ymax),np.float64);
Qye=np.zeros((Xmax,Ymax),np.float64);
Qze=np.zeros((Xmax,Ymax),np.float64);
uxi=np.zeros((Xmax,Ymax),np.float64);
uyi=np.zeros((Xmax,Ymax),np.float64);
uzi=np.zeros((Xmax,Ymax),np.float64);
uxe=np.zeros((Xmax,Ymax),np.float64);
uye=np.zeros((Xmax,Ymax),np.float64);
uze=np.zeros((Xmax,Ymax),np.float64);


#Initial Background field
for x in range(Xmax) :
    for y in range(Ymax) :
        By[x,y]=np.tanh((x-Xmax/2)/(Lx/dx));
        Bz[x,y]=alpha*1;

#perturbation
Ax=np.zeros((Xmax,Ymax),np.float64);
Ay=np.zeros((Xmax,Ymax),np.float64);
Az=np.zeros((Xmax,Ymax),np.float64);
for x in range(Xmax) :
    for y in range(Ymax) :
        Az[x,y]=-epsilon*np.exp(-(x-Xmax/2)**2/2/(Lx/dx)**2\
                                -(y-Ymax/2)**2/2/(Ly/dy)**2);
pertx,perty,pertz=f2d.curl_2D_periodic(Ax,Ay,Az,dx,dy);
Bx += pertx;
By += perty;
Bz += pertz;

#Initial Particle Flows + perturbation
curl_B_x, curl_B_y, curl_B_z = f2d.curl_2D_periodic(Bx,By,Bz,dx,dy);
curl_B = (curl_B_x, curl_B_y, curl_B_z);

(uxi, uyi, uzi)=tuple([1/(1+mass_ratio)*x for x in curl_B]);
(uxe, uye, uze)=tuple([-mass_ratio/(1+mass_ratio)*x for x in curl_B]);

del curl_B, curl_B_x, curl_B_y, curl_B_z;

#Initial Q + Pertubation
curl_ui_x, curl_ui_y, curl_ui_z = f2d.curl_2D_periodic(uxi,uyi,uzi,dx,dy);
(Qxi,Qyi,Qzi)=(mass_ratio*curl_ui_x+Bx,mass_ratio*curl_ui_y+By,\
                mass_ratio*curl_ui_z+Bz);

curl_ue_x, curl_ue_y, curl_ue_z = f2d.curl_2D_periodic(uxe,uye,uze,dx,dy);
(Qxe,Qye,Qze)=(curl_ue_x-Bx,curl_ue_y-By,curl_ue_z-Bz);

del curl_ui_x, curl_ui_y, curl_ui_z, curl_ue_x, curl_ue_y, curl_ue_z;


#test variables
max=[];

#solve main equations
t0 = time.time();
for t in range(Tmax):
    #store variables first
    if t % int(1/dt) == 0:
        with open('Qxi.txt','ab') as qxi_file:
            np.savetxt(qxi_file,Qxi);
        with open('Qyi.txt','ab') as qyi_file:
            np.savetxt(qyi_file,Qyi);
        with open('Qzi.txt','ab') as qzi_file:
            np.savetxt(qzi_file,Qzi);
        with open('Qxe.txt','ab') as qxe_file:
            np.savetxt(qxe_file,Qxe);
        with open('Qye.txt','ab') as qye_file:
            np.savetxt(qye_file,Qye);
        with open('Qze.txt','ab') as qze_file:
            np.savetxt(qze_file,Qze);
        with open('uxi.txt','ab') as uxi_file:
            np.savetxt(uxi_file,uxi);
        with open('uyi.txt','ab') as uyi_file:
            np.savetxt(uyi_file,uyi);
        with open('uzi.txt','ab') as uzi_file:
            np.savetxt(uzi_file,uzi);
        with open('uxe.txt','ab') as uxe_file:
            np.savetxt(uxe_file,uxe);
        with open('uye.txt','ab') as uye_file:
            np.savetxt(uye_file,uye);
        with open('uze.txt','ab') as uze_file:
            np.savetxt(uze_file,uze);
        with open('Bx.txt','ab') as Bx_file:
            np.savetxt(Bx_file,Bx);
        with open('By.txt','ab') as By_file:
            np.savetxt(By_file,By);
        with open('Bz.txt','ab') as Bz_file:
            np.savetxt(Bz_file,Bz);
        #measure time
        t1=time.time();
        time_elapsed=t1-t0;
        print(t*dt)
        print('Estimated Time Left: ' + str(time_elapsed*(float(Tmax-t))*dt))
        t0=time.time();
    #temp are temporary variables
    #advance Q_i
    tempx,tempy,tempz = f2d.cross(uxi,uyi,uzi,Qxi,Qyi,Qzi);
    tempx,tempy,tempz = f2d.curl_2D_periodic(tempx,tempy,tempz,dx,dy);
    (Qxi,Qyi,Qzi) = (Qxi+dt*tempx,Qyi+dt*tempy,Qzi+dt*tempz);
    
    #advance Q_e
    tempx,tempy,tempz = f2d.cross(uxe,uye,uze,Qxe,Qye,Qze);
    tempx,tempy,tempz = f2d.curl_2D_periodic(tempx,tempy,tempz,dx,dy);
    (Qxe,Qye,Qze) = (Qxe+dt*tempx,Qye+dt*tempy,Qze+dt*tempz);
    
    #find B
    tempx,tempy,tempz = Qxe-Qxi/mass_ratio,Qye-Qyi/mass_ratio,Qze-Qzi/mass_ratio;
    Bx,By,Bz = f2d.helmholtz_int_w_mass_ratio_2D(tempx,tempy,tempz,dx,dy,mass_ratio,int_error);
    
    #find u_i
    tempx,tempy,tempz = Bx-Qxi,By-Qyi,Bz-Qzi;
    tempx,tempy,tempz = f2d.curl_2D_periodic(tempx,tempy,tempz,dx,dy);
    tempx,tempy,tempz = tempx/mass_ratio,tempy/mass_ratio,tempz/mass_ratio;
    uxi,uyi,uzi = f2d.laplacian_int_2D(tempx,tempy,tempz,dx,dy,int_error);
    
    #find u_e
    tempx,tempy,tempz = -Bx-Qxe,-By-Qye,-Bz-Qze;
    tempx,tempy,tempz = f2d.curl_2D_periodic(tempx,tempy,tempz,dx,dy);
    uxe,uye,uze = f2d.laplacian_int_2D(tempx,tempy,tempz,dx,dy,int_error);