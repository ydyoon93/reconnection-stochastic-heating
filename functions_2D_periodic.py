# -*- coding: utf-8 -*-
"""
Created on Fri Mar 13 02:32:04 2018

@author: ydyoo

Define all 2D functions
"""
import numpy as np

#define cross product a x b; returns 3 dimensional tuple; cx,cy,cz=cross()
def cross(ax,ay,az,bx,by,bz):
    cx = np.multiply(ay,bz)-np.multiply(az,by);
    cy = np.multiply(az,bx)-np.multiply(ax,bz);
    cz = np.multiply(ax,by)-np.multiply(ay,bx);
    return cx, cy, cz

#define curl for 2D (d/dz = 0) 
#with neumann boundary conditions (derivatives at the boundaries are constant)
def curl_2D_periodic(ax,ay,az,dx,dy):
    xmax=ax.shape[0];
    ymax=ax.shape[1];
    cx = np.zeros((xmax,ymax),np.float64);
    cy = np.zeros((xmax,ymax),np.float64);
    cz = np.zeros((xmax,ymax),np.float64);
    x = slice(1,xmax-1);
    x1 = slice(2,xmax);
    x_1=slice(0,xmax-2);
    y = slice(1,ymax-1);
    y1 = slice(2,ymax);
    y_1 = slice(0,ymax-2);
    cx[x,y] = (az[x,y1]-az[x,y_1])/2/dy;
    cy[x,y] = -(az[x1,y]-az[x_1,y])/2/dx;
    cz[x,y] = (ay[x1,y]-ay[x_1,y])/2/dx \
            -(ax[x,y1]-ax[x,y_1])/2/dy;
    periodic_2D(cx);
    periodic_2D(cy);
    periodic_2D(cz);            
    return cx, cy, cz

#takes in a 2D vector field and applies periodic boundary condition
#symmetric in x, antisymmetric in y
def periodic_2D(ax):
    xmax=ax.shape[0];
    ymax=ax.shape[1];
    ax[0,:]=ax[2,:];
    ax[xmax-1,:]=ax[-3,:];
    ax[:,0]=ax[:,-2];
    ax[:,ymax-1]=ax[:,2];
    
#perform relaxation integration on laplacian(B)-B=Q
#where we know Q and are trying to find B
def helmholtz_int_2D(Qx,Qy,Qz,dx,dy,error):
    xmax=int(Qx.shape[0]);
    ymax=int(Qx.shape[1]);
    x = slice(1,xmax-1);
    x1 = slice(2,xmax);
    x_1=slice(0,xmax-2);
    y = slice(1,ymax-1);
    y1 = slice(2,ymax);
    y_1 = slice(0,ymax-2);
    halfx = slice(0,int(xmax/2));
    halfy = slice(0,int(ymax/2));
    Bx = np.zeros((xmax,ymax),np.float64);
    By = np.zeros((xmax,ymax),np.float64);
    Bz = np.zeros((xmax,ymax),np.float64);
    measured_error = [1];
    while measured_error[-1] > error:
        old_magnitude=np.sqrt(sum(sum(np.multiply(Bx[halfx,halfy],Bx[halfx,halfy])\
                              +np.multiply(By[halfx,halfy],By[halfx,halfy])\
                              +np.multiply(Bz[halfx,halfy],Bz[halfx,halfy]))));
        Bx[x,y]=(Qx[x,y]-(Bx[x1,y]+Bx[x_1,y])/dx**2-(Bx[x,y1]+Bx[x,y_1])/dy**2)\
                        /(-1-2/dx**2-2/dy**2);
        By[x,y]=(Qy[x,y]-(By[x1,y]+By[x_1,y])/dx**2-(By[x,y1]+By[x,y_1])/dy**2)\
                        /(-1-2/dx**2-2/dy**2);
        Bz[x,y]=(Qz[x,y]-(Bz[x1,y]+Bz[x_1,y])/dx**2-(Bz[x,y1]+Bz[x,y_1])/dy**2)\
                        /(-1-2/dx**2-2/dy**2);
        periodic_2D(Bx);
        periodic_2D(By);
        periodic_2D(Bz);
        new_magnitude=np.sqrt(sum(sum(np.multiply(Bx[halfx,halfy],Bx[halfx,halfy])\
                              +np.multiply(By[halfx,halfy],By[halfx,halfy])\
                              +np.multiply(Bz[halfx,halfy],Bz[halfx,halfy]))));
        measured_error.append(np.abs((new_magnitude-old_magnitude)/new_magnitude));
        if measured_error[-1] > measured_error[-2]:
            print('Oscillating Error')
            break;
    return Bx,By,Bz

#perform relaxation integration on laplacian(B)-(1+1/mass_ratio)*B=Q
#where we know Q and are trying to find B
def helmholtz_int_w_mass_ratio_2D(Qx,Qy,Qz,dx,dy,mass_ratio,error):
    xmax=int(Qx.shape[0]);
    ymax=int(Qx.shape[1]);
    x = slice(1,xmax-1);
    x1 = slice(2,xmax);
    x_1=slice(0,xmax-2);
    y = slice(1,ymax-1);
    y1 = slice(2,ymax);
    y_1 = slice(0,ymax-2);
    halfx = slice(0,int(xmax/2));
    halfy = slice(0,int(ymax/2));
    Bx = np.zeros((xmax,ymax),np.float64);
    By = np.zeros((xmax,ymax),np.float64);
    Bz = np.zeros((xmax,ymax),np.float64);
    measured_error = [1];
    while measured_error[-1] > error:
        old_magnitude=np.sqrt(sum(sum(np.multiply(Bx[halfx,halfy],Bx[halfx,halfy])\
                              +np.multiply(By[halfx,halfy],By[halfx,halfy])\
                              +np.multiply(Bz[halfx,halfy],Bz[halfx,halfy]))));
        Bx[x,y]=(Qx[x,y]-(Bx[x1,y]+Bx[x_1,y])/dx**2-(Bx[x,y1]+Bx[x,y_1])/dy**2)\
                        /(-(1+1/mass_ratio)-2/dx**2-2/dy**2);
        By[x,y]=(Qy[x,y]-(By[x1,y]+By[x_1,y])/dx**2-(By[x,y1]+By[x,y_1])/dy**2)\
                        /(-(1+1/mass_ratio)-2/dx**2-2/dy**2);
        Bz[x,y]=(Qz[x,y]-(Bz[x1,y]+Bz[x_1,y])/dx**2-(Bz[x,y1]+Bz[x,y_1])/dy**2)\
                        /(-(1+1/mass_ratio)-2/dx**2-2/dy**2);
        periodic_2D(Bx);
        periodic_2D(By);
        periodic_2D(Bz);
        new_magnitude=np.sqrt(sum(sum(np.multiply(Bx[halfx,halfy],Bx[halfx,halfy])\
                              +np.multiply(By[halfx,halfy],By[halfx,halfy])\
                              +np.multiply(Bz[halfx,halfy],Bz[halfx,halfy]))));
        measured_error.append(np.abs((new_magnitude-old_magnitude)/new_magnitude));
        if measured_error[-1] > measured_error[-2]:
            print('Oscillating Error')
            break;
    return Bx,By,Bz
        
#perform relaxation integration on laplacian(B)=Q
#where we know Q and are trying to find B
def laplacian_int_2D(Qx,Qy,Qz,dx,dy,error):
    xmax=int(Qx.shape[0]);
    ymax=int(Qx.shape[1]);
    x = slice(1,xmax-1);
    x1 = slice(2,xmax);
    x_1=slice(0,xmax-2);
    y = slice(1,ymax-1);
    y1 = slice(2,ymax);
    y_1 = slice(0,ymax-2);
    halfx = slice(0,int(xmax/2));
    halfy = slice(0,int(ymax/2));
    Bx = np.zeros((xmax,ymax),np.float64);
    By = np.zeros((xmax,ymax),np.float64);
    Bz = np.zeros((xmax,ymax),np.float64);
    measured_error = [1];
    while measured_error[-1] > error:
        old_magnitude=np.sqrt(sum(sum(np.multiply(Bx[halfx,halfy],Bx[halfx,halfy])\
                              +np.multiply(By[halfx,halfy],By[halfx,halfy])\
                              +np.multiply(Bz[halfx,halfy],Bz[halfx,halfy]))));
        Bx[x,y]=(Qx[x,y]-(Bx[x1,y]+Bx[x_1,y])/dx**2-(Bx[x,y1]+Bx[x,y_1])/dy**2)\
                        /(-2/dx**2-2/dy**2);
        By[x,y]=(Qy[x,y]-(By[x1,y]+By[x_1,y])/dx**2-(By[x,y1]+By[x,y_1])/dy**2)\
                        /(-2/dx**2-2/dy**2);
        Bz[x,y]=(Qz[x,y]-(Bz[x1,y]+Bz[x_1,y])/dx**2-(Bz[x,y1]+Bz[x,y_1])/dy**2)\
                        /(-2/dx**2-2/dy**2);
        periodic_2D(Bx);
        periodic_2D(By);
        periodic_2D(Bz);
        new_magnitude=np.sqrt(sum(sum(np.multiply(Bx[halfx,halfy],Bx[halfx,halfy])\
                              +np.multiply(By[halfx,halfy],By[halfx,halfy])\
                              +np.multiply(Bz[halfx,halfy],Bz[halfx,halfy]))));
        measured_error.append(np.abs((new_magnitude-old_magnitude)/new_magnitude));
        if measured_error[-1] > measured_error[-2]:
            print('Oscillating Error')
            break;
    return Bx,By,Bz

#perform relaxation integration on laplacian(B)=Q
#where we know Q and are trying to find B
#where B and Q are scalar quantities
def laplacian_scalar_int_2D(Q,dx,dy,error):
    xmax=int(Q.shape[0]);
    ymax=int(Q.shape[1]);
    x = slice(1,xmax-1);
    x1 = slice(2,xmax);
    x_1=slice(0,xmax-2);
    y = slice(1,ymax-1);
    y1 = slice(2,ymax);
    y_1 = slice(0,ymax-2);
    halfx = slice(0,int(xmax/2));
    halfy = slice(0,int(ymax/2));
    B = np.zeros((xmax,ymax),np.float64);
    measured_error = [1];
    while measured_error[-1] > error:
        old_magnitude=np.sqrt(sum(sum(np.multiply(B[halfx,halfy],B[halfx,halfy]))));
        B[x,y]=(Q[x,y]-(B[x1,y]+B[x_1,y])/dx**2-(B[x,y1]+B[x,y_1])/dy**2)\
                        /(-2/dx**2-2/dy**2);
        periodic_2D(B);
        new_magnitude=np.sqrt(sum(sum(np.multiply(B[halfx,halfy],B[halfx,halfy]))));
        measured_error.append(np.abs((new_magnitude-old_magnitude)/new_magnitude));
        if measured_error[-1] > measured_error[-2]:
            print('Oscillating Error')
            break;
    return B

#given a scalar field, compute the gradient vectors
def gradient_2D(field_scalar,dx,dy):
    xmax=int(field_scalar.shape[0]);
    ymax=int(field_scalar.shape[1]);
    x = slice(1,xmax-1);
    x1 = slice(2,xmax);
    x_1=slice(0,xmax-2);
    y = slice(1,ymax-1);
    y1 = slice(2,ymax);
    y_1 = slice(0,ymax-2);
    
    bx = np.zeros((xmax,ymax),np.float64);
    by = np.zeros((xmax,ymax),np.float64);
    bz = np.zeros((xmax,ymax),np.float64);
    
    bx[x,y] = (field_scalar[x1,y]-field_scalar[x_1,y])/2/dx;
    by[x,y] = (field_scalar[x,y1]-field_scalar[x,y_1])/2/dy;
    periodic_2D(bx);
    periodic_2D(by);
    periodic_2D(bz);
    
    return bx,by,bz

#given a vector field, calculate the divergence
def divergence_periodic_2D(ax,ay,az,dx,dy):
    xmax=int(ax.shape[0]);
    ymax=int(ax.shape[1]);
    x = slice(1,xmax-1);
    x1 = slice(2,xmax);
    x_1=slice(0,xmax-2);
    y = slice(1,ymax-1);
    y1 = slice(2,ymax);
    y_1 = slice(0,ymax-2);
    
    div_field = np.zeros((xmax,ymax),np.float64);
    
    div_field[x,y] = (ax[x1,y]-ax[x_1,y])/2/dx+(ay[x,y1]-ay[x,y_1])/2/dy;
    periodic_2D(div_field);
    
    return div_field


            