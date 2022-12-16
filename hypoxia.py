#!/usr/bin/env python
# coding: utf-8


import plotly.graph_objects as go
import matplotlib.pyplot as plt
import numpy as np


lamb=3e4
mu=3e4
nu=3e4
rho1=1/24;
rho2=1/48; 
d1=1*10**(-12)*3600;
d2=1*10**(-11)*3600; 
v0=-3/4;
k=np.sqrt(lamb**2+mu**2+nu**2);
v1=-1/(d2*k**2);
C=(d1*k**2*v1**2+rho1-rho2+v1)*rho1/(rho1**2-2*rho2*rho1+rho2**2);

g0=-rho2/rho1*(C*rho1-C*rho2-rho1)

print("lamb= ",str(lamb))
print("mu= ",str(mu))
print("nu= ",str(nu))
print("rho1= ",str(rho1))
print("rho2= ",str(rho2))
print("d1= ",str(d1))
print("d2= ",str(d2))
print("v0= ",str(v0))
print("g0= ",str(g0))
print("k= ",str(k))
print("v1= ",str(v1))
print("C=",str(C))
def xi(x,y,z,t):
    return lamb*x+mu*y+nu*z+t;

def v(x,y,z,t):
    return C+v0*np.exp(-v1*xi(x,y,z,t));

def u(x,y,z,t):
    return -C*rho2/rho1-v0*np.exp(-v1*xi(x,y,z,t))
r=0.0001;
a, b = np.mgrid[0:2*np.pi:20j, 0:np.pi:20j] #np.pi/2
x = r*np.cos(a)*np.sin(b)
y = r*np.sin(a)*np.sin(b)
z = r*np.cos(b)
time=0;

z_level=0.0
znew=z_level*np.ones(x.shape)  
fignew = go.Figure(data=[go.Surface(x=x, y=y, z=z, surfacecolor=u(x,y,z,time),colorscale='RdYlBu', 
                                    reversescale=True, opacity=0.5,colorbar=dict(len=0.5,title='u[cell/cm^3/h]')),
            
                        ])


fignew.update_layout(margin=dict(l=20, r=20, t=20, b=20),
    height = 800, width =800,yaxis_range=[-0.5,0.5],  scene={"aspectratio": {"x": 1, "y": 1, "z": 1},                                                        
                    'xaxis_title':'x [cm]',
                    'yaxis_title':'y [cm]',
                    'zaxis_title':'z [cm]'},
                     scene_camera = dict(    eye=dict(x=1.8, y=-1.8, z=0.5))
          )

          
fignew.show()


figneww = go.Figure(data=[go.Surface(x=x, y=y, z=z, surfacecolor=v(x,y,z,time),
                                     colorscale='RdYlBu', reversescale=True, 
                                     opacity=0.5,colorbar=dict(len=0.5,title='v[cell/cm^3/h]'))
                                              ])


figneww.update_layout(margin=dict(l=20, r=20, t=20, b=20),
    height = 800, width =800,yaxis_range=[-0.5,0.5],  scene={"aspectratio": {"x": 1, "y": 1, "z": 1},                                                        
                    'xaxis_title':'x [cm]',
                    'yaxis_title':'y [cm]',
                    'zaxis_title':'z [cm]'},
                      scene_camera = dict(eye=dict(x=1.8, y=-1.8, z=0.5))
                    
)

          
figneww.show()

####################


n=2;
fig, axs = plt.subplots(4,2,dpi=300,figsize=(5*n,5*n))
radius=0.0001
var=np.linspace(0,1)

axs[0][0].plot(var,u(var*radius,0,0,0),color="C3")
axs[0][0].plot(var,u(2*var*radius,0,0,0),color="C3",linestyle='dashed')
axs[0][0].plot(var,u(1/2*var*radius,0,0,0),color="C3",linestyle='dotted')
axs[0][0].set_ylabel('$u\quad[cell/\mu m^3/h]$')
axs[0][0].set_xlabel('$x\quad[\mu m]$')
axs[0][1].plot(var,v(var*radius,0,0,0))
axs[0][1].plot(var,v(2*var*radius,0,0,0),color="C0",linestyle='dashed')
axs[0][1].plot(var,v(1/2*var*radius,0,0,0),color="C0",linestyle='dotted')
axs[0][1].set_ylabel('$v\quad[cell/\mu m^3/h]$')
axs[0][1].set_xlabel('$x\quad[\mu cm]$')
axs[1][0].plot(var,u(0,var*radius,0,0),color="C3")
axs[1][0].plot(var,u(0,2*var*radius,0,0),color="C3",linestyle='dashed')
axs[1][0].plot(var,u(0,1/2*var*radius,0,0),color="C3",linestyle='dotted')
axs[1][0].set_ylabel('$u\quad[cell/\mu m^3/h]$')
axs[1][0].set_xlabel('$y\quad[\mu cm]$')
axs[1][1].plot(var,v(0,var*radius,0,0))
axs[1][1].plot(var,v(0,2*var*radius,0,0),color="C0",linestyle='dashed')
axs[1][1].plot(var,v(0,1/2*var*radius,0,0),color="C0",linestyle='dotted')
axs[1][1].set_ylabel('$v\quad[cell/\mu m^3/h]$')
axs[1][1].set_xlabel('$y\quad[\mu cm]$')
axs[2][0].plot(var,u(0,0,var*radius,0),color="C3")
axs[2][0].plot(var,u(0,0,2*var*radius,0),color="C3",linestyle='dashed')
axs[2][0].plot(var,u(0,0,1/2*var*radius,0),color="C3",linestyle='dotted')
axs[2][0].set_ylabel('$u\quad[cell/\mu m^3/h]$')
axs[2][0].set_xlabel('$z\quad[\mu cm]$')
axs[2][1].plot(var,v(0,0,var*radius,0))
axs[2][1].plot(var,v(0,0,2*var*radius,0),color="C0",linestyle='dashed')
axs[2][1].plot(var,v(0,0,1/2*var*radius,0),color="C0",linestyle='dotted')
axs[2][1].set_ylabel('$v\quad[cell/\mu m^3/h]$')
axs[2][1].set_xlabel('$z\quad[\mu cm]$')


var=np.linspace(0,0.75)
axs[3][0].plot(var,u(0,0,0,var*24),color="C3")
axs[3][0].plot(var,u(0,0,0,2*var*24),color="C3",linestyle='dashed')
axs[3][0].plot(var,u(0,0,0,1/2*var*24),color="C3",linestyle='dotted')
axs[3][0].set_ylabel('$u\quad[cell/cm^3/days]$')
axs[3][0].set_xlabel('$t\quad[days]$')
axs[3][1].plot(var,v(0,0,0,var*24))
axs[3][1].plot(var,v(0,0,0,2*var*24),color="C0",linestyle='dashed')
axs[3][1].plot(var,v(0,0,0,1/2*var*24),color="C0",linestyle='dotted')
axs[3][1].set_ylabel('$v\quad[cell/cm^3/days]$')
axs[3][1].set_xlabel('$t\quad[days]$')

plt.tight_layout( pad=1.08, h_pad=None, w_pad=None, rect=None)

###############




from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.path import Path
import matplotlib.patches as patches
import matplotlib as mpl

fig    = plt.figure(figsize=(10,12),dpi=300)
ax     = fig.gca(projection='3d')

r_ = np.linspace(0,0.0001,20)
theta_ = np.linspace(0,2*np.pi,20)

r, theta = np.meshgrid(r_, theta_)
x = r*np.sin(theta)
y = r*np.cos(theta)
radius=0.0001;
ratio=24;
X, Y   = np.meshgrid(x, y)
Z1 = u(x,y,0,0)
Z2 = u(x,y,0,0.5*ratio)
Z3 = u(x,y,0,1*ratio)
Z4 = u(x,y,0,1.5*ratio)
vmin=0.1
vmax=0.5

levels=np.linspace(Z1.min(), Z1.max(), 10)
p=ax.contourf(x,y,Z1, levels=levels, zdir='z', offset=0, cmap=plt.get_cmap('RdYlBu').reversed(),vmin = vmin, vmax = vmax)

levels=np.linspace(Z2.min(), Z2.max(), 10)
ax.contourf(x,y,Z2, levels=levels, zdir='z', offset=0.5, cmap=plt.get_cmap('RdYlBu').reversed(),vmin = vmin, vmax = vmax)

levels=np.linspace(Z3.min(), Z3.max(), 10)
ax.contourf(x,y,Z3, levels=levels, zdir='z', offset=1, cmap=plt.get_cmap('RdYlBu').reversed(),vmin = vmin, vmax = vmax)
levels=np.linspace(Z4.min(), Z4.max(), 10)
ax.contourf(x,y,Z4, levels=levels, zdir='z', offset=1.5, cmap=plt.get_cmap('RdYlBu').reversed(),vmin = vmin, vmax = vmax)
#
ax.set_xlim3d(-radius, radius)
ax.set_zlim3d(0,1.5)
ax.set_ylim3d(-radius, radius)
ax.view_init(-140, 60)
ax.set_xticklabels([-1,"","","",0,"","","",1])
ax.set_yticklabels([-1,"","","",0,"","","",1])

cmap = plt.get_cmap('RdYlBu').reversed()
norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
cbar=fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap),ax=ax)

ax.set_xlabel('$x\quad[\mu m]$',labelpad=10)
ax.set_ylabel('$y\quad[\mu m]$',labelpad=10)
ax.set_zlabel('$t\quad[days]$',labelpad=10)

cbar.set_label('$u(x,y,0,t)\quad[cells/\mu m^3/days]$')

#############



from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.path import Path
import matplotlib.patches as patches
import matplotlib as mpl

fig    = plt.figure(figsize=(10,12),dpi=300)
ax     = fig.gca(projection='3d')

r_ = np.linspace(0,0.0001,20)
theta_ = np.linspace(0,2*np.pi,20)

r, theta = np.meshgrid(r_, theta_)
x = r*np.sin(theta)
y = r*np.cos(theta)
radius=0.0001;
ratio=24;
X, Y   = np.meshgrid(x, y)
Z1 = v(x,y,0,0)
Z2 = v(x,y,0,0.5*ratio)
Z3 = v(x,y,0,1*ratio)
Z4 = v(x,y,0,1.5*ratio)
vmin=0
vmax=0.5

levels=np.linspace(Z1.min(), Z1.max(), 10)
p=ax.contourf(x,y,Z1, levels=levels, zdir='z', offset=0, cmap=plt.get_cmap('RdYlBu').reversed(),vmin = vmin, vmax = vmax)

levels=np.linspace(Z2.min(), Z2.max(), 10)
ax.contourf(x,y,Z2, levels=levels, zdir='z', offset=0.5, cmap=plt.get_cmap('RdYlBu').reversed(),vmin = vmin, vmax = vmax)

levels=np.linspace(Z3.min(), Z3.max(), 10)
ax.contourf(x,y,Z3, levels=levels, zdir='z', offset=1, cmap=plt.get_cmap('RdYlBu').reversed(),vmin = vmin, vmax = vmax)
levels=np.linspace(Z4.min(), Z4.max(), 10)
ax.contourf(x,y,Z4, levels=levels, zdir='z', offset=1.5, cmap=plt.get_cmap('RdYlBu').reversed(),vmin = vmin, vmax = vmax)
#
ax.set_xlim3d(-radius, radius)
ax.set_zlim3d(0,1.5)
ax.set_ylim3d(-radius, radius)
ax.view_init(-140, 60)
ax.set_xticklabels([-1,"","","",0,"","","",1])
ax.set_yticklabels([-1,"","","",0,"","","",1])

cmap = plt.get_cmap('RdYlBu').reversed()
norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
cbar=fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap),ax=ax)

ax.set_xlabel('$x\quad[\mu m]$',labelpad=10)
ax.set_ylabel('$y\quad[\mu m]$',labelpad=10)
ax.set_zlabel('$t\quad[days]$',labelpad=10)

cbar.set_label('$v(x,y,0,t)\quad[cells/\mu m^3/days]$')

