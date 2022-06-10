# -*- coding: utf-8 -*-
"""
Created on Fri Jun 10 10:24:32 2022

@author: Raphael Andrieux
"""

#pendule double:
%matplotlib qt

import numpy as np
import sympy as smp 
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from scipy.integrate import odeint
from scipy.integrate import solve_ivp
import scipy as  sp
from sympy import init_printing
from sympy.printing.mathml import print_mathml

plt.style.use("dark_background")
init_printing()
t,g=smp.symbols("t g")
m1,m2=smp.symbols("m1 m2")
L1,L2=smp.symbols("L1 L2")
the1,the2=smp.symbols(r'\theta_1,\theta_2',cls=smp.Function)
the1=the1(t)
the2=the2(t)
the1_d=smp.diff(the1,t)
the2_d=smp.diff(the2,t)
the1_dd=smp.diff(the1_d,t)
the2_dd=smp.diff(the2_d,t)

x1=L1*smp.sin(the1)
y1=-L1*smp.cos(the1)
x2=L1*smp.sin(the1)+L2*smp.sin(the2)
y2=-L1*smp.cos(the1)-L2*smp.cos(the2)
#definition des energie cinétique
C1=1/2*m1*(smp.diff(x1,t)**2+smp.diff(y1,t)**2)
C2=1/2*m1*(smp.diff(x2,t)**2+smp.diff(y2,t)**2)
C=C1+C2
#energie potentielle:
P1=m1*g*y1
P2=m2*g*y2
P=P1+P2
#lagrangienne
L=C-P

#les equations du mouvement de lagrange
LE1=smp.diff(L,the1)-smp.diff(smp.diff(L,the1_d),t).simplify()
LE2=smp.diff(L,the2)-smp.diff(smp.diff(L,the2_d),t).simplify()
#ici je veux solver pour la deuxieme derivative theta 2 en fonction du temps.
sols=smp.solve([LE1,LE2],(the1_dd,the2_dd),simplify=False, rational=False)

#ont a ainsi deux ODE du deuxieme ordre, maintenant il nous suffit de poser z=theta1', et u=theta2'
dzdt_f=smp.lambdify((t,g,m1,m2,L1,L2,the1,the2,the1_d,the2_d),sols[the1_dd])
dudt_f=smp.lambdify((t,g,m1,m2,L1,L2,the1,the2,the1_d,the2_d),sols[the2_dd])
dthe1dt_f=smp.lambdify(the1_d,the1_d)
dthe2dt_f=smp.lambdify(the2_d,the2_d)
def dSdt(S,t,g,m1,m2,L1,L2):
    the1,z,the2,u=S
    return [dthe1dt_f(z),dzdt_f(t,g,m1,m2,L1,L2,the1,the2,z,u),dthe2dt_f(u),dudt_f(t,g,m1,m2,L1,L2,the1,the2,z,u)]
t=np.linspace(0,40,1000)
g=9.81
m1=1
m2=1
L1=2
L2=1
bob_radius = 0.08
theta0=1
theta1ini=-3
sol=odeint(dSdt,y0=[theta0,-3,theta1ini,5],t=t,args=(g,m1,m2,L1,L2))
fpersec=len(t[t<1])#juste utile de connaitre les frames
theta1=sol.T[0]
theta2=sol.T[2]
def get_cords(t,the1,the2,L1,L2):
    return (L1*np.sin(the1),-L1*np.cos(the1), L1*np.sin(the1)+L2*np.sin(the2),-L1*np.cos(the1)-L2*np.cos(the2))
x1,y1,x2,y2=get_cords(t,sol.T[0],sol.T[2],L1,L2)
def animate(i):
    ligne1.set_data([0, x1[i], x2[i]], [0, y1[i], y2[i]])
    circle.set_center([x2[i],y2[i]])
    circle1.set_center([x1[i],y1[i]])
    
fig, ax = plt.subplots(1,1, figsize=(8,8))
circle = ax.add_patch(plt.Circle([np.cos(theta0),np.sin(theta0)], bob_radius,
                      fc='r', zorder=3))
circle1 = ax.add_patch(plt.Circle([np.cos(theta1ini),np.sin(theta1ini)], bob_radius,
                      fc='r', zorder=3))

ax.get_xaxis().set_ticks([])
ax.get_yaxis().set_ticks([])#enleve les ticks
ligne1,=plt.plot([],[],"-",lw=3,markersize=8)
ax.set_ylim(-4,4)
ax.set_xlim(-4,4)
ani=animation.FuncAnimation(fig,animate,frames=1000,interval=50)
ani.save("pen.gif",writer="pillow",fps=25)