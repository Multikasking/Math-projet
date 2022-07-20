%matplotlib qt
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from scipy.integrate import odeint
from scipy.integrate import solve_ivp
import scipy as  sp
plt.style.use("dark_background")
#equation differentiel du deuxieme ordre tq theta''=-g/l*sin(theta)
#psi=theta', on a alors psi'=-b/m*theta'-g/l*sin(theta)
def dthetadt(t,S):
    t,omega=S
    g,l=9.81,1
    b=1
    m=100
    return [omega,-b/m*omega -g/l*np.sin(t)]
L=1
g=9.81
t0=np.radians(90)
psi0=(2.000001*g/L)**(0.5)
S0=[t0,psi0]
t=np.linspace(0,20,1000)
sol=odeint(dthetadt,y0=S0,t=t,tfirst=True)
plt.plot(t,sol.T[0],label="w")

plt.plot(t,sol.T[1],label="acceleration angulaire")
leg=plt.legend(loc="upper right")
plt.show()
plt.plot(sol.T[0],sol.T[1])
plt.show()
fig=plt.figure()
ax=fig.add_subplot(aspect="equal")
ax.set_xlim(-L*1.2,L*1.2)
ax.set_ylim(-L*1.2, L*1.2)
x,y=[],[]
theta0=t0
x1,y1=L*np.sin(sol.T[0]),-L*np.cos(sol.T[0])
bob_radius = 0.08
line, = ax.plot([0, np.cos(theta0)], [0, np.sin(theta0)], lw=3, c='w')
circle = ax.add_patch(plt.Circle([np.cos(theta0),np.sin(theta0)], bob_radius,
                      fc='r', zorder=3))
def animate(i):
    x.append(x1[i])
    y.append(y1[i])
    circle.set_center([x[i],y[i]])
    line.set_data([0, x[i]], [0, y[i]])
    #plt.plot(x,y, scaley=True, scalex=True, color="blue",linestyle="-")
anim=animation.FuncAnimation(fig=fig,func=animate,interval=25)
anim.save("pen.gif",writer="pillow",fps=25)



    

