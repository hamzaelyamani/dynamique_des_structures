import numpy as np
import matplotlib.pyplot as plt

"""
Données en entrée
"""

m = 2.0
c = 30.0
k = 1000.0
x0 = 0.05
v0 = 0.0
dt = 0.01

"""
Paramètres intermédiaires
"""

w = np.sqrt(k/m)
D = (c/m)**2 - (4*k/m)

"""
Solution theorique
"""

r1 = 0.0
r2 = 0.0
r0 = 0.0
Rr = 0.0
Ir = 0.0
A = 0.0
B = 0.0
i = 0
t = np.zeros(201)
x = np.zeros(201)

if (D>0.0):
    r1 = -c/(2.0*m) + np.sqrt(D)/2.0
    r2 = -c/(2.0*m) - np.sqrt(D)/2.0
    A = ((x0*r2)-v0)/(r2-r1)
    B = (v0-(x0*r1))/(r2-r1)
    for i in range(0, 201):
        t[i]=dt*i
        x[i]=A*np.exp(r1*t[i])+B*np.exp(r2*t[i])
else:
    if(D==0.0):
        r0 = -c/(2*m)
        A = v0-x0*r0
        B = x0
        for i in range(0, 201):
            t[i]=dt*i
            x[i]=(A*t[i]+B)*np.exp(r0*t[i])
    else:
        Rr = -c/(2.0*m)
        Ir = np.sqrt(-D)/2.0
        A = x0
        B = (v0-x0*Rr)/(Ir)
        for i in range(0, 201):
            t[i]=dt*i
            x[i]=np.exp(Rr*t[i])*(A*np.cos(Ir*t[i])+B*np.sin(Ir*t[i]))

plt.plot(t,x)
plt.ylabel('Déplacement (m)')
plt.xlabel('Temps (s)')
plt.show()
