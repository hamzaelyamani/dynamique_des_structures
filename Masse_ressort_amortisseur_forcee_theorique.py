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
f = 2.0
Amp = 0.1

"""
Paramètres intermédiaires
"""

w = np.sqrt(k/m)
D = (c/m)**2 - (4*k/m)
O = 2*np.pi*f
mu = c/(2*w*m)
Dy = 1.0/(np.sqrt((1-(O/w)**2)**2+(2*mu*O/w)**2))
A1 = Dy*Amp
Phi = np.arctan((2*mu*O/w)/(1-(O/w)**2))

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
y = np.zeros(201)
z = np.zeros(201)

if (D>0.0):
    r1 = -c/(2.0*m) + np.sqrt(D)/2.0
    r2 = -c/(2.0*m) - np.sqrt(D)/2.0
    A = (((x0-A1*np.cos(Phi))*r2)-(v0-A1*O*np.sin(Phi)))/(r2-r1)
    B = ((v0-A1*O*np.sin(Phi))-((x0-A1*np.cos(Phi))*r1))/(r2-r1)
    for i in range(0, 201):
        t[i]=dt*i
        x[i]=A*np.exp(r1*t[i])+B*np.exp(r2*t[i])+A1*np.cos(O*t[i]-Phi)
        y[i]=A1*np.cos(O*t[i]-Phi)
        z[i]=A*np.exp(r1*t[i])+B*np.exp(r2*t[i])
else:
    if(D==0.0):
        r0 = -c/(2*m)
        A = (v0-A1*O*np.sin(Phi))-(x0-A1*np.cos(Phi))*r0
        B = (x0-A1*np.cos(Phi))
        for i in range(0, 201):
            t[i]=dt*i
            x[i]=(A*t[i]+B)*np.exp(r0*t[i])+A1*np.cos(O*t[i]-Phi)
            y[i]=A1*np.cos(O*t[i]-Phi)
            z[i]=(A*t[i]+B)*np.exp(r0*t[i])
    else:
        Rr = -c/(2.0*m)
        Ir = np.sqrt(-D)/2.0
        A = (x0-A1*np.cos(Phi))
        B = ((v0-A1*O*np.sin(Phi))-(x0-A1*np.cos(Phi))*Rr)/(Ir)
        for i in range(0, 201):
            t[i]=dt*i
            x[i]=np.exp(Rr*t[i])*(A*np.cos(Ir*t[i])+B*np.sin(Ir*t[i]))+A1*np.cos(O*t[i]-Phi)
            y[i]=A1*np.cos(O*t[i]-Phi)
            z[i]=np.exp(Rr*t[i])*(A*np.cos(Ir*t[i])+B*np.sin(Ir*t[i]))

plt.plot(t,x)
plt.plot(t,y)
plt.plot(t,z)
plt.ylabel('Déplacement (m)')
plt.xlabel('Temps (s)')
plt.show()
