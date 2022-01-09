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
Solution numérique
"""

r1 = 0.0
r2 = 0.0
r0 = 0.0
Rr = 0.0
Ir = 0.0
i = 0
t = np.zeros(201)
y1 = np.ones(201)*x0
y2 = np.ones(201)*v0
z = np.zeros(201)
z1 = np.zeros(201)

Rr = -c/(2.0*m)
Ir = np.sqrt(-D)/2.0
A = (x0-A1*np.cos(Phi))
B = ((v0-A1*O*np.sin(Phi))-(x0-A1*np.cos(Phi))*Rr)/(Ir)
z[0] = A1*np.cos(O*t[0]-Phi)
z1[0] = np.exp(Rr*t[0])*(A*np.cos(Ir*t[0])+B*np.sin(Ir*t[0]))

for i in range(200):
    t[i+1]=dt*(i+1)
    j11 = y2[i]
    j12 = -(k/m)*y1[i]-(c/m)*y2[i]+(k*Amp/m)*np.cos(O*t[i])
    j21 = y2[i]+j12*dt/2.0
    j22 = -(k/m)*(y1[i]+j11*dt/2.0)-(c/m)*(y2[i]+j12*dt/2.0)+(k*Amp/m)*np.cos(O*(t[i]+dt/2.0))
    j31 = y2[i]+j22*dt/2.0
    j32 = -(k/m)*(y1[i]+j21*dt/2.0)-(c/m)*(y2[i]+j22*dt/2.0)+(k*Amp/m)*np.cos(O*(t[i]+dt/2.0))
    j41 = y2[i]+j32*dt
    j42 = -(k/m)*(y1[i]+j31*dt)-(c/m)*(y2[i]+j32*dt)+(k*Amp/m)*np.cos(O*(t[i]+dt))
    y1[i+1] = y1[i]+(dt/6.0)*(j11+2.0*j21+2.0*j31+j41)
    y2[i+1] = y2[i]+(dt/6.0)*(j12+2.0*j22+2.0*j32+j42)
    z[i+1] = A1*np.cos(O*t[i+1]-Phi)
    z1[i+1]=np.exp(Rr*t[i+1])*(A*np.cos(Ir*t[i+1])+B*np.sin(Ir*t[i+1]))

plt.plot(t,y1)
plt.plot(t,z)
plt.plot(t,z1)
plt.ylabel('Déplacement (m)')
plt.xlabel('Temps (s)')
plt.show()
