import numpy as np
import matplotlib.pyplot as plt

"""
Données en entrée
"""

m = 2.0
k = 1000.0
x0 = 0.1
v0 = 0.5
dt = 0.01

"""
Paramètres intermédiaires
"""

w = np.sqrt(k/m)

"""
Solution theorique
"""

i=0
t = np.zeros(201)
x = np.zeros(201)
for i in range(201):
    t[i]=dt*i
    x[i]=x0*np.cos(w*i*dt)+(v0/w)*np.sin(w*i*dt)

plt.plot(t,x)
plt.ylabel('Déplacement (m)')
plt.xlabel('Temps (s)')
plt.show()