import numpy as np
import matplotlib.pyplot as plt


def processus_vitesse(Du,tau_v,delta_t,u0):
    normale=np.random.normal(0,1,1) 
    u1 = u0*(1-delta_t/tau_v)+np.sqrt(Du*delta_t)*normale
    return u1

def position(x0,u0,delta_t):
    x1 = x0 + u0*delta_t
    return x1

def variance(L):
    moy = np.mean(L)
    var = np.mean([(x-moy)**2 for x in L])
    return var

nu = 10**(-6)
kB = 1.38064852*10**(-23)
T = 293

D_p = 10*10**(-6)
rho_p = 1000
mp = 4/3*np.pi*(D_p/2)**3*rho_p

tau_v = mp/(3*np.pi*nu*D_p)
sigma = np.sqrt(3*kB*T/mp) 
Du = (2*sigma**2)/tau_v

var_th = sigma**2
Ec_th = 3/2*kB*T
print("var_th =", var_th)
print("Ec_th =", Ec_th)

print(tau_v)
print(Du)
print(sigma)

delta_t = tau_v/10
nb_iteration = int(100*tau_v/delta_t)
nb_realisation = 100

L_var = [0 for k in range(nb_realisation)]
L_pos_fin = [[0,0,0] for k in range(nb_realisation)]

# Simulation

for real in range(nb_realisation):
    u0 = np.array([0.00])
    ux1 = [u0 for k in range(nb_iteration)]
    uy1 = [u0 for k in range(nb_iteration)]
    uz1 = [u0 for k in range(nb_iteration)]

    pos0 = 0
    x1 = [pos0 for k in range(nb_iteration)]
    y1 = [pos0 for k in range(nb_iteration)]
    z1 = [pos0 for k in range(nb_iteration)]

    for k in range(nb_iteration-1):
        ux1[k+1]= processus_vitesse(Du,tau_v,delta_t,ux1[k])
        uy1[k+1]=processus_vitesse(Du,tau_v,delta_t,uy1[k])
        uz1[k+1]=processus_vitesse(Du,tau_v,delta_t,uz1[k])

    for k in range(nb_iteration-1):
        x1[k+1]=position(x1[k],ux1[k],delta_t)
        y1[k+1]=position(y1[k],uy1[k],delta_t)
        z1[k+1]=position(z1[k],uz1[k],delta_t)

    L_var[real] = np.var([ux1,uy1,uz1])
    L_pos_fin[real] = [x1[-1],y1[-1],z1[-1]]

# Statistiques

var = np.var([ux1,uy1,uz1])
Ec_pr = 0.5*mp*var
print("var_pr =", var)
print("Ec_pr =", Ec_pr)

L_dist = [np.sqrt(pos[0]**2+pos[1]**2+pos[2]**2) for pos in L_pos_fin]

# Courbes
    
plt.figure(1)
plt.plot([k*delta_t for k in range(nb_iteration)],ux1)

plt.figure(3)
plt.hist(L_var, bins=50)
plt.title("Histogramme des variances")
plt.xlabel("Variance")
plt.ylabel("Nombre de réalisations")

#plt.figure(4)
#plt.plot([k for i in range(nb_realisation)],L_dist,'+')

plt.show()