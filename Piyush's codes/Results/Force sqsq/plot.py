import numpy as np
import matplotlib.pyplot as plt
import sys

sys.path.append("/scratch/users/ppanchal/mpltools/mpltools")
from mpltools import annotation

data = np.loadtxt('smooth_cutoff.txt')
fbx = data[:,1]
fby = data[:,2]
fvx = data[:,3]
fvy = data[:,4]
databem1 = np.loadtxt('sc1_0_0')
databem2 = np.loadtxt('sc2_0_0')
fbemx = databem1[:200,4]
fbemy = databem2[:200,4] 
fbem = np.sqrt(fbemx**2+fbemy**2)
fsgx = databem1[:200,2]
fsgy = databem2[:200,2] 
fsg = np.sqrt(fsgx**2+fsgy**2)
hbem = 28/databem1[:118,0]

hmean = data[:,7]

fb = np.sqrt(fbx**2+fby**2)
fv = np.sqrt(fvx**2+fvy**2)

# Using last value from the egg shell method as converged
#f = fv[-1]
# Using value from BEM computations
f = np.sqrt(0.113596622341585**2+0.113596622341585**2)

fig = plt.figure()
ax = fig.add_subplot(111)

plt.loglog(hmean,abs(fb-f)/f,'.-',color='blue',zorder = 2)
plt.loglog(hmean,abs(fv-f)/f,'.-',color='green',zorder = 2)
plt.loglog(hbem,abs(fbem-f)/f,'.-',color='orange',zorder = 2)
plt.loglog(hbem,abs(fsg-f)/f,'.-',color='red',zorder = 2)
#plt.legend(['Boundary formula FEM','Volume formula'])
#plt.legend(['Boundary formula FEM','Volume formula','Boundary formula BEM'])
plt.legend(['Boundary formula FEM','Volume formula','Boundary formula BEM', 'BEM shape derivative'])
plt.xlabel('h')
plt.ylabel('Relative error')

# Fitting lines to these plots
fit_bnd = np.polyfit(np.log(hmean),np.log(abs(fb-f)/abs(f)),1)
fit_vol = np.polyfit(np.log(hmean),np.log(abs(fv-f)/abs(f)),1)
fit_bem = np.polyfit(np.log(hbem),np.log(abs(fbem-f)/f),1)
fit_sg = np.polyfit(np.log(hbem),np.log(abs(fsg-f)/f),1)
y_bnd = fit_bnd[0] * np.log(hmean) + fit_bnd[1]
y_vol = fit_vol[0] * np.log(hmean) + fit_vol[1]
y_bem = fit_bem[0] * np.log(hbem) + fit_bem[1]
y_sg = fit_sg[0] * np.log(hbem) + fit_sg[1]

bestfitcolor = "silver"
plt.loglog(hmean,np.exp(y_bnd),'--',color=bestfitcolor,zorder = 1)
plt.loglog(hmean,np.exp(y_vol),'--',color=bestfitcolor,zorder = 1)
plt.loglog(hbem,np.exp(y_bem),'--',color=bestfitcolor,zorder = 1)
plt.loglog(hbem,np.exp(y_sg),'--',color=bestfitcolor,zorder = 1)

# Printing the triangle
def plot_tri(h1,h2,v1,m,plt):
	v2 = v1 * (h2/h1)**(m)
	plt.loglog(np.array([h1,h2]),np.array([v1,v1]),'k-')
	plt.loglog(np.array([h2,h2]),np.array([v1,v2]),'k-')
	plt.loglog(np.array([h1,h2]),np.array([v1,v2]),'k-')
	
# Bdry formula
annotation.slope_marker((0.03, 2e-1), (0.29, 1),ax = ax)
# Volume formula
annotation.slope_marker((0.03, 1e-3), (1.10, 1),ax = ax)
# Bem bdry formula
annotation.slope_marker((0.03, 3e-2), (0.31, 1),ax = ax)
# Bem sg formula
annotation.slope_marker((0.03, 1.2e-4), (1.36, 1),ax = ax)


print("boundary rate", fit_bnd[0])
print("eggshell rate", fit_vol[0])
print("BEM rate", fit_bem[0])
print("BEM rate", fit_sg[0])

plt.savefig('eggshell_vs_bnd.eps')
plt.show()
