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

databem1 = np.loadtxt('sqkt_bdry_1_0_0')
databem2 = np.loadtxt('sqkt_bdry_2_0_0')
# Slicing the bem data!
#temp1 = databem1[-3:,:]
#temp2 = databem2[-3:,:]
#temp11 = databem1[:-3:10,:]
#temp22 = databem2[:-3:10,:]
temp1 = databem1[-3:,:]
temp2 = databem2[-3:,:]
temp11 = databem1[:-3,:]
temp22 = databem2[:-3,:]
test1 = np.append(temp11,temp1,axis=0)
test2 = np.append(temp22,temp2,axis=0)
databem1 = test1
databem2 = test2

nbem = 200
fbemx = databem1[:,4]
fbemy = databem2[:,4] 
fbem = np.sqrt(fbemx**2+fbemy**2)
hbem = 36/databem1[:,0]*2
fsgx = databem1[:,2]
fsgy = databem2[:,2] 
fsg = np.sqrt(fsgx**2+fsgy**2)

hmean = data[:,9]

fb = np.sqrt(fbx**2+fby**2)
fv = np.sqrt(fvx**2+fvy**2)

# Using last value from the egg shell method as converged
#f = fv[-1]
# Using value from BEM computations
#f = np.sqrt(0.258587323659794**2+0.177566366971898**2)

# 7290 panels
f = np.sqrt(0.258587358913174**2+0.177566402431972**2)

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
fit_bnd = np.polyfit(np.log(hmean),np.log(abs(fb-f)/f),1)
fit_vol = np.polyfit(np.log(hmean),np.log(abs(fv-f)/f),1)
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

annotation.slope_marker((2e-2, 2.5e-2), (0.43, 1),ax = ax)
annotation.slope_marker((2e-2, 3e-6), (2.28, 1),ax = ax)
annotation.slope_marker((3.9e-2, 2e-3), (1.85, 1),invert = True, ax = ax)
annotation.slope_marker((2e-2, 2e-8), (3.03, 1),ax = ax)

print("boundary rate", fit_bnd[0])
print("eggshell rate", fit_vol[0])
print("BEM rate", fit_bem[0])
print("BEM rate", fit_sg[0])

plt.savefig('eggshell_vs_bnd.eps')
plt.show()
