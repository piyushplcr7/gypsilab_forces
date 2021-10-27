import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt('smooth_cutoff.txt')
fbx = data[:,1]
fby = data[:,2]
fvx = data[:,3]
fvy = data[:,4]

hmean = data[:,9]

fb = np.sqrt(fbx**2+fby**2)
fv = np.sqrt(fvx**2+fvy**2)

# Using last value from the egg shell method as converged
#f = fv[-1]
# Using value from BEM computations
f = np.sqrt(0.258587323659794**2+0.177566366971898**2)

plt.loglog(hmean,abs(fb-f),'.-')
plt.loglog(hmean,abs(fv-f),'.-')
plt.legend(['Boundary formula','Egg-shell formula'])
plt.xlabel('h')
plt.ylabel('Absolute error')

# Fitting lines to these plots
fit_bnd = np.polyfit(np.log(hmean),np.log(abs(fb-f)),1)
fit_vol = np.polyfit(np.log(hmean),np.log(abs(fv-f)),1)
y_bnd = fit_bnd[0] * np.log(hmean) + fit_bnd[1]
y_vol = fit_vol[0] * np.log(hmean) + fit_vol[1]

plt.loglog(hmean,np.exp(y_bnd),'--')
plt.loglog(hmean,np.exp(y_vol),'--')

print("boundary rate", fit_bnd[0])
print("eggshell rate", fit_vol[0])

plt.savefig('eggshell_vs_bnd.eps')
plt.show()
