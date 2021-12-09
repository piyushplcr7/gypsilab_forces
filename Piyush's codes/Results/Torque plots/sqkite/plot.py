import numpy as np
import matplotlib.pyplot as plt
import sys

sys.path.append("/scratch/users/ppanchal/mpltools/mpltools")
from mpltools import annotation

data = np.loadtxt('sqkg3_0_0')
step = 1
tbem = data[::step,2]
tbdry = data[::step,4]
print("tbem shape", tbem.shape)
hbem = 18/data[::step,0]
print("hbem shape", hbem.shape)

# Using last value from the egg shell method as converged
#f = fv[-1]
# Using value from BEM computations
t = -0.00647118123921061 # 5400 panels

fig = plt.figure()
ax = fig.add_subplot(111)

plt.loglog(hbem,abs(tbem-t)/abs(t),'.-',color='blue',zorder = 2)
plt.loglog(hbem,abs(tbdry-t)/abs(t),'.-',color='green',zorder = 2)
plt.legend(['Bem shape gradient','Boundary shape gradient'])
plt.xlabel('h')
plt.ylabel('Relative error')

# Fitting lines to these plots
fit_bnd = np.polyfit(np.log(hbem),np.log(abs(tbdry-t)/abs(t)),1)
fit_bem = np.polyfit(np.log(hbem),np.log(abs(tbem-t)/abs(t)),1)
y_bem = fit_bem[0] * np.log(hbem) + fit_bem[1]
y_bnd = fit_bnd[0] * np.log(hbem) + fit_bnd[1]

bestfitcolor = "silver"
plt.loglog(hbem,np.exp(y_bnd),'--',color=bestfitcolor,zorder = 1)
plt.loglog(hbem,np.exp(y_bem),'--',color=bestfitcolor,zorder = 1)

	
annotation.slope_marker((0.06, 1e-3), (2.27, 1),ax = ax)

annotation.slope_marker((6e-2,5e-2), (1.80, 1),ax = ax)


print("boundary rate", fit_bnd[0])
print("BEM rate", fit_bem[0])

plt.savefig('bem_vs_bnd_torque.eps')
plt.show()
